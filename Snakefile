# SNIP-R Snakemake workflow.
#
# Provides three target rules:
#   prep_flanks       Part B step 1: regions BED -> FASTA / coords for CRISPick
#   pair_grnas        Part B step 2: CRISPick designs -> paired sgRNAs
#   screen_analysis   Part C: counts matrix + config -> Z-scores + plots
#
# Configuration is read from `config.yaml` in the run directory, e.g.:
#
#   # config.yaml
#   regions: example_data/Finaldf.bed.txt
#   reference: example_data/chr12.fa            # optional if skip_fasta: true
#   crispick_left: CRISPick/Left-sgrna-designs.txt
#   crispick_right: CRISPick/Right-sgrna-designs.txt
#   counts: example_data/SupplementaryTable4_RawScreenData.xlsx
#   screen_config: config/example_partC_config.yaml
#   out_dir: output
#   skip_fasta: false                           # set true to skip the bedtools FASTA-extraction step
#
# Usage:
#   snakemake -j1 prep_flanks
#   # then run CRISPick manually, save outputs in CRISPick/...
#   snakemake -j1 pair_grnas screen_analysis
#
# Inside the Docker image, the wrapper scripts are at /opt/snipr/scripts/.
# When invoked outside Docker, set SNIPR_SCRIPTS to that directory or place
# the scripts/ folder beside the Snakefile.

import os

configfile: "config.yaml"

SCRIPTS = os.environ.get("SNIPR_SCRIPTS", "scripts")
OUT     = config.get("out_dir", "output")

rule all:
    input:
        f"{OUT}/paired_grnas.tsv",
        f"{OUT}/partC_done.flag"

# ---- Part B step 1 --------------------------------------------------------
rule prep_flanks:
    """Prepare flanking FASTA + coordinate files for manual CRISPick run.

    If `skip_fasta: true` is set in config.yaml, the bedtoolsr/bedtools step
    is skipped: only coordinate files are written, and `reference` does not
    need to be provided. Empty placeholder .fa files are created so the rule's
    declared outputs all exist.
    """
    input:
        regions = config.get("regions", "")
    output:
        conversion = f"{OUT}/ROItoFlank.tsv",
        left_fa    = f"{OUT}/left.fa",
        right_fa   = f"{OUT}/right.fa"
    params:
        shrink_cutoff   = config.get("shrink_cutoff", 500),
        shrink_amount   = config.get("shrink_amount", 50),
        flank_size      = config.get("flank_size", 550),
        reference_arg   = (f"--reference {config['reference']}"
                           if config.get("reference") else ""),
        skip_fasta_arg  = "--skip-fasta" if config.get("skip_fasta", False) else ""
    shell:
        """
        Rscript {SCRIPTS}/partB_prep_flanks.R \
            --regions {input.regions} \
            {params.reference_arg} \
            --out-dir {OUT} \
            --shrink-cutoff {params.shrink_cutoff} \
            --shrink-amount {params.shrink_amount} \
            --flank-size {params.flank_size} \
            {params.skip_fasta_arg}
        # If --skip-fasta was used, the R script doesn't produce .fa files;
        # create empty placeholders so Snakemake's output-existence check
        # is satisfied. A real run produces non-empty .fa files.
        [ -e {output.left_fa} ]  || : > {output.left_fa}
        [ -e {output.right_fa} ] || : > {output.right_fa}
        """

# ---- Part B step 2 --------------------------------------------------------
rule pair_grnas:
    """Score CRISPick outputs and pair sgRNAs across left/right flanks."""
    input:
        left       = config.get("crispick_left", ""),
        right      = config.get("crispick_right", ""),
        conversion = f"{OUT}/ROItoFlank.tsv"
    output:
        pairs = f"{OUT}/paired_grnas.tsv"
    params:
        top_n     = config.get("top_n_per_flank", 5),
        on_w      = config.get("on_weight",  0.667),
        off_w     = config.get("off_weight", 0.333),
        pool      = config.get("pairs_pool", 25),
        per_roi   = config.get("pairs_per_roi", 10),
        max_reuse = config.get("max_grna_reuse", 3)
    shell:
        """
        Rscript {SCRIPTS}/partB_pair_grnas.R \
            --left-crispick {input.left} \
            --right-crispick {input.right} \
            --conversion-table {input.conversion} \
            --out {output.pairs} \
            --top-n-per-flank {params.top_n} \
            --on-weight {params.on_w} \
            --off-weight {params.off_w} \
            --pairs-pool {params.pool} \
            --pairs-per-roi {params.per_roi} \
            --max-grna-reuse {params.max_reuse}
        """

# ---- Part C ---------------------------------------------------------------
rule screen_analysis:
    """Run Part C screen analysis: per-pair / per-region / pan-donor + plots."""
    input:
        counts        = config.get("counts", ""),
        screen_config = config.get("screen_config", "")
    output:
        flag = f"{OUT}/partC_done.flag"
    shell:
        """
        Rscript {SCRIPTS}/partC_screen_analysis.R \
            --counts {input.counts} \
            --config {input.screen_config} \
            --out-dir {OUT}
        touch {output.flag}
        """
