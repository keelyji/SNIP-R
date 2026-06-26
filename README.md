# SNIP-R

SNIP-R (*Systematic Non-coding region Interrogation by Paired sgRNAs*) is a reproducible pipeline for designing paired-sgRNA deletion libraries and analyzing cis-regulatory element CRISPR screens in primary T cells. This repository contains the code used in the manuscript:

*Scalable hit-and-run platform for enhancer deletion reveals state-specific IFNG regulation in primary human T cells.*

## Overview

This repository is organized into two components:

1. **Reusable command-line pipeline** for preparing target-region flanks, pairing CRISPick sgRNA designs, and analyzing SNIP-R screen count matrices.
2. **Paper-reproduction notebooks and supporting files** for reproducing the IFNG screen-region selection, published analyses, ATAC-seq processing notes, Hi-C visualization, and supporting BED files.

## Workflow

For a new SNIP-R screen, the recommended workflow is:

| Step     | Script or folder                                | Purpose                                                                                                                                                    |
| -------- | ----------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Part A   | `notebooks/PartA_SNIP-R_ScreenRegionSelection/` | Select candidate target regions and generate a BED file. This step is study- and locus-specific, so it is provided as a reproducibility notebook/template. |
| Part B-1 | `scripts/partB_prep_flanks.R`                   | Convert target regions into left/right flanking windows for CRISPick.                                                                                      |
| Manual   | CRISPick web tool                               | Design sgRNAs for left and right flanking windows.                                                                                                         |
| Part B-2 | `scripts/partB_pair_grnas.R`                    | Pair and filter CRISPick sgRNA designs into SNIP-R paired-sgRNA libraries.                                                                                 |
| Part C   | `scripts/partC_screen_analysis.R`               | Analyze SNIP-R screen count matrices using a YAML configuration file.                                                                                      |

Parts B-1, B-2, and C are standalone command-line scripts and can be run with Docker or Snakemake. Part A is provided as an R Markdown notebook because candidate-region selection depends on the locus, cell type, chromatin-accessibility datasets, and biological question.

## Quickstart with Docker

Build the Docker image and run the bundled demo:

```bash
git clone https://github.com/keelyji/SNIP-R.git
cd SNIP-R

docker build -t snipr:latest .
docker run --rm -v "$PWD":/work -w /work snipr:latest bash tests/run_demo.sh
```

The demo uses bundled IFNG example data and runs the reusable SNIP-R pipeline through flank preparation, sgRNA-pairing, and screen analysis.

## Repository layout

```text
SNIP-R/
├── README.md
├── INSTALL.md
├── Dockerfile
├── Snakefile
├── LICENSE
├── scripts/
│   ├── partB_prep_flanks.R
│   ├── partB_pair_grnas.R
│   └── partC_screen_analysis.R
├── config/
│   └── example_partC_config.yaml
├── example_data/
│   ├── target_regions/
│   │   └── PartA_Output_Finaldf.bed.txt
│   ├── CRISPick/
│   │   ├── Left-sgrna-designs.txt
│   │   └── Right-sgrna-designs.txt
│   └── screen_counts/
│       └── SupplementaryTable4_RawScreenData.xlsx
├── notebooks/
│   ├── README.md
│   ├── PartA_SNIP-R_ScreenRegionSelection/
│   ├── PartB_SNIP-R_sgRNApairing/
│   ├── PartC_SNIP-R_ScreenAnalysis/
│   ├── PartD_other_code/
│   └── PartE_other_datasets/
└── tests/
    └── run_demo.sh
```

## Part A — Selecting candidate target regions

Part A documents how candidate IFNG regulatory regions were selected for the published SNIP-R screen. The notebook integrates chromatin-accessibility datasets, filters candidate regions within the IFNG locus, removes unsuitable targets, and produces the target-region BED file used as input for sgRNA design.

The final Part A output used for the IFNG demo is included at:

```text
example_data/target_regions/PartA_Output_Finaldf.bed.txt
```

Users applying SNIP-R to a new locus or cell type can replace this file with any BED file of target regions containing at least three columns:

```text
chr    start    end
```

An optional fourth column may be included as the target-region name.

## Part B-1 — Prepare flanking regions for CRISPick

Input:

* BED file of target regions.
* Reference genome FASTA, such as `hg19.fa`.

Example:

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest \
  Rscript /opt/snipr/scripts/partB_prep_flanks.R \
    --regions example_data/target_regions/PartA_Output_Finaldf.bed.txt \
    --reference /path/to/hg19.fa \
    --out-dir output/flanks/
```

Outputs:

* `left.fa`
* `right.fa`
* `left_coords.tsv`
* `right_coords.tsv`
* `ROItoFlank.tsv`

The left and right FASTA or coordinate files are used as inputs for CRISPick.

## Manual step — Run CRISPick

CRISPick is maintained by the Broad Institute. Open the CRISPick web tool and use the left and right flank files generated by Part B-1 as two separate design jobs.

For the IFNG screen, the following settings were used:

* Reference genome: Human GRCh37 / hg19
* Mechanism: CRISPRko
* Enzyme: SpyoCas9, Chen tracrRNA, RuleSet3
* Quota: 20

Save the CRISPick outputs as:

```text
CRISPick/Left-sgrna-designs.txt
CRISPick/Right-sgrna-designs.txt
```

Example CRISPick outputs are bundled in:

```text
example_data/CRISPick/
```

## Part B-2 — Pair and filter sgRNAs

Example:

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest \
  Rscript /opt/snipr/scripts/partB_pair_grnas.R \
    --left-crispick example_data/CRISPick/Left-sgrna-designs.txt \
    --right-crispick example_data/CRISPick/Right-sgrna-designs.txt \
    --conversion-table output/flanks/ROItoFlank.tsv \
    --out output/paired_grnas.tsv
```

This generates a table of paired sgRNAs, including left and right sgRNA sequences, on-target efficacy scores, paired scores, reuse filtering, and BsmBI-site flags for cloning compatibility.

## Part C — SNIP-R screen analysis

Screen analysis is controlled by a YAML configuration file describing donors, conditions, controls, baseline columns, derived buckets, plotting labels, and pan-donor analysis settings.

An annotated example configuration is provided at:

```text
config/example_partC_config.yaml
```

Example:

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest \
  Rscript /opt/snipr/scripts/partC_screen_analysis.R \
    --counts example_data/screen_counts/SupplementaryTable4_RawScreenData.xlsx \
    --config config/example_partC_config.yaml \
    --out-dir output/screen_analysis/
```

Main outputs include:

* per-pair Z-score tables
* per-region Z-score tables
* pan-donor hit tables
* volcano plots
* donor concordance plots
* distribution plots for selected control and target categories

## Snakemake workflow

A Snakemake workflow is included for chaining the command-line steps.

Example configuration:

```yaml
regions: example_data/target_regions/PartA_Output_Finaldf.bed.txt
reference: /path/to/hg19.fa
crispick_left: example_data/CRISPick/Left-sgrna-designs.txt
crispick_right: example_data/CRISPick/Right-sgrna-designs.txt
counts: example_data/screen_counts/SupplementaryTable4_RawScreenData.xlsx
screen_config: config/example_partC_config.yaml
out_dir: output
```

Run:

```bash
snakemake -j1 prep_flanks
# Run CRISPick manually if using new target regions.
snakemake -j1 pair_grnas screen_analysis
```

Snakemake is also available inside the Docker image:

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest snakemake -j1 all
```

## Running without Docker

For local installation instructions, including R 4.4.x setup, system dependencies, bedtools, and required R packages, see:

```text
INSTALL.md
```

## Reproducing the manuscript analyses

The `notebooks/` directory contains the R Markdown notebooks and supporting files used to reproduce the manuscript analyses:

```text
notebooks/
├── PartA_SNIP-R_ScreenRegionSelection/
├── PartB_SNIP-R_sgRNApairing/
├── PartC_SNIP-R_ScreenAnalysis/
├── PartD_other_code/
└── PartE_other_datasets/
```

Some raw public datasets, including dbGaP-controlled ATAC-seq data and Hi-C data from GEO, must be downloaded separately using the accessions listed in the manuscript and notebook documentation.

## Notes on CRISPick versions

CRISPick output formats may change over time. The bundled pairing script was written for the CRISPick output schema used for this repository. If CRISPick returns a different column structure, users may need to rename columns or open a GitHub issue.

## Contact

For paper-related questions, please contact the corresponding author. For code or reproducibility issues, please open a GitHub issue.

## License

MIT — see `LICENSE`.
