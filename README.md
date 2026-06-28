# SNIP-R

A reproducible pipeline for designing and analyzing cis-regulatory element CRISPR screens in primary T cells, as used in the paper *"Scalable platform for enhancer deletion reveals state-specific IFNG regulation in primary human T cells."*

The pipeline has three executable steps:

| Step | Script                           | Purpose                                                                            |
|------|----------------------------------|------------------------------------------------------------------------------------|
| B-1  | `partB_prep_flanks.R`            | Take a BED of target regions, emit flanking FASTA + coordinates for CRISPick      |
| B-2  | `partB_pair_grnas.R`             | Take CRISPick design outputs, produce scored & filtered left/right sgRNA pairs    |
| C    | `partC_screen_analysis.R`        | Take a counts matrix + a YAML describing your screen, produce Z-scored hit tables and plots |

Each script is a standalone command-line tool that accepts input/output paths and parameters via flags. All three are packaged in a Docker image; the same scripts can also be chained with the included Snakemake workflow.

The `notebooks/` directory contains the original R Markdown documents that produced the published figures (Parts A–E). They are kept for transparency but the CLI scripts and Docker image are the recommended interface for running the pipeline on new data.

---

## Quickstart (Docker)

The fastest way to confirm the pipeline runs on your machine:

```bash
git clone https://github.com/keelyji/SNIP-R.git
cd SNIP-R
docker build -t snipr:latest .
docker run --rm -v "$PWD":/work -w /work snipr:latest bash tests/run_demo.sh
```

This downloads the repo, builds the image (~10 min the first time), and executes the demo against the bundled IFNG screen data. The demo runs all three pipeline steps and checks that the IFNG promoter is recovered as a top hit in the pan-donor analysis. Total runtime after the image is built: ~2 minutes.

---

## Running on your own data

Part A, which documents candidate target-region selection for the published IFNG screen, is available in `notebooks/PartA_SNIP-R_ScreenRegionSelection/`. This step is provided as an R Markdown notebook for transparency and reproducibility. For new screens, users can supply their own BED file of target regions and begin with Step B-1 below.

### Step B-1 — Prepare flanking regions for CRISPick

Inputs:
- A BED file of target regions (columns: `chr`, `start`, `end`, optional `name`).
- A reference genome FASTA (e.g. `hg19.fa`).

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest \
  Rscript /opt/snipr/scripts/partB_prep_flanks.R \
    --regions   my_regions.bed \
    --reference /path/to/hg19.fa \
    --out-dir   output/
```

Outputs (in `output/`):
- `left.fa`, `right.fa` — FASTA sequences to paste into CRISPick
- `left_coords.tsv`, `right_coords.tsv` — coordinate form for CRISPick
- `ROItoFlank.tsv` — required for step B-2

Tunable parameters:

| Flag | Default | Meaning |
|------|---------|---------|
| `--shrink-cutoff` | `500` | Regions larger than this (bp) are shrunk inward |
| `--shrink-amount` | `50`  | Amount (bp) to shrink each side of large regions |
| `--flank-size`    | `550` | Length (bp) of left/right windows |
| `--skip-fasta`    | off   | Skip FASTA extraction if `bedtools` is unavailable |

### Manual step — Run CRISPick

CRISPick is a web tool maintained by the Broad Institute. Open <https://portals.broadinstitute.org/gppx/crispick/public>, set:

- Reference: Human GRCh37 (hg19)
- Mechanism: CRISPRko
- Enzyme: SpyoCas9 (Chen tracrRNA, RuleSet3)
- Quota: 20

Submit `left.fa` and `right.fa` as two separate jobs. Save the design outputs as `CRISPick/Left-sgrna-designs.txt` and `CRISPick/Right-sgrna-designs.txt`.

> **Note on CRISPick versions.** The CRISPick output schema changes occasionally. The `partB_pair_grnas.R` script handles the Feb 2026 schema. If you see "missing column" errors, your CRISPick output is a different version — open an issue.

### Step B-2 — Pair and filter sgRNAs

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest \
  Rscript /opt/snipr/scripts/partB_pair_grnas.R \
    --left-crispick     CRISPick/Left-sgrna-designs.txt \
    --right-crispick    CRISPick/Right-sgrna-designs.txt \
    --conversion-table  output/ROItoFlank.tsv \
    --out               output/paired_grnas.tsv
```

Output is a TSV with one row per sgRNA pair, columns for left/right sgRNA sequences, on-target efficacy scores, summed scores, and flags indicating whether either sgRNA contains a BsmBI Golden Gate cut site (which would interfere with cloning).

Tunable parameters:

| Flag | Default | Meaning |
|------|---------|---------|
| `--top-n-per-flank` | `5`     | Top N sgRNAs retained per flank by weighted score |
| `--on-weight`       | `0.667` | Weight on on-target rank in the weighted score |
| `--off-weight`      | `0.333` | Weight on off-target (CFD) rank |
| `--pairs-pool`      | `25`    | Top-K left×right pairs considered per ROI before reuse filtering |
| `--pairs-per-roi`   | `10`    | Final number of pairs per ROI |
| `--max-grna-reuse`  | `3`     | Max times any single sgRNA appears in an ROI's pairs |
| `--bsmbi-sites`     | `CGTCTC,GAGACG` | Motifs flagged as Golden Gate cut sites |

### Step C — Screen analysis

This step is fully config-driven. Write a YAML file describing your sample design (donors, conditions, controls, baseline column, etc.); the example at `config/example_partC_config.yaml` is fully annotated and reproduces the published IFNG analysis.

The config covers, at minimum:

- `donors`: list of `{prefix, label}` — one entry per donor in your matrix. Each donor may optionally have a `samples_exclude` list naming sample columns to skip (e.g. for a condition that failed QC in one donor; the demo config excludes donor 1's chronic samples for exactly this reason)
- `samples`: sample column names *without* the donor prefix
- `id_columns`: non-numeric identifier columns
- `input_filter_columns`: columns whose zero values cause a row to be dropped
- `baseline_column`: column subtracted (in log space) from every other column
- `zscore_control`: SimpleName value(s) used as the Z-score null distribution
- `zscore_conditions`: list of conditions to Z-score
- `computed_columns`: derived columns built by summing components (e.g. `Acute_Low = Acute_Neg + Acute_M1`)
- `center_filter`: parameters for trimming per-region distributions
- `pandonor`: conditions for cross-donor analysis, plus highlight categories, palette, volcano-label thresholds (`label_p_max`, `label_z_min`), `volcano_exclude` (drop control/positive-control SimpleNames from volcano plots), and `display_names` (rename volcano plots/files for readability)
- `multitrack_distribution`: list of regions (by SimpleName) to plot as a multi-row distribution figure

Run:

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest \
  Rscript /opt/snipr/scripts/partC_screen_analysis.R \
    --counts  my_counts.xlsx \
    --config  my_screen.yaml \
    --out-dir output/
```

Outputs (per condition / per donor where applicable):
- `per_pair_zscores_<donor>.tsv` — Z-scored gRNA pairs
- `per_pair_zscores_centered_<donor>.tsv` — center-N filtered version
- `per_region_zscores_<donor>.tsv` — per-region means
- `pandonor_<condition>.tsv` — pan-donor mean Z + t-test p-value per region
- `volcano_<name>.pdf` — pan-donor volcano plot. The filename uses `pandonor.display_names` if set (e.g. `volcano_Acute.pdf`, `volcano_Chronic.pdf`); otherwise it falls back to the canonical condition name (e.g. `volcano_Acute_Low.pdf`)
- `concordance_<refcondition>.pdf` — all donor-pair concordance plots
- `distribution_<refcondition>_<category>.pdf` — one per category in `output.plot_distribution_categories` (e.g. `NonTargeting`, `ClosedChromatinRegion`, `Promoter`), each showing per-donor Z-score distributions with the named category highlighted
- `multitrack_distribution_<refcondition>.pdf` — one row per region in `multitrack_distribution.regions`, with that region's per-donor Z-scores highlighted against the full distribution

---

## Snakemake (optional)

If you'd rather chain the steps:

```yaml
# config.yaml
regions:        my_regions.bed
reference:      /ref/hg19.fa
crispick_left:  CRISPick/Left-sgrna-designs.txt
crispick_right: CRISPick/Right-sgrna-designs.txt
counts:         my_counts.xlsx
screen_config:  my_screen.yaml
out_dir:        output
```

```bash
snakemake -j1 prep_flanks
# ...run CRISPick manually, save designs...
snakemake -j1 pair_grnas screen_analysis
```

Snakemake is bundled in the Docker image:

```bash
docker run --rm -v "$PWD":/work -w /work snipr:latest snakemake -j1 all
```

---

## Running without Docker

If you prefer to work in your existing R or RStudio environment, see [`INSTALL.md`](INSTALL.md) for step-by-step instructions covering R 4.4.x setup, system dependencies (bedtools), R package installation, and how to run each script from the command line or interactively from RStudio.

---

## Reproducing the paper figures

The notebooks under `notebooks/` reproduce the published figures exactly. Each notebook bundles its own input data, except for ATAC-seq raw data (dbGaP study accession `phs002510.v1.p1`) and Hi-C data (GEO `GSE126117`), which must be obtained separately. See `notebooks/README.md` for details.

The CLI pipeline produces equivalent quantitative results from the same inputs.

---

## Repository layout

```
.
├── README.md
├── INSTALL.md                       (running locally without Docker)
├── Dockerfile
├── Snakefile
├── LICENSE                          (MIT)
├── scripts/
│   ├── partB_prep_flanks.R
│   ├── partB_pair_grnas.R
│   ├── partC_screen_analysis.R
├── config/
│   └── example_partC_config.yaml    (annotated, reproduces the paper's analysis)
├── example_data/                    (bundled inputs for the demo)
│   ├── Finaldf.bed.txt
│   ├── CRISPick/
│   └── SupplementaryTable4_RawScreenData.xlsx
├── notebooks/                       (original R Markdown — paper figures)
│   ├── PartA_SNIP-R_ScreenRegionSelection/
│   ├── PartB_SNIP-R_sgRNApairing/
│   ├── PartC_SNIP-R_ScreenAnalysis/
│   ├── PartD_other_code/
│   └── PartE_other_datasets/
└── tests/
    └── run_demo.sh                  (end-to-end smoke test)
```

---

## Contact

For paper-related questions, please contact the corresponding author. For code/reproducibility issues, open a GitHub issue at <https://github.com/keelyji/SNIP-R/issues>.

## License

MIT — see [`LICENSE`](LICENSE).
