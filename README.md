# SNIP-R

A reproducible pipeline for designing and analyzing cis-regulatory element CRISPR screens in primary T cells, as used in the paper *"Scalable hit-and-run platform for enhancer deletion reveals state-specific IFNG regulation in primary human T cells."*

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

<<<<<<< HEAD
```bash
git clone https://github.com/keelyji/SNIP-R.git
cd SNIP-R
docker build -t snipr:latest .
docker run --rm -v "$PWD":/work -w /work snipr:latest bash tests/run_demo.sh
=======

## Part B: To design sgRNA pairs for SNIP-R screen
1. All code for sgRNA pairing are performed in R and R studio, while the gRNA design is performed in Broad Institute's CRISPick. Please note that input and output of CRISPick may change overtime as it is updated by Broad Institute, so the result output might vary and code may need adjustment according to the CRISPick version used. We updated the sgRNA pairing code to reflect the latest CRISPick update as of February 2026 for this repository. In general, the following CRISPick setting was used: Human GRCh37 (hg19), CRISPRko, SpyoCas9 (Chen tracrRNA, RuleSet3), CRISPickQuota: 20. Guide designs were provided by the CRISPick web tool of the GPP at the Broad Institute (Sanson et al. 2018, Doench et al. 2016)

- 01_SNIP-R_Pairing_V2.Rmd file contains all the code necessary for this step. The exact R package and version used in the paper are noted in the respective html file.
2. Necessary dataset are bed file of target regions (e.g., output from Part B).
3. Expected time: demo dataset (1 hour). Expected output for demo: PartB_Output.txt

## Part C: To analyse the result of SNIP-R screen
1. All code for target region selection are performed in R and R studio. Key packages used in this steps are: readr, tidyverse, readxl, ggrepel, and patchwork.
- 02_ScreenAnalysis.Rmd file contains all the code necessary for this step.The exact R package and version used in the paper are noted in the respective html file.
2. Necessary dataset are:
- Screen sequencing data in count matrix form. We provided a counts matrix called "SupplementaryTable4_RawScreenData.xlsx" available in our paper and in the SNIP-R screen analysis folder located in this repository.
3. Analysis flowchart and details on the analysis are noted in our paper.
4. We also provided code to replicate main figures in our paper.
5. Expected time: demo dataset (1 hour). Expected output for demo: 02_ScreenAnalysis.html.

## Part D: Other codes
- Contains ATAC-seq script used for the paper
- Contains R markdown file for Hi-C visualization with Plotgardener. All dataset (except hi-c data) is available in PartE other datasets. Requirement: Hi-C dataset analyzed from Bediaga et al. (GEO: GSE126117). Expected output: HiCplot_expectedOutput.pdf.
  
## Part E: Other datasets
- Contains bed files of all regulatory elements screened IFNG SNIP-R screen and the identified regulatory elements from the screen

# Contact
Please contact corresponding author for any questions, comments, or concerns regarding the paper in general. For issues with code/reproducibility, please open a GitHub issue.


# Compute Environment and Version

```
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 15.7.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] bedtoolsr_2.30.0-6 patchwork_1.3.0    ggrepel_0.9.6      lubridate_1.9.3   
 [5] forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2       
 [9] tidyr_1.3.1        tibble_3.3.1       ggplot2_4.0.1      tidyverse_2.0.0   
[13] readxl_1.4.3       readr_2.1.6       

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       compiler_4.4.1     Rcpp_1.0.13        tidyselect_1.2.1  
 [5] scales_1.4.0       yaml_2.3.10        fastmap_1.2.0      R6_2.5.1          
 [9] generics_0.1.3     knitr_1.48         pillar_1.9.0       RColorBrewer_1.1-3
[13] tzdb_0.4.0         rlang_1.1.4        utf8_1.2.4         stringi_1.8.4     
[17] xfun_0.47          S7_0.2.0           timechange_0.3.0   cli_3.6.3         
[21] withr_3.0.1        magrittr_2.0.3     digest_0.6.37      grid_4.4.1        
[25] rstudioapi_0.16.0  hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5       
[29] evaluate_1.0.0     glue_1.8.0         farver_2.1.2       cellranger_1.1.0  
[33] fansi_1.0.6        rmarkdown_2.28     tools_4.4.1        pkgconfig_2.0.3   
[37] htmltools_0.5.8.1 
>>>>>>> b31344848305db5cb1b90242956dab1e42f80b1d
```

This downloads the repo, builds the image (~10 min the first time), and executes the demo against the bundled IFNG screen data. The demo runs all three pipeline steps and checks that the IFNG promoter is recovered as a top hit in the pan-donor analysis. Total runtime after the image is built: ~2 minutes.

---

## Running on your own data

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
