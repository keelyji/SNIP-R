# Running SNIP-R locally (without Docker)

This guide covers running the pipeline directly from R or RStudio. The Docker route (see `README.md` quickstart) is still the simplest path if you don't already have R set up — these instructions are for users who prefer working in their existing R environment.

## 1. Install system prerequisites

### R 4.4.1

The pipeline was developed and tested on R 4.4.1. Other 4.4.x and 4.5.x versions should work; older versions may fail because of newer dplyr/tidyselect syntax.

- **macOS / Windows**: download from <https://cran.r-project.org/> (CRAN binaries for 4.4.1 are archived at <https://cran.r-project.org/bin/macosx/big-sur-arm64/base/> for ARM Macs, and <https://cran.r-project.org/bin/windows/base/old/4.4.1/> for Windows).
- **Linux**:
  ```bash
  # Ubuntu/Debian
  sudo apt-get update && sudo apt-get install -y r-base
  ```

Confirm: `R --version` should show 4.4.x.

### RStudio (optional)

If you prefer a GUI, install RStudio Desktop from <https://posit.co/download/rstudio-desktop/>. It uses whatever R installation it finds.

### bedtools (needed only for Part B FASTA extraction)

Required for `partB_prep_flanks.R` unless you pass `--skip-fasta`.

- **macOS** (Homebrew): `brew install bedtools`
- **Linux**: `sudo apt-get install -y bedtools`
- **Windows**: bedtools is not natively supported. Either use WSL (Windows Subsystem for Linux) with the Linux instructions, or run `partB_prep_flanks.R` with `--skip-fasta` and obtain the FASTAs separately (e.g. via UCSC's online "Get DNA" tool against your coordinate files).

Confirm: `bedtools --version`.

## 2. Install R packages

From a terminal:

```bash
Rscript -e 'install.packages(c(
  "optparse", "yaml", "readr", "readxl", "dplyr", "tidyr",
  "ggplot2", "ggrepel", "patchwork", "bedtoolsr"),
  repos = "https://cloud.r-project.org")'
```

Or, from inside R / RStudio:

```r
install.packages(c(
  "optparse", "yaml", "readr", "readxl", "dplyr", "tidyr",
  "ggplot2", "ggrepel", "patchwork", "bedtoolsr"))
```

Notes:

- `bedtoolsr` (from PhanstielLab) wraps the `bedtools` command-line binary. If it isn't on CRAN for your platform, install from GitHub:
  ```r
  install.packages("remotes")
  remotes::install_github("PhanstielLab/bedtoolsr")
  ```
- If a package fails to install on macOS due to missing system libraries, install Xcode Command Line Tools (`xcode-select --install`) and retry.
- On Linux you may need `sudo apt-get install -y libxml2-dev libssl-dev libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev` before the R packages can compile.

## 3. Get the pipeline

```bash
git clone https://github.com/keelyji/SNIP-R.git
cd SNIP-R
```

## 4. Run the demo

From a terminal:

```bash
bash tests/run_demo.sh
```

You should see `=== ALL CHECKS PASSED ===` at the end. Total runtime: ~30 seconds. Outputs appear under a temp directory printed in the last line.

## 5. Running the scripts on your own data

### From the command line

```bash
# Part B step 1 — prepare flanks
Rscript scripts/partB_prep_flanks.R \
  --regions   my_regions.bed \
  --reference /path/to/hg19.fa \
  --out-dir   output/

# (Manually run CRISPick — see README.md)

# Part B step 2 — pair sgRNAs
Rscript scripts/partB_pair_grnas.R \
  --left-crispick    CRISPick/Left-sgrna-designs.txt \
  --right-crispick   CRISPick/Right-sgrna-designs.txt \
  --conversion-table output/ROItoFlank.tsv \
  --out              output/paired_grnas.tsv

# Part C — screen analysis
Rscript scripts/partC_screen_analysis.R \
  --counts  my_counts.xlsx \
  --config  my_screen.yaml \
  --out-dir output/
```

Append `--help` to any script for the full list of options.

### From inside RStudio

The three CLI scripts are plain R files — you can also run them interactively in RStudio. Two approaches:

**Approach A: source with arguments (simplest).** Open `scripts/partC_screen_analysis.R` in RStudio. The `optparse` block at the top reads `commandArgs(trailingOnly = TRUE)`, so set them before sourcing:

```r
# In the RStudio console, with the SNIP-R/ folder as your working directory:
setwd("~/path/to/SNIP-R")

commandArgs <- function(...) c(
  "--counts",  "example_data/SupplementaryTable4_RawScreenData.xlsx",
  "--config",  "config/example_partC_config.yaml",
  "--out-dir", "output/"
)
source("scripts/partC_screen_analysis.R")
```

The script will run end-to-end inside your R session, leaving all intermediate objects (`norm_frames`, `pd_results`, etc.) available in the global environment for interactive inspection.

**Approach B: use RStudio's terminal.** RStudio includes a Terminal tab (Tools → Terminal → New Terminal). Use the command-line invocations from section 5 above. Outputs land in your project directory and can be opened with `File → Open File`.

**Approach C: convert any script to a notebook.** The CLI scripts are intentionally linear (no functions hiding the pipeline). You can copy the body of `partC_screen_analysis.R` into a fresh `.Rmd`, replace the `opt <- parse_args(...)` block with hard-coded values, and run chunk by chunk for development. The original paper notebooks under `notebooks/` are organized this way and may be a better starting point if your goal is exploration rather than batch analysis.

## 6. Adapting the Part C config to your screen

The single most important file when running on your own data is the Part C YAML config. `config/example_partC_config.yaml` is fully commented and reproduces the IFNG screen analysis from the paper. Copy it and edit:

```bash
cp config/example_partC_config.yaml config/my_screen.yaml
# edit my_screen.yaml in your favourite editor
```

Key fields to adjust:

| Field | What to change |
|-------|----------------|
| `input.counts_file` | Path to your counts matrix (xlsx/csv/tsv) |
| `donors` | List of `{prefix, label}` — one per donor in your matrix (any number is fine) |
| `donors[].samples_exclude` | *(optional)* For any donor, list samples that should be skipped (e.g. a condition that failed QC for that donor). The demo config uses this to exclude donor 1's chronic-stim samples; the pan-donor analysis then automatically uses only donors 2 and 3 for chronic conditions. |
| `samples` | Column names *without* the donor prefix that should be processed |
| `id_columns` | Non-numeric identifier columns to carry through |
| `input_filter_columns` | Columns whose zero values should drop a row |
| `baseline_column` | Column subtracted (in log space) from every other column |
| `zscore_control` | SimpleName value identifying control sgRNAs (default: `NonTargeting`) |
| `zscore_conditions` | Conditions to compute Z-scores for |
| `computed_columns` | Optional summed bucket columns (e.g. `Acute_Low = Acute_Neg + Acute_M1`) |
| `pandonor.conditions` | Conditions to include in the pan-donor analysis |
| `pandonor.display_names` | *(optional)* Map condition → short name for volcano plot titles and filenames (e.g. `Acute_Low → Acute`) |
| `pandonor.volcano_exclude` | *(optional)* SimpleName values to drop from volcano plots (e.g. positive controls like `sgIFNG`) |
| `pandonor.label_p_max`, `label_z_min` | Volcano label thresholds — a region is labeled if `pval < label_p_max` AND `|meanZ| > label_z_min` |
| `pandonor.palette` / `pandonor.highlight_categories` | Plot colors and per-category SimpleName lists |
| `multitrack_distribution.regions` | List of SimpleNames to draw as a multi-row distribution plot |

The script validates the config against your counts matrix and will fail loudly with a list of missing columns if names don't match.

## 7. Troubleshooting

| Symptom | Likely cause / fix |
|---------|--------------------|
| `bedtoolsr is not installed` | `install.packages("bedtoolsr")`, or run `partB_prep_flanks.R` with `--skip-fasta` and supply FASTAs another way |
| `bedtools: command not found` (raised by bedtoolsr) | Install the bedtools binary (section 1 above) — bedtoolsr is just a wrapper |
| `Left CRISPick file missing column: <name>` | CRISPick output schema changed. The script expects: `Input`, `sgRNA.Sequence`, `On.Target.Efficacy.Score`, `On.Target.Rank`, `Aggregate.CFD.Score`. Open an issue with the version of CRISPick you used. |
| `Counts matrix is missing required columns` | The column names in your counts file don't match `<donor.prefix>_<sample>` — update either the file or the YAML config. If the columns are deliberately missing for a known QC reason (e.g. a chronic-stim arm that failed for one donor), add those sample names to that donor's `samples_exclude` list instead of editing the file. |
| `Control SD for condition 'X' is zero/NA` | All control sgRNAs have the same value for condition X — check that your `zscore_control` SimpleName actually matches control rows in your data |
| Plots look strange / labels overlap | Pass a different palette in `pandonor.palette`, or run the script and then post-process the saved PDFs |
