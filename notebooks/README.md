# SNIP-R notebooks

These R Markdown documents are the original code used to generate the figures and tables in the paper. They are preserved here for transparency. **For running the pipeline on your own data, use the command-line scripts under `../scripts/` (see top-level `README.md`).**

## Setup for the Part B notebook

`PartB_SNIP-R_sgRNApairing/01_SNIP-R_Pairing_V2.Rmd` needs an uncompressed `chr12.fa` in the same directory. The repo ships `chr12.fa.gz` (~41 MB compressed). Extract once:

    cd notebooks/PartB_SNIP-R_sgRNApairing/
    gunzip -k chr12.fa.gz   # keeps chr12.fa.gz; produces chr12.fa (~130 MB)

The CLI pipeline (`../scripts/partB_prep_flanks.R`) does not need this — it takes a reference FASTA path via `--reference`, so use whatever uncompressed `.fa` file you already have for your genome.

## Contents

| Folder | Notebook | What it does |
|--------|----------|--------------|
| `PartA_SNIP-R_ScreenRegionSelection/` | `00_ScreenRegionSelection_V2.Rmd` | Selects screening regions from ATAC-seq data and filters |
| `PartB_SNIP-R_sgRNApairing/`          | `01_SNIP-R_Pairing_V2.Rmd`        | Notebook version of the sgRNA pairing workflow |
| `PartC_SNIP-R_ScreenAnalysis/`        | `02_ScreenAnalysis.Rmd`           | Notebook version of the screen analysis (produces the paper figures) |
| `PartD_other_code/`                   | `HiCplot_figure1.Rmd`, `ATACpipeline.txt` | Hi-C plotting and ATAC-seq pipeline notes |
| `PartE_other_datasets/`               | (data only)                       | Bed files of all regulatory elements screened and identified |

Each Rmd has a rendered `.html` alongside it with the exact session info and package versions used.

## External datasets required

- Part A: ATAC-seq from Yates et al. 2021 — dbGaP study accession `phs002510.v1.p1`
- Part D Hi-C plotting: GEO `GSE126117` (Bediaga et al.)

The other inputs needed by each notebook are bundled within the respective subfolder.
