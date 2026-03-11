# SNIP-R

This repository provides code for designing and analyzing cis-regulatory element CRISPR screens in primary T cells with the Systematic Non-coding element Interrogation by Paired sgRNAs (SNIP-R) system introduced in the paper "Scalable hit-and-run platform for enhancer deletion reveals state-specific IFNG regulation in primary human T cells." The code for SNIP-R is mainly divided into three sections: A) SNIP-R screen region selection, B) SNIP-R sgRNA pairing, and C) SNIP-R screen analysis.

The associated data for SNIP-R screen region selection requires ATAC-seq dataset from Yates et al, 2021 (dbGaP study accession: phs002510.v1.p1). Other datasets used for SNIP-R screen region selection are provided in this repository. For SNIP-R screen analysis code, raw CRISPR screening data from this paper is needed and provided in Supplementary Table 4.

For reference, SNIP-R identified enhancers from acute and chronic screen in the paper are provided as a bed track called "AllRegulatoryElement_AcuteChronic_SNIP-Rscreens.hg19.bed.txt"

Code and dataset to reproduce main figures from the paper are located in the Part C SNIP-R screen analysis folder.

# Contents
- **Part A SNIP-R screen region selection:** directory containing the data and scripts for reproducing the screen region selection
- **Part B SNIP-R sgRNA pairing:** directory containing the data and scripts for sgRNA pairing for dual-sgRNA mediated deletion
- **Part C SNIP-R screen analysis:** directory containing the data and scripts for reproducing the screen analyses in the figures
- **Part D other codes:** ATAC-seq code used for the paper and Hi-C visualization with Plotgardener.
- **Part E other datasets:** Bed files of all regulatory elements screened and identified by IFNG SNIP-R screens 

# Instructions
1. Download the github content or clone the github content (expected time: one minute).
2. Follow below instructions, which are separated into different parts. Note, specific package version is available in the html output file in respective folder.
3. Downloading R and R studio is necessary to run the codes for SNIP-R region selection, design, and analysis.

## Part A: To select regions to target with SNIP-R
1. All code for target region selection are performed in R and R studio. Key packages used in this steps are: bedtoolsr from PhanstielLab, tidyverse, and readr.
- 00_ScreenRegionSelection_V2.Rmd file contains all the code necessary for this step. The exact R package and version used in the paper are noted in the respective html file.
3. Necessary dataset (Note: we included example datasets used in our paper as reference to replicate our result. As our paper focus on IFNG, located in human chromosome 12, most files only contain chromosome 12 information).
- ATAC-seq counts matrix and the associated peak universe. In our paper, we focused on in vivo CD8+ T cell ATAC-seq peaks (dbGaP study accession (phs002510.v1.p1) and filtered with in vitro CD8+ T cell ATAC-seq datasets (as we wanted to screen for in vivo relevant peaks in vitro culture setting). ATAC-seq peak universe is used to identify open chromatin regions to target, while counts matrix is used to filter for differentially accessible regions. We included the peak universe from phs002510.v1.p1 called "Yates_all_human_Tcell_ATAC.hg19.merged.bed", "Yates_nondiffOCR_chr12.bed", and "InVitro_OCR_chr12.bed" in this repository.
- General ATAC-seq pipeline to create count matrix and peak universe from fastq files are provided in "ATACpipeline.txt" file provided in this repository.
- UCSC protein coding gene bed file downloaded from UCSC. We included "UCSC_proteincodinggenes_hg19.bed.zip" in this repository.
- human genome chromosome size file. For instance, download hg19.chrom.sizes from https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/. We included "human.hg19.genome file" in this repository.
- bed file containing chr start stop for the genomic locus interested in screening. We included "IFNGlocus.bed.txt" in this repository. 
- All gene locus bed file downloaded from NCBI. We included "NCBI_hg19_refseq_chr12.bed" for our repository.
4. Expected time: demo dataset (30 minute)



## Part B: To design sgRNA pairs for SNIP-R screen
1. All code for sgRNA pairing are performed in R and R studio, while the gRNA design is performed in Broad Institute's CRISPick. Please note that input and output of CRISPick may change overtime as it is updated by Broad Institute, so the result output might vary and code may need adjustment according to the CRISPick version used. We updated the sgRNA pairing code to reflect the latest CRISPick update as of February 2026 for this repository. In general, the following CRISPick setting was used: Human GRCh37 (hg19), CRISPRko, SpyoCas9 (Chen tracrRNA, RuleSet3), CRISPickQuota: 20. Guide designs were provided by the CRISPick web tool of the GPP at the Broad Institute (Sanson et al. 2018, Doench et al. 2016)

- 01_SNIP-R_Pairing_V2.Rmd file contains all the code necessary for this step. The exact R package and version used in the paper are noted in the respective html file.
2. Necessary dataset are bed file of target regions (e.g., output from Part B).
3. Expected time: demo dataset (1 hour)

## Part C: To analyse the result of SNIP-R screen
1. All code for target region selection are performed in R and R studio. Key packages used in this steps are: readr, tidyverse, readxl, ggrepel, and patchwork.
- 02_ScreenAnalysis.Rmd file contains all the code necessary for this step.The exact R package and version used in the paper are noted in the respective html file.
2. Necessary dataset are:
- Screen sequencing data in count matrix form. We provided a counts matrix called "SupplementaryTable4_RawScreenData.xlsx" available in our paper and in the SNIP-R screen analysis folder located in this repository.
3. Analysis flowchart and details on the analysis are noted in our paper.
4. We also provided code to replicate main figures in our paper.
5. Expected time: demo dataset (1 hour)

## Part D: Other codes
- Contains ATAC-seq script used for the paper
- Contains R markdown file for Hi-C visualization with Plotgardener. All files to reproduce figure 1 Hi-C plots are located in Part E: data folder

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
```

