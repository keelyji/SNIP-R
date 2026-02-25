# SNIP-R

This repository provides code for designing and analyzing cis-regulatory element CRISPR screens in primary T cells with the Systematic Non-coding element Interrogation by Paired sgRNAs (SNIP-R) system introduced in the paper "Scalable hit-and-run platform for enhancer deletion reveals state-specific IFNG regulation in primary human T cells." The code for SNIP-R is mainly divided into three sections: A) SNIP-R screen region selection, B) SNIP-R sgRNA pairing, and C) SNIP-R screen analysis.

The associated data for SNIP-R screen region selection requires ATAC-seq dataset from Yates et al, 2021 (dbGaP study accession: phs002510.v1.p1). Other datasets used for SNIP-R screen region selection are provided in this repository. For SNIP-R screen analysis code, raw CRISPR screening data from this paper is needed and provided in supplementary table 4.

For reference, SNIP-R identified enhancers from acute and chronic screen in the paper are provided as a bed track called "AllRegulatoryElement_AcuteChronic_SNIP-Rscreens.hg19.bed.txt"

# Contents
- SNIP-R screen region selection: directory containing the data and scripts for reproducing the screen region selection
- SNIP-R sgRNA pairing: directory containing the data and scripts for sgRNA pairing for dual-sgRNA mediated deletion
- SNIP-R screen analysis: directory containing the data and scripts for reproducing the screen analyses in the figures

# Instructions
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



## Part B: To design sgRNA pairs for SNIP-R screen
1. All code for sgRNA pairing are performed in R and R studio, while the gRNA design is performed in Broad Institute's CRISPick. Please note that input and output of CRISPick may change overtime as it is updated by Broad Institute, so the result output might vary and code may need adjustment according to the CRISPick version used. We updated the sgRNA pairing code to reflect the latest CRISPick update as of February 2026 for this repository. In general, the following CRISPick setting was used: Human GRCh37 (hg19), CRISPRko, SpyoCas9 (Chen tracrRNA, RuleSet3), CRISPickQuota: 20.
- 01_SNIP-R_Pairing_V2.Rmd file contains all the code necessary for this step. The exact R package and version used in the paper are noted in the respective html file.
2. Necessary dataset are bed file of target regions (e.g., output from Part B).


## Part C: To analyse the result of SNIP-R screen
1. All code for target region selection are performed in R and R studio. Key packages used in this steps are: readr, tidyverse, readxl, ggrepel, and patchwork.
- 02_ScreenAnalysis.Rmd file contains all the code necessary for this step.The exact R package and version used in the paper are noted in the respective html file.
2. Necessary dataset are:
- Screen sequencing data in count matrix form. We provided a counts matrix called "SupplementaryTable4_RawScreenData.xlsx" available in our paper and in the SNIP-R screen analysis folder located in this repository.
3. Analysis flowchart and details on the analysis are noted in our paper.
4. We also provided code to replicate main figures in our paper.

# Contact
Please contact corresponding author for any questions, comments, or concerns regarding the paper in general. For issues with code/reproducibility, please open a GitHub issue.


