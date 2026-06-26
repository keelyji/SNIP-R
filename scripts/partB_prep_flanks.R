#!/usr/bin/env Rscript
#
# partB_prep_flanks.R
#
# SNIP-R Part B, step 1 of 2: prepare flanking regions for CRISPick sgRNA design.
#
# Reads a BED file of target regions and a reference FASTA, optionally shrinks
# large regions inward, builds left/right flanking windows, and writes:
#   - left.fa             FASTA of left flanks  (paste into CRISPick)
#   - right.fa            FASTA of right flanks (paste into CRISPick)
#   - left_coords.tsv     Coordinates of left flanks  (alt CRISPick input)
#   - right_coords.tsv    Coordinates of right flanks (alt CRISPick input)
#   - ROItoFlank.tsv      Mapping of ROI <-> left flank <-> right flank
#                         (required by partB_pair_grnas.R)
#
# Usage:
#   Rscript partB_prep_flanks.R \
#     --regions Finaldf.bed.txt \
#     --reference chr12.fa \
#     --out-dir output/
#
# CRISPick settings used in the paper (Feb 2026 run):
#   Human GRCh37 (hg19), CRISPRko, SpyoCas9 (Chen tracrRNA, RuleSet3),
#   CRISPick quota: 20.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ---- args ------------------------------------------------------------------

option_list <- list(
  make_option(c("-r", "--regions"), type = "character",
              help = "BED file of target regions (chr, start, end, [name])"),
  make_option(c("-f", "--reference"), type = "character",
              help = "Reference genome FASTA (must have .fai index in same dir, or one will be made)"),
  make_option(c("-o", "--out-dir"), type = "character", default = "output",
              help = "Output directory [default: %default]"),
  make_option("--shrink-cutoff", type = "integer", default = 500,
              help = "Regions larger than this (bp) are shrunk inward [default: %default]"),
  make_option("--shrink-amount", type = "integer", default = 50,
              help = "Amount (bp) to shrink each side of large regions [default: %default]"),
  make_option("--flank-size", type = "integer", default = 550,
              help = "Length (bp) of left/right flank windows [default: %default]"),
  make_option("--skip-fasta", action = "store_true", default = FALSE,
              help = "Skip FASTA extraction (only emit coordinate files). Useful if bedtools is not installed.")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- validate --------------------------------------------------------------

stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)
stop_if(is.null(opt$regions),  "--regions is required")
stop_if(is.null(opt$reference) && !opt$`skip-fasta`,
        "--reference is required unless --skip-fasta is set")
stop_if(!file.exists(opt$regions), paste("Regions file not found:", opt$regions))
if (!is.null(opt$reference)) {
  stop_if(!file.exists(opt$reference), paste("Reference FASTA not found:", opt$reference))
}

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

message("[partB_prep_flanks] regions:      ", opt$regions)
message("[partB_prep_flanks] reference:    ", if (is.null(opt$reference)) "(skipped)" else opt$reference)
message("[partB_prep_flanks] shrink:       >", opt$`shrink-cutoff`, "bp by ", opt$`shrink-amount`, "bp/side")
message("[partB_prep_flanks] flank size:   ", opt$`flank-size`, "bp")
message("[partB_prep_flanks] out-dir:      ", opt$`out-dir`)

# ---- read regions ----------------------------------------------------------

df <- read_tsv(opt$regions, col_names = FALSE, show_col_types = FALSE,
               comment = "#", progress = FALSE)
if (ncol(df) < 3) stop("BED file must have at least 3 columns (chr, start, end)", call. = FALSE)
# Pad to 4 columns: chr, start, end, name
if (ncol(df) == 3) df$X4 <- paste0("region_", seq_len(nrow(df)))
df <- df[, 1:4]
colnames(df) <- c("Chr", "Start", "End", "Name")
df$Size <- df$End - df$Start

message("[partB_prep_flanks] regions read: ", nrow(df))

# ---- shrink large regions inward ------------------------------------------

df_small <- df %>% dplyr::filter(.data$Size <= opt$`shrink-cutoff`)
df_large <- df %>% dplyr::filter(.data$Size >  opt$`shrink-cutoff`) %>%
  mutate(Start = .data$Start + opt$`shrink-amount`,
         End   = .data$End   - opt$`shrink-amount`)
df2 <- bind_rows(df_small, df_large) %>%
  mutate(Size = .data$End - .data$Start)
message("[partB_prep_flanks]   shrunk:     ", nrow(df_large))
message("[partB_prep_flanks]   left as-is: ", nrow(df_small))

# ---- build flanking windows ------------------------------------------------

flanked <- df2 %>%
  mutate(Left_Start  = .data$Start - opt$`flank-size`,
         Left_End    = .data$Start,
         Right_Start = .data$End,
         Right_End   = .data$End + opt$`flank-size`)

# Sanity: no negative coordinates
if (any(flanked$Left_Start < 0)) {
  warning("Some left flanks have negative coordinates and were clipped to 0.")
  flanked$Left_Start <- pmax(0L, flanked$Left_Start)
}

left_bed  <- flanked %>% select("Chr", "Left_Start",  "Left_End")
right_bed <- flanked %>% select("Chr", "Right_Start", "Right_End")

# Build the ROI <-> flank mapping table used downstream
conversion <- flanked %>%
  transmute(
    ROI        = paste(.data$Chr, .data$Start,       .data$End,       sep = "-"),
    LeftFlank  = paste(.data$Chr, .data$Left_Start,  .data$Left_End,  sep = "-"),
    RightFlank = paste(.data$Chr, .data$Right_Start, .data$Right_End, sep = "-"),
    Type = .data$Name,
    Size = .data$Size
  )

conv_path <- file.path(opt$`out-dir`, "ROItoFlank.tsv")
write_tsv(conversion, conv_path)
message("[partB_prep_flanks] wrote ", conv_path)

# Also write coordinate-only files (alt CRISPick input)
left_coords_path  <- file.path(opt$`out-dir`, "left_coords.tsv")
right_coords_path <- file.path(opt$`out-dir`, "right_coords.tsv")
write_tsv(left_bed,  left_coords_path,  col_names = FALSE)
write_tsv(right_bed, right_coords_path, col_names = FALSE)
message("[partB_prep_flanks] wrote ", left_coords_path)
message("[partB_prep_flanks] wrote ", right_coords_path)

# ---- FASTA extraction ------------------------------------------------------

if (!opt$`skip-fasta`) {
  if (!requireNamespace("bedtoolsr", quietly = TRUE)) {
    stop("bedtoolsr is not installed. Install it, or rerun with --skip-fasta.",
         call. = FALSE)
  }
  fasta_in <- normalizePath(opt$reference)

  # bedtoolsr::bt.getfasta returns a data.frame whose single column alternates
  # ">name" and sequence lines, as bedtools getfasta would emit on stdout.
  extract_fasta <- function(bed_df, label) {
    out <- bedtoolsr::bt.getfasta(fi = fasta_in, bed = as.data.frame(bed_df))
    if (is.data.frame(out)) out <- out[[1]]
    out
  }

  left_fa  <- extract_fasta(left_bed,  "left")
  right_fa <- extract_fasta(right_bed, "right")

  left_fa_path  <- file.path(opt$`out-dir`, "left.fa")
  right_fa_path <- file.path(opt$`out-dir`, "right.fa")
  writeLines(left_fa,  left_fa_path)
  writeLines(right_fa, right_fa_path)
  message("[partB_prep_flanks] wrote ", left_fa_path)
  message("[partB_prep_flanks] wrote ", right_fa_path)
}

# ---- next-step instructions ------------------------------------------------

cat("\n",
    "------------------------------------------------------------\n",
    " NEXT STEPS (manual CRISPick run)\n",
    "------------------------------------------------------------\n",
    " 1. Open https://portals.broadinstitute.org/gppx/crispick/public\n",
    " 2. Settings:\n",
    "      Reference   : Human GRCh37 (hg19)\n",
    "      Mechanism   : CRISPRko\n",
    "      Enzyme      : SpyoCas9 (Chen tracrRNA, RuleSet3)\n",
    "      Quota       : 20\n",
    " 3. Submit LEFT  flanks: paste contents of  ",
        file.path(opt$`out-dir`, "left.fa"),  "\n",
    "    Save the resulting design file as     CRISPick/Left-sgrna-designs.txt\n",
    " 4. Submit RIGHT flanks: paste contents of  ",
        file.path(opt$`out-dir`, "right.fa"), "\n",
    "    Save the resulting design file as     CRISPick/Right-sgrna-designs.txt\n",
    " 5. Run partB_pair_grnas.R, supplying both CRISPick outputs\n",
    "    and the ROItoFlank.tsv from this step.\n",
    "------------------------------------------------------------\n",
    sep = "")
