#!/usr/bin/env Rscript
#
# partB_pair_grnas.R
#
# SNIP-R Part B, step 2 of 2: score CRISPick outputs, pair left/right sgRNAs
# per ROI, filter to a reusable subset, and flag Golden Gate cut sites.
#
# Inputs:
#   --left-crispick     Left  flank CRISPick design file (.txt)
#   --right-crispick    Right flank CRISPick design file (.txt)
#   --conversion-table  ROItoFlank.tsv from partB_prep_flanks.R
#
# Output (TSV) columns (one row per gRNA pair):
#   Left_ROI, Left_flank, Left_sgRNA.Sequence, Left_On.Target.Efficacy.Score,
#   Right_ROI, Right_flank, Right_sgRNA.Sequence, Right_On.Target.Efficacy.Score,
#   Summed_Score, has_BsmBI_left, has_BsmBI_right
#
# Usage:
#   Rscript partB_pair_grnas.R \
#     --left-crispick  CRISPick/Left-sgrna-designs.txt \
#     --right-crispick CRISPick/Right-sgrna-designs.txt \
#     --conversion-table output/ROItoFlank.tsv \
#     --out output/paired_grnas.tsv

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
})

# ---- args ------------------------------------------------------------------

option_list <- list(
  make_option("--left-crispick",     type = "character",
              help = "CRISPick design file for LEFT flanks"),
  make_option("--right-crispick",    type = "character",
              help = "CRISPick design file for RIGHT flanks"),
  make_option("--conversion-table",  type = "character",
              help = "ROItoFlank.tsv produced by partB_prep_flanks.R"),
  make_option(c("-o", "--out"),      type = "character", default = "paired_grnas.tsv",
              help = "Output TSV path [default: %default]"),
  make_option("--top-n-per-flank",   type = "integer", default = 5,
              help = "Top N sgRNAs to retain per flank by weighted score [default: %default]"),
  make_option("--on-weight",         type = "double",  default = 0.667,
              help = "Weight on on-target rank [default: %default]"),
  make_option("--off-weight",        type = "double",  default = 0.333,
              help = "Weight on off-target rank [default: %default]"),
  make_option("--pairs-pool",        type = "integer", default = 25,
              help = "Top-K left*right pairs (by summed on-target score) considered per ROI before reuse filtering [default: %default]"),
  make_option("--pairs-per-roi",     type = "integer", default = 10,
              help = "Final number of pairs to emit per ROI [default: %default]"),
  make_option("--max-grna-reuse",    type = "integer", default = 3,
              help = "Max times any single sgRNA may be reused within an ROI's pairs [default: %default]"),
  make_option("--bsmbi-sites",       type = "character", default = "CGTCTC,GAGACG",
              help = "Comma-separated motifs flagged as Golden Gate cut sites [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)
stop_if(is.null(opt$`left-crispick`),    "--left-crispick is required")
stop_if(is.null(opt$`right-crispick`),   "--right-crispick is required")
stop_if(is.null(opt$`conversion-table`), "--conversion-table is required")
for (f in c(opt$`left-crispick`, opt$`right-crispick`, opt$`conversion-table`)) {
  stop_if(!file.exists(f), paste("File not found:", f))
}

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
bsmbi_motifs <- strsplit(opt$`bsmbi-sites`, ",", fixed = TRUE)[[1]]

# ---- read inputs -----------------------------------------------------------

Left  <- read_tsv(opt$`left-crispick`,  show_col_types = FALSE, progress = FALSE)
Right <- read_tsv(opt$`right-crispick`, show_col_types = FALSE, progress = FALSE)
# Normalize column names to syntactic identifiers (e.g. "sgRNA Sequence"
# -> "sgRNA.Sequence") for consistency with the published CRISPick fields.
colnames(Left)  <- make.names(colnames(Left))
colnames(Right) <- make.names(colnames(Right))
conv  <- read_tsv(opt$`conversion-table`, show_col_types = FALSE, progress = FALSE)

needed <- c("Input", "sgRNA.Sequence", "On.Target.Efficacy.Score",
            "On.Target.Rank", "Aggregate.CFD.Score")
for (col in needed) {
  if (!col %in% colnames(Left))  stop("Left CRISPick file missing column: ",  col, call. = FALSE)
  if (!col %in% colnames(Right)) stop("Right CRISPick file missing column: ", col, call. = FALSE)
}

message("[partB_pair_grnas]  left  CRISPick rows: ", nrow(Left))
message("[partB_pair_grnas]  right CRISPick rows: ", nrow(Right))
message("[partB_pair_grnas]  conversion rows:     ", nrow(conv))

# ---- compute off-target rank from CFD score within each flank -------------
# Newer CRISPick outputs no longer include Off.Target.Rank directly, so we
# derive it from Aggregate.CFD.Score (lower CFD = better; min_rank gives
# rank 1 to the lowest CFD).

rank_and_score <- function(crispick_df) {
  crispick_df %>%
    group_by(.data$Input) %>%
    mutate(Off.Target.Rank = min_rank(.data$Aggregate.CFD.Score)) %>%
    ungroup() %>%
    mutate(Weighted_Score =
             opt$`on-weight`  * .data$On.Target.Rank +
             opt$`off-weight` * .data$Off.Target.Rank)
}

Left_scored  <- rank_and_score(Left)
Right_scored <- rank_and_score(Right)

# Top N sgRNAs per flank by weighted score (lower = better)
keep_top <- function(scored_df, n) {
  scored_df %>%
    group_by(.data$Input) %>%
    slice_min(order_by = .data$Weighted_Score, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    select("Input", "sgRNA.Sequence", "On.Target.Efficacy.Score",
           "On.Target.Rank", "Off.Target.Rank", "Weighted_Score")
}

Left_top  <- keep_top(Left_scored,  opt$`top-n-per-flank`)
Right_top <- keep_top(Right_scored, opt$`top-n-per-flank`)
top <- bind_rows(Left_top, Right_top) %>%
  mutate(Input = gsub(":", "-", .data$Input, fixed = TRUE)) %>%
  rename(flank = "Input")

# ---- attach ROI via conversion table --------------------------------------

conv_long <- bind_rows(
  conv %>% transmute(ROI = .data$ROI, flank = .data$LeftFlank,
                     Type = .data$Type, Size = .data$Size, side = "left"),
  conv %>% transmute(ROI = .data$ROI, flank = .data$RightFlank,
                     Type = .data$Type, Size = .data$Size, side = "right")
)

joined <- inner_join(top, conv_long, by = "flank") %>%
  arrange(.data$ROI) %>%
  select("ROI", "flank", "side", "sgRNA.Sequence",
         "On.Target.Efficacy.Score")

n_rois <- length(unique(joined$ROI))
message("[partB_pair_grnas]  unique ROIs with guides: ", n_rois)

# ---- pair left x right per ROI by summed on-target score ------------------

pair_one_roi <- function(roi_label, df, pool_size) {
  rows <- df %>% dplyr::filter(.data$ROI == roi_label)
  left  <- rows %>% dplyr::filter(.data$side == "left")
  right <- rows %>% dplyr::filter(.data$side == "right")
  if (nrow(left) == 0 || nrow(right) == 0) return(NULL)

  L <- left  %>% rename_with(~ paste0("Left_",  .x))
  R <- right %>% rename_with(~ paste0("Right_", .x))
  combs <- tidyr::expand_grid(L, R) %>%
    mutate(Summed_Score = .data$Left_On.Target.Efficacy.Score +
                          .data$Right_On.Target.Efficacy.Score) %>%
    arrange(desc(.data$Summed_Score)) %>%
    slice(seq_len(pool_size))
  combs
}

paired_pool <- bind_rows(
  lapply(unique(joined$ROI), pair_one_roi, df = joined,
         pool_size = opt$`pairs-pool`)
)

# ---- enforce sgRNA reuse cap per ROI --------------------------------------

filter_reuse <- function(pool_df, roi_label, n, max_reuse) {
  rows <- pool_df %>% dplyr::filter(.data$Left_ROI == roi_label)
  if (!nrow(rows)) return(NULL)
  picked <- list()
  left_use <- character(0); right_use <- character(0)
  for (i in seq_len(nrow(rows))) {
    pair <- rows[i, ]
    l <- pair$Left_sgRNA.Sequence
    r <- pair$Right_sgRNA.Sequence
    if (sum(left_use  == l) < max_reuse &&
        sum(right_use == r) < max_reuse) {
      picked[[length(picked) + 1L]] <- pair
      left_use  <- c(left_use,  l)
      right_use <- c(right_use, r)
      if (length(picked) == n) break
    }
  }
  if (!length(picked)) return(NULL)
  dplyr::bind_rows(picked)
}

final <- bind_rows(
  lapply(unique(paired_pool$Left_ROI), filter_reuse,
         pool_df    = paired_pool,
         n          = opt$`pairs-per-roi`,
         max_reuse  = opt$`max-grna-reuse`)
)

# ---- BsmBI flagging --------------------------------------------------------

flag_motifs <- function(seq, motifs) {
  vapply(seq, function(s) any(vapply(motifs, grepl, logical(1), x = s, fixed = TRUE)),
         logical(1))
}
final$has_BsmBI_left  <- flag_motifs(final$Left_sgRNA.Sequence,  bsmbi_motifs)
final$has_BsmBI_right <- flag_motifs(final$Right_sgRNA.Sequence, bsmbi_motifs)

n_flag <- sum(final$has_BsmBI_left | final$has_BsmBI_right)
if (n_flag) {
  message("[partB_pair_grnas]  WARNING: ", n_flag,
          " pairs contain a Golden Gate (BsmBI) motif (", opt$`bsmbi-sites`, ")",
          " - inspect before ordering.")
}

# ---- write -----------------------------------------------------------------

write_tsv(final, opt$out)
message("[partB_pair_grnas]  wrote ", opt$out, "  (", nrow(final), " pairs across ",
        length(unique(final$Left_ROI)), " ROIs)")
