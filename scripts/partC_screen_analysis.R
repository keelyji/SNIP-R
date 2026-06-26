#!/usr/bin/env Rscript
#
# partC_screen_analysis.R
#
# SNIP-R Part C: end-to-end screen analysis from a counts matrix and a
# YAML config describing the sample design. Produces per-pair Z-scores,
# per-region summaries, pan-donor mean Z + t-test tables, and optional
# concordance / volcano / distribution plots.
#
# Usage:
#   Rscript partC_screen_analysis.R \
#     --counts SupplementaryTable4_RawScreenData.xlsx \
#     --config config/example_partC_config.yaml \
#     --out-dir output/
#
# See config/example_partC_config.yaml for the full schema.

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(tools)
})

# ---- args ------------------------------------------------------------------

option_list <- list(
  make_option(c("-c", "--counts"),  type = "character",
              help = "Counts matrix (.xlsx/.xls/.tsv/.csv). Overrides input.counts_file in config."),
  make_option(c("-y", "--config"),  type = "character",
              help = "YAML config describing the sample design"),
  make_option(c("-o", "--out-dir"), type = "character", default = "output",
              help = "Output directory [default: %default]"),
  make_option("--sheet",            type = "character", default = NULL,
              help = "Sheet name for Excel input (overrides config)")
)
opt <- parse_args(OptionParser(option_list = option_list))

stop_if <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)
stop_if(is.null(opt$config), "--config is required")
stop_if(!file.exists(opt$config), paste("Config not found:", opt$config))

cfg <- yaml::read_yaml(opt$config)

# CLI overrides config
if (!is.null(opt$counts)) cfg$input$counts_file <- opt$counts
if (!is.null(opt$sheet))  cfg$input$sheet <- opt$sheet

stop_if(is.null(cfg$input$counts_file), "input.counts_file must be set (config or --counts)")
stop_if(!file.exists(cfg$input$counts_file),
        paste("Counts file not found:", cfg$input$counts_file))

dir.create(opt$`out-dir`, recursive = TRUE, showWarnings = FALSE)

message("[partC]  counts:  ", cfg$input$counts_file)
message("[partC]  config:  ", opt$config)
message("[partC]  out-dir: ", opt$`out-dir`)

# ---- helpers ---------------------------------------------------------------

read_counts <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (is.null(sheet)) readxl::read_excel(path) else readxl::read_excel(path, sheet = sheet)
  } else if (ext == "csv") {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  }
}

as_char_vec <- function(x) {
  if (is.null(x)) character(0) else as.character(unlist(x))
}

# `%||%` for older R / no rlang dep
`%||%` <- function(a, b) if (is.null(a)) b else a

# Validate every expected donor x sample column exists in the input matrix.
# Per-donor samples_exclude entries (set via the donors[].samples_exclude config
# key) are skipped — e.g. for screens where a specific condition failed QC in
# one donor and that donor's data for those conditions should not be loaded.
check_columns <- function(df, donors, samples) {
  expected <- unlist(lapply(donors, function(d) {
    excl <- as_char_vec(d$samples_exclude)
    paste0(d$prefix, "_", setdiff(samples, excl))
  }))
  missing  <- setdiff(expected, colnames(df))
  if (length(missing)) {
    stop("Counts matrix is missing required columns:\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
  }
}

# Build a per-donor sub-frame: id columns + the donor's samples, with the
# donor prefix stripped so all downstream code is donor-agnostic. Sample
# columns listed in `exclude` are read as NA so they propagate through the
# pipeline without contributing to any statistics.
make_donor_frame <- function(df, id_cols, prefix, samples, exclude = character(0)) {
  effective_samples <- setdiff(samples, exclude)
  donor_cols <- paste0(prefix, "_", effective_samples)
  sub <- df[, c(id_cols, donor_cols), drop = FALSE]
  colnames(sub) <- c(id_cols, effective_samples)
  for (s in exclude) sub[[s]] <- NA_real_
  # Reorder to canonical sample order
  sub[, c(id_cols, samples), drop = FALSE]
}

# ---- read & validate -------------------------------------------------------

df <- read_counts(cfg$input$counts_file, cfg$input$sheet)
id_cols      <- as_char_vec(cfg$id_columns)
samples      <- as_char_vec(cfg$samples)
donor_specs  <- cfg$donors
donor_labels <- vapply(donor_specs, function(d) d$label, character(1))

# Ensure id columns exist
missing_id <- setdiff(id_cols, colnames(df))
if (length(missing_id))
  stop("ID columns missing from counts matrix: ", paste(missing_id, collapse = ", "),
       call. = FALSE)

check_columns(df, donor_specs, samples)

message("[partC]  donors:  ", paste(donor_labels, collapse = ", "))
message("[partC]  samples per donor: ", length(samples))

# ---- per-donor preprocessing ----------------------------------------------
#
# Order of operations per donor:
#   (1) split out donor's columns, strip prefix
#   (2) compute derived columns (sums of components) on raw counts
#   (3) log2-CPM normalize all numeric columns: log2(x/sum(x)*1e6 + 1)
#   (4) [joint across donors] drop rows where any filter column == 0
#   (5) subtract baseline_column (log) from every numeric column
#   (6) Z-score each `zscore_condition` against the `zscore_control`
#       distribution within that donor

computed_specs   <- cfg$computed_columns %||% list()
filter_cols      <- as_char_vec(cfg$input_filter_columns)
baseline_col     <- cfg$baseline_column
zscore_conds     <- as_char_vec(cfg$zscore_conditions)
zscore_ctrl_vals <- as_char_vec(cfg$zscore_control)

donor_frames_raw <- lapply(donor_specs, function(d)
  make_donor_frame(df, id_cols, d$prefix, samples,
                   exclude = as_char_vec(d$samples_exclude)))
names(donor_frames_raw) <- donor_labels

# Report excluded donor x sample combinations for transparency
for (d in donor_specs) {
  excl <- as_char_vec(d$samples_exclude)
  if (length(excl))
    message("[partC]  excluding samples for donor ", d$label, ": ",
            paste(excl, collapse = ", "))
}

# (2) computed columns on raw counts. If any component column is entirely NA
# (e.g. because excluded for this donor), the computed column is also NA so
# the missingness propagates correctly.
add_computed <- function(donor_df, specs) {
  for (s in specs) {
    name <- s$name
    comp <- as_char_vec(s$components)
    missing <- setdiff(comp, colnames(donor_df))
    if (length(missing))
      stop("computed_columns[", name, "] components missing: ",
           paste(missing, collapse = ", "), call. = FALSE)
    if (any(vapply(comp, function(col) all(is.na(donor_df[[col]])), logical(1)))) {
      donor_df[[name]] <- NA_real_
    } else {
      donor_df[[name]] <- rowSums(donor_df[, comp, drop = FALSE], na.rm = TRUE)
    }
  }
  donor_df
}
donor_frames_raw <- lapply(donor_frames_raw, add_computed, specs = computed_specs)

all_sample_cols <- c(samples, vapply(computed_specs, function(s) s$name, character(1)))

# (3) log2-CPM normalize numeric columns. All-NA columns (excluded for this
# donor) are passed through unchanged.
log_norm <- function(donor_df, num_cols) {
  for (col in num_cols) {
    if (all(is.na(donor_df[[col]]))) next
    x <- donor_df[[col]]
    donor_df[[col]] <- log2(x / sum(x, na.rm = TRUE) * 1e6 + 1)
  }
  donor_df
}
donor_frames_log <- lapply(donor_frames_raw, log_norm, num_cols = all_sample_cols)

# (4) joint filter: keep rows where every filter_col is non-zero in every donor
# that has data for that column. For a donor where the column is fully excluded
# (all NA), that donor is treated as "passing" — we don't drop rows based on a
# column the donor doesn't have.
if (length(filter_cols)) {
  keep <- Reduce(`&`, lapply(donor_frames_log, function(d) {
    Reduce(`&`, lapply(filter_cols, function(col) {
      x <- d[[col]] != 0
      ifelse(is.na(x), TRUE, x)
    }))
  }))
  n_before <- nrow(donor_frames_log[[1]])
  donor_frames_log <- lapply(donor_frames_log, function(d) d[keep, , drop = FALSE])
  message("[partC]  filter kept ", sum(keep), " / ", n_before, " rows")
}

# (5) subtract baseline column. All-NA columns are left as NA.
if (!is.null(baseline_col)) {
  if (!baseline_col %in% all_sample_cols)
    stop("baseline_column '", baseline_col, "' is not among samples or computed_columns",
         call. = FALSE)
  donor_frames_log <- lapply(donor_frames_log, function(d) {
    base <- d[[baseline_col]]
    for (col in all_sample_cols) {
      if (all(is.na(d[[col]]))) next
      d[[col]] <- d[[col]] - base
    }
    d
  })
  message("[partC]  subtracted baseline column: ", baseline_col)
}

# (6) Z-score conditions against control distribution: for each condition,
#     Norm_<cond> = (<cond> - mean(<cond> in controls)) / sd(<cond> in controls)
# Conditions whose data is fully NA for this donor (excluded samples) get a
# Norm_ column of NAs rather than failing, so downstream pan-donor logic can
# automatically use only the donors that have data.
zscore_donor <- function(donor_df, conds, ctrl_vals) {
  ctrl <- donor_df %>% dplyr::filter(.data$SimpleName %in% ctrl_vals)
  if (!nrow(ctrl))
    stop("No control rows (SimpleName in {",
         paste(ctrl_vals, collapse = ", "), "}) for Z-scoring", call. = FALSE)
  for (c in conds) {
    if (!c %in% colnames(donor_df))
      stop("Z-score condition '", c, "' not found in donor frame", call. = FALSE)
    if (all(is.na(donor_df[[c]]))) {
      donor_df[[paste0("Norm_", c)]] <- NA_real_
      next
    }
    mu <- mean(ctrl[[c]], na.rm = TRUE)
    sd <- stats::sd(ctrl[[c]], na.rm = TRUE)
    if (is.na(sd) || sd == 0)
      stop("Control SD for condition '", c, "' is zero/NA -- cannot Z-score", call. = FALSE)
    donor_df[[paste0("Norm_", c)]] <- (donor_df[[c]] - mu) / sd
  }
  donor_df
}

norm_frames <- lapply(donor_frames_log, zscore_donor,
                      conds = zscore_conds, ctrl_vals = zscore_ctrl_vals)

# ---- write per-pair Z-scores ----------------------------------------------

if (isTRUE(cfg$output$write_per_pair)) {
  for (i in seq_along(norm_frames)) {
    out <- norm_frames[[i]] %>%
      select(all_of(id_cols), starts_with("Norm_"))
    path <- file.path(opt$`out-dir`,
                      paste0("per_pair_zscores_", donor_labels[i], ".tsv"))
    write_tsv(out, path)
    message("[partC]  wrote ", path)
  }
}

# ---- per-region (per-donor) means -----------------------------------------

per_region_by_donor <- lapply(norm_frames, function(d) {
  d %>%
    group_by(.data$SimpleName) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
})

if (isTRUE(cfg$output$write_per_region)) {
  for (i in seq_along(per_region_by_donor)) {
    out <- per_region_by_donor[[i]] %>%
      select("SimpleName", starts_with("Norm_"))
    path <- file.path(opt$`out-dir`,
                      paste0("per_region_zscores_", donor_labels[i], ".tsv"))
    write_tsv(out, path)
    message("[partC]  wrote ", path)
  }
}

# ---- center-N filtering ---------------------------------------------------
# For each region, keep the N gRNA pairs whose value of the reference
# condition is closest to the region's mean. Controls are kept with a
# different per-name count.

cf            <- cfg$center_filter %||% list()
ref_cond      <- cf$reference_condition
n_per_roi     <- cf$n_per_roi    %||% 8L
n_per_control <- cf$n_per_control %||% 24L
ctrl_names    <- as_char_vec(cf$control_names)

centered_per_donor <- NULL
if (!is.null(ref_cond)) {
  ref_col <- paste0("Norm_", ref_cond)
  centered_per_donor <- lapply(norm_frames, function(d) {
    if (!ref_col %in% colnames(d))
      stop("center_filter.reference_condition '", ref_cond,
           "' has no Norm_ column (was it Z-scored?)", call. = FALSE)
    non_ctrl <- d %>% dplyr::filter(!.data$SimpleName %in% ctrl_names) %>%
      group_by(.data$SimpleName) %>%
      mutate(AvgB4 = mean(.data[[ref_col]], na.rm = TRUE),
             D     = abs(.data[[ref_col]] - .data$AvgB4)) %>%
      slice_min(order_by = .data$D, n = n_per_roi, with_ties = FALSE) %>%
      ungroup()
    ctrl <- d %>% dplyr::filter(.data$SimpleName %in% ctrl_names) %>%
      group_by(.data$SimpleName) %>%
      mutate(AvgB4 = mean(.data[[ref_col]], na.rm = TRUE),
             D     = abs(.data[[ref_col]] - .data$AvgB4)) %>%
      slice_min(order_by = .data$D, n = n_per_control, with_ties = FALSE) %>%
      ungroup()
    bind_rows(ctrl, non_ctrl)
  })
  if (isTRUE(cfg$output$write_per_pair)) {
    for (i in seq_along(centered_per_donor)) {
      path <- file.path(opt$`out-dir`,
                        paste0("per_pair_zscores_centered_",
                               donor_labels[i], ".tsv"))
      out <- centered_per_donor[[i]] %>%
        select(all_of(id_cols), starts_with("Norm_"))
      write_tsv(out, path)
      message("[partC]  wrote ", path)
    }
  }
}

# ---- pan-donor analysis ---------------------------------------------------

pd_conds <- as_char_vec(cfg$pandonor$conditions)
pd_ctrl  <- as_char_vec(cfg$pandonor$control %||% cfg$zscore_control)

if (length(pd_conds)) {
  # Join donor frames on identifiers, suffixing Norm_ columns by donor label.
  join_cols <- intersect(c(id_cols, "SimpleName"), colnames(norm_frames[[1]]))
  per_donor_norm <- lapply(seq_along(norm_frames), function(i) {
    d <- norm_frames[[i]] %>% select(all_of(join_cols), starts_with("Norm_"))
    norm_cols <- grep("^Norm_", colnames(d), value = TRUE)
    colnames(d)[match(norm_cols, colnames(d))] <-
      paste0(donor_labels[i], "_", norm_cols)
    d
  })
  joined <- Reduce(function(x, y) full_join(x, y, by = join_cols), per_donor_norm)

  pd_results <- list()
  for (cond in pd_conds) {
    norm_col <- paste0("Norm_", cond)
    donor_cols <- paste0(donor_labels, "_", norm_col)
    missing <- setdiff(donor_cols, colnames(joined))
    if (length(missing)) {
      warning("Skipping pandonor condition '", cond,
              "' -- missing columns: ", paste(missing, collapse = ", "))
      next
    }
    pdz_df <- joined %>%
      mutate(PDZ = rowMeans(across(all_of(donor_cols)), na.rm = TRUE)) %>%
      select(all_of(join_cols), "PDZ")

    ctrl_pdz <- pdz_df %>%
      dplyr::filter(.data$SimpleName %in% pd_ctrl) %>%
      pull(.data$PDZ)
    pdz_df$is_ctrl <- pdz_df$SimpleName %in% pd_ctrl

    region_stats <- pdz_df %>%
      dplyr::filter(!.data$is_ctrl) %>%
      group_by(.data$SimpleName) %>%
      summarise(
        n_pairs = dplyr::n(),
        meanZ   = mean(.data$PDZ, na.rm = TRUE),
        tval    = {
          x <- na.omit(.data$PDZ)
          if (length(x) < 2 || length(na.omit(ctrl_pdz)) < 2) NA_real_
          else unname(stats::t.test(x, ctrl_pdz)$statistic)
        },
        pval    = {
          x <- na.omit(.data$PDZ)
          if (length(x) < 2 || length(na.omit(ctrl_pdz)) < 2) NA_real_
          else stats::t.test(x, ctrl_pdz)$p.value
        },
        .groups = "drop"
      )

    # Append a single row for the pooled control distribution (for plotting)
    if (length(ctrl_pdz)) {
      region_stats <- bind_rows(
        region_stats,
        tibble(SimpleName = paste(pd_ctrl, collapse = "+"),
               n_pairs    = length(ctrl_pdz),
               meanZ      = mean(ctrl_pdz, na.rm = TRUE),
               tval       = 0,
               pval       = 1)
      )
    }
    pd_results[[cond]] <- region_stats

    if (isTRUE(cfg$output$write_pandonor)) {
      path <- file.path(opt$`out-dir`,
                        paste0("pandonor_", cond, ".tsv"))
      write_tsv(region_stats, path)
      message("[partC]  wrote ", path)
    }
  }
} else {
  pd_results <- list()
}

# ---- plots -----------------------------------------------------------------

plot_palette <- cfg$pandonor$palette %||% list()
highlight    <- cfg$pandonor$highlight_categories %||% list()

# Build a Color column from highlight rules (exact match on SimpleName).
add_color_column <- function(df, default = "OCR") {
  df$Color <- default
  for (cat in names(highlight)) {
    vals <- as_char_vec(highlight[[cat]])
    df$Color[df$SimpleName %in% vals] <- cat
  }
  df
}

# ---- volcano plots --------------------------------------------------------

if (isTRUE(cfg$output$plot_volcano) && length(pd_results)) {
  volcano_exclude_vals <- as_char_vec(cfg$pandonor$volcano_exclude)
  display_names <- cfg$pandonor$display_names %||% list()

  for (cond in names(pd_results)) {
    display <- display_names[[cond]] %||% cond
    pd <- pd_results[[cond]] %>%
      dplyr::filter(!is.na(.data$pval),
                    !.data$SimpleName %in% volcano_exclude_vals) %>%
      add_color_column()

    pal_vec <- unlist(plot_palette)
    pal_vec <- pal_vec[intersect(names(pal_vec), unique(pd$Color))]

    # Label any hit (depleted or enriched) above the significance + effect-size
    # thresholds. Configurable per-screen if you need to tighten/loosen.
    label_p_max  <- cfg$pandonor$label_p_max  %||% 0.01
    label_z_min  <- cfg$pandonor$label_z_min  %||% 0.5
    labelers <- pd %>%
      dplyr::filter(.data$pval < label_p_max,
                    abs(.data$meanZ) > label_z_min)

    p <- ggplot(pd, aes(x = .data$meanZ, y = -log10(.data$pval),
                        color = .data$Color)) +
      geom_point(alpha = 0.8) +
      { if (length(pal_vec)) scale_color_manual(values = pal_vec) else NULL } +
      ggrepel::geom_text_repel(data = labelers,
                               aes(label = .data$SimpleName),
                               color = "black", max.overlaps = Inf) +
      theme_bw() +
      labs(title = paste0("Pan-donor volcano: ", display),
           x = "Pan-donor mean Z", y = "-log10(p)")
    path <- file.path(opt$`out-dir`, paste0("volcano_", display, ".pdf"))
    ggsave(path, p, width = 6, height = 4)
    message("[partC]  wrote ", path)
  }
}

# ---- concordance plot (all donor pairs, using first pandonor condition) ---

if (isTRUE(cfg$output$plot_concordance) &&
    length(donor_labels) >= 2 && length(pd_conds) >= 1 &&
    length(per_region_by_donor)) {
  ref_pd_cond <- pd_conds[1]
  ref_pd_col  <- paste0("Norm_", ref_pd_cond)

  per_region_long <- bind_rows(lapply(seq_along(per_region_by_donor), function(i) {
    d <- per_region_by_donor[[i]]
    if (!ref_pd_col %in% colnames(d)) return(NULL)
    tibble(SimpleName = d$SimpleName,
           donor      = donor_labels[i],
           value      = d[[ref_pd_col]])
  }))
  if (nrow(per_region_long)) {
    wide <- per_region_long %>%
      pivot_wider(names_from = "donor", values_from = "value") %>%
      add_color_column()

    pal_vec <- unlist(plot_palette)
    pal_vec <- pal_vec[intersect(names(pal_vec), unique(wide$Color))]

    plots <- list()
    pairs <- utils::combn(donor_labels, 2, simplify = FALSE)
    for (p in pairs) {
      a <- p[1]; b <- p[2]
      if (!(a %in% colnames(wide)) || !(b %in% colnames(wide))) next
      sub <- wide[!is.na(wide[[a]]) & !is.na(wide[[b]]), , drop = FALSE]
      if (nrow(sub) < 3) next
      stat <- suppressWarnings(stats::cor.test(sub[[a]], sub[[b]], method = "pearson"))
      gp <- ggplot(sub, aes(x = .data[[a]], y = .data[[b]], color = .data$Color)) +
        geom_point() +
        { if (length(pal_vec)) scale_color_manual(values = pal_vec) else NULL } +
        geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
        annotate("text", x = -Inf, y = Inf,
                 hjust = -0.1, vjust = 1.2,
                 label = sprintf("r=%.3f  p=%.2g",
                                 unname(stat$estimate), stat$p.value)) +
        theme_bw() +
        labs(title = paste0("Concordance: ", a, " vs ", b, " (", ref_pd_cond, ")"))
      plots[[paste(a, b, sep = "_v_")]] <- gp
    }
    if (length(plots)) {
      combined <- patchwork::wrap_plots(plots, ncol = min(3, length(plots)))
      path <- file.path(opt$`out-dir`,
                        paste0("concordance_", ref_pd_cond, ".pdf"))
      ggsave(path, combined,
             width = 5 * min(3, length(plots)),
             height = 4 * ceiling(length(plots) / 3))
      message("[partC]  wrote ", path)
    }
  }
}

# ---- distribution plots (per donor x per highlight category) --------------
# Produces one PDF per category in `output.plot_distribution_categories`
# (each PDF showing all donors stacked vertically, with that category's
# gRNA pairs highlighted in red against gray for everything else).
# This reproduces the three-panel paper figure (NTC / CCR / Promoter).

if (isTRUE(cfg$output$plot_distribution) && !is.null(centered_per_donor)) {
  ref_col   <- paste0("Norm_", ref_cond)
  hl_table  <- cfg$pandonor$highlight_categories %||% list()
  # Categories to plot: explicit list in config, or every highlight category,
  # or fall back to the single legacy "Promoter" category if nothing is set.
  dist_cats <- as_char_vec(cfg$output$plot_distribution_categories %||%
                           names(hl_table) %||% "Promoter")

  for (cat_name in dist_cats) {
    hl_simple <- as_char_vec(hl_table[[cat_name]])
    if (!length(hl_simple)) {
      warning("plot_distribution_categories: '", cat_name,
              "' is not in pandonor.highlight_categories; skipping")
      next
    }
    panels <- list()
    for (i in seq_along(centered_per_donor)) {
      d <- centered_per_donor[[i]]
      if (!ref_col %in% colnames(d)) next
      d$color <- "#D4D4D4"                            # light gray
      d$color[d$SimpleName %in% hl_simple] <- "#B22244"  # raspberry red
      d$alpha <- 0.10
      d$alpha[d$SimpleName %in% hl_simple] <- 1.0
      d <- d[order(d$alpha), ]                        # red drawn on top
      p <- ggplot(d, aes(x = .data[[ref_col]], y = 1)) +
        geom_vline(xintercept = d[[ref_col]],
                   color = d$color, alpha = d$alpha, linewidth = 1) +
        ylim(0, 0.5) +
        scale_x_continuous(limits = c(-4.5, 4.5), n.breaks = 6) +
        theme_minimal() +
        theme(axis.text.y  = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid   = element_blank()) +
        labs(title = paste0("Donor ", i),
             x = paste0("Z-score for each sgRNA pair (", ref_cond, ")"),
             y = NULL)
      panels[[donor_labels[i]]] <- p
    }
    if (length(panels)) {
      combined <- patchwork::wrap_plots(panels, ncol = 1) +
        patchwork::plot_annotation(
          title = cat_name,
          theme = theme(plot.title = element_text(color = "#B22244",
                                                  face = "bold", size = 14)))
      path <- file.path(opt$`out-dir`,
                        paste0("distribution_", ref_cond, "_", cat_name, ".pdf"))
      ggsave(path, combined, width = 6, height = 1.4 * length(panels))
      message("[partC]  wrote ", path)
    }
  }
}

# ---- multi-track distribution plot ----------------------------------------
# One row per requested region. In each row:
#   - GRAY lines: per-region mean Z (across all regions × all donors) -- the
#     background distribution of perturbations.
#   - RED  lines: per-donor mean Z for THAT row's region (so 1 line per donor
#     per row, clustered if the region is reproducible across donors).
# Configured by `multitrack_distribution.regions` (a list of SimpleName values).

mt <- cfg$multitrack_distribution %||% list()
mt_regions <- as_char_vec(mt$regions)

if (isTRUE(cfg$output$plot_multitrack_distribution) && length(mt_regions)) {
  mt_ref_cond <- mt$reference_condition %||% ref_cond
  mt_ref_col  <- paste0("Norm_", mt_ref_cond)

  # Long-format per-region table: SimpleName / donor / value
  all_per_region <- bind_rows(lapply(seq_along(per_region_by_donor), function(i) {
    d <- per_region_by_donor[[i]]
    if (!mt_ref_col %in% colnames(d)) return(NULL)
    tibble(SimpleName = as.character(d$SimpleName),
           donor      = donor_labels[i],
           value      = d[[mt_ref_col]])
  }))

  if (nrow(all_per_region)) {
    # Validate every requested region exists in the data
    missing <- setdiff(mt_regions, unique(all_per_region$SimpleName))
    if (length(missing)) {
      warning("multitrack_distribution.regions not found in data (skipping): ",
              paste(missing, collapse = ", "))
      mt_regions <- intersect(mt_regions, unique(all_per_region$SimpleName))
    }
  }

  if (length(mt_regions)) {
    panels <- lapply(mt_regions, function(sn) {
      d <- all_per_region
      d$color <- "#D4D4D4"
      d$color[d$SimpleName == sn] <- "#B22244"
      d$alpha <- 0.10
      d$alpha[d$SimpleName == sn] <- 1.0
      d <- d[order(d$alpha), ]    # red drawn on top
      ggplot(d, aes(x = .data$value, y = 1)) +
        geom_vline(xintercept = d$value,
                   color = d$color, alpha = d$alpha, linewidth = 1) +
        ylim(0, 0.5) +
        theme_minimal() +
        theme(axis.text.x     = element_blank(),
              axis.ticks.x    = element_blank(),
              axis.text.y     = element_blank(),
              axis.ticks.y    = element_blank(),
              panel.grid      = element_blank(),
              plot.title      = element_text(size = 10, face = "plain")) +
        labs(title = sn, x = NULL, y = NULL)
    })
    # Re-enable x-axis ticks/labels on the bottom panel only
    panels[[length(panels)]] <- panels[[length(panels)]] +
      theme(axis.text.x  = element_text(size = 9),
            axis.ticks.x = element_line()) +
      labs(x = paste0("Per-region mean Z-score (", mt_ref_cond, ")"))

    combined <- patchwork::wrap_plots(panels, ncol = 1) +
      patchwork::plot_annotation(
        title = paste0("Perturbation distribution in ", mt_ref_cond,
                       " across ", length(donor_labels), " donors"),
        theme = theme(plot.title = element_text(size = 12)))

    path <- file.path(opt$`out-dir`,
                      paste0("multitrack_distribution_", mt_ref_cond, ".pdf"))
    ggsave(path, combined, width = 6, height = 0.7 * length(mt_regions) + 1)
    message("[partC]  wrote ", path)
  }
}

message("[partC]  done.")
