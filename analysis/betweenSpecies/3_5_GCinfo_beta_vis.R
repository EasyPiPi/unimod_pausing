library(tidyverse)
library(plyranges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(genomation)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(STADyUM)
library(patchwork)
library(UpSetR)

# ---- Paths ---------------------------------------------------------------

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", "."),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

lrt_dir    <- file.path(.paths$outputs, "publish/betweenSpecies/LRT")
seq_dir    <- file.path(.paths$outputs, "publish/betweenSpecies/sequence")
result_dir <- file.path(.paths$outputs, "publish/betweenSpecies")

# ---- Helper: read one GC info CSV ----------------------------------------

read_gc_info <- function(path, species) {
  df <- read.table(path, sep = ",")
  colnames(df) <- c(species, paste0(species, "_gc_content"), paste0(species, "_gc_skew"))
  df
}

# ---- Read GC info files --------------------------------------------------

gc <- list(
  human_cd4   = read_gc_info(file.path(seq_dir, "human_cd4_tss_sort_1000bp_extend_200_gc_info.csv"),   "human"),
  human_cd14  = read_gc_info(file.path(seq_dir, "human_cd14_tss_sort_1000bp_extend_200_gc_info.csv"),  "human"),
  rhesus_cd4  = read_gc_info(file.path(seq_dir, "rhesus_cd4_tss_sort_1000bp_extend_200_gc_info.csv"),  "rhesus"),
  rhesus_cd14 = read_gc_info(file.path(seq_dir, "rhesus_cd14_tss_sort_1000bp_extend_200_gc_info.csv"), "rhesus"),
  baboon_cd4  = read_gc_info(file.path(seq_dir, "baboon_cd4_tss_sort_1000bp_extend_200_gc_info.csv"),  "baboon"),
  baboon_cd14 = read_gc_info(file.path(seq_dir, "baboon_cd14_tss_sort_1000bp_extend_200_gc_info.csv"), "baboon")
)

# ---- Core processing function --------------------------------------------

prep_lrt_effect_df <- function(
  lrt_path,
  gc_info_1, gc_info_2,
  key_1, key_2,
  sp1_prefix, sp2_prefix,
  cell_type_label,
  tss_col    = "tssType",
  comparison_label = NULL,
  beta_col   = "betaCategory",
  drop_level = "Up",
  keep_cols  = c("geneId", "gc_content_lfc", "delta_gc_skew", "betaCategory")
) {
  keep_cols <- c(keep_cols, tss_col)

  if (!is.null(drop_level) && !(drop_level %in% c("Up", "Down"))) {
    stop("drop_level must be NULL, 'Up', or 'Down'.")
  }
  conflicted::conflict_prefer("select", "dplyr")

  df <- read.table(lrt_path, sep = "\t", header = TRUE)

  required_cols <- c("geneId", "betaCategory", "tssType", key_1, key_2)
  missing_lrt <- base::setdiff(required_cols, colnames(df))
  if (length(missing_lrt) > 0) {
    stop("LRT file missing required columns: ", paste(missing_lrt, collapse = ", "),
         "\n  File: ", lrt_path)
  }

  df <- df %>%
    left_join(gc_info_1, by = setNames(key_1, key_1)) %>%
    left_join(gc_info_2, by = setNames(key_2, key_2))

  gc1_content <- paste0(sp1_prefix, "_gc_content")
  gc2_content <- paste0(sp2_prefix, "_gc_content")
  gc1_skew    <- paste0(sp1_prefix, "_gc_skew")
  gc2_skew    <- paste0(sp2_prefix, "_gc_skew")

  missing_cols <- base::setdiff(c(gc1_content, gc2_content, gc1_skew, gc2_skew), colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required GC columns after join: ", paste(missing_cols, collapse = ", "))
  }

  if (is.null(comparison_label)) comparison_label <- paste0(sp1_prefix, "/", sp2_prefix)

  out <- df %>%
    mutate(
      gc_content_lfc = log2(.data[[gc1_content]] / .data[[gc2_content]]),
      delta_gc_skew  = .data[[gc1_skew]] - .data[[gc2_skew]]
    ) %>%
    dplyr::select(all_of(keep_cols)) %>%
    dplyr::rename(tss_type = all_of(tss_col)) %>%
    dplyr::mutate(cell_type = cell_type_label, comparison = comparison_label)

  if (!is.null(drop_level)) {
    out <- out[out[[beta_col]] != drop_level, , drop = FALSE]
    changed_level <- base::setdiff(c("Up", "Down"), drop_level)[1]
  } else {
    changed_level <- NULL
  }

  out <- out %>%
    mutate(
      !!beta_col := case_when(
        !is.null(changed_level) & .data[[beta_col]] == changed_level ~ "Changed",
        .data[[beta_col]] == "Others"                                 ~ "Unchanged",
        TRUE                                                          ~ .data[[beta_col]]
      )
    )

  non_conserved <- out[out[["tss_type"]] == "non-Conserved TSS", , drop = FALSE]
  conserved     <- out[out[["tss_type"]] == "Conserved TSS",     , drop = FALSE]
  non_conserved[["tss_type"]] <- NULL
  conserved[["tss_type"]]     <- NULL

  list(filtered = out, non_conserved = non_conserved, conserved = conserved)
}

# ---- Comparison configs --------------------------------------------------
#
# drop_level = "Up" for all: keeps the "Down" direction as Changed.
# For rhesus/baboon, the new files already have rhesus as sp1, so no sign
# flip is needed — the direction and gc_content_lfc are already correct.

comparisons <- list(
  list(lrt = "human_rhesus_cd4_lrt.txt",   sp1 = "human",  sp2 = "rhesus", ct_key = "cd4",  cell_type = "CD4+ T cells",    label = "Human/Rhesus"),
  list(lrt = "human_rhesus_cd14_lrt.txt",  sp1 = "human",  sp2 = "rhesus", ct_key = "cd14", cell_type = "CD14+ Monocytes", label = "Human/Rhesus"),
  list(lrt = "human_baboon_cd4_lrt.txt",   sp1 = "human",  sp2 = "baboon", ct_key = "cd4",  cell_type = "CD4+ T cells",    label = "Human/Baboon"),
  list(lrt = "human_baboon_cd14_lrt.txt",  sp1 = "human",  sp2 = "baboon", ct_key = "cd14", cell_type = "CD14+ Monocytes", label = "Human/Baboon"),
  list(lrt = "rhesus_baboon_cd4_lrt.txt",  sp1 = "rhesus", sp2 = "baboon", ct_key = "cd4",  cell_type = "CD4+ T cells",    label = "Rhesus/Baboon"),
  list(lrt = "rhesus_baboon_cd14_lrt.txt", sp1 = "rhesus", sp2 = "baboon", ct_key = "cd14", cell_type = "CD14+ Monocytes", label = "Rhesus/Baboon")
)

results <- lapply(comparisons, function(cfg) {
  prep_lrt_effect_df(
    lrt_path         = file.path(lrt_dir, cfg$lrt),
    gc_info_1        = gc[[paste0(cfg$sp1, "_", cfg$ct_key)]],
    gc_info_2        = gc[[paste0(cfg$sp2, "_", cfg$ct_key)]],
    key_1            = cfg$sp1,
    key_2            = cfg$sp2,
    sp1_prefix       = cfg$sp1,
    sp2_prefix       = cfg$sp2,
    cell_type_label  = cfg$cell_type,
    comparison_label = cfg$label,
    drop_level       = "Up"
  )
})

df_all <- bind_rows(lapply(results, `[[`, "filtered"))

# ---- Plotting helpers ----------------------------------------------------

celltype_color     <- c("#6d6e71", "#d1d3d4")
names(celltype_color) <- c("CD4+ T cells", "CD14+ Monocytes")
comparison_levels  <- c("Human/Rhesus", "Rhesus/Baboon", "Human/Baboon")

make_effect_df <- function(data, metric_col) {
  data %>%
    group_by(tss_type, comparison, cell_type) %>%
    summarise(
      mean_changed   = mean(.data[[metric_col]][betaCategory == "Changed"],   na.rm = TRUE),
      mean_unchanged = mean(.data[[metric_col]][betaCategory == "Unchanged"], na.rm = TRUE),
      n_changed      = sum(betaCategory == "Changed"),
      n_unchanged    = sum(betaCategory == "Unchanged"),
      sd_changed     = sd(.data[[metric_col]][betaCategory == "Changed"],   na.rm = TRUE),
      sd_unchanged   = sd(.data[[metric_col]][betaCategory == "Unchanged"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      effect_size = mean_changed - mean_unchanged,
      se          = sqrt(sd_changed^2 / n_changed + sd_unchanged^2 / n_unchanged),
      ci_low      = effect_size - 1.96 * se,
      ci_high     = effect_size + 1.96 * se,
      comparison  = factor(comparison, levels = comparison_levels)
    )
}

make_effect_plot <- function(effect_df, y_limits, y_label) {
  ggplot(effect_df, aes(x = comparison, y = effect_size, color = cell_type)) +
    facet_wrap(~ tss_type) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(width = 0.4), size = 2) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.2, position = position_dodge(width = 0.4)) +
    scale_color_manual(values = celltype_color) +
    theme_cowplot() +
    scale_y_continuous(limits = y_limits) +
    coord_flip() +
    ylab(y_label) +
    xlab("") +
    theme(
      legend.title     = element_blank(),
      strip.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
      legend.position  = "top",
      strip.text       = element_text(size = 10),
      axis.text        = element_text(size = 12)
    )
}

# ---- GC content plot -----------------------------------------------------

p_content <- make_effect_plot(
  make_effect_df(df_all, "gc_content_lfc"),
  y_limits = c(-0.1, 0.5),
  y_label  = expression(Delta * "log"[2] * "FC(GC content)")
)
ggsave(file.path(result_dir, "all_beta_gc_content.pdf"), plot = p_content, width = 8, height = 3.5)

# ---- GC skew plot --------------------------------------------------------

p_skew <- make_effect_plot(
  make_effect_df(df_all, "delta_gc_skew"),
  y_limits = c(-0.1, 0.1),
  y_label  = expression(Delta * "GC skew")
)
ggsave(file.path(result_dir, "all_beta_gc_skew.pdf"), plot = p_skew, width = 8, height = 3.5)

### THE END
