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
  tss_col          = "tssType",
  comparison_label = NULL,
  keep_cols        = c("geneId", "gc_content_lfc", "delta_gc_skew", "pause_change")
) {
  conflicted::conflict_prefer("select", "dplyr")

  df <- read.table(lrt_path, sep = "\t", header = TRUE)

  required_cols <- c("geneId", "pause_change", "tssType", key_1, key_2)
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

  missing_gc <- base::setdiff(c(gc1_content, gc2_content, gc1_skew, gc2_skew), colnames(df))
  if (length(missing_gc) > 0) {
    stop("Missing required GC columns after join: ", paste(missing_gc, collapse = ", "))
  }

  if (is.null(comparison_label)) comparison_label <- paste0(sp1_prefix, "/", sp2_prefix)

  out <- df %>%
    mutate(
      gc_content_lfc = log2(.data[[gc1_content]] / .data[[gc2_content]]),
      delta_gc_skew  = .data[[gc1_skew]] - .data[[gc2_skew]]
    ) %>%
    dplyr::select(all_of(c(keep_cols, tss_col))) %>%
    dplyr::rename(tss_type = all_of(tss_col)) %>%
    dplyr::mutate(cell_type = cell_type_label, comparison = comparison_label)

  out$pause_change[out$pause_change %in% c("Broad_to_Broad", "Sharp_to_Sharp")] <- "Unchanged"

  out
}

# ---- Comparison configs --------------------------------------------------
#
# All six comparisons use the same tssType column and drop_level is not
# applicable here. The rhesus_baboon files are already oriented as
# Rhesus/Baboon — no sign flip or pause_change swap needed.

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
    comparison_label = cfg$label
  )
})

df_all <- bind_rows(results)

# ---- Plotting helpers ----------------------------------------------------

chang_color <- c("#FC8D62", "#66C2A5")
names(chang_color) <- c("Broad_to_Sharp", "Sharp_to_Broad")
comparison_levels <- c("Human/Rhesus", "Rhesus/Baboon", "Human/Baboon")

make_effect_df <- function(data, metric_col) {
  data %>%
    dplyr::filter(pause_change %in% c("Broad_to_Sharp", "Sharp_to_Broad", "Unchanged")) %>%
    group_by(tss_type, comparison, cell_type, pause_change) %>%
    summarise(
      mean_gc = mean(.data[[metric_col]], na.rm = TRUE),
      sd_gc   = sd(.data[[metric_col]],   na.rm = TRUE),
      n       = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = pause_change, values_from = c(mean_gc, sd_gc, n)) %>%
    mutate(
      effect_BS  = mean_gc_Broad_to_Sharp - mean_gc_Unchanged,
      se_BS      = sqrt(sd_gc_Broad_to_Sharp^2 / n_Broad_to_Sharp + sd_gc_Unchanged^2 / n_Unchanged),
      effect_SB  = mean_gc_Sharp_to_Broad  - mean_gc_Unchanged,
      se_SB      = sqrt(sd_gc_Sharp_to_Broad^2  / n_Sharp_to_Broad  + sd_gc_Unchanged^2 / n_Unchanged),
      ci_low_BS  = effect_BS - 1.96 * se_BS,
      ci_high_BS = effect_BS + 1.96 * se_BS,
      ci_low_SB  = effect_SB - 1.96 * se_SB,
      ci_high_SB = effect_SB + 1.96 * se_SB
    )
}

make_plot_df <- function(effect_df) {
  bind_rows(
    effect_df %>% transmute(tss_type, comparison, cell_type,
                            contrast    = "Broad_to_Sharp",
                            effect_size = effect_BS, se = se_BS,
                            ci_low = ci_low_BS, ci_high = ci_high_BS),
    effect_df %>% transmute(tss_type, comparison, cell_type,
                            contrast    = "Sharp_to_Broad",
                            effect_size = effect_SB, se = se_SB,
                            ci_low = ci_low_SB, ci_high = ci_high_SB)
  ) %>%
    mutate(
      contrast   = factor(contrast,   levels = c("Broad_to_Sharp", "Sharp_to_Broad")),
      comparison = factor(comparison, levels = comparison_levels)
    )
}

make_effect_plot <- function(plot_df, y_label, y_limits = NULL, y_breaks = waiver()) {
  ggplot(plot_df, aes(x = comparison, y = effect_size, color = contrast)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.2, position = position_dodge(width = 0.4)) +
    facet_grid(rows = vars(cell_type), cols = vars(tss_type)) +
    coord_flip() +
    labs(x = "", y = y_label) +
    scale_color_manual(values = chang_color, name = expression(sigma ~ change)) +
    theme_cowplot() +
    theme(
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
      legend.position  = "top",
      strip.placement  = "outside",
      strip.text.y     = element_text(angle = 270),
      strip.text       = element_text(size = 10)
    ) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks)
}

# ---- GC content plot -----------------------------------------------------

p_content <- make_effect_plot(
  make_plot_df(make_effect_df(df_all, "gc_content_lfc")),
  y_label  = expression(Delta * "log"[2] * "FC(GC content)"),
  y_limits = c(-0.3, 0.3)
)
ggsave(file.path(result_dir, "all_sigma_gc_content.pdf"), plot = p_content, width = 8, height = 4)

# ---- GC skew plot --------------------------------------------------------

p_skew <- make_effect_plot(
  make_plot_df(make_effect_df(df_all, "delta_gc_skew")),
  y_label  = expression(Delta * "GC skew"),
  y_limits = c(-0.1, 0.1),
  y_breaks = c(-0.05, 0, 0.05)
)
ggsave(file.path(result_dir, "all_sigma_gc_skew.pdf"), plot = p_skew, width = 8, height = 4)

### THE END
