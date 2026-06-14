# ============================================================
# 3_1_nucleosomeOccupancy_cd4_new.r
#
# Compare nucleosome occupancy (Micro-C) across species for CD4 cells,
# stratified by beta and chi LRT category.
#
# Output:
#   <result_dir>/all_cd4_beta_comparison.pdf
#   <result_dir>/all_cd4_chi_comparison.pdf
#
# Self-contained: no external R scripts are sourced.
# Dependencies: ggpubr, tidyverse, cowplot.
# ============================================================

library(ggpubr)
library(tidyverse)
library(cowplot)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))

# ── Path configuration ────────────────────────────────────────────────────────
lrt_root   <- file.path(.paths$outputs, 'publish/betweenSpecies/LRT/')
nuc_root   <- file.path(.paths$outputs, 'publish/singleCellType/2_nucleosomePositioning/')
result_dir <- file.path(.paths$outputs, 'publish/betweenSpecies/Nucleosome/')
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)


# ── Plot style ────────────────────────────────────────────────────────────────

beta_color_value <- c(
  "Others" = "gray",
  "Up"     = "#E41A1C",
  "Down"   = "#377EB8"
)

# Copied from helper_function_plot.R; no external dependencies.
my_boxplot_style <- function() {
  list(
    geom_boxplot(outlier.shape = NA, linewidth = 0.3),
    theme_cowplot(),
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 14),
      axis.title.x    = element_text(size = 18),
      axis.text       = element_text(size = 14),
      strip.text      = element_text(size = 14),
      axis.line       = element_line(size = 0.3),
      axis.ticks      = element_line(size = 0.3)
    )
  )
}


# ── Comparison configurations ─────────────────────────────────────────────────
#
# species1 = denominator in nuc_log2FC; species2 = numerator.
# nuc_log2FC = log2(species2_nuc_mean_200 / species1_nuc_mean_200)
# This matches the LRT lfc direction (log2(beta_species2 / beta_species1)).
#
comparisons <- list(

  list(
    name    = "human_rhesus_cd4",
    species1 = "rhesus",
    species2 = "human",
    label   = "Human/Rhesus"
  ),

  list(
    name    = "human_baboon_cd4",
    species1 = "baboon",
    species2 = "human",
    label   = "Human/Baboon"
  ),

  list(
    name    = "rhesus_baboon_cd4",
    species1 = "baboon",
    species2 = "rhesus",
    label   = "Rhesus/Baboon"
  )

)


# ── Validation helpers ────────────────────────────────────────────────────────

check_lrt_cols <- function(df, name) {
  required <- c("tssType", "betaCategory", "chi_lfc_group")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "LRT table '%s' is missing columns: %s\n  Re-run 2_LRT_between_species.R to regenerate.",
      name, paste(missing, collapse = ", ")
    ))
  }
}

check_nuc_obj <- function(nuc_obj, species) {
  if (!"rate_df" %in% names(nuc_obj))
    stop(sprintf("Nucleosome RDS for '%s' has no 'rate_df' element.", species))
  required <- c("geneId", "cd4_nuc_mean_200")
  missing  <- setdiff(required, colnames(nuc_obj$rate_df))
  if (length(missing) > 0)
    stop(sprintf("Nucleosome rate_df for '%s' is missing columns: %s",
                 species, paste(missing, collapse = ", ")))
}


# ── Load nucleosome data (one RDS per species, shared across comparisons) ─────

species_names <- unique(unlist(lapply(comparisons, function(x) c(x$species1, x$species2))))
nuc_data <- setNames(
  lapply(species_names, function(sp) {
    path <- paste0(nuc_root, sp, '/cd4_ns_matrix.RDS')
    obj  <- readRDS(path)
    check_nuc_obj(obj, sp)
    obj
  }),
  species_names
)


# ── Per-comparison processing ─────────────────────────────────────────────────

process_comparison <- function(cfg) {
  # Load and validate LRT table.
  lrt <- read.table(paste0(lrt_root, cfg$name, '_lrt.txt'),
                    sep = '\t', header = TRUE)
  check_lrt_cols(lrt, cfg$name)

  # Also check that species-specific gene ID columns are present.
  id_cols <- c(cfg$species1, cfg$species2)
  missing_ids <- setdiff(id_cols, colnames(lrt))
  if (length(missing_ids) > 0)
    stop(sprintf("LRT table '%s' is missing gene ID columns: %s",
                 cfg$name, paste(missing_ids, collapse = ", ")))

  # Join species1 nucleosome data -> suffix .x (denominator).
  nuc1_df <- nuc_data[[cfg$species1]]$rate_df[, c("geneId", "cd4_nuc_mean_200")]
  lrt <- left_join(lrt, nuc1_df, by = setNames("geneId", cfg$species1))

  # Join species2 nucleosome data -> suffix .y (numerator).
  nuc2_df <- nuc_data[[cfg$species2]]$rate_df[, c("geneId", "cd4_nuc_mean_200")]
  lrt <- left_join(lrt, nuc2_df, by = setNames("geneId", cfg$species2))

  # Keep only rows where both nuc values are present and positive.
  n_before <- nrow(lrt)
  lrt <- lrt[
    !is.na(lrt$cd4_nuc_mean_200.x) & lrt$cd4_nuc_mean_200.x > 0 &
    !is.na(lrt$cd4_nuc_mean_200.y) & lrt$cd4_nuc_mean_200.y > 0, ]

  # nuc_log2FC = log2(species2 / species1), consistent with LRT lfc direction.
  lrt$nuc_log2FC <- log2(lrt$cd4_nuc_mean_200.y / lrt$cd4_nuc_mean_200.x)
  n_after <- nrow(lrt)

  message(sprintf("%s: %d genes kept after filtering (%d removed with missing/non-positive nuc values)",
                  cfg$label, n_after, n_before - n_after))

  # Return slim table for plotting.
  lrt %>%
    dplyr::transmute(
      betaCategory  = factor(betaCategory, levels = c("Down", "Others", "Up")),
      chi_lfc_group = factor(chi_lfc_group,
                             levels = c("Decrease", "Unchanged", "Increase")),
      tss_type      = tssType,
      nuc_log2FC    = nuc_log2FC,
      comparison    = cfg$label
    )
}


# ── Combine all comparisons ───────────────────────────────────────────────────

all_lrt_nucleosome <- bind_rows(lapply(comparisons, process_comparison))

all_lrt_nucleosome$comparison <- factor(
  all_lrt_nucleosome$comparison,
  levels = c("Human/Rhesus", "Human/Baboon", "Rhesus/Baboon")
)


# ── Figure 1: beta category vs nucleosome log2FC ─────────────────────────────

p_beta <- all_lrt_nucleosome %>%
  ggplot(aes(x = betaCategory, y = nuc_log2FC, fill = betaCategory)) +
  facet_grid(tss_type ~ comparison) +
  stat_compare_means(label.x = 1, label.y = 2.5, label = 'p') +
  ylab(expression(Log[2] * "FC(Micro-C)")) +
  labs(x = '', fill = expression(beta * " change")) +
  coord_cartesian(ylim = c(-2, 3)) +
  my_boxplot_style() +
  scale_fill_manual(values = beta_color_value) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    strip.background  = element_blank(),
    legend.position   = 'bottom',
    legend.justification = "center"
  )

ggsave(paste0(result_dir, 'all_cd4_beta_comparison.pdf'), p_beta,
       width = 7, height = 5)


# ── Figure 2: chi lfc group vs nucleosome log2FC ─────────────────────────────

p_chi <- all_lrt_nucleosome %>%
  ggplot(aes(x = chi_lfc_group, y = nuc_log2FC, fill = chi_lfc_group)) +
  facet_grid(tss_type ~ comparison) +
  stat_compare_means(label.x = 1, label.y = 2.5, label = 'p') +
  ylab(expression(Log[2] * "FC(Micro-C)")) +
  labs(x = '', fill = expression(chi * " change")) +
  coord_cartesian(ylim = c(-2, 3)) +
  my_boxplot_style() +
  scale_fill_brewer(palette = "Greens") +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    strip.background  = element_blank(),
    legend.position   = 'bottom',
    legend.justification = "center"
  )

ggsave(paste0(result_dir, 'all_cd4_chi_comparison.pdf'), p_chi,
       width = 7, height = 7)
