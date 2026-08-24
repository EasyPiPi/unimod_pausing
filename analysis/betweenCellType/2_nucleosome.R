# =============================================================================
# Load packages
# =============================================================================
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(cowplot)

# =============================================================================
# Set paths
# =============================================================================
root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/betweenCellType/nucleosomePositioning")

# =============================================================================
# Species configuration
# =============================================================================
# To add a new species (e.g., baboon), add a new entry below with its paths,
# scale factor, and plot parameters. Set has_entropy = FALSE until a $raw
# nucleosome matrix is confirmed available.

species_configs <- list(

  human = list(
    species            = "human",
    cd4_nucleosome_in  = file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/human/cd4_ns_matrix.RDS"),
    cd14_nucleosome_in = file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/human/cd14_ns_matrix.RDS"),
    lrt_result_in      = file.path(.paths$outputs, "publish/betweenCellType/LRT/human/lrt_result.txt"),
    out_subdir         = "human",
    # Normalize CD14 read depth to CD4 library size
    nuc_scale_factor   = 577010542 / 530304146,
    p1_ylim            = c(-1, 2),
    p1_label_y          = 1.7,
    has_entropy         = TRUE,
    entropy_matrix_slot = "raw",
    entropy_cols        = 1001:1201
  ),

  rhesus = list(
    species            = "rhesus",
    cd4_nucleosome_in  = file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/rhesus/cd4_ns_matrix.RDS"),
    cd14_nucleosome_in = file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/rhesus/cd14_ns_matrix.RDS"),
    lrt_result_in      = file.path(.paths$outputs, "publish/betweenCellType/LRT/rhesus/lrt_result.txt"),
    out_subdir         = "rhesus",
    # Rhesus libraries are already size-matched; no scaling needed
    nuc_scale_factor   = 1,
    p1_ylim            = c(-2, 3),
    p1_label_y          = 2.5,
    has_entropy         = TRUE,
    entropy_matrix_slot = "raw",
    entropy_cols        = 1001:1201
  )

)

# =============================================================================
# Helper functions
# =============================================================================

# Match LRT genes to their nucleosome signal matrices, scaling CD14 counts.
get_match_matrix <- function(lrt_df, cd4_rate_df, cd14_rate_df,
                             cd4_matrix, cd14_matrix, scale_factor) {
  lrt_df$index_cd4  <- match(lrt_df$geneId, cd4_rate_df$geneId)
  lrt_df$index_cd14 <- match(lrt_df$geneId, cd14_rate_df$geneId)
  list(
    cd4_lrt_matrix  = cd4_matrix[lrt_df$index_cd4, ],
    cd14_lrt_matrix = cd14_matrix[lrt_df$index_cd14, ] * scale_factor
  )
}

# Normalized entropy: higher value = more diffuse nucleosome signal.
compute_entropy_norm <- function(signal, eps = 1e-8) {
  signal <- pmax(signal, 0)
  if (sum(signal) == 0) return(NA_real_)
  p     <- signal / sum(signal)
  H     <- -sum(p * log(p + eps))
  H_max <- log(length(p))
  H / H_max
}

# Shared boxplot theme used across figures.
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

# =============================================================================
# Analysis functions
# =============================================================================

# Load nucleosome RDS files and LRT table; attach mean signal and compute
# Micro-C log2FC (CD14 vs CD4) for every gene.
prepare_lrt_data <- function(cfg) {
  cd4_nucleosome  <- readRDS(cfg$cd4_nucleosome_in)
  cd14_nucleosome <- readRDS(cfg$cd14_nucleosome_in)
  lrt_df          <- read.table(cfg$lrt_result_in, sep = "\t", header = TRUE)

  lrt_df <- left_join(lrt_df,
                      cd4_nucleosome$rate_df[, c("geneId", "cd4_nuc_mean_200")],
                      by = "geneId")
  lrt_df <- left_join(lrt_df,
                      cd14_nucleosome$rate_df[, c("geneId", "cd14_nuc_mean_200")],
                      by = "geneId")

  lrt_df$nuc_log2FC <- log2(lrt_df$cd14_nuc_mean_200 * cfg$nuc_scale_factor /
                              lrt_df$cd4_nuc_mean_200)
  lrt_df <- lrt_df[!is.infinite(lrt_df$nuc_log2FC), ]
  lrt_df <- lrt_df[!is.na(lrt_df$nuc_log2FC), ]

  lrt_df <- lrt_df %>%
    mutate(
      chi_lfc_group = factor(chi_lfc_group, levels = c("Decrease", "Unchanged", "Increase"))
    ) %>%
    droplevels()

  list(lrt_df = lrt_df, cd4_nucleosome = cd4_nucleosome, cd14_nucleosome = cd14_nucleosome)
}

# Plot Micro-C log2FC grouped by expression-change category (chi group).
plot_chi_boxplot <- function(lrt_df, cfg, result_dir) {
  out_path <- file.path(result_dir, cfg$out_subdir, "chi_boxplot_allgenes.pdf")

  p <- lrt_df %>%
    ggplot(aes(x = chi_lfc_group, y = nuc_log2FC, fill = chi_lfc_group)) +
    stat_compare_means(label.x = 0.8, label.y = 1.2, label = "p") +
    ylab(expression(Log[2] * "FC(Micro-C)")) +
    xlab("") +
    coord_cartesian(ylim = c(-0.5, 1.4)) +
    my_boxplot_style() +
    scale_fill_brewer(palette = "Greens") +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1),
          axis.title.x = element_blank())

  ggsave(out_path, p, width = 2.5, height = 4)
  invisible(p)
}

# Plot Micro-C log2FC by beta-change group, faceted by expression group.
# y-axis limits and label position are species-specific (from cfg).
plot_beta_by_chi_boxplot <- function(lrt_df, cfg, result_dir) {
  out_path <- file.path(result_dir, cfg$out_subdir, "beta_groupby_chi_boxplot_allgenes.pdf")

  # Recode factor levels to full labels for facet display
  lrt_df$chi_lfc_group <- recode(
    lrt_df$chi_lfc_group,
    "Decrease"  = "Decreased expression",
    "Unchanged" = "Unchanged expression",
    "Increase"  = "Increased expression"
  )

  p <- lrt_df %>%
    ggplot(aes(x = betaCategory, y = nuc_log2FC, fill = betaCategory)) +
    geom_violin(width = 0.9, alpha = 0.5, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", linewidth = 0.3) +
    facet_wrap(~ chi_lfc_group, ncol = 3) +
    stat_compare_means(label.x = 1, label.y = cfg$p1_label_y, label = "p") +
    ylab(expression(Log[2] * "FC (Micro-C signal)")) +
    coord_cartesian(ylim = cfg$p1_ylim) +
    scale_fill_manual(values = c("Others" = "#BFBFBF",
                                 "Up"     = "#E45756",
                                 "Down"   = "#4C78A8")) +
    labs(fill = expression(beta * " change")) +
    xlab("") +
    cowplot::theme_cowplot() +
    theme(
      strip.background     = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      axis.text.x          = element_blank(),
      axis.ticks.x         = element_blank(),
      legend.position      = "bottom",
      legend.justification = "center"
    )

  ggsave(out_path, p, width = 8, height = 3)
  invisible(p)
}

# Entropy-based +1 nucleosome positioning analysis.
# Matrix slot and column range are taken from cfg (entropy_matrix_slot, entropy_cols).
run_entropy_analysis <- function(lrt_df, cd4_nucleosome, cd14_nucleosome, cfg, result_dir) {
  out_path <- file.path(result_dir, cfg$out_subdir, "+1_nucleosome_positioning_raw.pdf")

  cd4_mat  <- cd4_nucleosome[[cfg$entropy_matrix_slot]][, cfg$entropy_cols]
  cd14_mat <- cd14_nucleosome[[cfg$entropy_matrix_slot]][, cfg$entropy_cols]

  lrt_match_matrix <- get_match_matrix(
    lrt_df,
    cd4_nucleosome$rate_df,
    cd14_nucleosome$rate_df,
    cd4_mat,
    cd14_mat,
    cfg$nuc_scale_factor
  )

  cd4_entropy  <- apply(lrt_match_matrix$cd4_lrt_matrix,  1, compute_entropy_norm)
  cd14_entropy <- apply(lrt_match_matrix$cd14_lrt_matrix, 1, compute_entropy_norm)

  # Positive log2FC = sharper (more focused) nucleosome positioning in CD14
  lrt_df$nuc_positioning_lfc <- log2((1 - cd14_entropy) / (1 - cd4_entropy))

  lrt_df$pausing_switch <- paste0(lrt_df$CD4_sdGroup, "_to_", lrt_df$CD14_sdGroup)
  # Broad_to_Sharp is excluded: the transition is not informative for positioning
  lrt_df_filtered <- lrt_df[lrt_df$pausing_switch != "Broad_to_Sharp", ]

  p <- lrt_df_filtered %>%
    ggplot(aes(x = pausing_switch, y = nuc_positioning_lfc, fill = pausing_switch)) +
    stat_compare_means(label.x = 0.8, label.y = 0.6, label = "p") +
    ylab(expression(Log[2] * "FC(1-Entropy)")) +
    xlab("") +
    coord_cartesian(ylim = c(-0.5, 0.8)) +
    my_boxplot_style() +
    scale_fill_manual(
      name   = expression(sigma ~ "change"),
      values = c("Broad_to_Broad" = "grey",
                 "Sharp_to_Broad" = "#9E77CF",
                 "Sharp_to_Sharp" = "grey")
    ) +
    theme(axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          legend.position = "right")

  ggsave(out_path, p, width = 5, height = 4)
  invisible(p)
}

# Run the full analysis pipeline for one species config.
run_species_analysis <- function(cfg, result_dir) {
  message("Running nucleosome positioning analysis for: ", cfg$species)

  data            <- prepare_lrt_data(cfg)
  lrt_df          <- data$lrt_df
  cd4_nucleosome  <- data$cd4_nucleosome
  cd14_nucleosome <- data$cd14_nucleosome

  plot_chi_boxplot(lrt_df, cfg, result_dir)
  plot_beta_by_chi_boxplot(lrt_df, cfg, result_dir)

  if (cfg$has_entropy) {
    run_entropy_analysis(lrt_df, cd4_nucleosome, cd14_nucleosome, cfg, result_dir)
  }

  invisible(NULL)
}

# =============================================================================
# Run analysis for all species
# =============================================================================
for (cfg in species_configs) {
  run_species_analysis(cfg, result_dir)
}
