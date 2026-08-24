# ============================================================
# Visualize transcription rates for human, rhesus, and baboon.
#
# Reads the STADyUM rate RDS files produced by the processing
# script and generates main-figure and supplementary PDFs.
# ============================================================


# ============================================================
# LIBRARIES
# ============================================================

library(STADyUM)
library(tidyverse)
library(ggplot2)
library(ggpointdensity)
library(cowplot)
library(patchwork)


# ============================================================
# CONFIGURATION
# ============================================================

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", "."),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/")

# One row per sample; file paths point to the RDS objects saved by the
# processing script.
rate_info <- tibble(
  species  = c("Human",        "Human",           "Rhesus",       "Rhesus",          "Baboon",       "Baboon"),
  celltype = c("CD4+ T cell",  "CD14+ monocyte",  "CD4+ T cell",  "CD14+ monocyte",  "CD4+ T cell",  "CD14+ monocyte"),
  file     = c(
    file.path(result_dir, "human/cd4_rate.RDS"),
    file.path(result_dir, "human/cd14_rate.RDS"),
    file.path(result_dir, "rhesus/cd4_rate.RDS"),
    file.path(result_dir, "rhesus/cd14_rate.RDS"),
    file.path(result_dir, "baboon/cd4_rate.RDS"),
    file.path(result_dir, "baboon/cd14_rate.RDS")
  )
)


# ============================================================
# HELPER FUNCTIONS
# Copied from helper_function_plot.R — only the two functions
# used by this script are included.
# ============================================================

# Scatter plot of log10(beta) vs log10(chi) with Spearman rho annotation.
plot_beta_vs_chi <- function(
    rate_df,
    label_x,
    label_y,
    xlim      = NULL,
    ylim      = NULL,
    cor_size  = 5
) {
  rho <- cor(
    log10(rate_df$betaAdp),
    log10(rate_df$chi),
    method = "spearman",
    use    = "complete.obs"
  )

  p <- ggplot(rate_df, aes(x = log10(betaAdp), y = log10(chi))) +
    geom_pointdensity(adjust = 0.3, size = 0.6) +
    scale_color_viridis_c() +
    theme_cowplot() +
    labs(
      x     = expression(log[10] ~ beta),
      y     = expression(log[10] ~ chi),
      color = "Density"
    ) +
    annotate(
      "text",
      x     = label_x,
      y     = label_y,
      label = paste0("rho = ", round(rho, 2)),
      hjust = 0,
      size  = cor_size
    ) +
    theme(
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 14),
      strip.text       = element_text(size = 14),
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm"),
      legend.spacing   = unit(0.2, "cm")
    )

  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  }
  p
}

# Scatter plot of fkMean (mu) vs fkSD (sigma) with Spearman rho annotation.
plot_fkMean_vs_fkSD <- function(
    rate_df,
    cor_size = 5,
    label_x  = NULL,
    label_y  = NULL
) {
  if (is.null(label_x)) label_x <- quantile(rate_df$fkMean, 0.02, na.rm = TRUE)
  if (is.null(label_y)) label_y <- quantile(rate_df$fkSD,   0.98, na.rm = TRUE)

  rho <- cor(rate_df$fkMean, rate_df$fkSD, method = "spearman", use = "complete.obs")

  ggplot(rate_df, aes(x = fkMean, y = fkSD)) +
    geom_pointdensity(adjust = 0.3, size = 0.6) +
    scale_color_viridis_c() +
    theme_cowplot() +
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.2)) +
    labs(
      x     = expression(mu),
      y     = expression(sigma),
      color = "Density"
    ) +
    annotate(
      "text",
      x     = -Inf,
      y     = Inf,
      label = paste0("rho = ", round(rho, 2)),
      hjust = -0.1,
      vjust = 1.2,
      size  = cor_size
    ) +
    coord_cartesian(clip = "off") +
    theme(
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 14),
      strip.text       = element_text(size = 14),
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 8),
      legend.key.size  = unit(0.4, "cm"),
      legend.spacing   = unit(0.2, "cm"),
      plot.margin      = margin(t = 8, r = 5, b = 5, l = 5)
    )
}


# ============================================================
# CUSTOM PLOTTING FUNCTIONS
# Defined in the original visualization script.
# ============================================================

# Add a centered title to any ggplot.
add_plot_title <- function(p, title, hjust = 0.5, face = "plain") {
  p +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = hjust, face = face))
}

# 2D density plot of log10(beta) vs sigma with Spearman rho annotation.
# A horizontal dashed line marks the Sharp/Broad boundary if sdGroup is present.
plot_beta_vs_sigma_density <- function(
    rate_df,
    xlim      = c(-5, -2),
    ylim      = c(0, 70),
    bins      = 6,
    alpha     = 0.8,
    text_size = 14,
    cor_size  = 5,
    label_x   = NULL,
    label_y   = NULL
) {
  if (is.null(label_x)) label_x <- xlim[1] + 0.01 * diff(xlim)
  if (is.null(label_y)) label_y <- ylim[2] - 0.01 * diff(ylim)

  # Compute Sharp/Broad boundary as midpoint between the two groups.
  sd_threshold <- NULL
  if (!"sdGroup" %in% names(rate_df)) {
    warning("plot_beta_vs_sigma_density: 'sdGroup' column not found; skipping boundary line.")
  } else {
    sharp_max <- max(rate_df$fkSD[rate_df$sdGroup == "Sharp"], na.rm = TRUE)
    broad_min <- min(rate_df$fkSD[rate_df$sdGroup == "Broad"], na.rm = TRUE)
    if (!is.finite(sharp_max) || !is.finite(broad_min)) {
      warning("plot_beta_vs_sigma_density: Sharp or Broad group is absent or all-NA; skipping boundary line.")
    } else {
      sd_threshold <- (sharp_max + broad_min) / 2
    }
  }

  rho <- cor(
    log10(rate_df$betaAdp),
    rate_df$fkSD,
    method = "spearman",
    use    = "complete.obs"
  )

  p <- ggplot(rate_df, aes(x = log10(betaAdp), y = fkSD)) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom  = "polygon",
      bins  = bins,
      alpha = alpha
    ) +
    annotate(
      "text",
      x     = label_x,
      y     = label_y,
      label = paste0("rho = ", round(rho, 2)),
      hjust = 0,
      size  = cor_size
    ) +
    scale_fill_viridis_c(name = "Density") +
    labs(
      x = expression(log[10] ~ beta),
      y = expression(sigma)
    ) +
    cowplot::theme_cowplot() +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme(text = element_text(size = text_size))

  if (!is.null(sd_threshold)) {
    p <- p + geom_hline(yintercept = sd_threshold, linetype = "dashed")
  }
  p
}


# ============================================================
# DATA LOADING
# ============================================================

rate_info <- rate_info %>%
  mutate(
    rate_obj = map(file, readRDS),
    rate_df  = map(rate_obj, rates),
    title    = paste(species, celltype)
  )


# ============================================================
# MAIN FIGURE PLOTS
# ============================================================

# --- Main: beta vs chi (Human CD4 only) ---
p_main <- plot_beta_vs_chi(
  rate_info$rate_df[[1]],
  label_x = -3.5, label_y = -3.8, c(-6, -1), c(-4, 1)
)
ggsave(file.path(result_dir, "human/cd4_beta_chi_scatter_new.pdf"),
       p_main, width = 4, height = 3.5)

# --- Main: beta vs sigma density (Human CD4 only) ---
p_main <- plot_beta_vs_sigma_density(rate_info$rate_df[[1]])
ggsave(file.path(result_dir, "human/beta_fk_new.pdf"),
       p_main, width = 4.5, height = 4)


# ============================================================
# SUPPLEMENTARY FIGURE PLOTS
# ============================================================

# --- Supplementary: beta vs chi (all except Human CD4) ---
# Label positions and x-axis range differ per sample, so each is explicit.
p2 <- plot_beta_vs_chi(rate_info$rate_df[[2]], label_x = -3.5, label_y = -3.5, c(-6, -1), c(-4, 1)) %>% add_plot_title(rate_info$title[[2]])
p3 <- plot_beta_vs_chi(rate_info$rate_df[[3]], label_x = -3.5, label_y = -3.6, c(-6, -1), c(-4, 1)) %>% add_plot_title(rate_info$title[[3]])
p4 <- plot_beta_vs_chi(rate_info$rate_df[[4]], label_x = -3.5, label_y = -3.7, c(-6, -1), c(-4, 1)) %>% add_plot_title(rate_info$title[[4]])
p5 <- plot_beta_vs_chi(rate_info$rate_df[[5]], label_x = -3,   label_y = -3.7, c(-5,  0), c(-4, 1)) %>% add_plot_title(rate_info$title[[5]])
p6 <- plot_beta_vs_chi(rate_info$rate_df[[6]], label_x = -2.2, label_y = -3.7, c(-5,  0), c(-4, 1)) %>% add_plot_title(rate_info$title[[6]])

supp_beta_chi <- (p2 | p3 | p4) / (p5 | p6 | plot_spacer())
ggsave(file.path(result_dir, "beta_chi_scatter_all.pdf"),
       supp_beta_chi, width = 10, height = 5)

# --- Supplementary: mu vs sigma (all 6 samples) ---
mu_sigma_plots <- map2(
  rate_info$rate_df, rate_info$title,
  ~add_plot_title(plot_fkMean_vs_fkSD(.x), .y)
)
supp_mu_sigma <- (mu_sigma_plots[[1]] | mu_sigma_plots[[2]] | mu_sigma_plots[[3]]) /
                 (mu_sigma_plots[[4]] | mu_sigma_plots[[5]] | mu_sigma_plots[[6]])
ggsave(file.path(result_dir, "mu_sigma_scatter_all.pdf"),
       supp_mu_sigma, width = 10, height = 5)

# --- Supplementary: beta vs sigma density (all except Human CD4) ---
beta_sigma_supp <- map2(
  rate_info$rate_df[-1], rate_info$title[-1],
  ~add_plot_title(plot_beta_vs_sigma_density(.x), .y)
)
supp_beta_sigma <- (beta_sigma_supp[[1]] | beta_sigma_supp[[2]] | beta_sigma_supp[[3]]) /
                   (beta_sigma_supp[[4]] | beta_sigma_supp[[5]] | plot_spacer())
ggsave(file.path(result_dir, "beta_sigma_2d_all.pdf"),
       supp_beta_sigma, width = 12, height = 7)
