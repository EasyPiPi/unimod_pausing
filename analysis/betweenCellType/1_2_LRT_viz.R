library(tidyverse)
library(cowplot)
library(ggpubr)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", "."),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/betweenCellType/LRT/")

# ==============================================================
# PLOTTING FUNCTIONS
# ==============================================================

##
plot_beta_quantile_heatmap <- function(df,
                                       cd4_col,
                                       cd14_col,
                                       fill_label = "Gene count",
                                       show_text = TRUE,
                                       text_size = 5,
                                       show_legend = FALSE) {
    
    mat <- table(df[[cd4_col]], df[[cd14_col]])
    plot_df <- as.data.frame(mat)
    
    p <- ggplot(plot_df, aes(x = Var1, y = Var2, fill = Freq)) +
        geom_tile()
    
    if (show_text) {
        p <- p + geom_text(aes(label = Freq), size = text_size)
    }
    
    p +
        scale_fill_viridis_c() +
        labs(
            x = expression(beta[CD4]),
            y = expression(beta[CD14]),
            fill = fill_label
        ) +
        coord_equal() +
        cowplot::theme_cowplot() +
        theme(legend.position = if (show_legend) "right" else "none")
}


# lrt_beta_scatter: log10(beta1) vs log10(beta2), colored by betaCategory
plot_beta_scatter <- function(beta_tbl,
                              label_x = -6,
                              label_y = -1,
                              x_lim   = c(-6, -1),
                              y_lim   = c(-6, -1)) {
  ggplot(beta_tbl, aes(x = log10(beta1), y = log10(beta2))) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", color = "gray60", linewidth = 0.8) +
    geom_point(aes(color = betaCategory), alpha = 0.3, size = 0.5) +
    scale_color_manual(values = c(
      "Others" = "gray",
      "Up"     = "#E41A1C",
      "Down"   = "#377EB8"
    )) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(
      x     = expression(log[10] * "(" * beta[CD4] * ")"),
      y     = expression(log[10] * "(" * beta[CD14] * ")"),
      color = expression(beta * " change")
    ) +
    stat_cor(
      aes(label = gsub("R", "rho", after_stat(r.label))),
      label.x = label_x,
      label.y = label_y
    ) +
    theme_cowplot() +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12),
      axis.line  = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3)
    )
}



# deltaBeta_deltaSigma: density scatter of lfc vs deltaSD with Spearman rho
# filter_lfc: if TRUE, restricts to abs(lfc) < 5 before computing rho and plotting (human only)
plot_delta_beta_sigma <- function(beta_tbl, filter_lfc = FALSE) {
  if (filter_lfc) {
    beta_tbl <- beta_tbl[abs(beta_tbl$lfc) < 5, ]
  }

  rho <- cor(beta_tbl$lfc, beta_tbl$deltaSD,
             method = "spearman", use = "complete.obs")

  beta_tbl %>%
    ggplot(aes(x = lfc, y = deltaSD)) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom  = "polygon",
      bins  = 8,
      alpha = 0.8
    ) +
    scale_fill_viridis_c(option = "C", name = "Density") +
    annotate(
      "text",
      x     = Inf, y = Inf,
      label = paste0("Spearman rho = ", round(rho, 3)),
      hjust = 1.1, vjust = 1.5,
      size  = 4
    ) +
    theme_cowplot() +
    labs(
      x = expression(Delta * beta),
      y = expression(Delta * sigma)
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    coord_cartesian(xlim = c(-5, 5)) +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12),
      axis.line  = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3)
    )
}


# fksd_density: 2D density of CD4_fkSD vs CD14_fkSD with quadrant labels
plot_fksd_density <- function(lrt, cutoff_cd4, cutoff_cd14) {
  xmax <- max(lrt$CD4_fkSD, na.rm = TRUE)
  ymax <- max(lrt$CD14_fkSD, na.rm = TRUE)

  label_df <- lrt %>%
    count(pause_change) %>%
    mutate(
      label = case_when(
        pause_change == "Sharp_to_Sharp" ~ paste0("Stable Sharp\nn = ", n),
        pause_change == "Broad_to_Broad" ~ paste0("Stable Broad\nn = ", n),
        pause_change == "Sharp_to_Broad" ~ paste0("Sharp → Broad\nn = ", n),
        pause_change == "Broad_to_Sharp" ~ paste0("Broad → Sharp\nn = ", n),
        TRUE ~ paste0(pause_change, "\nn = ", n)
      ),
      x = case_when(
        pause_change %in% c("Sharp_to_Sharp", "Sharp_to_Broad") ~ cutoff_cd4 * 0.45,
        pause_change %in% c("Broad_to_Sharp", "Broad_to_Broad") ~ cutoff_cd4 + (xmax - cutoff_cd4) * 0.55
      ),
      y = case_when(
        pause_change %in% c("Sharp_to_Sharp", "Broad_to_Sharp") ~ cutoff_cd14 * 0.45,
        pause_change %in% c("Sharp_to_Broad", "Broad_to_Broad") ~ cutoff_cd14 + (ymax - cutoff_cd14) * 0.55
      )
    )

  ggplot(lrt, aes(x = CD4_fkSD, y = CD14_fkSD)) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom  = "polygon",
      bins  = 8,
      alpha = 0.8
    ) +
    scale_fill_viridis_c(name = "Density") +
    geom_vline(xintercept = cutoff_cd4,  linetype = "dashed", color = "gray30", linewidth = 0.5) +
    geom_hline(yintercept = cutoff_cd14, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    geom_label(
      data        = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size        = 4,
      label.size  = 0.2,
      label.color = "gray70",
      fill        = "white",
      alpha       = 0.85
    ) +
    labs(
      x = expression(sigma[CD4]),
      y = expression(sigma[CD14])
    ) +
    theme_cowplot() +
    theme(
      axis.title      = element_text(size = 18),
      axis.text       = element_text(size = 14),
      axis.line       = element_line(size = 0.3),
      axis.ticks      = element_line(size = 0.3),
      legend.position = "right"
    )
}

# Cutoff = midpoint between max(Sharp) and min(Broad) for a given fkSD column
get_sdgroup_cutoff <- function(df, value_col, group_col) {
  sharp_max <- max(df[[value_col]][df[[group_col]] == "Sharp"], na.rm = TRUE)
  broad_min <- min(df[[value_col]][df[[group_col]] == "Broad"], na.rm = TRUE)
  (sharp_max + broad_min) / 2
}

# ==============================================================
# SPECIES CONFIGURATION
# ==============================================================

species_config <- list(

  human = list(
    result_subdir    = 'human',
    beta_scatter     = list(label_x = -6,   label_y = -1,
                            x_lim   = c(-6,   -1), y_lim = c(-6,   -1)),
    delta_beta_sigma = list(filter_lfc = TRUE)    # filter abs(lfc) < 5 before rho + plot
  ),

  rhesus = list(
    result_subdir    = 'rhesus',
    beta_scatter     = list(label_x = -5.5, label_y = -1,
                            x_lim   = c(-5.5, -1), y_lim = c(-5.5, -1)),
    delta_beta_sigma = list(filter_lfc = FALSE)
  ),

  baboon = list(
    result_subdir    = 'baboon',
    beta_scatter     = list(label_x = -6,   label_y = 0,
                            x_lim   = c(-6,   0),  y_lim = c(-6,   0)),
    delta_beta_sigma = list(filter_lfc = FALSE)
  )

)

# ==============================================================
# MAIN PLOTTING LOOP
# ==============================================================

for (sp in names(species_config)) {
  cfg     <- species_config[[sp]]
  out_dir <- paste0(result_dir, cfg$result_subdir, '/')
  message("Plotting: ", sp)

  lrt_tbl <- read.table(paste0(out_dir, 'lrt_result.txt'),
                        sep = '\t', header = TRUE, stringsAsFactors = FALSE)

  # Restore factor levels consumed by color mappings
  lrt_tbl$betaCategory <- factor(lrt_tbl$betaCategory,
                                 levels = c("Down", "Others", "Up"))

  # --- Plot 0: beta_quantile heatmap ---
  p0 <-plot_beta_quantile_heatmap(lrt_tbl, 'CD4_betaGroup', 'CD14_betaGroup')
  ggsave(paste0(out_dir, 'beta_quantile.pdf'), p0, width=3,height=3)

  # --- Plot 1: lrt_beta_scatter ---
  bs <- cfg$beta_scatter
  p1 <- plot_beta_scatter(lrt_tbl,
                          label_x = bs$label_x,
                          label_y = bs$label_y,
                          x_lim   = bs$x_lim,
                          y_lim   = bs$y_lim)
  ggsave(paste0(out_dir, 'lrt_beta_scatter_test.pdf'), p1, width = 4, height = 3)

  # --- Plot 2: deltaBeta_deltaSigma ---
  p2 <- plot_delta_beta_sigma(lrt_tbl,
                              filter_lfc = cfg$delta_beta_sigma$filter_lfc)
  ggsave(paste0(out_dir, 'deltaBeta_deltaSigma_test.pdf'), p2, width = 4.5, height = 3)

  # --- Plot 3: fksd_density ---
  fksd_tbl    <- lrt_tbl[!is.na(lrt_tbl$CD4_sdGroup) & !is.na(lrt_tbl$CD14_sdGroup), ]
  cutoff_cd4  <- get_sdgroup_cutoff(fksd_tbl, "CD4_fkSD",  "CD4_sdGroup")
  cutoff_cd14 <- get_sdgroup_cutoff(fksd_tbl, "CD14_fkSD", "CD14_sdGroup")
  message("[", sp, "] CD4 cutoff = ",  cutoff_cd4)
  message("[", sp, "] CD14 cutoff = ", cutoff_cd14)
  p3 <- plot_fksd_density(fksd_tbl,
                          cutoff_cd4  = cutoff_cd4,
                          cutoff_cd14 = cutoff_cd14)
  ggsave(paste0(out_dir, 'fksd_density_test.pdf'), p3, width = 5, height = 4)
}
