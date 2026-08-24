# ============================================================
# 2_3_LRT_stat.R
#
# Visualization of between-species LRT results.
# Reads the standardized _lrt.txt files produced by
# 2_LRT_between_species.R and generates:
#
#   Main-figure plots (human_rhesus_cd4):
#     beta_scatter.pdf
#     deltaBeta_deltaSigma_density.pdf
#     fksd_density.pdf
#
#   Supplementary multi-panel figures (5 other comparisons):
#     beta_lrt_supp.pdf
#     all_sigma_lrt.pdf
#     all_sigma_beta.pdf
#
# Dependencies: R packages only (ggpubr, tidyverse, cowplot, patchwork).
#   plot_beta_scatter() is loaded from helper_function_plot.R.
# ============================================================

library(ggpubr)
library(tidyverse)
library(cowplot)
library(patchwork)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

# ── Path configuration ────────────────────────────────────────────────────────
lrt_root   <- file.path(.paths$outputs, 'publish/betweenSpecies/LRT/')
result_dir <- lrt_root

plot_beta_scatter <- function(beta_tbl, x_label, y_label,
                              label_x = -6, label_y = -1,
                              x_lim = c(-6, -1), y_lim = c(-6, -1)) {
  cor_df <- beta_tbl %>%
    dplyr::group_by(tssType) %>%
    dplyr::summarise(
      rho = cor(log10(beta1), log10(beta2), method = "spearman", use = "complete.obs"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(label = paste0("rho = ", round(rho, 3)))

  ggplot(beta_tbl, aes(x = log10(beta1), y = log10(beta2))) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "gray60", linewidth = 0.8) +
    geom_point(aes(color = betaCategory), alpha = 0.3, size = 0.5) +
    scale_color_manual(values = c("Others" = "gray", "Up" = "#E41A1C",
                                  "Down" = "#377EB8")) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(x = x_label, y = y_label, color = expression(beta * " change")) +
    facet_wrap(~ tssType) +
    geom_text(data = cor_df, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.2, size = 3, inherit.aes = FALSE) +
    cowplot::theme_cowplot() +
    theme(
      axis.title = element_text(size = 14), axis.text = element_text(size = 12),
      axis.line  = element_line(size = 0.3), axis.ticks = element_line(size = 0.3),
      strip.background = element_rect(fill = "white", colour = "black",
                                      linewidth = 0.4),
      strip.text = element_text(size = 12)
    )
}


# ── Supplementary comparison configurations ───────────────────────────────────
#
# human_rhesus_cd4 is excluded here because its plots go in the main figure.
# rhesus_baboon comparisons already have baboon as species1 in the LRT output —
# no direction flip is needed.

comparisons <- list(

  list(
    name        = "human_baboon_cd4",
    prefix1     = "baboonCD4",
    prefix2     = "humanCD4",
    label1      = "Baboon",
    label2      = "Human",
    show_legend = FALSE
  ),

  list(
    name        = "rhesus_baboon_cd4",
    prefix1     = "baboonCD4",
    prefix2     = "rhesusCD4",
    label1      = "Baboon",
    label2      = "Rhesus",
    show_legend = TRUE
  ),

  list(
    name        = "human_rhesus_cd14",
    prefix1     = "rhesusCD14",
    prefix2     = "humanCD14",
    label1      = "Rhesus",
    label2      = "Human",
    show_legend = FALSE
  ),

  list(
    name        = "human_baboon_cd14",
    prefix1     = "baboonCD14",
    prefix2     = "humanCD14",
    label1      = "Baboon",
    label2      = "Human",
    show_legend = FALSE
  ),

  list(
    name        = "rhesus_baboon_cd14",
    prefix1     = "baboonCD14",
    prefix2     = "rhesusCD14",
    label1      = "Baboon",
    label2      = "Rhesus",
    show_legend = FALSE
  )

)


# ── Helper functions ──────────────────────────────────────────────────────────

# Validate that required columns are present in a loaded LRT table.
# Pass extra_cols to check comparison-specific columns (e.g. sigma columns).
check_required_cols <- function(df, name, extra_cols = character(0)) {
  required <- c("tssType", "lfc", "deltaSD", "pause_change",
                "betaCategory", "beta1", "beta2",
                extra_cols)
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "LRT table '%s' is missing required columns: %s\n  Re-run 2_LRT_between_species.R to regenerate.",
      name, paste(missing, collapse = ", ")
    ))
  }
}

# Compute the midpoint between the largest Sharp fkSD and the smallest Broad fkSD.
# Used to draw the sdGroup boundary line in sigma density plots.
get_sigma_cutoff <- function(df, fkSD_col, sdGroup_col) {
  max_sharp <- max(df[[fkSD_col]][df[[sdGroup_col]] == "Sharp"], na.rm = TRUE)
  min_broad <- min(df[[fkSD_col]][df[[sdGroup_col]] == "Broad"], na.rm = TRUE)
  (max_sharp + min_broad) / 2
}

# Build a beta or sigma axis label expression from a species name string.
# e.g. make_label("beta", "Baboon") -> beta[Baboon]
make_label <- function(greek, species_name) {
  bquote(.(as.name(greek))[.(as.name(species_name))])
}


plot_fkSD_density <- function(lrt,
                              species1_col,
                              species2_col,
                              cutoff_species1,
                              cutoff_species2,
                              species1_label = "Species 1",
                              species2_label = "Species 2") {

  xmax <- max(lrt[[species1_col]], na.rm = TRUE)
  ymax <- max(lrt[[species2_col]], na.rm = TRUE)

  label_df <- lrt %>%
    count(tssType, pause_change) %>%
    mutate(
      label = case_when(
        pause_change == "Sharp_to_Sharp" ~ paste0("Stable Sharp\nn = ", n),
        pause_change == "Broad_to_Broad" ~ paste0("Stable Broad\nn = ", n),
        pause_change == "Sharp_to_Broad" ~ paste0("Sharp → Broad\nn = ", n),
        pause_change == "Broad_to_Sharp" ~ paste0("Broad → Sharp\nn = ", n),
        TRUE ~ paste0(pause_change, "\nn = ", n)
      ),
      x = case_when(
        pause_change %in% c("Sharp_to_Sharp", "Sharp_to_Broad") ~
          cutoff_species1 * 0.45,
        pause_change %in% c("Broad_to_Sharp", "Broad_to_Broad") ~
          cutoff_species1 + (xmax - cutoff_species1) * 0.55
      ),
      y = case_when(
        pause_change %in% c("Sharp_to_Sharp", "Broad_to_Sharp") ~
          cutoff_species2 * 0.45,
        pause_change %in% c("Sharp_to_Broad", "Broad_to_Broad") ~
          cutoff_species2 + (ymax - cutoff_species2) * 0.55
      )
    )

  ggplot(lrt, aes(x = .data[[species1_col]], y = .data[[species2_col]])) +
    facet_wrap(~ tssType) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      bins = 8,
      alpha = 0.8
    ) +
    scale_fill_viridis_c(name = "Density") +
    geom_vline(
      xintercept = cutoff_species1,
      linetype = "dashed", color = "gray30", linewidth = 0.5
    ) +
    geom_hline(
      yintercept = cutoff_species2,
      linetype = "dashed", color = "gray30", linewidth = 0.5
    ) +
    geom_label(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 2,
      label.size = 0.2,
      label.color = "gray70",
      fill = "white",
      alpha = 0.85
    ) +
    labs(
      x = bquote(sigma[.(species1_label)]),
      y = bquote(sigma[.(species2_label)])
    ) +
    cowplot::theme_cowplot() +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10),
      strip.text = element_text(size = 10),
      axis.line  = element_line(linewidth = 0.3),
      strip.background = element_rect(fill = "white", colour = "black",
                                      linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    )
}


plot_beta_sigma_density <- function(df,
                                    xlim = c(-3, 3),
                                    ylim = c(-0.6, 0.6),
                                    bins = 6,
                                    rho_digits = 3,
                                    point_size = 4) {

  rho_df <- df %>%
    group_by(tssType) %>%
    summarise(
      rho = cor(lfc, deltaSD, method = "spearman", use = "complete.obs"),
      .groups = "drop"
    )

  df %>%
    ggplot(aes(x = lfc, y = deltaSD)) +
    facet_wrap(~ tssType) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      bins = bins,
      alpha = 0.8
    ) +
    scale_fill_viridis_c(option = "C", name = "Density") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_text(
      data = rho_df,
      aes(label = paste0("rho = ", round(rho, rho_digits))),
      x = Inf, y = -Inf,
      hjust = 1.1, vjust = -0.5,
      size = point_size,
      inherit.aes = FALSE
    ) +
    cowplot::theme_cowplot() +
    labs(
      x = expression(Delta * beta),
      y = expression(Delta * sigma),
      fill = "Density"
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 12),
      strip.text = element_text(size = 12),
      axis.line  = element_line(linewidth = 0.3),
      strip.background = element_rect(fill = "white", colour = "black",
                                      linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10)
    )
}


# ── Load and validate supplementary LRT tables ───────────────────────────────

lrt_data <- lapply(comparisons, function(cfg) {
  path <- paste0(lrt_root, cfg$name, '_lrt.txt')
  df   <- read.table(path, sep = '\t', header = TRUE)
  check_required_cols(df, cfg$name)
  df
})
names(lrt_data) <- sapply(comparisons, `[[`, "name")


# ── Supplementary Figure 1: beta_lrt_supp.pdf ────────────────────────────────
# beta1 vs beta2 scatter plots, faceted by tssType.

plots_beta <- Map(function(cfg, df) {
  p <- plot_beta_scatter(
    df,
    x_label = make_label("beta", cfg$label1),
    y_label  = make_label("beta", cfg$label2)
  )
  if (!cfg$show_legend) p <- p + theme(legend.position = "none")
  p
}, comparisons, lrt_data)

all_beta_lrt <- (plots_beta[[1]] | plots_beta[[2]] | patchwork::plot_spacer()) /
                (plots_beta[[3]] | plots_beta[[4]] | plots_beta[[5]])

ggsave(paste0(result_dir, 'beta_lrt_supp.pdf'), all_beta_lrt,
       width = 13, height = 6)


# ── Supplementary Figure 2: all_sigma_lrt.pdf ────────────────────────────────
# fkSD (sigma) density: species1 vs species2, faceted by tssType.
# Sigma cutoffs are derived automatically from sdGroup assignments.

plots_sigma <- Map(function(cfg, df) {
  fkSD1_col    <- paste0(cfg$prefix1, "_fkSD")
  fkSD2_col    <- paste0(cfg$prefix2, "_fkSD")
  sdGroup1_col <- paste0(cfg$prefix1, "_sdGroup")
  sdGroup2_col <- paste0(cfg$prefix2, "_sdGroup")

  cutoff1 <- get_sigma_cutoff(df, fkSD1_col, sdGroup1_col)
  cutoff2 <- get_sigma_cutoff(df, fkSD2_col, sdGroup2_col)

  p <- plot_fkSD_density(
    lrt             = df,
    species1_col    = fkSD1_col,
    species2_col    = fkSD2_col,
    cutoff_species1 = cutoff1,
    cutoff_species2 = cutoff2,
    species1_label  = cfg$label1,
    species2_label  = cfg$label2
  )
  if (!cfg$show_legend) p <- p + theme(legend.position = "none")
  p
}, comparisons, lrt_data)

all_sigma_lrt <- (plots_sigma[[1]] | plots_sigma[[2]] | patchwork::plot_spacer()) /
                 (plots_sigma[[3]] | plots_sigma[[4]] | plots_sigma[[5]])

ggsave(paste0(result_dir, 'all_sigma_lrt.pdf'), all_sigma_lrt,
       width = 13, height = 6)


# ── Supplementary Figure 3: all_sigma_beta.pdf ───────────────────────────────
# delta-beta vs delta-sigma density, faceted by tssType.

plots_beta_sigma <- Map(function(cfg, df) {
  df_plot <- df[!is.na(df$lfc) & !is.na(df$deltaSD), ]
  p <- plot_beta_sigma_density(df_plot)
  if (!cfg$show_legend) p <- p + theme(legend.position = "none")
  p
}, comparisons, lrt_data)

all_sigma_beta <- (plots_beta_sigma[[1]] | plots_beta_sigma[[2]] | patchwork::plot_spacer()) /
                  (plots_beta_sigma[[3]] | plots_beta_sigma[[4]] | plots_beta_sigma[[5]])

ggsave(paste0(result_dir, 'all_sigma_beta.pdf'), all_sigma_beta,
       width = 13, height = 7)


# ── Main-figure plots: human_rhesus_cd4 ──────────────────────────────────────
# These three plots are for the main figure (not the supplementary panels).
# Output is written to the human_rhesus_cd4 comparison sub-directory.

hr_cd4_dir <- paste0(lrt_root, 'human_rhesus_cd4_')

hr_cd4 <- read.table(paste0(hr_cd4_dir, 'lrt.txt'), sep = '\t', header = TRUE)

check_required_cols(
  hr_cd4, "human_rhesus_cd4",
  extra_cols = c("rhesusCD4_fkSD", "humanCD4_fkSD",
                 "rhesusCD4_sdGroup", "humanCD4_sdGroup")
)

# Sigma cutoffs derived from sdGroup assignments (no hard-coded values).
cutoff_rhesus <- get_sigma_cutoff(hr_cd4, "rhesusCD4_fkSD", "rhesusCD4_sdGroup")
cutoff_human  <- get_sigma_cutoff(hr_cd4, "humanCD4_fkSD",  "humanCD4_sdGroup")

# 1. Beta scatter plot.
p_beta <- plot_beta_scatter(
  hr_cd4,
  x_label = expression(beta[Rhesus]),
  y_label  = expression(beta[Human])
)
ggsave(paste0(hr_cd4_dir, 'beta_scatter.pdf'), p_beta, width = 6.5, height = 3)

# 2. Delta-beta vs delta-sigma density plot.
hr_cd4_filtered <- hr_cd4[!is.na(hr_cd4$lfc) & !is.na(hr_cd4$deltaSD), ]

p_bs <- plot_beta_sigma_density(hr_cd4_filtered)
ggsave(paste0(hr_cd4_dir, 'deltaBeta_deltaSigma_density.pdf'), p_bs,
       width = 6, height = 3)

# 3. Sigma density plot (rhesus fkSD vs human fkSD).
p_fksd <- plot_fkSD_density(
  lrt             = hr_cd4,
  species1_col    = "rhesusCD4_fkSD",
  species2_col    = "humanCD4_fkSD",
  cutoff_species1 = cutoff_rhesus,
  cutoff_species2 = cutoff_human,
  species1_label  = "Rhesus",
  species2_label  = "Human"
)
ggsave(paste0(hr_cd4_dir, 'fksd_density.pdf'), p_fksd, width = 8, height = 3.5)
