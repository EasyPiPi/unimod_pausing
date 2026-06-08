# ==============================================================
# LIBRARIES
# ==============================================================

library(tidyverse)
library(cowplot)
library(patchwork)
library(RColorBrewer)

# ==============================================================
# GLOBAL PATHS AND CONFIGURATION
# ==============================================================

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))

# Shared directory for both input RDS files and output PDFs
ns_dir <- file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning")

# Color palette for sigma (sd) group
sd_color <- c("#1b9e77", "#762a83")

# Reusable theme strips for multi-panel layouts
remove_x      <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
remove_y      <- theme(axis.title.y = element_blank())
remove_legend <- theme(legend.position = "none")

# ==============================================================
# INPUT FILE TABLE
# ==============================================================

# ns_matrix RDS files produced by 2_2_1_nucleosome_position.R
ns_files <- tibble(
    species  = c("Human",       "Human",           "Rhesus",      "Rhesus",          "Baboon",      "Baboon"),
    celltype = c("CD4+ T cell", "CD14+ monocyte",  "CD4+ T cell", "CD14+ monocyte",  "CD4+ T cell", "CD14+ monocyte"),
    file = file.path(ns_dir, c(
        "human/cd4_ns_matrix.RDS",
        "human/cd14_ns_matrix.RDS",
        "rhesus/cd4_ns_matrix.RDS",
        "rhesus/cd14_ns_matrix.RDS",
        "baboon/cd4_ns_matrix.RDS",
        "baboon/cd14_ns_matrix.RDS"
    ))
)

# ==============================================================
# PLOTTING FUNCTIONS
# ==============================================================

# Returns a line-profile plot faceted by motif class, coloring lines by group.
plot_beta_group_line_faceted_with_promoter <- function(count_matrix,
                                                       rate_df,
                                                       group_col     = "beta_group",
                                                       class_col     = "motif_class",
                                                       range         = NULL,
                                                       colors        = NULL,
                                                       x_label       = c(0, 200),
                                                       shift         = 1001,
                                                       legend_title  = expression(beta ~ group),
                                                       y_label       = "Mean Micro-C signal") {
    stopifnot(nrow(count_matrix) == nrow(rate_df))

    df        <- as.data.frame(count_matrix)
    df$group  <- rate_df[[group_col]]
    df$class  <- factor(rate_df[[class_col]])

    group_means_df <- df %>%
        group_by(class, group) %>%
        summarise(across(starts_with("V"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

    df_long <- group_means_df %>%
        pivot_longer(cols = starts_with("V"), names_to = "position_str", values_to = "value") %>%
        mutate(position = as.numeric(str_remove(position_str, "V")) - shift)

    if (!is.null(range)) {
        df_long <- df_long %>% filter(position >= range[1], position <= range[2])
    }
    df_long <- df_long %>% filter(!class %in% c("Unknown", "Other") & !is.na(class))

    if (is.null(colors)) {
        colors <- RColorBrewer::brewer.pal(6, "Blues")[2:6]
    }

    ggplot(df_long, aes(x = position, y = value, color = group)) +
        geom_line(linewidth = 0.5) +
        facet_wrap(~ class, ncol = 6) +
        labs(x = "Distance from TSS (bp)", y = y_label, color = legend_title) +
        scale_color_manual(values = colors) +
        scale_x_continuous(breaks = x_label) +
        scale_y_continuous(n.breaks = 4) +
        theme_cowplot() +
        theme(
            panel.grid       = element_blank(),
            strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
            legend.position  = "right",
            axis.ticks       = element_line(color = "black", size = 0.3),
            axis.line        = element_line(color = "black", size = 0.3),
            legend.box       = "horizontal",
            axis.text        = element_text(size = 11),
            strip.text       = element_text(size = 9)
        )
}

# Wrapper applying standard range/x-label defaults used throughout this script.
make_microc_profile_plot <- function(ns_matrix, group_col, class_col, ...) {
    plot_beta_group_line_faceted_with_promoter(
        ns_matrix$smooth,
        ns_matrix$rate_df,
        range     = c(0, 1000),
        x_label   = c(0, 500),
        group_col = group_col,
        class_col = class_col,
        ...
    )
}

# Adds a centered title to an existing ggplot.
add_plot_title <- function(p, title, hjust = 0.5, face = "plain") {
    p + ggtitle(title) +
        theme(plot.title = element_text(hjust = hjust, face = face, size = 12))
}

# ==============================================================
# MAIN FIGURE
# ==============================================================

# Load Human CD4 and relabel expression groups
human_cd4 <- readRDS(ns_files$file[1])
human_cd4$rate_df$chiGroup <- factor(
    human_cd4$rate_df$chiGroup,
    levels = levels(human_cd4$rate_df$chiGroup),
    labels = c("Low expression", "Medium expression", "High expression")
)

# Beta groups stratified by expression level (chi group)
beta_by_chi <- make_microc_profile_plot(human_cd4, group_col = "betaGroup", class_col = "chiGroup")
ggsave(file.path(ns_dir, "human/cd4_Micro-C_beta_gb_chi.pdf"), beta_by_chi, width = 8, height = 3)

# Sigma (sd) groups stratified by expression level
sigma_by_chi <- make_microc_profile_plot(
    human_cd4,
    group_col    = "sdGroup",
    class_col    = "chiGroup",
    colors       = sd_color,
    legend_title = expression(sigma ~ group)
)
ggsave(file.path(ns_dir, "human/cd4_Micro-C_fkSD_gb_chi.pdf"), sigma_by_chi, width = 8, height = 3)

# ==============================================================
# SUPPLEMENTARY FIGURE DATA PREPARATION
# ==============================================================

# Builds beta_by_chi, sigma_by_chi, and sigma_by_beta plots for one sample.
make_ns_plots <- function(file, species, celltype) {
    ns_matrix <- readRDS(file)

    ns_matrix$rate_df$chiGroup <- factor(
        ns_matrix$rate_df$chiGroup,
        levels = levels(ns_matrix$rate_df$chiGroup),
        labels = c("Low expression", "Medium expression", "High expression")
    )

    label <- paste(species, celltype)

    list(
        beta_by_chi = make_microc_profile_plot(ns_matrix, "betaGroup", "chiGroup") %>%
            add_plot_title(label),

        sigma_by_chi = make_microc_profile_plot(
            ns_matrix, "sdGroup", "chiGroup",
            colors = sd_color, legend_title = expression(sigma ~ group)
        ) %>% add_plot_title(label),

        sigma_by_beta = make_microc_profile_plot(
            ns_matrix, "sdGroup", "betaGroup",
            colors = sd_color, legend_title = expression(sigma ~ group)
        ) %>% add_plot_title(label)
    )
}

# Generate plots for all 6 species × cell type combinations
plot_list <- pmap(ns_files, ~ make_ns_plots(file = ..3, species = ..1, celltype = ..2))

# ==============================================================
# SUPPLEMENTARY FIGURE LAYOUT HELPERS
# ==============================================================

# 5-panel layout (2+2+1): omits Human CD4 (already in main figure, index 1).
# Axis labels appear only on the bottom-most panel of each column.
make_supp_microC_plot <- function(plot_list, plot_name) {
    p1 <- plot_list[[2]][[plot_name]] + remove_x + remove_y + remove_legend
    p2 <- plot_list[[3]][[plot_name]] + remove_x + remove_y + remove_legend
    p3 <- plot_list[[4]][[plot_name]] + remove_x + remove_legend
    p4 <- plot_list[[5]][[plot_name]] + remove_y  + remove_legend
    p5 <- plot_list[[6]][[plot_name]] + remove_y

    (p1 | p2) / (p3 | p4) / (p5 | patchwork::plot_spacer())
}

# 6-panel layout (2+2+2): all species × cell type combinations.
make_supp_microC_plot_6 <- function(plot_list, plot_name, plot_idx = 1:6) {
    p <- map(plot_idx, ~ plot_list[[.x]][[plot_name]])

    p[[1]] <- p[[1]] + remove_x + remove_y + remove_legend
    p[[2]] <- p[[2]] + remove_x + remove_y + remove_legend
    p[[3]] <- p[[3]] + remove_x + remove_legend
    p[[4]] <- p[[4]] + remove_x + remove_y + remove_legend
    p[[5]] <- p[[5]] + remove_y + remove_legend
    p[[6]] <- p[[6]] + remove_y + remove_legend

    (p[[1]] | p[[2]]) / (p[[3]] | p[[4]]) / (p[[5]] | p[[6]])
}

# ==============================================================
# SAVE SUPPLEMENTARY FIGURES
# ==============================================================

# Supp: beta groups by expression level, all samples except Human CD4
supp_beta_chi_microC <- make_supp_microC_plot(plot_list, plot_name = "beta_by_chi")
ggsave(file.path(ns_dir, "beta_chi_Micro_all.pdf"), supp_beta_chi_microC, width = 15, height = 8)

# Supp: sigma groups by expression level, all samples except Human CD4
supp_sigma_chi_microC <- make_supp_microC_plot(plot_list, plot_name = "sigma_by_chi")
ggsave(file.path(ns_dir, "sigma_chi_Micro_all.pdf"), supp_sigma_chi_microC, width = 15, height = 6)

# Supp: sigma groups by beta group, all 6 samples
supp_sigma_beta_microC <- make_supp_microC_plot_6(plot_list, plot_name = "sigma_by_beta", plot_idx = 1:6)
ggsave(file.path(ns_dir, "sigma_beta_Micro_all.pdf"), supp_sigma_beta_microC, width = 15, height = 8)
