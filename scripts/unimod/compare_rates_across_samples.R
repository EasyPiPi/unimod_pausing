#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(corrr)
library(ComplexHeatmap)
# library(tidymodels)
library(tidyverse)
library(GGally)
library(ggpointdensity)
library(recipes)
library(ggpubr)
library(cowplot)

#### snakemake files ####
# root_dir <- "~/Desktop/github/unimod_human"
root_dir <- "~/projects/Snakemake_projects/unimod_human"
rate_dirs <- list.files(file.path(root_dir, "outputs/within_sample"))
result_dir <- file.path(root_dir, "outputs/between_samples/figures")

#### end of parsing arguments ####
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

rate_tbls <-
  map(rate_dirs,
      ~ readRDS(file.path("outputs/within_sample", .x, "pause_release/rate.RDS")))
names(rate_tbls) <- rate_dirs

rate_tbl <- bind_rows(rate_tbls, .id = "sample_id")

# subset some samples
rate_tbl_control <- rate_tbl %>% filter(sample_id %in% str_subset(rate_dirs, "control"))

# use rplot in correlate, no scatter
# rate_tbl %>%
#   select(sample_id, gene_id, chi) %>%
#   pivot_wider(names_from = sample_id, values_from = chi) %>%
#   select(-gene_id) %>%
#   correlate(method = "spearman", use = "pairwise.complete.obs", quiet = TRUE, diagonal = 1) %>%
#   rearrange() %>%  # rearrange by correlations
#   shave() %>%
#   rplot(print_cor = TRUE) +
#   theme(axis.text.x = element_text(angle = 15, hjust = 1))

# function for lower half of the correlation plot
# option 1: ggpointdensity
# lowerfun <- function(data, mapping){
#   ggplot(data = data, mapping = mapping)+
#     ggpointdensity::geom_pointdensity() +
#     scale_color_viridis_c(alpha = 1, direction = 1, option = "D") +
#     geom_abline(intercept = 0, slope = 1, color="grey50",
#                 linetype="dashed", size = 1, alpha = 0.5)}

# option 2: geom_point with alpha
lowerfun <- function(data, mapping, axis_lim){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.1, size = 0.1) +
    geom_abline(intercept = 0, slope = 1, color="grey50",
                linetype="dashed", size = 1, alpha = 0.5) +
    coord_cartesian(xlim = axis_lim, ylim = axis_lim)}

scatter_matrix <- function(rate_tbl, col, axis_lim = c(-16, 2), size = 9) {
  # clean up sample names for better visualization in correlation matrix
  rate_tbl$sample_id <-
    sapply(str_split(rate_tbl$sample_id, pattern = "-"),
           function(x) paste0(x[3:4], collapse = "-"))
  rate_tbl %>%
    dplyr::select(sample_id, gene_id, {{col}}) %>%
    pivot_wider(names_from = sample_id, values_from = {{col}}) %>%
    dplyr::select(-gene_id) %>%
    log2() %>%
    ggpairs(upper = list(continuous = wrap("cor", method = "spearman",
                                           digits = 2, size = size)),
            lower = list(continuous = wrap(lowerfun, axis_lim = axis_lim))) +
    cowplot::theme_cowplot()
}

p <- scatter_matrix(rate_tbl = rate_tbl, col = "chi", size = 5, axis_lim = c(-12, 2))
ggsave(filename = file.path(result_dir, "chi_correlations_between_samples.png"),
       plot = p, width = 12, height = 10)

p <- scatter_matrix(rate_tbl = rate_tbl %>% filter(str_detect(sample_id, "control")), col = "chi")
ggsave(filename = file.path(result_dir, "chi_correlations_between_control_samples.png"),
       plot = p, width = 12, height = 10)

p <- scatter_matrix(rate_tbl = rate_tbl, col = "beta", axis_lim = c(-20, -3), size = 5)
ggsave(filename = file.path(result_dir, "beta_correlations_between_samples.png"),
       plot = p, width = 12, height = 10)

p <- scatter_matrix(rate_tbl = rate_tbl %>% filter(str_detect(sample_id, "control")),
                    col = "beta", axis_lim = c(-20, -3))
ggsave(filename = file.path(result_dir, "beta_correlations_between_control_samples.png"),
       plot = p, width = 12, height = 10)

# get heat shock and Celastrol sample
rate_tbl_subset <-
  rate_tbl %>%
  filter(sample_id %in% c("PROseq-K562-dukler-control-SE",
                          "PROseq-K562-vihervaara-control-SE"))

# scatterplot with density
pointdensity <- function(col_name, title, axis_lim = NULL) {
  p <- rate_tbl_subset %>%
    dplyr::select(sample_id, gene_id, {{col_name}}) %>%
    mutate({{col_name}} := log2({{col_name}})) %>%
    pivot_wider(names_from = sample_id, values_from = {{col_name}}) %>%
    ggplot(mapping = aes(x = `PROseq-K562-dukler-control-SE`, y = `PROseq-K562-vihervaara-control-SE`)) +
    geom_pointdensity() +
    viridis::scale_color_viridis(alpha = 1, direction = 1, option = "D") +
    # geom_abline(intercept = 0, slope = 1, color="grey50",
    #             linetype="dashed", size = 1, alpha = 0.5) +
    labs(title = title,
         x = expression(italic("Dukler et al.  ")*log[2]*"(Rates)"),
         y = expression(italic("Vihervaara et al.  ")*log[2]*"(Rates)")) +
    cowplot::theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5, size = 15))
  if (is.null(axis_lim)) {
    return(p + coord_equal() + ggpubr::stat_cor())
  } else {
    return(p +
             coord_cartesian(xlim = axis_lim, ylim = axis_lim) +
             ggpubr::stat_cor(label.x = axis_lim[1] + 1, label.y = axis_lim[2] - 1))
  }
}

p <- pointdensity(chi, title = "Chi estimates")
ggsave(filename = file.path(result_dir, "pointdensity_initiation_dukler_vs_vihervaara.png"),
       plot = p, width = 6, height = 5)

p <- pointdensity(beta, title = "Pause Release Rates", axis_lim = c(-20, -3))
ggsave(filename = file.path(result_dir, "pointdensity_pause_dukler_vs_vihervaara.png"),
       plot = p, width = 6, height = 5)

# plots for rate estimates
process_df <- function(sample_name, sel_col, log_convert = TRUE) {
  df <- rate_tbl %>% 
    filter(str_detect(sample_id, sample_name)) %>%
    select(sample_id, gene_id, {{sel_col}}) %>%
    pivot_wider(id_cols = gene_id, names_from = sample_id, values_from = {{sel_col}}) %>%
    na.omit() %>%
    pivot_longer(cols = contains(sample_name)) 
  if (log_convert) {
    df <- df %>% mutate(value = log2(value))
  }
  return(df)
}

violin_rate <- function(df) {
  p <- df %>%
    ggviolin(x = "name", y = "value", fill = "name",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(label.x = 1.4, label.y.npc = "top") +
    labs(x = "", y = expression(log[2]*beta)) +
    theme(legend.position = "none")
  return(p)
}

boxplot_rate <- function(df) {
  df %>%
    ggboxplot(x = "name", y = "value", fill = "name",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
              outlier.shape = NA) +
    stat_compare_means(label.x = 1.4, label.y = -8.5) +
    labs(x = "", y = expression(log[2]*beta)) +
    coord_cartesian(ylim = c(-17, -7.5)) +
    theme(legend.position = "none")
}

boxplot_site <- function(df, y_axis_lab) {
  df %>%
    ggboxplot(x = "name", y = "value", fill = "name",
              palette = c("#00AFBB", "#E7B800", "#FC4E07")) +
    stat_compare_means(label.x = 1.4, label.y.npc = "top") +
    labs(x = "", y = y_axis_lab) +
    theme(legend.position = "none")
  
}

# Celastrol
sample_name <- "dukler"
df <- process_df(sample_name, beta)
p <- violin_rate(df)
ggsave(file.path(result_dir, paste0(sample_name, "_violin_beta_distribution.png")), plot = p,
       width = 8, height = 6)
p <- boxplot_rate(df) 
ggsave(file.path(result_dir, paste0(sample_name, "_boxplot_beta_distribution.png")), plot = p,
       width = 8, height = 6)

df <- process_df(sample_name, fk_mean, log_convert = FALSE)
p <- boxplot_site(df,  "Position of Pause Sites")
ggsave(file.path(result_dir, paste0(sample_name, "_boxplot_fk_mean.png")), plot = p,
       width = 8, height = 6)

df <- process_df(sample_name, fk_var, log_convert = FALSE)
p <- boxplot_site(df,  "Variance of Pause Sites")
ggsave(file.path(result_dir, paste0(sample_name, "_boxplot_fk_variance.png")), plot = p,
       width = 8, height = 6)

# Heat shock
sample_name <- "vihervaara"
df <- process_df(sample_name, beta)
p <- violin_rate(df)
ggsave(file.path(result_dir, paste0(sample_name, "_violin_beta_distribution.png")), plot = p,
       width = 8, height = 6)
p <- boxplot_rate(df) 
ggsave(file.path(result_dir, paste0(sample_name, "_boxplot_beta_distribution.png")), plot = p,
       width = 8, height = 6)

df <- process_df(sample_name, fk_mean, log_convert = FALSE)
p <- boxplot_site(df,  "Position of Pause Sites")
ggsave(file.path(result_dir, paste0(sample_name, "_boxplot_fk_mean.png")), plot = p,
       width = 8, height = 6)

df <- process_df(sample_name, fk_var, log_convert = FALSE)
p <- boxplot_site(df,  "Variance of Pause Sites")
ggsave(file.path(result_dir, paste0(sample_name, "_boxplot_fk_variance.png")), plot = p,
       width = 8, height = 6)

# fk mean vs. var
theme_set(cowplot::theme_cowplot())
p <- rate_tbl %>%
  filter(str_detect(sample_id, "vihervaara")) %>%
  mutate(fk_std = fk_var ^ 0.5,
         Condition = ifelse(str_detect(sample_id, "control"), "NHS", "HS")) %>%
  ggplot(aes(x = fk_mean, y = fk_std, color = Condition)) +
  geom_density_2d() +
  labs(x = expression(Mean~of~k), y = expression(SD~of~k)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))
  
ggsave(filename = file.path(result_dir, "fk_mean_vs_std_vihervaara.png"), plot = p,
       width = 6, height = 4)
ggsave(filename = file.path(result_dir, "fk_mean_vs_std_vihervaara.pdf"), plot = p,
       width = 6, height = 4)

p <- rate_tbl %>%
  filter(str_detect(sample_id, "dukler")) %>%
  mutate(fk_std = fk_var ^ 0.5,
         Condition = ifelse(str_detect(sample_id, "control"), "Control", "Treated")) %>%
  ggplot(aes(x = fk_mean, y = fk_std, color = Condition)) +
  geom_density_2d() +
  labs(x = expression(Mean~of~k), y = expression(SD~of~k)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

ggsave(filename = file.path(result_dir, "fk_mean_vs_std_dukler.png"), plot = p,
       width = 6, height = 4)
ggsave(filename = file.path(result_dir, "fk_mean_vs_std_dukler.pdf"), plot = p,
       width = 6, height = 4)

#### figures for steric hindrance ####
# # tables for matching omega
# ovp_tbls <- tibble(sample = rate_dirs,
#        tbl = map(rate_dirs,
#     ~ read_csv(file.path("outputs/within_sample", .x, "steric_hindrance/gressel_overlap.csv"),
#                col_types = cols(gene_id = col_character())))) %>%
#   filter(str_detect(sample, "dukler|vihervaara"))
# 
# ovp_tbls <- ovp_tbls %>%
#   mutate(n = map_dbl(tbl, NROW),
#          median_omega = map_dbl(tbl, ~ .x %>% pull(initiation.rate.per.cell) %>% median(na.rm = TRUE)))

# tables with matching a lower initiation rate
st_tbls <- 
  tibble(sample = rate_dirs,
         tbl = map(rate_dirs,
                   ~ read_csv(file.path("outputs/within_sample", .x, "steric_hindrance/rate.csv"),
                              col_types = cols(gene_id = col_character()))))

st_tbls <- st_tbls %>%
  filter(str_detect(sample, "dukler|vihervaara"))

# tables with matching a higher initiation rate
st_mt_tbls <- 
  tibble(sample = rate_dirs,
         tbl = map(rate_dirs,
                   ~ read_csv(file.path("outputs/within_sample", .x,
                                        "steric_hindrance/rate_matched_gressel.csv"),
                              col_types = cols(gene_id = col_character()))))

st_mt_tbls <- st_mt_tbls %>%
  filter(str_detect(sample, "dukler|vihervaara"))

# pull out genes with robust expressions, i.e. high chi
st_tbls <- st_tbls %>%
  mutate(high_chi = map(tbl,
                        ~ .x %>% slice_max(order_by = chi, prop = 0.8) %>%
                          pull(gene_id)))

# get overlapped genes between control and treated samples
get_overlop_gene <- function(st_tbls, dataset) {
  x <- st_tbls %>% filter(str_detect(sample, {{dataset}})) %>% pull(high_chi)
  return(intersect(x[[1]], x[[2]]))
}

dukler_gn <- get_overlop_gene(st_tbls, "dukler")
vihervaara_gn <- get_overlop_gene(st_tbls, "vihervaara")

process_tbl <- function(st_tbls, dataset, gene_id) {
  st_tbls %>%
    filter(str_detect(sample, {{dataset}})) %>%
    select(sample, tbl) %>%
    unnest(tbl) %>%
    filter(gene_id %in% {{gene_id}})
}

dukler_tbl <-
  process_tbl(st_tbls, "dukler", dukler_gn) %>%
  mutate(sample = ifelse(str_detect(sample, "control"), "Control", "Treated"),
         sample = factor(sample, levels = c("Control", "Treated")))
  
vihervaara_tbl <-
  process_tbl(st_tbls, "vihervaara", vihervaara_gn) %>%
  mutate(sample = ifelse(str_detect(sample, "control"), "NHS", "HS"),
         sample = factor(sample, levels = c("NHS", "HS")))

vihervaara_figs <- list()
dukler_figs <- list()

## phi estimates between untreated and treated samples ##
format_digits <- function(x, digits) format(round(x, digits = digits), nsmall = digits)

hist_phi <- function(tbl) {
  ano_sum <- tbl %>%
    group_by(sample) %>%
    summarise(mean_phi = mean(phi),
              median_phi = median(phi))
  
  annotations <- data.frame(
    xpos = c(0.25, 0.25, 0.6, 0.6), ypos = rep(-Inf, 4),
    annotateText = c(paste0("Mean: ", format_digits(ano_sum$mean_phi[[1]], digits = 2)),
                     paste0("Median: ", format_digits(ano_sum$median_phi[[1]], digits = 2)),
                     paste0("Mean: ", format_digits(ano_sum$mean_phi[[2]], digits = 2)),
                     paste0("Median: ", format_digits(ano_sum$median_phi[[2]], digits = 2))),
    hjustvar = rep(0, 4),
    vjustvar = c(-11, -13, -5, -7))
  
  tbl %>%
    ggplot() +
    geom_histogram(aes(x = phi, fill = sample), bins = 50, color="gray",
                   position="identity", alpha = 0.5) + 
    scale_fill_manual(values=c("#377eb8", "#e41a1c")) +
    geom_text(data = annotations,
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar,
                  label = annotateText), size = 5) +
    # scale_y_continuous(labels = scales::percent) +
    labs(x = expression("Estimated "*phi), y = "Frequency", fill = "") +
    theme(legend.position = c(.9, .95),
          legend.justification = c("right", "top"))
}

vihervaara_figs[["p1"]] <- hist_phi(vihervaara_tbl)
dukler_figs[["p1"]] <- hist_phi(dukler_tbl)

## phi estimates by low or high omega scaling ##
st_tbls$phi_frac <-
  map2_dbl(st_tbls$tbl, list(dukler_gn, dukler_gn, vihervaara_gn, vihervaara_gn),
           ~ .x %>%
               filter(gene_id %in% .y) %>%
               summarise(phi_frac = mean(phi > 0.95)) %>% 
               pull(phi_frac))
st_tbls$group <- "L" 

st_mt_tbls$phi_frac <-
  map2_dbl(st_mt_tbls$tbl, list(dukler_gn, dukler_gn, vihervaara_gn, vihervaara_gn),
           ~ .x %>%
             filter(gene_id %in% .y) %>%
             summarise(phi_frac = mean(phi > 0.95)) %>% 
             pull(phi_frac))
st_mt_tbls$group <- "H" 

frac_df <- bind_rows(st_tbls %>% select(sample, phi_frac, group),
                     st_mt_tbls %>% select(sample, phi_frac, group)) %>%
  mutate(group = factor(group, level = c("L", "H"))) %>%
  arrange(sample, group)

vihervaara_frac_tbl <- frac_df %>%
  filter(str_detect(sample, "vihervaara")) %>%
  mutate(sample = ifelse(str_detect(sample, "control"), "NHS", "HS"),
         sample = factor(sample, levels = c("NHS", "HS")))
dukler_frac_tbl <- frac_df %>%
  filter(str_detect(sample, "dukler")) %>%
  mutate(sample = ifelse(str_detect(sample, "control"), "Control", "Treated"),
         sample = factor(sample, levels = c("Control", "Treated")))

phi_bar <- function(tbl) {
  tbl %>%
    ggplot(aes(x = group, y = phi_frac, fill = sample)) +
    geom_col(alpha = 0.5, width = 0.6, position = position_dodge(width=0.8)) +
    scale_fill_manual(values=c("#377eb8", "#e41a1c")) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Percentage", fill = "") +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = c(.4, .95),
          legend.justification = c("right", "top"))
}

vihervaara_figs[["p2"]]  <- phi_bar(vihervaara_frac_tbl)
dukler_figs[["p2"]] <- phi_bar(dukler_frac_tbl)

## initiation rates alpha vs. omega ##
init_hist <- function(tbl) {
  x <- summary(tbl$omega_zeta)
  y <- summary(tbl$alpha_zeta)
  
  annotation_x <- data.frame(
    xpos = c(0.4, 0.4), ypos =  c(0.15, 0.135),
    annotateText = c(paste0(""),
                     paste0("Median: ", format_digits(x[["Median"]], digits = 2))
                     # paste0("Mean: ", round(x[["Mean"]], digits = 2)),
    ),
    hjustvar = c(0, 0),
    vjustvar = c(0, 0))
  
  annotation_y <- data.frame(
    xpos = c(0.8, 0.8), ypos =  c(0.045, 0.03),
    annotateText = c(paste0(""),
                     paste0("Median: ", format_digits(y[["Median"]], digits = 2))
                     # paste0("Mean: ", round(y[["Mean"]], digits = 2)),
    ),
    hjustvar = c(0, 0),
    vjustvar = c(0, 0))
  
  x_lab <- seq(0, 2, 0.5)
  names(x_lab) <- x_lab  
  x_lab[length(x_lab)] <- "≥2.0"
  
  tbl %>%
    select(omega_zeta, alpha_zeta) %>% 
    pivot_longer(cols = c(omega_zeta, alpha_zeta)) %>%
    mutate(value = ifelse(value > 2, 2, value)) %>% 
    ggplot(aes(x = value)) +
    # geom_histogram(aes(y =), bins = bin_n) +
    geom_histogram(aes(y =  (..count..)/sum(..count..), fill = name), position="identity",
                   color = "gray", binwidth = 0.1, alpha=0.3) +
    geom_text(data = annotation_x,
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar,
                  label = annotateText), size = 5) +
    geom_text(data = annotation_y,
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar,
                  label = annotateText), size = 5) +
    geom_text(aes(x = 1.2, y = 0.095, label = "Steric Hindrance"), fontface = 'italic', size = 4) +
    geom_segment(aes(x = 1.5, y = 0.08, xend = 0.5, yend = 0.08),
                 lineend = "round", linejoin = "round", color = "gray30",
                 arrow = arrow(length = unit(0.2, "cm")), size = 1.5) +
    geom_vline(xintercept = 1.95, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 2)) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks=seq(0, 2, 0.5), labels = x_lab) +
    scale_fill_discrete(labels=c(expression(alpha*zeta), expression(omega*zeta))) +
    labs(x = "Initiation Rate (events/min.)", y = "Percentage", fill = "") +
    theme(legend.position = c(.9, .95),
          legend.justification = c("right", "top"))
}

vihervaara_figs[["p3"]] <- init_hist(vihervaara_tbl %>% filter(sample == "NHS"))
vihervaara_figs[["p4"]] <- init_hist(vihervaara_tbl %>% filter(sample == "HS"))

dukler_figs[["p3"]] <- init_hist(dukler_tbl %>% filter(sample == "Control"))
dukler_figs[["p4"]] <- init_hist(dukler_tbl %>% filter(sample == "Treated"))

composite_plot <- function(fig_ls) {
  g1 <- plot_grid(fig_ls[["p1"]], fig_ls[["p2"]],
                  rel_widths = c(2, 1), align = "hv", axis = "l")
  g2 <- plot_grid(fig_ls[["p3"]], fig_ls[["p4"]],
                  rel_widths = c(1, 1), align = "hv", axis = "l")
  plot_grid(g1, g2, nrow = 2, align = "hv", axis = "l")
}

p <- composite_plot(vihervaara_figs)
ggsave(file.path(result_dir, "vihervaara_steric_hindrance.png"),
       width = 10, height = 8, device = png)

p <- composite_plot(dukler_figs)
ggsave(file.path(result_dir, "dukler_steric_hindrance.png"),
       width = 10, height = 8, device = png)

# vihervaara_tbl %>%
#   group_by(sample) %>%
#   summarise(mean_omega_zeta = mean(omega_zeta))

# # PCA analysis
# col_names = c("assay", "cell_line", "reference", "group", "read_type")
# 
# pca_and_visualization <- function(df, col, pca_title) {
#   # perform PCA analysis
#   dr <- df %>%
#     # filter(sample_id != "PROseq-K562-chivu-treated-PE") %>%
#     na.omit() %>%
#     dplyr::select(sample_id, gene_id, {{col}}) %>%
#     pivot_wider(names_from=gene_id, values_from={{col}})
# 
#   dr <- dr %>%
#     # replace na with 0s, only a few of them in each sample
#     dplyr::mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0))) %>%
#     recipes::recipe(~ .) %>%
#     # center the data
#     step_center(all_numeric()) %>%
#     # center the data
#     step_scale(all_numeric()) %>%
#     # pca on all numeric variables
#     step_pca(all_numeric()) %>%
#     prep()
# 
#   # percentage of variance could be explained by PCs
#   # std deviation
#   sdev <- dr$steps[[3]]$res$sdev
#   #sdev <- dr$steps[[1]]$res$sdev
#   percent_variation <- sdev^2 / sum(sdev^2)
# 
#   dr_juice <- dr %>% juice() %>%
#     separate(sample_id, into = col_names)
# 
#   # visualize the first four PCs
#   p1 <- dr_juice %>%
#     ggplot(aes(PC1, PC2)) +
#     geom_point(aes(color = group, shape = reference), alpha = 0.7, size = 2) +
#     scale_color_brewer(palette="Set1") +
#     # geom_text(aes(label = id), hjust = -0.2, vjust = 0, size = 1) +
#     labs(x = paste0("PC1 (", round(percent_variation[[1]] * 100, digits = 2), "% variance explained)"),
#          y = paste0("PC2 (", round(percent_variation[[2]] * 100, digits = 2), "% variance explained)"),
#          title = pca_title) +
#     theme_classic() +
#     theme(plot.title = element_text(hjust = 0.5))
# 
#   # p2 <- dr_juice %>%
#   #   ggplot(aes(PC2, PC3)) +
#   #   geom_point(aes(color = species, shape = cell_type), alpha = 0.7, size = 2) +
#   #   scale_color_brewer(palette="Set1") +
#   #   # geom_text(aes(label = id), hjust = -0.2, vjust = 0, size = 1) +
#   #   labs(x = paste0("PC2 (", round(percent_variation[[2]] * 100, digits = 2), "% variance explained)"),
#   #        y = paste0("PC3 (", round(percent_variation[[3]] * 100, digits = 2), "% variance explained)"),
#   #        title = pca_title) +
#   #   theme_classic() +
#   #   theme(plot.title = element_text(hjust = 0.5))
# 
#   p2 <- dr_juice %>%
#     ggplot(aes(PC3, PC4)) +
#     geom_point(aes(color = group, shape = reference), alpha = 0.7, size = 2) +
#     scale_color_brewer(palette="Set1") +
#     # geom_text(aes(label = id), hjust = -0.2, vjust = 0, size = 1) +
#     labs(x = paste0("PC3 (", round(percent_variation[[3]] * 100, digits = 2), "% variance explained)"),
#          y = paste0("PC4 (", round(percent_variation[[4]] * 100, digits = 2), "% variance explained)"),
#          title = pca_title) +
#     theme_classic() +
#     theme(plot.title = element_text(hjust = 0.5))
# 
#   return(list("dr" = dr, "p1" = p1, "p2" = p2))
# }
# 
# message("Numbers of NA values in each sample: ")
# rate_tbl %>% group_by(sample_id) %>% summarise_all(~ sum(is.na(.x)))
# 
# rate_tbl_wo_chivu <- rate_tbl %>% filter(!str_detect(sample_id, "chivu"))
# 
# pca_chi <-
#   pca_and_visualization(df = rate_tbl_wo_chivu, col = "chi",
#                         pca_title = "Principal component analysis for alpha")
# ggsave(file.path(result_dir, "alpha_pca_between_samples.png"), plot = pca_alpha$p1,
#        width = 6, height = 4)
# 
# pca_beta <-
#   pca_and_visualization(df = rate_tbl_wo_chivu, col = "beta",
#                         pca_title = "Principal component analysis for beta")
# ggsave(file.path(result_dir, "beta_pca_between_samples.png"), plot = pca_beta$p1,
#        width = 6, height = 4)
