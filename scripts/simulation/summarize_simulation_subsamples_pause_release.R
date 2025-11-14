#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
meta_pr_in <- snakemake@input[["meta_pr"]]
meta_st_in <- snakemake@input[["meta_st"]]

table_dir <- snakemake@params[["table"]]
figure_dir <- snakemake@params[["figure"]]

#### load packages ####
library(tidyverse)
library(RColorBrewer)

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_human"
# 
# meta_pr_in <- file.path(root_dir, "metadata/simulation_params_pause_release.csv")
# meta_st_in <- file.path(root_dir, "metadata/simulation_params_steric_hindrance.csv")
# 
# table_dir <- file.path(root_dir, "outputs/simulation/tables/rate/pause_release")
# figure_dir <- file.path(root_dir, "outputs/simulation/figures/pause_release_spacing50_varied_zeta")
# #### end of parsing arguments ####

walk(c(table_dir, figure_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

#### set up parameters ####
phi_cap <- 1
kmin <- 1
kmax <- 200
k <- 50
zeta <- 2000
#### end of setting up parameters ####

# read in simulation parameters in metadata frame
rate_tbls <- bind_rows(map(list(meta_pr_in, meta_st_in), read_csv, show_col_types = FALSE))
colnames(rate_tbls) <-  c("id", "k", "ksd", "m", "a", "b", "g", "z", "t", "n", "s", "h", "l")

# use a subset of combinations for the pause release section
rate_tbls <- rate_tbls %>%
  filter(k == 50, h == 17, ksd %in% c(0, 5, 15, 25))

# read in potential and effective initiation rates based on simulations
rnap_tbl <- bind_rows(map(file.path(table_dir, paste0(rate_tbls$id, ".csv")),
              read_csv, show_col_types = FALSE))
# scale it back to the unit of 1 min, 5e-5 min per slice, 100 slices
rnap_tbl <- rnap_tbl %>% mutate(across(.fn = ~ .x * 100))
rate_tbls <- rate_tbls %>% bind_cols(rnap_tbl)

# rate_tbls <- rate_tbls %>%
#   filter(k == 50, h == 37, ksd %in% c(0, 5, 15, 25))

# read in rate estimates
rate_tbls$tbl <- map(file.path(table_dir, paste0(rate_tbls$id, ".RDS")), readRDS)

# generate table to visualize results
plot_tbl <- rate_tbls %>%
  select(id, ksd, a, b, tbl, potential, actual) %>%
  mutate(ksd = case_when(
    ksd == 0 ~ "spike",
    ksd == 5 ~ "narrow",
    ksd == 15 ~ "median",
    ksd == 25 ~ "broad",
    TRUE ~ as.character(ksd)),
         ksd = factor(ksd, levels = c("spike", "narrow", "median", "broad"))) %>%
  dplyr::rename("initiation" = a, "pause" = b) %>%
  unnest(cols = tbl) %>%
  select(id, ksd, initiation, pause, chi, beta_org, beta_adp, beta_max_rc,
         R, R_pause, polII_prop, alpha_empirical, Xk, Yk, fk, potential, actual) %>%
  mutate(t = map_dbl(Yk, ~ sum(.x)))

## scaling rates for better visualization ##
# # Optional: alpha could be calculated by using t if adapted model has been introduced
# # calculate mean of t for genes with highest initiation and lowest pause release rates
# mean_t <- plot_tbl %>%
#   filter(initiation == max(.$initiation), pause == min(.$pause)) %>%
#   summarise(t = mean(t, na.rm = TRUE))
# 
# plot_tbl <- plot_tbl %>%
#   mutate(# get effetive initiation rate omega
#     omega = chi / mean_t$t)

# Method 1: multiply rates with elongation rates, beta could be directly scaled
# while alpha has to be scaled first by the above method, otherwise alpha go with 
# method 2
plot_tbl <- plot_tbl %>%
  mutate(across(.cols = c(beta_org, beta_adp, beta_max_rc), ~ .x * zeta, .names = "{.col}_zeta"))
         # across(.cols = c(beta_org, beta_adp, alpha, omega), ~ .x * zeta))

# Method 2: match one of the estimates to known rates to get scaling factor
# chi
scale_factor_chi <- plot_tbl %>%
  # pick the one without pausing
  filter(initiation == 1, pause == 1000) %>%
  group_by(ksd) %>%
  summarise(chi = median(chi))

# # beta could be scaled by the same logic, but we have scaled it using method 1
# scale_factor_beta <- plot_tbl %>%
#   filter(initiation == 1, pause == 1) %>%
#   group_by(ksd) %>%
#   summarise(beta_org = median(beta_org),
#             beta_adp = median(beta_adp))

plot_tbl_gb <- plot_tbl %>%
  group_by(ksd) %>%
  nest() %>%
  # left_join(scale_factor_beta, by = "ksd") %>%
  left_join(scale_factor_chi, by = "ksd")

# plot_tbl_gb <- plot_tbl_gb %>%
#   mutate(data =
#            pmap(list(data, chi, beta_org, beta_adp),
#          function(df, a, b, c) {
#            df %>% mutate(alpha = chi / a,
#                          alpha_emprical = alpha_empirical / a,
#                          beta_org_raw = beta_org,
#                          beta_adp_raw = beta_adp,
#                          # use adapted model to scale beta
#                          beta_org = beta_org / c,
#                          beta_adp = beta_adp / c)
#          }))

plot_tbl_gb <- plot_tbl_gb %>%
  mutate(data =
           pmap(list(data, chi),
                function(df, a) {
                  df %>% mutate(omega = chi / (a * zeta),
                                alpha_zeta_empirical = alpha_empirical / a)
                }))

plot_tbl <- plot_tbl_gb %>%
  select(data, ksd) %>%
  unnest(cols = data)

# estimate phi and potential initiation rate, alpha
plot_tbl <- plot_tbl %>%
  mutate(
    ## multi-state Markov Model to calculate phi ##
    # this is a version to calculate phi without EM
    phi = omega / beta_adp + omega * (k - 1),
    # cap some edge cases with phi slightly higher than 1
    phi_cap = if_else(phi > phi_cap, phi_cap, phi),
    alpha = omega / (1 - phi_cap),
    omega_zeta = omega * zeta,
    alpha_zeta = alpha * zeta)

# filter out some cases
plot_tbl <- plot_tbl %>% 
  # filter out unrealistic extreme values (i.e., 0.01, 1000)
 filter(!pause %in% c(0.01, 1000), !initiation %in% c(0.01, 1000))

#### plots to explore the full parameter landscape ####
# plot estimated alpha and beta
# to adjust the color
# https://stackoverflow.com/questions/27055094/can-i-adjust-the-lower-limit-of-scale-color-brewer
# labels <- c(0.01, 0.1, 1, 10, 100, 1000)
# yaxis_labels <- c(0.1, 1, 10)
labels <- yaxis_labels <- c(0.1, 1, 10, 50, 100)

plot_tbl <- plot_tbl %>%
  mutate(initiation = factor(initiation, levels = labels),
         pause = factor(pause, levels = labels))

suffix <- "_subsample_cells"

# handy function to avoid repeated code
boxplot_rate <- function(plot_tbl, x, y, fill, color, xlab, ylab, facet,
                         yaxis_labels = NULL, colorlab = NULL) {
  
  if (facet == "type1") {
    p <- plot_tbl %>% ggplot(aes_string(x = x, y = y, color = fill)) +
      geom_boxplot() +
      scale_color_brewer(palette = color) +
      labs(x = xlab, y = ylab, color = colorlab) +
      cowplot::theme_cowplot(font_family = "Helvetica") +
      facet_wrap(ksd ~ ., scales = "free", nrow = 2)
    
    if (!is.null(yaxis_labels)) {
      p <- p + geom_hline(yintercept = yaxis_labels,
                          linetype="dashed", color = "gray")
      }
    } else if (facet == "type2") {
    p <- plot_tbl %>% ggplot(aes_string(x = x, y = y, color = fill)) +
        geom_boxplot() +
        scale_color_brewer(palette = color) +
        labs(x = xlab, y = ylab, color = colorlab) +
        facet_wrap(c("ksd", x), scales = "free", nrow = 2)
    }
  
  return(p)
}

# Chi and Omega
p <- boxplot_rate(plot_tbl, x = "initiation", y = "chi", facet = "type2",
             fill = "pause", color = "Set1",
             xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
             ylab = expression("Estimated "*chi), yaxis_labels = yaxis_labels,
             colorlab = expression("True "*beta*zeta)) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("chi_freey", suffix, ".png")),
       plot = p, width = 16, height = 8)

p <- boxplot_rate(plot_tbl, x = "initiation", y = "alpha_zeta", facet = "type2",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*alpha*zeta), yaxis_labels = yaxis_labels,
                  colorlab = expression("True "*beta*zeta)) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("alpha_freey", suffix, ".png")),
       plot = p, width = 16, height = 8)

# underestimate initiation rate due to steric hindrance
p <- boxplot_rate(plot_tbl, x = "initiation", y = "omega_zeta", facet = "type2",
             fill = "pause", color = "Set1",
             xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
             ylab = expression("Estimated "*omega*zeta), yaxis_labels = yaxis_labels,
             colorlab = expression("True "*beta*zeta)) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("omega_freey", suffix, ".png")),
       plot = p, width = 16, height = 8)

p <- boxplot_rate(plot_tbl, x = "initiation", y = "alpha_zeta_empirical", facet = "type2",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Empirical "*alpha*zeta), yaxis_labels = yaxis_labels,
                  colorlab = expression("True "*beta*zeta)) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("empirical_alpha_freey", suffix, ".png")), plot = p,
       width = 16, height = 8)

# phi and R
p <- boxplot_rate(plot_tbl, x = "initiation", y = "phi", facet = "type1",
                  fill = "pause", color = "Greens",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*phi),
                  colorlab = expression("True "*beta*zeta))

ggsave(file.path(figure_dir, paste0("estimated_phi", suffix, ".png")), plot = p,
       width = 16, height = 8)

p <- boxplot_rate(plot_tbl, x = "initiation", y = "phi_cap", facet = "type1",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Capped "*phi),
                  colorlab = expression("True "*beta*zeta))

ggsave(file.path(figure_dir, paste0("estimated_phi_cap", suffix, ".png")), plot = p,
       width = 16, height = 8)

p <- boxplot_rate(plot_tbl, x = "initiation", y = "polII_prop", facet = "type1",
                  fill = "pause", color = "Greens",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Empirical "*phi),
                  colorlab = expression("True "*beta*zeta))

ggsave(file.path(figure_dir, paste0("empirical_pause_site_occupied_by_polII", suffix, ".png")), plot = p,
       width = 16, height = 8)

p <- boxplot_rate(plot_tbl, x = "initiation", y = "R", facet = "type1",
                  fill = "pause", color = "Purples",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Average number of RNAP cell"^{-1}*" gene"^{-1}),
                  colorlab = expression("True "*beta*zeta))

ggsave(file.path(figure_dir, paste0("polII_per_transcript_per_cell", suffix, ".png")), plot = p,
       width = 16, height = 8)

# compare empirical vs. estimated phi 
p <- plot_tbl %>%
  ggplot(aes(polII_prop, phi, color = pause)) +
  geom_point() +
  labs(x = expression("Empirical "*phi),
       y = expression("Estimated "*phi),
       title = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
       color = expression("True "*beta*zeta)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  facet_grid(ksd ~ initiation) +
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("phi_empirical_vs_markov_model", suffix, ".png")), plot = p,
       width = 12, height = 8)

# Pause release rate
p <- boxplot_rate(plot_tbl, x = "pause", y = "beta_org_zeta", facet = "type1",
                  fill = "initiation", color = "Reds",
                  xlab = expression("True "*beta*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*beta*zeta),
                  colorlab = expression("True "*alpha*zeta))

ggsave(file.path(figure_dir, paste0("beta_original_model", suffix, ".png")), plot = p,
       width = 16, height = 8)

p <- boxplot_rate(plot_tbl, x = "pause", y = "beta_org_zeta", facet = "type2",
                  fill = "initiation", color = "Set1",
                  xlab = expression("True "*beta*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*beta*zeta),
                  colorlab = expression("True "*alpha*zeta)) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("beta_original_model_freey", suffix, ".png")), plot = p,
       width = 16, height = 8)

p <- boxplot_rate(plot_tbl, x = "pause", y = "beta_adp_zeta", facet = "type1",
                  fill = "initiation", color = "Reds",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*beta*zeta),
                  colorlab = expression("True "*alpha*zeta))

ggsave(file.path(figure_dir, paste0("beta_adapted_model", suffix, ".png")), plot = p,
       width = 16, height = 8)

p <- plot_tbl %>%
  boxplot_rate(x = "pause", y = "beta_adp_zeta", facet = "type2",
               fill = "initiation", color = "Set1",
               xlab = expression("True "*beta*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
               ylab = expression("Estimated "*beta*zeta),
               colorlab = expression("True "*alpha*zeta)) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("beta_adapted_model_freey", suffix, ".png")), plot = p,
       width = 16, height = 8)

# compare beta estimates from the original and adapted models
p <- plot_tbl %>%
  mutate(beta_org_zeta = log10(beta_org_zeta),
         beta_adp_zeta = log10(beta_adp_zeta)) %>%
  ggplot(aes(beta_org_zeta, beta_adp_zeta, color = pause)) +
  geom_point() +
  coord_cartesian(xlim = c(-2, 6), ylim = c(-2, 6)) +
  labs(x = expression(log[10]*beta*zeta*" (original model)"),
       y = expression(log[10]*beta*zeta*" (adapted model)"), color = expression("True "*beta*zeta),
       title = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  facet_grid(ksd ~ initiation) +
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("beta_original_vs_adapted_model", suffix, ".png")), plot = p,
       width = 10, height = 8)

#### figures for main text ####
# theme_set(cowplot::theme_cowplot(font_family = "Helvetica"))
theme_set(cowplot::theme_cowplot(font_family = "sans"))

# figure 2 rates estimates by original model, omega is effective initiation rate
# we call it alpha here since adapted model hasn't been introduced yet
p1 <- plot_tbl %>%
  filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10), ksd == "broad") %>%
  ggplot() +
  geom_boxplot(aes(x=initiation, y=omega_zeta, color = pause), lwd=1.5, outlier.shape = NA) +
  # geom_point(aes(x=initiation, y=actual, group = pause), position=position_dodge(width=0.75), color = "red") +
  labs(x = expression("True "*alpha*zeta), y = expression("Estimated "*alpha*zeta), color = expression("True "*beta*zeta)) +
  # geom_abline(intercept = c(0.1, 1, 10), slope = 0, linetype = "dashed", color = "gray", lwd = 1.5) +
  facet_wrap(. ~ initiation, scales = "free", nrow = 1) +
  geom_hline(data = subset(plot_tbl, initiation == 0.1), aes(yintercept = 0.1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 1), aes(yintercept = 1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 10), aes(yintercept = 10),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  # set xlim or ylim for each facet
  # https://stackoverflow.com/questions/63550588/ggplot2coord-cartesian-on-facets
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 0.125)),
      scale_y_continuous(limits = c(0, 1.25)),
      scale_y_continuous(limits = c(0, 12.5))
    )
  ) +
  theme(strip.text.x = element_blank(),
        # legend.position = c(0.83, 0.3),
        text = element_text(size = 25),
        axis.text = element_text(size = 20))

x <- plot_tbl %>%
  filter(pause == 1, initiation %in% c(0.1, 1, 10), ksd %in% c("spike", "broad")) 

x1 <- x %>%
  filter(ksd == "spike") %>%
  select(ksd, initiation, beta_org_zeta) %>%
  rename(beta = beta_org_zeta)

x2 <- x %>%
  filter(ksd == "broad") %>%
  mutate(ksd = "broad_ave") %>% 
  select(ksd, initiation, beta_org_zeta) %>%
  rename(beta = beta_org_zeta)

x3 <- x %>%
  filter(ksd == "broad") %>%
  mutate(ksd = "broad_max") %>% 
  select(ksd, initiation, beta_max_rc_zeta) %>%
  rename(beta = beta_max_rc_zeta)

# sum of Xk
x4 <- x %>%
  filter(ksd == "broad") %>%
  mutate(beta = chi / map_dbl(Xk, sum) * zeta) %>%
  select(ksd, initiation, beta) 

x <- bind_rows(x1, x2, x4) %>%
  mutate(ksd = factor(ksd, levels = c("spike", "broad_ave", "broad")))

p2 <- x %>%
  ggplot(aes(x=ksd, y=beta, color = initiation)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +
  labs(x = "", y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed", color = "gray", lwd = 1.5) +
  facet_wrap(. ~ ksd, scales = "free", nrow = 1) +
  scale_x_discrete(labels=c("broad_ave" = "Variable k (ave.)", 
                            "broad" = "Variable k (sum)", 
                            "spike" = "Fixed k")) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 2)),
      scale_y_continuous(limits = c(0, 400)),
      scale_y_continuous(limits = c(0, 2))
  )) +
  theme(strip.text.x = element_blank(),
        # axis.text.x = element_blank(),
        # legend.position = c(0.83, 0.8),
        text = element_text(size = 25),
        axis.text = element_text(size = 20))

p3 <- x4 %>%
  ggplot(aes(x=ksd, y=beta, color = initiation)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +
  labs(x = "", y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed", color = "gray", lwd = 1.5) +
  scale_x_discrete(labels=c("broad" = "Variable k (sum)")) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 20))

ggsave(file.path(figure_dir, paste0("fig2_original_model_sum_of_Xk", suffix, ".png")), plot = p3,
       width = 5, height = 5)

p <- cowplot::plot_grid(p1, NULL, p2, ncol = 1, align = "hv", axis = "tblr",
                        rel_heights = c(10, 1, 10)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

ggsave(file.path(figure_dir, paste0("fig2_original_model_freex", suffix, ".png")), plot = p,
       width = 12, height = 10)

# figure 3 beta estimates based on adapted model
p <- plot_tbl %>%
  filter(pause == 1, initiation %in% c(0.1, 1, 10), ksd %in% c("spike", "broad")) %>%
  ggplot(aes(x=ksd, y=beta_adp_zeta, color = initiation)) +
  geom_boxplot(lwd = 1.5,  outlier.shape = NA) +
  labs(x = "", y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed", color = "gray", lwd = 1.5) +
  coord_cartesian(ylim = c(0, 2.5)) +
  facet_wrap(. ~ ksd, scales = "free_x", nrow = 1) +
  scale_x_discrete(labels=c("broad" = "Variable k", "spike" = "Fixed k")) +
  theme(strip.text.x = element_blank(),
        # axis.text.x = element_blank(),
        # legend.position = c(0.83, 0.8),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

ggsave(file.path(figure_dir, paste0("fig3_beta_adapted_model_freex", suffix, ".png")), plot = p,
       width = 8, height = 4)

# figure 3 Xk and fk 
plot_Xk <- function(plot_tbl, ksd_name, ksd) {
  df <- plot_tbl %>%
    filter(initiation == 10, pause == 0.1, ksd == ksd_name) %>%
    slice_head(n = 1) %>%
    slice_tail(n = 1) %>% 
    unnest(cols = c(Xk, fk))
  
  scale_yaxis <- sum(df$Xk)
  
  # plot first 100bp for clarity
  kmax_p <- min(100, kmax)
  
  p <- df[kmin:kmax_p,] %>%
    mutate(Xk = Xk / scale_yaxis) %>% ggplot() +
    geom_col(aes(x = seq(kmin:kmax_p), y = Xk)) +
    geom_line(aes(x = seq(kmin:kmax_p), y = fk), color = "red") +
    stat_function(fun = dnorm, args = list(k, ksd), color = "blue", linetype = "dashed") +
    scale_y_continuous(
      "Xk",
      labels = ~ . * scale_yaxis,
      sec.axis = sec_axis(~ ., name = "fk"),
      breaks = seq(0, 20, 2) / scale_yaxis
    ) +
    labs(x = "") +
    cowplot::theme_cowplot()
  
  return(p)
}

p <- plot_Xk(plot_tbl, "narrow", 5)
ggsave(file.path(figure_dir, paste0("fig3_simulated_read_count_narrow_pause_peak", suffix, ".png")), plot = p,
       width = 6, height = 2)

p <- plot_Xk(plot_tbl, "median", 15)
ggsave(file.path(figure_dir, paste0("fig3_simulated_read_count_median_pause_peak", suffix, ".png")), plot = p,
       width = 6, height = 2)

p <- plot_Xk(plot_tbl, "broad", 25)
ggsave(file.path(figure_dir, paste0("fig3_simulated_read_count_broad_pause_peak", suffix, ".png")), plot = p,
       width = 6, height = 2)

# figure 4
p <- plot_tbl %>%
  filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10)) %>%
  ggplot(aes(x=initiation, y=alpha_zeta, color = pause)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +
  labs(x = expression("True "*alpha*zeta), y = expression("Estimated "*alpha*zeta), color = expression("True "*beta*zeta)) +
  facet_wrap(ksd ~ initiation, scales = "free", nrow = 2) +
  geom_hline(data = subset(plot_tbl, initiation == 0.1), aes(yintercept = 0.1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 1), aes(yintercept = 1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 10), aes(yintercept = 10),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  ggh4x::facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(0, 2)),
      scale_y_continuous(limits = c(0, 20))
    ), 4)
  ) +
  theme(# strip.text.x = element_blank(),
    # legend.position = c(0.83, 0.3),
    text = element_text(size = 25),
    axis.text = element_text(size = 25))

ggsave(file.path(figure_dir, paste0("alpha_adapted_model_freey", suffix, ".png")), plot = p,
       width = 24, height = 10)

#### supplement ####
# for fig. 2
p1 <- plot_tbl %>%
  filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10), ksd == "spike") %>%
  ggplot() +
  geom_boxplot(aes(x=initiation, y=polII_prop, color = pause), lwd=1.5, outlier.shape = NA) +
  labs(x = expression("True "*alpha*zeta), y = "landing-pad occupancy", color = expression("True "*beta*zeta)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20))

p2 <- plot_tbl %>% 
  filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10), ksd == "spike") %>%
  ggplot() +
  geom_point(aes(x = polII_prop, y = actual / potential, color = pause)) +
  geom_abline(intercept = 1, slope = -1, color = "gray", linetype = "dashed", lwd=1.5) +
  labs(x = "landing-pad occupancy", y = expression("Effective initiation rate"/"True "*alpha*zeta),
       color = expression("True "*beta*zeta)) +
  theme(strip.text.x = element_blank(),
        # legend.position = c(0.83, 0.3),
        text = element_text(size = 20),
        axis.text = element_text(size = 20))

p3 <- plot_tbl %>% 
  filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10), ksd == "spike") %>%
  ggplot() +
  geom_point(aes(x = actual, y = omega_zeta, color = pause, shape = initiation), size = 2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed", lwd=1.5) +
  labs(x = "Effective initiation rate", y = expression("Estimated "*alpha*zeta),
       color = expression("True "*beta*zeta), shape = expression("True "*alpha*zeta)) +
  facet_wrap(. ~ pause, scales = "free") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(0, 0.15)),
      scale_x_continuous(limits = c(0, 1)),
      scale_x_continuous(limits = c(0, 5))
    ),
    y = list(
      scale_y_continuous(limits = c(0, 0.15)),
      scale_y_continuous(limits = c(0, 1)),
      scale_y_continuous(limits = c(0, 5))
    )
  ) +
  theme(strip.text.x = element_blank(),
        # legend.position = c(0.83, 0.3),
        text = element_text(size = 20),
        axis.text = element_text(size = 20))

p <- cowplot::plot_grid(p1, p3, p2, ncol = 1) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) 

ggsave(file.path(figure_dir, paste0("fig2_supplement", suffix, ".png")), plot = p,
       width = 15, height = 15)

# for fig.2B and fig. 3A
x <- plot_tbl %>%
  filter(pause %in% c(0.1, 1, 10), initiation %in% c(0.1, 1, 10), ksd %in% c("spike", "broad")) 

x1 <- x %>%
  filter(ksd == "spike") %>%
  select(ksd, initiation, beta_org_zeta, pause) %>%
  rename(beta = beta_org_zeta)

x2 <- x %>%
  filter(ksd == "broad") %>%
  mutate(ksd = "broad_ave") %>% 
  select(ksd, initiation, beta_org_zeta, pause) %>%
  rename(beta = beta_org_zeta)

x3 <- x %>%
  filter(ksd == "broad") %>%
  mutate(ksd = "broad_max") %>% 
  select(ksd, initiation, beta_max_rc_zeta, pause) %>%
  rename(beta = beta_max_rc_zeta)

# sum of Xk
x4 <- x %>%
  filter(ksd == "broad") %>%
  mutate(ksd = "broad_sum") %>% 
  mutate(beta = chi / map_dbl(Xk, sum) * zeta) %>%
  select(ksd, initiation, beta, pause) 

x <- bind_rows(x1, x2, x3, x4) %>%
  mutate(ksd = factor(ksd, levels = c("spike", "broad_ave", "broad_max", "broad_sum")))

p <- x %>%
  ggplot(aes(x=ksd, y=beta, color = initiation)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +
  labs(x = "", y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
  facet_wrap(pause ~ ksd, scales = "free", nrow = 3, drop = TRUE) +
  scale_x_discrete(labels=c("broad_ave" = "Variable k (ave.)", "spike" = "Fixed k",
                            "broad_sum" = "Variable k (sum)", "broad_max" = "Variable k (max)")) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(0, 40)),
      scale_y_continuous(limits = c(0, 4)),
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(0, 2)),
      scale_y_continuous(limits = c(0, 400)),
      scale_y_continuous(limits = c(0, 40)),
      scale_y_continuous(limits = c(0, 2)),
      scale_y_continuous(limits = c(0, 20)),
      scale_y_continuous(limits = c(0, 3000)),
      scale_y_continuous(limits = c(0, 200)),
      scale_y_continuous(limits = c(0, 20))
      )
  ) +
  geom_hline(data = subset(x, initiation == 0.1 & pause == 0.1 & ksd %in% c("spike", "broad_ave", "broad_max", "broad_sum")),
             aes(yintercept = 0.1), linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(x, initiation == 0.1 & pause == 1 & ksd %in% c("spike", "broad_ave", "broad_max", "broad_sum")),
             aes(yintercept = 1), linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(x, initiation == 0.1 & pause == 10 & ksd %in% c("spike", "broad_ave", "broad_max", "broad_sum")),
             aes(yintercept = 10), linetype = "dashed", color = "gray", lwd = 1.5) +
  theme(strip.text.x = element_blank(),
        # axis.text.x = element_blank(),
        # legend.position = c(0.83, 0.8),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

ggsave(file.path(figure_dir, paste0("fig2_supplement_beta_original_model", suffix, ".png")), plot = p,
       width = 20, height = 12)
ggsave(file.path(figure_dir, paste0("fig2_supplement_beta_original_model", suffix, ".pdf")), plot = p,
       width = 20, height = 12)

p <- plot_tbl %>%
  filter(pause %in% c(0.1, 1, 10), initiation %in% c(0.1, 1, 10), ksd %in% c("spike", "broad")) %>%
  ggplot(aes(x=ksd, y=beta_adp_zeta, color = initiation)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +
  labs(x = "", y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
  facet_wrap(pause ~ ksd, scales = "free", nrow = 3) +
  scale_x_discrete(labels=c("broad" = "Variable k", "spike" = "Fixed k")) +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(0, 2.5)),
      scale_y_continuous(limits = c(0, 2.5)),
      scale_y_continuous(limits = c(0, 30)),
      scale_y_continuous(limits = c(0, 30))
    )
  ) +
  geom_hline(data = subset(plot_tbl, initiation == 0.1 & pause == 0.1 & ksd %in% c("spike", "broad")),
             aes(yintercept = 0.1), linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 0.1 & pause == 1 & ksd %in% c("spike", "broad")),
             aes(yintercept = 1), linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 0.1 & pause == 10 & ksd %in% c("spike", "broad")),
             aes(yintercept = 10), linetype = "dashed", color = "gray", lwd = 1.5) +
  theme(strip.text.x = element_blank(),
        # axis.text.x = element_blank(),
        # legend.position = c(0.83, 0.8),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

ggsave(file.path(figure_dir, paste0("fig2_supplement_beta_adapted_model", suffix, ".png")), plot = p,
       width = 10, height = 12)
ggsave(file.path(figure_dir, paste0("fig2_supplement_beta_adapted_model", suffix, ".pdf")), plot = p,
       width = 10, height = 12)
