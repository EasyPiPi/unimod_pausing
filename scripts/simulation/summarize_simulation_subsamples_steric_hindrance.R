#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
meta_in <- snakemake@input[["meta_st"]]

table_dir <- snakemake@params[["table"]]
figure_dir <- snakemake@params[["figure"]]

#### load packages ####
library(tidyverse)
library(RColorBrewer)
library(cowplot)

#### testing files ####
# root_dir <- "~/Desktop/github/unimod_human"
# 
# meta_in <- file.path(root_dir, "metadata/simulation_params_steric_hindrance.csv")
# 
# table_dir <- file.path(root_dir, "outputs/simulation/tables/rate/steric_hindrance")
# figure_dir <- file.path(root_dir, "outputs/simulation/figures/steric_hindrance_varied_zeta")

#### end of parsing arguments ####
walk(c(table_dir, figure_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

#### set up parameters ####
phi_cap <- 1
kmin <- 1
kmax <- 200

k <- c("50" = 50, "70" = 70)
zeta <- 2000
#### end of setting up parameters ####

# read in simulation parameters in metadata frame
rate_tbls <- read_csv(meta_in)
colnames(rate_tbls) <-  c("id", "k", "ksd", "m", "a", "b", "g", "z", "t", "n", "s", "h", "l")

# read in potential and effective initiation rates based on simulations
rnap_tbl <- bind_rows(map(file.path(table_dir, paste0(rate_tbls$id, ".csv")),
                          read_csv, show_col_types = FALSE))
# scale it back to the unit of 1 min, 5e-5 min per slice, 100 slices
rnap_tbl <- rnap_tbl %>% mutate(across(.fn = ~ .x * 100))
rate_tbls <- rate_tbls %>% bind_cols(rnap_tbl)

# read in rate estimates
rate_tbls$tbl <- map(file.path(table_dir, paste0(rate_tbls$id, ".RDS")), readRDS)

# generate table to visualize results
plot_tbl <- rate_tbls %>%
  select(id, k, a, b, s, h, potential, actual, tbl) %>%
  mutate(rnap_size = s + h) %>% 
  dplyr::rename("initiation" = a, "pause" = b) %>%
  unnest(cols = tbl) %>%
  select(id, k, rnap_size, rnap_n, initiation, pause, chi, beta_org, beta_adp,
         R, R_pause, polII_prop, alpha_empirical, Xk, Yk, fk,
         potential, actual, phi) %>%
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
  mutate(across(.cols = c(beta_org, beta_adp), ~ .x * zeta, .names = "{.col}_zeta"))
# across(.cols = c(beta_org, beta_adp, alpha, omega), ~ .x * zeta))

# Method 2: match one of the estimates to known rates to get scaling factor
# chi
scale_factor_chi <- plot_tbl %>%
  # pick the one without pausing
  filter(initiation == 1, pause == 1000) %>%
  group_by(k, rnap_size) %>%
  summarise(chi = mean(chi))

# # beta could be scaled by the same logic, but we have scaled it using method 1
# scale_factor_beta <- plot_tbl %>%
#   filter(initiation == 1, pause == 1) %>%
#   group_by(r) %>%
#   summarise(beta_org = median(beta_org),
#             beta_adp = median(beta_adp))

plot_tbl_gb <- plot_tbl %>%
  group_by(k, rnap_size) %>%
  nest() %>%
  # left_join(scale_factor_beta, by = "r") %>%
  inner_join(scale_factor_chi, by = c("k", "rnap_size"))

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
                  df %>% mutate(omega_zeta = chi / a,
                                alpha_zeta_empirical = alpha_empirical / a)
                }))

## Estimate phi after EM, now has been incorporated into EM ##
# # get parameters for calculating phi
# calculate_f_and_kmean <- function(s, k) {
#   # sd is set as 25 here
#   x <- round(rnorm(1e4, mean = k, sd = 25))
#   x <- x[x >= 17 & x <= 200]
#   f <- mean(x > s) 
#   f1 <- mean((x > s) & (x <= 2 * s))
#   f2 <- mean(x > 2 * s)
#   return(c("f" = f, "f1" = f1, "f2" = f2))
# }
# 
# plot_tbl_gb <- plot_tbl_gb %>%
#   bind_cols(bind_rows(map2(.$rnap_size,
#                            .$k, calculate_f_and_kmean)))
# 
# plot_tbl <- plot_tbl_gb %>%
#   select(data, k, rnap_size, f, f1, f2) %>%
#   unnest(cols = data)
# 
# pi.func <- function(pi, beta, omega, f) {
#   retval <- (omega^2 - beta*omega) + (2*beta*omega - beta^2 - beta^2*f)*pi +
#             (2*beta^2 - beta^3/omega)*pi^2 + beta^3/omega * pi^3
#   return(retval)
# }
# 
# compute.phi = function(beta, omega, f) {
#   pi.root <- list("root" = NA_integer_)
#   try(pi.root <- uniroot(pi.func, c(-1,1), beta, omega, f), silent = TRUE)	
#   return (1-pi.root$root)
# }
# 
# phi.func.3RNAP = function(phi, beta, omega, f1, f2) {
#   alpha = omega/(1-phi)
#   return((1-f1-f2)*alpha/(alpha+beta) + f1*alpha^2/(alpha^2+beta^2+alpha*beta) +
#            f2*alpha^3/(beta^2*alpha + beta^3 + alpha^2*beta + alpha^3) - phi)	
# }
# 
# compute.phi.3RNAP = function(beta, omega, f1, f2) {
#   phi.root <- list("root" = NA_integer_)
#   try(phi.root <- uniroot(phi.func.3RNAP, c(-1, 0.999999999), beta, omega, f1, f2))	
#   return (phi.root$root)
# }

# estimate phi and potential initiation rate, alpha
plot_tbl <- plot_tbl_gb %>%
  ungroup() %>% 
  select(k, rnap_size, data) %>%
  unnest(data) %>% 
  mutate(omega = omega_zeta / zeta) %>%
  mutate(
    ## multi-state Markov Model to calculate phi ##
    # compute phi after EM, now has been incorporated into EM
    # phi = pmap_dbl(list(beta_adp_raw, omega, f), compute.phi),
    # phi = pmap_dbl(list(beta_adp_raw, omega, f1, f2), compute.phi.3RNAP),
    # cap some edge cases with phi slightly higher than 1
    # phi_cap = if_else(phi > phi_cap, phi_cap, phi),
    alpha_zeta = omega_zeta / (1 - phi)) %>%
  mutate(k = factor(paste0("k", k), levels = c("k50", "k70")))

# make names nicer in plots
plot_tbl <- plot_tbl %>%
  mutate(rnap_size = factor(paste0("s=", rnap_size), levels = c("s=33", "s=50", "s=70")))

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
      facet_wrap(rnap_size ~ ., scales = "free", nrow = 3, drop = FALSE)
    
    if (!is.null(yaxis_labels)) {
      p <- p + geom_hline(yintercept = yaxis_labels,
                          linetype="dashed", color = "gray")
    }
  } else if (facet == "type2") {
    p <- plot_tbl %>% ggplot(aes_string(x = x, y = y, color = fill)) +
      geom_boxplot() +
      scale_color_brewer(palette = color) +
      labs(x = xlab, y = ylab, color = colorlab) +
      facet_wrap(c("rnap_size", x), scales = "free", nrow = 3, drop = FALSE)
  }
  
  return(p)
}

# Chi and Omega
plot_tbl_gb <- plot_tbl %>% group_by(k) %>% nest()

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "chi", facet = "type2",
                                fill = "pause", color = "Set1",
                                xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                                ylab = expression("Estimated "*chi), yaxis_labels = yaxis_labels,
                                colorlab = expression("True "*beta*zeta)) +
      cowplot::theme_cowplot())

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("chi_freey", suffix, ".png")),
       plot = p, width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "alpha_zeta", facet = "type2",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*alpha*zeta), yaxis_labels = yaxis_labels,
                  colorlab = expression("True "*beta*zeta)) +
  cowplot::theme_cowplot())

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("alpha_freey", suffix, ".png")),
       plot = p, width = 16, height = 16)

# underestimate initiation rate due to steric hindrance
p <-  map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "omega_zeta", facet = "type2",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*omega*zeta), yaxis_labels = yaxis_labels,
                  colorlab = expression("True "*beta*zeta)) +
  cowplot::theme_cowplot())

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("omega_freey", suffix, ".png")),
       plot = p, width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "alpha_zeta_empirical", facet = "type2",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Empirical "*alpha*zeta), yaxis_labels = yaxis_labels,
                  colorlab = expression("True "*beta*zeta)) +
  cowplot::theme_cowplot())

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("empirical_alpha_freey", suffix, ".png")), plot = p,
       width = 16, height = 16)

# phi and R
p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "phi", facet = "type1",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*phi),
                  colorlab = expression("True "*beta*zeta)))

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("estimated_phi", suffix, ".png")), plot = p,
       width = 16, height = 16)

# p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "phi_cap", facet = "type1",
#                   fill = "pause", color = "Set1",
#                   xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
#                   ylab = expression("Capped "*phi),
#                   colorlab = expression("True "*beta*zeta)))
# 
# p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)
# 
# ggsave(file.path(figure_dir, paste0("estimated_phi_cap", suffix, ".png")), plot = p,
#        width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "polII_prop", facet = "type1",
                  fill = "pause", color = "Set1",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Empirical "*phi),
                  colorlab = expression("True "*beta*zeta)))

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("empirical_pause_site_occupied_by_polII", suffix, ".png")), plot = p,
       width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "R_pause", facet = "type1",
                                          fill = "pause", color = "Set1",
                                          xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                                          ylab = expression("Average number of RNAP cell"^{-1}*" gene"^{-1}),
                                          colorlab = expression("True "*beta*zeta)))

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("polII_per_pause_peak_per_cell", suffix, ".png")), plot = p,
       width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "initiation", y = "R", facet = "type1",
                  fill = "pause", color = "Purples",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Average number of RNAP cell"^{-1}*" gene"^{-1}),
                  colorlab = expression("True "*beta*zeta)))

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("polII_per_transcript_per_cell", suffix, ".png")), plot = p,
       width = 16, height = 16)

# compare empirical vs. estimated phi 
p <- map(plot_tbl_gb$data, ~ .x %>%
           ggplot(aes(polII_prop, phi, color = pause)) +
           geom_point() +
           labs(x = expression("Empirical "*phi),
                y = expression("Estimated "*phi),
                title = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                color = expression("True "*beta*zeta)) +
           geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
           facet_grid(rnap_size ~ initiation) +
           cowplot::theme_cowplot() +
           theme(plot.title = element_text(hjust = 0.5)))

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("phi_empirical_vs_markov_model", suffix, ".png")), plot = p,
       width = 20, height = 16)

# Pause release rate
p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "pause", y = "beta_org_zeta", facet = "type1",
                  fill = "initiation", color = "Reds",
                  xlab = expression("True "*beta*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*beta*zeta),
                  colorlab = expression("True "*alpha*zeta)))

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("beta_original_model", suffix, ".png")), plot = p,
       width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "pause", y = "beta_org_zeta", facet = "type2",
                  fill = "initiation", color = "Set1",
                  xlab = expression("True "*beta*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*beta*zeta),
                  colorlab = expression("True "*alpha*zeta)) +
  cowplot::theme_cowplot())

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("beta_original_model_freey", suffix, ".png")), plot = p,
       width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ boxplot_rate(.x, x = "pause", y = "beta_adp_zeta", facet = "type1",
                  fill = "initiation", color = "Reds",
                  xlab = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
                  ylab = expression("Estimated "*beta*zeta),
                  colorlab = expression("True "*alpha*zeta)))

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("beta_adapted_model", suffix, ".png")), plot = p,
       width = 16, height = 16)

p <- map(plot_tbl_gb$data, ~ .x %>%
  boxplot_rate(x = "pause", y = "beta_adp_zeta", facet = "type2",
               fill = "initiation", color = "Set1",
               xlab = expression("True "*beta*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
               ylab = expression("Estimated "*beta*zeta),
               colorlab = expression("True "*alpha*zeta)) +
  cowplot::theme_cowplot())

p <- plot_grid(plotlist = p, ncol = 1, labels = plot_tbl_gb$k)

ggsave(file.path(figure_dir, paste0("beta_adapted_model_freey", suffix, ".png")), plot = p,
       width = 16, height = 16)

# number of RNAP located before pause site
x <- plot_tbl %>%
  ungroup() %>% select(id, rnap_n) %>% unnest(rnap_n) %>%
  group_by(id) %>% count(rnap_n)

rnap_n_tbl <- plot_tbl %>%
  select(id, initiation, pause, k, rnap_size) %>%
  unique() %>% 
  left_join(x) %>% 
  group_by(id) %>%
  mutate(pct = n / sum(n))

rnap_n_tbl_gb <- rnap_n_tbl %>% group_by(k, rnap_size) %>% nest()

y <- map(rnap_n_tbl_gb$data, ~ .x %>%
      ggplot(aes(x = rnap_n, y = pct)) +
      geom_col() +
      facet_wrap(pause ~ initiation) +
      labs(title = expression(alpha*zeta*" on cols, "*beta*zeta*" on rows"))) 

walk2(file.path(figure_dir,
               paste0(rnap_n_tbl_gb$k, "_", rnap_n_tbl_gb$rnap_size, "_rnap_n.png")),
     y, ggsave, width = 12, height = 8)

#### figures for main text ####
# theme_set(cowplot::theme_cowplot(font_family = "Helvetica"))
theme_set(cowplot::theme_cowplot(font_family = "sans"))

# figure 2 rates estimates by original model, omega is effective initiation rate
# we call it alpha here since adapted model hasn't been introduced yet
p1 <- 
  map(plot_tbl_gb$data, function(plot_tbl) {
    plot_tbl %>%
      filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10)) %>%
      ggplot() +
      geom_boxplot(aes(x=initiation, y=omega_zeta, color = pause), lwd=1.5, outlier.shape = NA) +
      # geom_point(aes(x=initiation, y=actual, group = pause), position=position_dodge(width=0.75), color = "red") +
      labs(x = expression("True "*alpha*zeta), y = expression("Estimated "*alpha*zeta), color = expression("True "*beta*zeta)) +
      facet_wrap(rnap_size ~ initiation, scales = "free", nrow = 3) +
      geom_hline(data = subset(plot_tbl, initiation == 0.1), aes(yintercept = 0.1),
                 linetype = "dashed", color = "gray", lwd = 1.5) +
      geom_hline(data = subset(plot_tbl, initiation == 1), aes(yintercept = 1),
                 linetype = "dashed", color = "gray", lwd = 1.5) +
      geom_hline(data = subset(plot_tbl, initiation == 10), aes(yintercept = 10),
                 linetype = "dashed", color = "gray", lwd = 1.5) +
      # set xlim or ylim for each facet
      # https://stackoverflow.com/questions/63550588/ggplot2coord-cartesian-on-facets
      ggh4x::facetted_pos_scales(
        y = rep(list(
          scale_y_continuous(limits = c(0, 0.125)),
          scale_y_continuous(limits = c(0, 1.25)),
          scale_y_continuous(limits = c(0, 12.5))
        ), 3)
      ) +
      theme(# strip.text.x = element_blank(),
            # legend.position = c(0.83, 0.3),
            text = element_text(size = 25),
            axis.text = element_text(size = 25))
  })
  

p <- cowplot::plot_grid(plotlist = p1, ncol = 2, labels = plot_tbl_gb$k,
                        label_size = 25, vjust = -1, hjust = 0) +
  theme(plot.margin = unit(c(2,0.5,0.5,0.5), "cm"))

ggsave(file.path(figure_dir, paste0("fig2_alpha_original_model_freex", suffix, ".png")),
       plot = p, width = 24, height = 10)

p2 <-
  map(plot_tbl_gb$data, function(plot_tbl) {
    plot_tbl %>%
      filter(pause %in% c(0.1, 1, 10), initiation %in% c(0.1, 1, 10)) %>%
      ggplot(aes(x=pause, y=beta_org_zeta, color = initiation)) +
      geom_boxplot(lwd=1.5, outlier.shape = NA) +
      labs(x =  expression("True "*beta*zeta), y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
      facet_wrap(rnap_size ~ pause, scales = "free", nrow = 3) +
      geom_hline(data = subset(plot_tbl, pause == 0.1), aes(yintercept = 0.1),
                 linetype = "dashed", color = "gray", lwd = 1.5) +
      geom_hline(data = subset(plot_tbl, pause == 1), aes(yintercept = 1),
                 linetype = "dashed", color = "gray", lwd = 1.5) +
      geom_hline(data = subset(plot_tbl, pause == 10), aes(yintercept = 10),
                 linetype = "dashed", color = "gray", lwd = 1.5) +
      ggh4x::facetted_pos_scales(
        y = rep(list(
          scale_y_continuous(limits = c(0, 30)),
          scale_y_continuous(limits = c(0, 400)),
          scale_y_continuous(limits = c(0, 3000))
        ), 3)
      ) +
      theme(# strip.text.x = element_blank(),
            # axis.text.x = element_blank(),
            # legend.position = c(0.83, 0.8),
            text = element_text(size = 25),
            axis.text = element_text(size = 25))
    })
  
p <- cowplot::plot_grid(plotlist = p2, ncol = 2, labels = plot_tbl_gb$k,
                        label_size = 25, vjust = -1, hjust = 0) +
  theme(plot.margin = unit(c(2,0.5,0.5,0.5), "cm"))

ggsave(file.path(figure_dir, paste0("fig2_beta_original_model_freex", suffix, ".png")), plot = p,
       width = 24, height = 10)

# figure 3 beta estimates based on adapted model
p <- 
  map(plot_tbl_gb$data, function(plot_tbl) {
    plot_tbl %>%
          filter(pause %in% c(0.1, 1, 10), initiation %in% c(0.1, 1, 10), rnap_size %in% c("s=33", "s=70")) %>%
          ggplot(aes(x=pause, y=beta_adp_zeta, color = initiation)) +
          geom_boxplot(lwd = 1.5,  outlier.shape = NA) +
          labs(x = expression("True "*beta*zeta), y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
          facet_wrap(rnap_size ~ pause, scales = "free", nrow = 2, ) +
          geom_hline(data = subset(plot_tbl, pause == 0.1 & rnap_size %in% c("s=33", "s=70")), aes(yintercept = 0.1),
                     linetype = "dashed", color = "gray", lwd = 1.5) +
          geom_hline(data = subset(plot_tbl, pause == 1 & rnap_size %in% c("s=33", "s=70")), aes(yintercept = 1),
                     linetype = "dashed", color = "gray", lwd = 1.5) +
          geom_hline(data = subset(plot_tbl, pause == 10 & rnap_size %in% c("s=33", "s=70")), aes(yintercept = 10),
                     linetype = "dashed", color = "gray", lwd = 1.5) +
          ggh4x::facetted_pos_scales(
            y = rep(list(
              scale_y_continuous(limits = c(0, 0.2)),
              scale_y_continuous(limits = c(0, 2)),
              scale_y_continuous(limits = c(0, 30))
            ), 3)
          ) +
          theme(# strip.text.x = element_blank(),
                # axis.text.x = element_blank(),
                # legend.position = c(0.83, 0.8),
                text = element_text(size = 25),
                axis.text = element_text(size = 25))
    })

# plot k50 separately for supplement figure
ggsave(file.path(figure_dir, paste0("beta_adapted_model_k50_varied_spacing", suffix, ".pdf")), plot = p[[1]],
       width = 12, height = 8)

p <- cowplot::plot_grid(plotlist = p, ncol = 2, labels = plot_tbl_gb$k,
                        label_size = 25, vjust = -1, hjust = 0) +
  theme(plot.margin = unit(c(2,0.5,0.5,0.5), "cm"))

ggsave(file.path(figure_dir, paste0("fig3_beta_adapted_model_freex", suffix, ".png")), plot = p,
       width = 24, height = 10)

# figure 4 alpha estimates based on adapted model
p <- map(plot_tbl_gb$data, function(plot_tbl) {
  plot_tbl %>%
  filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10)) %>%
  ggplot(aes(x=initiation, y=alpha_zeta, color = pause)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +
  labs(x = expression("True "*alpha*zeta), y = expression("Estimated "*alpha*zeta), color = expression("True "*beta*zeta)) +
  facet_wrap(rnap_size ~ initiation, scales = "free", nrow = 3) +
  geom_hline(data = subset(plot_tbl, initiation == 0.1), aes(yintercept = 0.1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 1), aes(yintercept = 1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, initiation == 10), aes(yintercept = 10),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  ggh4x::facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(0, 10)),
      scale_y_continuous(limits = c(0, 50))
    ), 3)
  ) +
  theme(# strip.text.x = element_blank(),
        # legend.position = c(0.83, 0.3),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))
  })

p <- cowplot::plot_grid(plotlist = p, ncol = 2, labels = plot_tbl_gb$k,
                        label_size = 25, vjust = -1, hjust = 0) +
  theme(plot.margin = unit(c(2,0.5,0.5,0.5), "cm"))

ggsave(file.path(figure_dir, paste0("alpha_adapted_model_freex_k50_k70", suffix, ".png")), plot = p,
       width = 24, height = 10)

# subset k50 and s50
df <- plot_tbl %>%
  filter(initiation %in% c(0.1, 1, 10), pause %in% c(0.1, 1, 10), k == "k50", rnap_size == "s=50")

# figure 3 pause release
p <- df %>% ggplot(aes(x=pause, y=beta_adp_zeta, color = initiation)) +
  geom_boxplot(lwd = 1.5,  outlier.shape = NA) +
  labs(x = expression("True "*beta*zeta), y = expression("Estimated "*beta*zeta), color = expression("True "*alpha*zeta)) +
  facet_wrap(. ~ pause, scales = "free", nrow = 1) +
  geom_hline(data = subset(plot_tbl, pause == 0.1 & rnap_size %in% c("s=50")), aes(yintercept = 0.1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, pause == 1 & rnap_size %in% c("s=50")), aes(yintercept = 1),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  geom_hline(data = subset(plot_tbl, pause == 10 & rnap_size %in% c("s=50")), aes(yintercept = 10),
             linetype = "dashed", color = "gray", lwd = 1.5) +
  ggh4x::facetted_pos_scales(
    y = rep(list(
      scale_y_continuous(limits = c(0, 0.2)),
      scale_y_continuous(limits = c(0, 2)),
      scale_y_continuous(limits = c(0, 20))
    ), 3)
  ) +
  theme(strip.text.x = element_blank(),
    text = element_text(size = 25),
    axis.text = element_text(size = 25))

ggsave(file.path(figure_dir, paste0("beta_adapted_model_k50_s50", suffix, ".pdf")), plot = p,
       width = 12, height = 4)

# figure 4
# compare empirical vs. estimated phi
# for plotting means
mean_phi <-
  df %>% group_by(initiation, pause) %>%
  summarise(mean_polII_prop = mean(polII_prop, na.rm = TRUE),
            mean_phi = mean(phi, na.rm = TRUE))
# for labeling strip
col_lab <- paste0("αζ = ", c(0.1, 1, 10))
names(col_lab) <- c(0.1, 1, 10)

p <- df %>%
  ggplot(aes(polII_prop, phi, color = pause)) +
  geom_point(alpha = 0.4) +
  geom_segment(data = mean_phi, aes(x = mean_polII_prop, xend = mean_polII_prop,
                                    y = mean_phi - 0.05, yend = mean_phi + 0.05)) +
  geom_segment(data = mean_phi, aes(x = mean_polII_prop - 0.05, xend = mean_polII_prop + 0.05,
                                    y = mean_phi, yend =  mean_phi)) +
  labs(x = expression("Empirical "*phi),
       y = expression("Estimated "*phi),
       # title = expression("True "*alpha*zeta*" (RNAP min"^{-1}*" cell"^{-1}*")"),
       color = expression("True "*beta*zeta)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  facet_wrap(initiation ~ ., scales = "free_y", nrow = 1,
             labeller = labeller(initiation = col_lab)) +
  ylim(-0.05, 1.05) +
  xlim(-0.05, 1.05) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.text = element_text(size = 20))

# how to render greek letters correctly
# https://stackoverflow.com/questions/66762896/saving-ggplots-with-greek-alphabets-in-pdf
ggsave(file.path(figure_dir, paste0("fig4_phi_adapted_model_freex", suffix, ".pdf")), plot = p,
       width = 14, height = 5, device = cairo_pdf)

ggsave(file.path(figure_dir, paste0("fig4_phi_adapted_model_freex", suffix, ".png")), plot = p,
       width = 14, height = 5)
     
# alpha estimates
# # chekc infinities first
# df %>% group_by(initiation, pause) %>% summarise(sum(is.infinite(alpha)))
# df %>% group_by(initiation, pause) %>% summarise(sum(is.na(alpha)))
# make some big values to keep those combinations still have their spots in figure
# x <- df %>%
#   filter((initiation == 1 & pause == 0.1) | (initiation == 10 & pause %in% c(0.1, 1))) %>%
#   mutate(alpha = 99)

p <- df %>%
  # filter(!(initiation == 1 & pause == 0.1)) %>%
  # filter(!(initiation == 10 & pause %in% c(0.1, 1))) %>%
  # bind_rows(x) %>%
  filter(initiation %in% c(0.1, 1, 10)) %>%
  # filter(initiation %in% c(0.1, 1)) %>% 
  ggplot(aes(x=initiation, y=alpha_zeta, color = pause)) +
  geom_boxplot(lwd=1.5, outlier.shape = NA) +
  labs(x = expression("True "*alpha*zeta),  y = expression("Estimated "*alpha*zeta), color = expression("True "*beta*zeta)) +
  facet_wrap(initiation ~ ., scales = "free", drop = TRUE) +
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
    ), 3)
  ) +
  theme(strip.text.x = element_blank(),
        text = element_text(size = 20),
        axis.text = element_text(size = 20))

ggsave(file.path(figure_dir, paste0("fig4_alpha_adapted_model_freex", suffix, ".pdf")), plot = p,
       width = 12, height = 4)

#### supplement ####
# p <- plot_tbl %>%
#   filter(pause %in% c(0.1, 1, 10), initiation %in% c(0.1, 1, 10)) %>%
#   boxplot_rate(x = "initiation", y = "phi", facet = "type1",
#                fill = "pause", color = "Set1",
#                xlab = expression("True "*alpha*zeta),
#                ylab = expression("Estimated "*phi),
#                colorlab = expression("True "*beta*zeta))
# 
# ggsave(file.path(figure_dir, paste0("estimated_phi_k50_varied_spacing", suffix, ".pdf")), plot = p,
#        width = 6, height = 5)

mean_phi <-
  plot_tbl %>%
  filter(pause %in% c(0.1, 1, 10), initiation %in% c(0.1, 1, 10), k == "k50") %>%
  group_by(initiation, pause, rnap_size) %>%
  summarise(mean_polII_prop = mean(polII_prop, na.rm = TRUE),
            mean_phi = mean(phi, na.rm = TRUE))

p <- plot_tbl %>%
  filter(pause %in% c(0.1, 1, 10), initiation %in% c(0.1, 1, 10), k == "k50") %>%
  ggplot(aes(polII_prop, phi, color = pause)) +
  geom_segment(data = mean_phi, aes(x = mean_polII_prop, xend = mean_polII_prop,
                                    y = mean_phi - 0.05, yend = mean_phi + 0.05)) +
  geom_segment(data = mean_phi, aes(x = mean_polII_prop - 0.05, xend = mean_polII_prop + 0.05,
                                    y = mean_phi, yend =  mean_phi)) +
  geom_point(alpha = 0.4) +
  labs(x = expression("Empirical "*phi),
       y = expression("Estimated "*phi),
       # title = expression("True "*alpha*zeta),
       color = expression("True "*beta*zeta)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  facet_grid(rnap_size ~ initiation,
             labeller = labeller(initiation = col_lab)) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("estimated_vs_empirical_phi_k50_varied_spacing", suffix, ".pdf")), plot = p,
       width = 8, height = 6, device = cairo_pdf)
