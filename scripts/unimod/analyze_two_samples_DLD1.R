#### note ####
# this script currently takes care of samples with or without spike-ins
# note that it uses sum of reads across samples

#### log file ####
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
rc1_in <- snakemake@input[["rc1"]]
rc2_in <- snakemake@input[["rc2"]]

spike_in <- snakemake@input[["spike_in"]]

helper_tc_in <- snakemake@params[["helper_tc"]]
helper_pr_in <- snakemake@params[["helper_pr"]]

result_dir <- snakemake@params[["result_dir"]]

chi_out <- snakemake@output[["chi"]]
beta_out <- snakemake@output[["beta"]]
fk_out <- snakemake@output[["fk"]]
#### end of parsing arguments ####

#### load packages ####
library(tidyverse)
library(ggpubr)
library(ggpointdensity)
library(viridis)
library(ggExtra)

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_pausing"

# helper_tc_in <- file.path(root_dir, "scripts/unimod/helper_function_em_two_condition.R")
# helper_pr_in <- file.path(root_dir, "scripts/unimod/helper_function_em_pause_release.R")

# rc1_in <-
#   file.path(
#     root_dir,
#     "outputs/within_sample/PROseq-DLD1-aoi-NELFC_NVP2_Ctrl-SE/pause_release/rate.RDS"
#   )
# rc2_in <-
#   file.path(
#     root_dir,
#     "outputs/within_sample/PROseq-DLD1-aoi-NELFC_NVP2-SE/pause_release/rate.RDS"
#   )

# # rc1_in <- file.path(root_dir, "outputs/within_sample/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE/pause_release/rate.RDS")
# # rc2_in <- file.path(root_dir, "outputs/within_sample/PROseq-DLD1-aoi-NELFC_Fp-SE/pause_release/rate.RDS")

# # rc1_in <- file.path(root_dir, "outputs/within_sample/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE/pause_release/rate.RDS")
# # rc2_in <- file.path(root_dir, "outputs/within_sample/PROseq-DLD1-aoi-NELFC_Auxin-SE/pause_release/rate.RDS")

# spike_in <- file.path(root_dir, "metadata/scaling_factor.csv")

# result_dir <-
#   file.path(
#     root_dir, "outputs/between_samples",
#     paste0("NELFC_NVP2_Ctrl", "_vs_", "NELFC_NVP2")
#   )

# chi_out <- file.path(result_dir, "chi.csv")
# beta_out <- file.path(result_dir, "beta.csv")
# fk_out <- file.path(result_dir, "fk.csv")
# #### end of parsing arguments ####

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
source(helper_tc_in)
source(helper_pr_in)

## set up parameters ##
# plotting parameters
theme_set(cowplot::theme_cowplot())
scatter_colors <- c("#6e8fb2", "gray", "#c16e71")
paired_colors <- c(
  "#EAB67A", "#F5D8B7",
  "#7DA494", "#B6C9C0",
  "#58539f", "#bbbbd6",
  "#d86967", "#eebabb"
)
bar_colors <- c("#00AFBB", "#E7B800", "#FC4E07")
# scatter_colors <- c("#E41A1C", "#377EB8", "gray")

# model parameters
k <- 50
kmin <- 1
kmax <- 200 # also used as k on the poisson case

rnap_size <- 50
zeta <- 2000

# criteria for significance, p.adj < 0.05 and 20% differences
sig_p <- 0.05
lfc1 <- log2(1.2 / 1)
lfc2 <- log2(0.8 / 1)

## read in data ##
# read in summary statistics from the varied pause site model
rc1 <- readRDS(rc1_in)
rc2 <- readRDS(rc2_in)

# union set of genes being analyzed
gn_union <- intersect(rc1$gene_id, rc2$gene_id)
rc1 <- rc1[match(gn_union, rc1$gene_id), ]
rc2 <- rc2[match(gn_union, rc2$gene_id), ]

# read in spike-in scaling factors then compute lambda1 and lambda2
scale_tbl <- read_csv(spike_in, show_col_types = FALSE)
scale_tbl <- scale_tbl[str_detect(rc2_in, scale_tbl$sample), ]

lambda1 <-
  scale_tbl$control_1 + ifelse(is.na(scale_tbl$control_2), 0, scale_tbl$control_2)
lambda2 <-
  scale_tbl$treated_1 + ifelse(is.na(scale_tbl$treated_2), 0, scale_tbl$treated_2)

## LRT for Chi ##
tao1 <- lambda1 / (lambda1 + lambda2)
tao2 <- 1 - tao1

chi_tbl <-
  tibble(
    gene_id = rc1$gene_id,
    chi1 = rc1$chi,
    chi2 = rc2$chi * lambda1 / lambda2,
    lfc = log2(chi2 / chi1)
  )

chi_tbl <- chi_tbl %>%
  bind_cols(bind_rows(map2(rc1$s, rc2$s, omega_lrt, tao1 = tao1, tao2 = tao2)))

chi_tbl <- chi_tbl %>%
  mutate(padj = p.adjust(p, method = "BH"))

chi_tbl <- chi_tbl %>%
  mutate(
    chi = (chi1 + chi2) / 2,
    logchi = log2(chi),
    category =
      case_when(
        (padj < 0.05) & (lfc > lfc1) ~ "Up",
        (padj < 0.05) & (lfc < lfc2) ~ "Down",
        TRUE ~ "Others"
      ),
    category = factor(category, levels = c("Up", "Down", "Others"))
  )

## LRT for beta ##
# need to jointly do EM one more time for H0, which assume betas are the same
# between conditions
# initialize fk with some reasonable values based on heuristic
fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
fk_int <- fk_int / sum(fk_int)
# try uniform distribution which gives similar results
# fk_int <- rep(1 / kmax, kmax)

# collect and construct stats
s1 <- rc1$s
s2 <- rc2$s
t1_h1 <- map_dbl(rc1$Yk, sum)
t2_h1 <- map_dbl(rc2$Yk, sum)
Xk1 <- rc1$Xk
Xk2 <- rc2$Xk
# gene body length, assumed to be same between conditions
M <- rc1$N

# some values inherit from EM for H1, could be further integrated into EM here
chi_hat <- (s1 + s2) / M
# beta_int <- chi_hat / (t1_h1 + t2_h1)
beta_int <- chi_hat / (map_dbl(rc1$Xk, sum) + map_dbl(rc2$Xk, sum))

scale_factor <- lambda1 / lambda2
# chi_hat for control and test sets
chi_hat1 <- rc1$chi
chi_hat2 <- rc2$chi
# max iterations and tolerance for EM
max_itr <- 500
tor <- 1e-6
# run EM for multiple combinations of parameters
em_res <- pmap(
  list(Xk1, Xk2, beta_int, chi_hat, chi_hat1, chi_hat2),
  function(x, y, z, k, m, n) {
    tryCatch(main_EM_h0(fk_int,
      Xk1 = x, Xk2 = y, kmin, kmax, beta_int = z,
      chi_hat = k, chi_hat1 = m, chi_hat2 = n,
      max_itr, tor
    ),
    error = function(err) {
      # handling the error, one of the cases is when
      # there is no read counts in the pause region
      list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
    }
    )
  }
)

h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])

beta_tbl <-
  tibble(
    gene_id = rc1$gene_id,
    beta1 = rc1$beta,
    beta2 = rc2$beta,
    lfc = log2(beta2 / beta1),
    fk_mean1 = rc1$fk_mean,
    fk_mean2 = rc2$fk_mean,
    fk_var1 = rc1$fk_var,
    fk_var2 = rc2$fk_var,
    # use eq (25) instead of (31) to compute T stats
    t_stats_beta = rc1$likelihood + rc2$likelihood - h0_likelihood
  )

# some genes with negative T stats, fix them
idx <- beta_tbl$t_stats_beta < 0

# use parameter estimates from h0 as initial values for EM in h1
h0_beta <- map_dbl(em_res, "beta")
h0_fk1 <- map(em_res, "fk1")
h0_fk2 <- map(em_res, "fk2")

# condition 1
em_hc <- pmap(
  list(h0_fk1[idx], Xk1[idx], h0_beta[idx], chi_hat1[idx]),
  function(x, y, z, k) {
    tryCatch(main_EM(
      fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
      beta_int = z, chi_hat = k,
      max_itr = max_itr, tor = tor
    ),
    error = function(err) {
      list(
        "beta" = NA, "Yk" = NA,
        "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA
      )
    }
    )
  }
)

# condition 2
em_ht <- pmap(
  list(h0_fk2[idx], Xk2[idx], h0_beta[idx], chi_hat2[idx]),
  function(x, y, z, k) {
    tryCatch(main_EM(
      fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
      beta_int = z, chi_hat = k,
      max_itr = max_itr, tor = tor
    ),
    error = function(err) {
      list(
        "beta" = NA, "Yk" = NA,
        "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA
      )
    }
    )
  }
)

h1_likelihood1 <- map_dbl(em_hc, ~ .x$likelihoods[[length(.x$likelihoods)]])
h1_likelihood2 <- map_dbl(em_ht, ~ .x$likelihoods[[length(.x$likelihoods)]])

beta_tbl_idx <-
  tibble(
    gene_id = names(em_hc),
    beta1 = map_dbl(em_hc, "beta"),
    beta2 = map_dbl(em_ht, "beta"),
    lfc = log2(beta2 / beta1),
    fk_mean1 = map_dbl(em_hc, "fk_mean"),
    fk_mean2 = map_dbl(em_ht, "fk_mean"),
    fk_var1 = map_dbl(em_hc, "fk_var"),
    fk_var2 = map_dbl(em_ht, "fk_var"),
    t_stats_beta = h1_likelihood1 + h1_likelihood2 - h0_likelihood[idx]
  )

beta_tbl <- bind_rows(beta_tbl[!idx, ], beta_tbl_idx)

beta_tbl <- beta_tbl %>%
  mutate(p_beta = pchisq(2 * t_stats_beta, df = 1, ncp = 0, lower.tail = F, log.p = FALSE))

beta_tbl <- beta_tbl %>% mutate(padj_beta = p.adjust(p_beta, method = "BH"))

beta_tbl <- beta_tbl %>%
  mutate(
    beta = (beta1 + beta2) / 2,
    logbeta_zeta = log2(beta * zeta),
    category =
      case_when(
        (padj_beta < 0.05) & (lfc > lfc1) ~ "Up",
        (padj_beta < 0.05) & (lfc < lfc2) ~ "Down",
        TRUE ~ "Others"
      ),
    category = factor(category, levels = c("Up", "Down", "Others"))
  )

beta_tbl <- beta_tbl %>%
  mutate(
    fk_std1 = fk_var1^0.5,
    fk_std2 = fk_var2^0.5
  )

## LRT for beta distribution##
beta_int1 <- rc1$beta
beta_int2 <- rc2$beta

# run EM for multiple combinations of parameters
em_res <- pmap(
  list(Xk1, Xk2, beta_int1, beta_int2, chi_hat, chi_hat1, chi_hat2),
  function(x, y, z1, z2, k, m, n) {
    tryCatch(main_EM_fk_h0(fk_int,
      Xk1 = x, Xk2 = y, kmin, kmax, beta_int1 = z1, beta_int2 = z2,
      chi_hat = k, chi_hat1 = m, chi_hat2 = n, scale_factor = scale_factor,
      max_itr, tor
    ),
    error = function(err) {
      # handling the error, one of the cases is when
      # there is no read counts in the pause region
      list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
    }
    )
  }
)

h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])
h0_fk <- map(em_res, ~ .x$fk)
h0_fk_mean <- map_dbl(em_res, ~ .x$fk_mean[[length(.x$fk_mean)]])
h0_fk_var <- map_dbl(em_res, ~ .x$fk_var[[length(.x$fk_var)]])

fk_tbl <-
  tibble(
    gene_id = rc1$gene_id,
    beta1 = rc1$beta,
    beta2 = rc2$beta,
    lfc = log2(beta2 / beta1),
    fk_mean1 = rc1$fk_mean,
    fk_mean2 = rc2$fk_mean,
    h0_fk = h0_fk,
    h0_fk_mean = h0_fk_mean,
    h0_fk_var = h0_fk_var,
    fk_var1 = rc1$fk_var,
    fk_var2 = rc2$fk_var,
    # use eq (25) instead of (31) to compute T stats
    t_stats_fk = rc1$likelihood + (rc2$likelihood * scale_factor) - h0_likelihood
  )

# some genes with negative T stats, fix them
idx <- fk_tbl$t_stats_fk < 0

# use parameter estimates from h0 as initial values for EM in h1
h0_fk <- map(em_res, "fk")

em_hc <- pmap(
  list(h0_fk[idx], Xk1[idx], beta_int1[idx], chi_hat1[idx]),
  function(x, y, z, k) {
    tryCatch(main_EM(
      fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
      beta_int = z, chi_hat = k,
      max_itr = max_itr, tor = tor
    ),
    error = function(err) {
      list(
        "beta" = NA, "Yk" = NA,
        "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA
      )
    }
    )
  }
)

em_ht <- pmap(
  list(h0_fk[idx], Xk2[idx], beta_int2[idx], chi_hat2[idx]),
  function(x, y, z, k) {
    tryCatch(main_EM(
      fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
      beta_int = z, chi_hat = k,
      max_itr = max_itr, tor = tor
    ),
    error = function(err) {
      list(
        "beta" = NA, "Yk" = NA,
        "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA
      )
    }
    )
  }
)

h1_likelihood1 <- map_dbl(em_hc, ~ .x$likelihoods[[length(.x$likelihoods)]])
h1_likelihood2 <- map_dbl(em_ht, ~ .x$likelihoods[[length(.x$likelihoods)]])

fk_tbl_idx <-
  tibble(
    gene_id = names(em_hc),
    beta1 = map_dbl(em_hc, "beta"),
    beta2 = map_dbl(em_ht, "beta"),
    lfc = log2(beta2 / beta1),
    fk_mean1 = map_dbl(em_hc, "fk_mean"),
    fk_mean2 = map_dbl(em_ht, "fk_mean"),
    h0_fk = map(em_res, "fk")[idx],
    h0_fk_mean = map_dbl(em_res, "fk_mean")[idx],
    h0_fk_var = map_dbl(em_res, "fk_var")[idx],
    fk_var1 = map_dbl(em_hc, "fk_var"),
    fk_var2 = map_dbl(em_ht, "fk_var"),
    t_stats_fk = h1_likelihood1 + h1_likelihood2 - h0_likelihood[idx]
  )

fk_tbl <- bind_rows(fk_tbl[!idx, ], fk_tbl_idx)

fk_tbl <- fk_tbl %>%
  mutate(p_fk = pchisq(2 * t_stats_fk, df = 2, ncp = 0, lower.tail = F, log.p = FALSE))

fk_tbl <- fk_tbl %>% mutate(padj_fk = p.adjust(p_fk, method = "BH"))

fk_tbl$delta_mean <- fk_tbl$fk_mean2 - fk_tbl$fk_mean1
fk_tbl$delta_sd <-
  (sqrt(fk_tbl$fk_var2) - sqrt(fk_tbl$fk_var1)) / (sqrt(fk_tbl$fk_var2) + sqrt(fk_tbl$fk_var1))

fk_tbl <- fk_tbl %>%
  mutate(
    significant =
      case_when(
        (padj_fk < 0.05) & ((abs(delta_sd) > 0.2) | (abs(delta_mean) > 20)) ~ "Yes",
        TRUE ~ "No"
      ),
    significant = factor(significant, levels = c("Yes", "No"))
  )

#### visualize results ####
p <- chi_tbl %>%
  select(chi1, chi2) %>%
  pivot_longer(cols = contains("chi")) %>%
  mutate(value = log2(value)) %>%
  ggboxplot(
    x = "name", y = "value", fill = "name",
    palette = paired_colors,
    outlier.shape = NA, notch = TRUE
  ) +
  stat_compare_means() +
  scale_x_discrete(labels = c("chi1" = "Control", "chi2" = "Treated")) +
  labs(x = "", y = expression(log[2] * chi)) +
  coord_cartesian(ylim = c(scale_tbl$chi_ymin, scale_tbl$chi_ymax)) +
  theme(legend.position = "none")

ggsave(file.path(result_dir, "chi_distribution.pdf"),
  plot = p,
  width = 4, height = 5
)

p <- chi_tbl %>%
  ggplot(aes(x = logchi, y = lfc, color = category)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = scatter_colors[c(3, 1, 2)]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ylim(-6, 6) +
  labs(
    y = expression(log[2] * "FC(Treated/Control)"),
    x = expression(log[2] * "mean(" * chi * ")"),
    color = "Category"
  ) +
  cowplot::theme_cowplot()

ggsave(file.path(result_dir, "chi_mean_vs_lfc.pdf"),
  plot = p,
  width = 5, height = 3
)

p <- beta_tbl %>%
  select(beta1, beta2) %>%
  pivot_longer(cols = contains("beta")) %>%
  mutate(value = log2(value * zeta)) %>%
  ggboxplot(
    x = "name", y = "value", fill = "name",
    palette = paired_colors[c(3, 4)],
    outlier.shape = NA, notch = TRUE
  ) +
  stat_compare_means(label.x.npc = "left", label.y = scale_tbl$beta_ymax) +
  scale_x_discrete(labels = c("beta1" = "Control", "beta2" = "Treated")) +
  labs(x = "", y = expression(log[2] * beta * zeta)) +
  coord_cartesian(ylim = c(scale_tbl$beta_ymin, scale_tbl$beta_ymax)) +
  theme(legend.position = "none")

ggsave(file.path(result_dir, "beta_distribution.pdf"),
  plot = p,
  width = 4, height = 5
)

p <- beta_tbl %>%
  ggplot(aes(x = logbeta_zeta, y = lfc, color = category)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = scatter_colors[c(3, 1, 2)]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ylim(-6, 6) +
  xlim(-10, 10) +
  labs(
    y = expression(log[2] * "FC(Treated/Control)"),
    x = expression(log[2] * "mean(" * beta * zeta * ")"),
    color = "Category"
  ) +
  cowplot::theme_cowplot()

ggsave(file.path(result_dir, "beta_mean_vs_lfc.pdf"),
  plot = p,
  width = 5, height = 3
)

# Changes of initiation and pause release rates
lfc_tbl <- chi_tbl %>%
  select(gene_id, lfc, category) %>%
  inner_join(beta_tbl %>% select(gene_id, lfc, category),
    by = "gene_id", suffix = c("_chi", "_beta")
  )

p <- lfc_tbl %>%
  ggscatter(
    x = "lfc_chi", y = "lfc_beta",
    add = "reg.line", alpha = 0.2,
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE, # Add confidence interval
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = "pearson", label.sep = "\n",
      label.x.npc = "left", label.y.npc = "top"
    )
  ) +
  ylim(-8, 4) +
  xlim(-6, 6)

ggsave(file.path(result_dir, "lfc_chi_vs_beta.png"),
  plot = p,
  width = 6, height = 5
)

# Number of genes with differential rates
lfc_summary <- lfc_tbl %>%
  select(contains("category")) %>%
  pivot_longer(everything()) %>%
  group_by(name, value) %>%
  summarise(count = n()) %>%
  mutate(name = str_remove(name, "category_"))

p <- lfc_summary %>%
  ggplot(aes(x = name, y = count, fill = value)) +
  geom_col(position = "dodge") +
  scale_x_discrete(labels = c(
    "beta" = expression(beta),
    "chi" = expression(chi)
  )) +
  scale_fill_manual(values = scatter_colors[c(3, 1, 2)]) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.25) +
  labs(x = "", y = "Number of Genes", fill = "Category")

ggsave(file.path(result_dir, "gene_number_with_differential_rates.pdf"),
  plot = p,
  width = 6, height = 4
)

# mean and variance of pause sites in different beta categories
p <- beta_tbl %>%
  ggplot() +
  geom_density_2d(aes(x = fk_mean1, y = fk_std1), color = paired_colors[5]) +
  geom_density_2d(aes(x = fk_mean2, y = fk_std2), color = paired_colors[6]) +
  facet_wrap(. ~ category) +
  labs(x = expression(mu), y = expression(sigma))

ggsave(file.path(result_dir, "mean_vs_std_of_k.png"),
  plot = p,
  width = 9, height = 4
)

# plot mean and std of pause sites side by side
p1 <- beta_tbl %>%
  ggplot(aes(x = fk_mean1, y = fk_std1)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_vline(xintercept = c(50, 100), linetype = "dashed", color = "gray") +
  coord_cartesian(xlim = c(0, 200), ylim = c(0, 70)) +
  labs(x = "Mean of k", y = "SD of k")

p2 <- beta_tbl %>%
  ggplot(aes(x = fk_mean2, y = fk_std2)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_vline(xintercept = c(50, 100), linetype = "dashed", color = "gray") +
  coord_cartesian(xlim = c(0, 200), ylim = c(0, 70)) +
  labs(x = "Mean of k", y = "SD of k")

p <- cowplot::plot_grid(p1, p2)

ggsave(
  filename = file.path(result_dir, "mean_vs_std_of_k_pointdensity.pdf"), plot = p,
  width = 14, height = 5
)
# plot changes in mean and std of pause sites
p <- fk_tbl %>%
  ggplot(aes(x = delta_mean, y = delta_sd, color = significant)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = scatter_colors[c(3, 2)]) +
  geom_vline(xintercept = -20, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  labs(
    x = "Differences in mean of k",
    y = "Normalized difference in sd of k",
    color = "Significant"
  ) +
  cowplot::theme_cowplot() +
  theme(
    legend.position = c(0.75, 0.15)
  )

p <- ggMarginal(p, type = "histogram", bins = 50, fill = "grey60", color = "black")

ggsave(
  filename = file.path(result_dir, "lrt_mean_vs_std_of_k_pointdensity.pdf"), plot = p,
  width = 6, height = 5
)

## write out results ##
write_csv(chi_tbl, file = chi_out)
write_csv(beta_tbl, file = beta_out)
fk_tbl <- fk_tbl %>% select(-h0_fk)
write_csv(fk_tbl, file = fk_out)