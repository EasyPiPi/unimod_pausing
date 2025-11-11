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

omega_out <- snakemake@output[["omega"]]
beta_out <- snakemake@output[["beta"]]
#### end of parsing arguments ####

#### load packages ####
library(tidyverse)
library(ggpubr)
library(ggpointdensity)
library(viridis)

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

# omega_out <- file.path(result_dir, "omega.csv")
# beta_out <- file.path(result_dir, "beta.csv")

# # #### end of parsing arguments ####
theme_set(cowplot::theme_cowplot())

# set up parameters
quantile_normalization <- "identity"

k <- 50
kmin <- 1
kmax <- 200 # also used as k on the poisson case

rnap_size <- 50
zeta <- 2000

# criteria for significance, p.adj < 0.05 and 2 fold differences
sig_p <- 0.05
lfc1 <- 0
lfc2 <- 0

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
source(helper_tc_in)
source(helper_pr_in)

# read in summary statistics from the varied pause site model
rc1 <- readRDS(rc1_in)
rc2 <- readRDS(rc2_in)

# Number of genes in each sample
NROW(rc1)
NROW(rc2)

# union set of genes being analyzed
gn_union <- intersect(rc1$gene_id, rc2$gene_id)
rc1 <- rc1[match(gn_union, rc1$gene_id), ]
rc2 <- rc2[match(gn_union, rc2$gene_id), ]

#### Poisson-based Likelihood Ratio Tests ####
# read in number of spike-in or total number of mappable reads
# use them as scaling factor
scale_tbl <- read_csv(spike_in, show_col_types = FALSE)
# subset the right table
scale_tbl <- scale_tbl[str_detect(rc2_in, scale_tbl$sample), ]
#
# based on formula (27) and (28), cancel out M and zeta since they are the same
lambda1 <- scale_tbl$control_1 + ifelse(is.na(scale_tbl$control_2), 0, scale_tbl$control_2)
lambda2 <- scale_tbl$treated_1 + ifelse(is.na(scale_tbl$treated_2), 0, scale_tbl$treated_2)

## LRT for omega ##
tao1 <- lambda1 / (lambda1 + lambda2)
tao2 <- 1 - tao1

omega_tbl <-
  tibble(
    gene_id = rc1$gene_id,
    chi1 = rc1$chi,
    chi2 = rc2$chi * lambda1 / lambda2,
    lfc = log2(chi2 / chi1)
  )

omega_tbl <- omega_tbl %>%
  bind_cols(bind_rows(map2(rc1$s, rc2$s, omega_lrt, tao1 = tao1, tao2 = tao2)))

omega_tbl <- omega_tbl %>%
  mutate(padj = p.adjust(p, method = "BH"))

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
    tryCatch(
      main_EM_h0(fk_int,
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
# calculate t1 and t2 for H0
# t1_h0 <- map_dbl(em_res, ~ sum(.x[["Yk1"]]), .default = NA)
# t2_h0 <- map_dbl(em_res, ~ sum(.x[["Yk2"]]), .default = NA)

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
    t_stats = rc1$likelihood + rc2$likelihood - h0_likelihood
  )

# some genes with negative T stats, fix them
idx <- beta_tbl$t_stats < 0

# use parameter estimates from h0 as initial values for EM in h1
h0_beta <- map_dbl(em_res, "beta")
h0_fk1 <- map(em_res, "fk1")
h0_fk2 <- map(em_res, "fk2")

em_hc <- pmap(
  list(h0_fk1[idx], Xk1[idx], h0_beta[idx], chi_hat1[idx]),
  function(x, y, z, k) {
    tryCatch(
      main_EM(
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
  list(h0_fk2[idx], Xk2[idx], h0_beta[idx], chi_hat2[idx]),
  function(x, y, z, k) {
    tryCatch(
      main_EM(
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
    t_stats = h1_likelihood1 + h1_likelihood2 - h0_likelihood[idx]
  )

beta_tbl <- bind_rows(beta_tbl[!idx, ], beta_tbl_idx)

beta_tbl <- beta_tbl %>%
  mutate(p = pchisq(2 * t_stats, df = 1, ncp = 0, lower.tail = F, log.p = FALSE))

beta_tbl <- beta_tbl %>% mutate(padj = p.adjust(p, method = "BH"))

# troubleshoot for genes with negative T stats
idx1 <- which(beta_tbl$t_stats < 0)
idx2 <- which(beta_tbl$t_stats > 0)

if (length(idx1) > 0) {
  compare_pause_sites <- function(idx, col_name) {
    fk_tbl <- bind_cols(
      abs(rc1[idx, "fk_mean"] - rc2[idx, "fk_mean"]),
      abs(rc1[idx, "fk_var"]^0.5 - rc2[idx, "fk_var"]^0.5)
    )
    fk_tbl[["category"]] <- col_name

    return(fk_tbl)
  }

  fk_comparison <- bind_rows(
    compare_pause_sites(idx1, "negative"),
    compare_pause_sites(idx2, "positive")
  )

  p <- fk_comparison %>%
    ggplot(aes(x = fk_mean, color = category)) +
    geom_density() +
    xlab(expression("|" * mu[1] - mu[2] * "|"))

  ggsave(file.path(result_dir, "fk_mean_differences.png"),
    plot = p,
    width = 6, height = 3
  )

  p <- fk_comparison %>%
    ggplot(aes(x = fk_var, color = category)) +
    geom_density() +
    xlab(expression("|" * sigma[1] - sigma[2] * "|"))

  ggsave(file.path(result_dir, "fk_var_differences.png"),
    plot = p,
    width = 6, height = 3
  )
}

#### visualize results ####
# violion plot
# p <- omega_tbl %>%
#   select(contains("chi")) %>%
#   pivot_longer(cols = contains("chi")) %>%
#   mutate(value = log2(value)) %>%
#    ggviolin(x = "name", y = "value", fill = "name",
#            palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#            add = "median_q1q3", add.params = list(fill = "white")) +
#   stat_compare_means(label.x = 1, label.y.npc = "top") +
#   scale_x_discrete(labels=c("chi1" = "Control", "chi2" = "Treated")) +
#   labs(x = "", y = expression(log[2]*chi)) +
#   coord_cartesian(ylim = c(-10, 2)) +
#   theme(legend.position = "none")

p <- omega_tbl %>%
  select(contains("chi")) %>%
  pivot_longer(cols = contains("chi")) %>%
  mutate(value = log2(value)) %>%
  ggboxplot(
    x = "name", y = "value", fill = "name",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    outlier.shape = NA, notch = TRUE
  ) +
  stat_compare_means() +
  scale_x_discrete(labels = c("chi1" = "Control", "chi2" = "Treated")) +
  labs(x = "", y = expression(log[2] * chi)) +
  coord_cartesian(ylim = c(scale_tbl$chi_ymin, scale_tbl$chi_ymax)) +
  theme(legend.position = "none")

ggsave(file.path(result_dir, "chi_distribution.png"),
  plot = p,
  width = 4, height = 5
)

# p <- alpha_tbl %>%
#   select(contains("alpha_zeta")) %>%
#   pivot_longer(cols = contains("alpha_zeta")) %>%
#   mutate(value = log2(value)) %>%
#   ggboxplot(x = "name", y = "value", fill = "name",
#             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#             outlier.shape = NA, notch = TRUE) +
#   stat_compare_means(label.x.npc = "left", label.y = scale_tbl$alpha_ymax) +
#   scale_x_discrete(labels=c("alpha_zeta1" = "Control", "alpha_zeta2" = "Treated")) +
#   labs(x = "", y = expression(log[2]*alpha*zeta)) +
#   coord_cartesian(ylim = c(scale_tbl$alpha_ymin, scale_tbl$alpha_ymax)) +
#   theme(legend.position = "none")
#
# ggsave(file.path(result_dir, "alpha_distribution.png"), plot = p,
#        width = 4, height = 5)

# p <- beta_tbl %>%
#   select(contains("beta")) %>%
#   pivot_longer(cols = contains("beta")) %>%
#   mutate(value = log2(value * zeta)) %>%
#   ggviolin(x = "name", y = "value", fill = "name",
#            palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#            add = "median_q1q3", add.params = list(fill = "white")) +
#   stat_compare_means(label.x = 1.4, label.y.npc = "top") +
#   scale_x_discrete(labels=c("beta1" = "Control", "beta2" = "Treated")) +
#   labs(x = "", y = expression(log[2]*beta*zeta)) +
#   theme(legend.position = "none")

p <- beta_tbl %>%
  select(contains("beta")) %>%
  pivot_longer(cols = contains("beta")) %>%
  mutate(value = log2(value * zeta)) %>%
  ggboxplot(
    x = "name", y = "value", fill = "name",
    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
    outlier.shape = NA, notch = TRUE
  ) +
  stat_compare_means(label.x.npc = "left", label.y = scale_tbl$beta_ymax) +
  scale_x_discrete(labels = c("beta1" = "Control", "beta2" = "Treated")) +
  labs(x = "", y = expression(log[2] * beta * zeta)) +
  coord_cartesian(ylim = c(scale_tbl$beta_ymin, scale_tbl$beta_ymax)) +
  theme(legend.position = "none")

ggsave(file.path(result_dir, "beta_distribution.png"),
  plot = p,
  width = 4, height = 5
)

# volcano plot
#### functions for visualizatiob ####
volcano_plot <- function(df, sig_p) {
  df$significant <- df$padj < sig_p
  df %>%
    na.omit() %>%
    ggplot(aes(x = lfc, y = -log10(padj), color = significant)) +
    geom_point() +
    scale_color_manual(values = c("grey", "red")) +
    labs(
      x = expression(log[2] * "FC(Treated/Control)"),
      y = expression(-log[10] * "(padj)")
    ) +
    cowplot::theme_cowplot()
}

p <- volcano_plot(omega_tbl, sig_p = sig_p)
ggsave(file.path(result_dir, "chi_lrt_volcano.png"),
  plot = p,
  width = 8, height = 6
)

p <- volcano_plot(beta_tbl, sig_p = sig_p)
ggsave(file.path(result_dir, "beta_lrt_volcano.png"),
  plot = p,
  width = 8, height = 6
)

# mean vs. LFC
omega_tbl <- omega_tbl %>%
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

p <- omega_tbl %>%
  ggplot(aes(x = logchi, y = lfc, color = category)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "gray")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ylim(-6, 6) +
  labs(
    y = expression(log[2] * "FC(Treated/Control)"),
    x = expression(log[2] * "mean(" * chi * ")"),
    color = "Category"
  ) +
  cowplot::theme_cowplot()

ggsave(file.path(result_dir, "chi_mean_vs_lfc.png"),
  plot = p,
  width = 6, height = 3
)

beta_tbl <- beta_tbl %>%
  mutate(
    beta = (beta1 + beta2) / 2,
    logbeta_zeta = log2(beta * zeta),
    category =
      case_when(
        (padj < 0.05) & (lfc > lfc1) ~ "Up",
        (padj < 0.05) & (lfc < lfc2) ~ "Down",
        TRUE ~ "Others"
      ),
    category = factor(category, levels = c("Up", "Down", "Others"))
  )

p <- beta_tbl %>%
  ggplot(aes(x = logbeta_zeta, y = lfc, color = category)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "gray")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ylim(-6, 6) +
  xlim(-10, 10) +
  labs(
    y = expression(log[2] * "FC(Treated/Control)"),
    x = expression(log[2] * "mean(" * beta * zeta * ")"),
    color = "Category"
  ) +
  cowplot::theme_cowplot()

ggsave(file.path(result_dir, "beta_mean_vs_lfc.png"),
  plot = p,
  width = 6, height = 3
)

# Changes of initiation and pause release rates
lfc_tbl <- omega_tbl %>%
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
    "alpha" = expression(alpha),
    "beta" = expression(beta),
    "chi" = expression(chi)
  )) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.25) +
  labs(x = "", y = "Number of Genes", fill = "Category")

ggsave(file.path(result_dir, "gene_number_with_differential_rates.png"),
  plot = p,
  width = 6, height = 4
)

# Number of genes analyzed
NROW(beta_tbl)
# Number of genes showing significant downregulation in β
sum(beta_tbl$padj < 0.05 & beta_tbl$lfc < 0)
# Number of genes showing significant upregulation in β
sum(beta_tbl$padj < 0.05 & beta_tbl$lfc > 0)

# mean and variance of pause sites in different beta categories
beta_tbl <- beta_tbl %>%
  mutate(
    fk_std1 = fk_var1^0.5,
    fk_std2 = fk_var2^0.5
  )

p <- beta_tbl %>%
  ggplot() +
  geom_density_2d(aes(x = fk_mean1, y = fk_std1), color = "blue") +
  geom_density_2d(aes(x = fk_mean2, y = fk_std2), color = "red") +
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
  filename = file.path(result_dir, "mean_vs_std_of_k_pointdensity.png"), plot = p,
  width = 14, height = 5
)

# output tables
write_csv(omega_tbl, file = omega_out)
write_csv(beta_tbl, file = beta_out)
# write_csv(alpha_tbl, file = alpha_out)
write_csv(lfc_tbl, file = file.path(result_dir, "lfc.csv"))
write_csv(lfc_summary %>% pivot_wider(names_from = name, values_from = count),
  file = file.path(result_dir, "lfc_summary.csv")
)
