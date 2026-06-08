# ============================================================
# Perturbation analysis for DLD1 cells (paired-condition LRT)
#
# Compares pause-release rates between a control and a
# perturbation condition using:
#   - beta LRT  (pause index change)
#   - fk LRT    (pause-site position distribution change)
#
# Set CONDITION in the INPUT CONFIGURATION section below,
# then source this script.
# ============================================================


# ============================================================
# LOAD PACKAGES
# ============================================================

library(tidyverse)
library(ggpubr)
library(ggpointdensity)
library(viridis)
library(ggExtra)


# ============================================================
# LOAD CONFIGURATION
# ============================================================

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))


# ============================================================
# SOURCE HELPER SCRIPTS
# ============================================================

source(file.path(root_dir, "codes", "publish", "perturbation",
                 "helper_function_em_two_condition.R"))
source(file.path(root_dir, "codes", "publish", "perturbation",
                 "helper_function_em_pause_release.R"))


# ============================================================
# INPUT CONFIGURATION
# Set CONDITION to one of: "NVP2", "Fp", "Auxin"
# ============================================================

CONDITION <- "Fp"

.data_root <- file.path(.paths$data, "perturbation")

condition_configs <- list(

  NVP2 = list(
    rc1_in    = file.path(.data_root, "PROseq-DLD1-aoi-NELFC_NVP2_Ctrl-SE",
                          "pause_release", "rate.RDS"),
    rc2_in    = file.path(.data_root, "PROseq-DLD1-aoi-NELFC_NVP2-SE",
                          "pause_release", "rate.RDS"),
    scale_key = "NVP2",
    out_subdir = "NVP2"
  ),

  Fp = list(
    rc1_in    = file.path(.data_root, "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE",
                          "pause_release", "rate.RDS"),
    rc2_in    = file.path(.data_root, "PROseq-DLD1-aoi-NELFC_Fp-SE",
                          "pause_release", "rate.RDS"),
    scale_key = "Fp",
    out_subdir = "Fp"
  ),

  Auxin = list(
    rc1_in    = file.path(.data_root, "PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE",
                          "pause_release", "rate.RDS"),
    rc2_in    = file.path(.data_root, "PROseq-DLD1-aoi-NELFC_Auxin-SE",
                          "pause_release", "rate.RDS"),
    scale_key = "Auxin",
    out_subdir = "Auxin"
  )

)

cfg <- condition_configs[[CONDITION]]

result_dir <- file.path(.paths$outputs, "publish", "perturbation", cfg$out_subdir)
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)


# ============================================================
# READ DATA
# Reads EM summary statistics from the pause-release model.
# Applies a max-count filter (> 8) and intersects genes
# present in both conditions before analysis.
# ============================================================

rc1 <- readRDS(cfg$rc1_in)
rc2 <- readRDS(cfg$rc2_in)

message("Genes in control: ", NROW(rc1))
message("Genes in treated: ", NROW(rc2))

rc1$max_count <- map_dbl(rc1$Xk, max)
rc1 <- rc1[rc1$max_count > 8, ]

rc2$max_count <- map_dbl(rc2$Xk, max)
rc2 <- rc2[rc2$max_count > 8, ]

gn_union <- intersect(rc1$gene_id, rc2$gene_id)
rc1 <- rc1[match(gn_union, rc1$gene_id), ]
rc2 <- rc2[match(gn_union, rc2$gene_id), ]

message("Genes after intersection and filtering: ", length(gn_union))

# Spike-in normalization: compute scale factor from spike-in read counts
scale_tbl  <- read_csv(file.path(.data_root, "scaling_factor.csv"),
                       show_col_types = FALSE)
scale_row  <- scale_tbl[scale_tbl$sample == cfg$scale_key, ]

lambda1      <- scale_row$control_1 + ifelse(is.na(scale_row$control_2), 0, scale_row$control_2)
lambda2      <- scale_row$treated_1 + ifelse(is.na(scale_row$treated_2), 0, scale_row$treated_2)
scale_factor <- lambda1 / lambda2


# ============================================================
# MODEL PARAMETERS
# ============================================================

kmin <- 1
kmax <- 200   # also used as k in the Poisson case
zeta <- 2000

# significance thresholds: adjusted p < 0.05 and ±20% fold change
sig_p     <- 0.05
beta_lfc1 <- log2(1.2 / 1)
beta_lfc2 <- log2(0.8 / 1)

max_itr <- 500
tor     <- 1e-6

# Initial fk distribution (truncated normal; uniform gives similar results)
fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
fk_int <- fk_int / sum(fk_int)

# Convenience extractions shared by both LRTs
s1      <- rc1$s
s2      <- rc2$s
Xk1     <- rc1$Xk
Xk2     <- rc2$Xk
M       <- rc1$N          # gene body length, assumed same across conditions
chi_hat  <- (s1 + s2) / M
chi_hat1 <- rc1$chi
chi_hat2 <- rc2$chi


# ============================================================
# BETA LRT
# Tests whether the pause index (beta) differs between
# conditions. H0: beta is shared; H1: beta is condition-specific.
# T statistic: L(H1_ctrl) + L(H1_treat) - L(H0).
# Degrees of freedom: 1.
# ============================================================

beta_int <- chi_hat / (map_dbl(rc1$Xk, sum) + map_dbl(rc2$Xk, sum))

em_res <- pmap(
  list(Xk1, Xk2, beta_int, chi_hat, chi_hat1, chi_hat2),
  function(x, y, z, k, m, n) {
    tryCatch(
      main_EM_h0(fk_int, Xk1 = x, Xk2 = y, kmin, kmax, beta_int = z,
                 chi_hat = k, chi_hat1 = m, chi_hat2 = n,
                 scale_factor, max_itr, tor),
      error = function(err) list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
    )
  }
)

h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])

beta_tbl <- tibble(
  gene_id  = rc1$gene_id,
  beta1    = rc1$beta,
  beta2    = rc2$beta,
  lfc      = log2(beta2 / beta1),
  fk_mean1 = rc1$fk_mean,
  fk_mean2 = rc2$fk_mean,
  fk_var1  = rc1$fk_var,
  fk_var2  = rc2$fk_var,
  t_stats  = rc1$likelihood + rc2$likelihood - h0_likelihood
)

# Refit genes with negative T stats using H0 parameter estimates as H1 starting values
idx <- beta_tbl$t_stats < 0

h0_beta <- map_dbl(em_res, "beta")
h0_fk1  <- map(em_res, "fk1")
h0_fk2  <- map(em_res, "fk2")

em_hc <- pmap(
  list(h0_fk1[idx], Xk1[idx], h0_beta[idx], chi_hat1[idx]),
  function(x, y, z, k) {
    tryCatch(
      main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
              beta_int = z, chi_hat = k, max_itr = max_itr, tor = tor),
      error = function(err)
        list("beta" = NA, "Yk" = NA, "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
    )
  }
)

em_ht <- pmap(
  list(h0_fk2[idx], Xk2[idx], h0_beta[idx], chi_hat2[idx]),
  function(x, y, z, k) {
    tryCatch(
      main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
              beta_int = z, chi_hat = k, max_itr = max_itr, tor = tor),
      error = function(err)
        list("beta" = NA, "Yk" = NA, "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
    )
  }
)

h1_likelihood1 <- map_dbl(em_hc, ~ .x$likelihoods[[length(.x$likelihoods)]])
h1_likelihood2 <- map_dbl(em_ht, ~ .x$likelihoods[[length(.x$likelihoods)]])

beta_tbl_idx <- tibble(
  gene_id  = names(em_hc),
  beta1    = map_dbl(em_hc, "beta"),
  beta2    = map_dbl(em_ht, "beta"),
  lfc      = log2(beta2 / beta1),
  fk_mean1 = map_dbl(em_hc, "fk_mean"),
  fk_mean2 = map_dbl(em_ht, "fk_mean"),
  fk_var1  = map_dbl(em_hc, "fk_var"),
  fk_var2  = map_dbl(em_ht, "fk_var"),
  t_stats  = h1_likelihood1 + h1_likelihood2 - h0_likelihood[idx]
)

beta_tbl <- bind_rows(beta_tbl[!idx, ], beta_tbl_idx)

beta_tbl <- beta_tbl %>%
  mutate(
    p    = pchisq(2 * t_stats, df = 1, ncp = 0, lower.tail = FALSE),
    padj = p.adjust(p, method = "BH"),
    beta = (beta1 + beta2) / 2,
    logbeta_zeta = log2(beta * zeta),
    category = case_when(
      (padj < sig_p) & (lfc > beta_lfc1) ~ "Up",
      (padj < sig_p) & (lfc < beta_lfc2) ~ "Down",
      TRUE ~ "Others"
    ),
    category  = factor(category, levels = c("Up", "Others", "Down")),
    padj_plot = pmax(padj, 1e-300)
  )


# ============================================================
# FK LRT
# Tests whether the pause-site position distribution (fk)
# differs between conditions.
# H0: fk is shared; H1: fk is condition-specific.
# Degrees of freedom: 2.
# ============================================================

beta_int1 <- rc1$beta
beta_int2 <- rc2$beta

em_res <- pmap(
  list(Xk1, Xk2, beta_int1, beta_int2, chi_hat, chi_hat1, chi_hat2),
  function(x, y, z1, z2, k, m, n) {
    tryCatch(
      main_EM_fk_h0(fk_int, Xk1 = x, Xk2 = y, kmin, kmax,
                    beta_int1 = z1, beta_int2 = z2,
                    chi_hat = k, chi_hat1 = m, chi_hat2 = n,
                    scale_factor = scale_factor, max_itr, tor),
      error = function(err) list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
    )
  }
)

h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])
h0_fk         <- map(em_res, ~ .x$fk)
h0_fk_mean    <- map_dbl(em_res, ~ .x$fk_mean[[length(.x$fk_mean)]])
h0_fk_var     <- map_dbl(em_res, ~ .x$fk_var[[length(.x$fk_var)]])

fk_tbl <- tibble(
  gene_id    = rc1$gene_id,
  beta1      = rc1$beta,
  beta2      = rc2$beta,
  lfc        = log2(beta2 / beta1),
  fk_mean1   = rc1$fk_mean,
  fk_mean2   = rc2$fk_mean,
  h0_fk      = h0_fk,
  h0_fk_mean = h0_fk_mean,
  h0_fk_var  = h0_fk_var,
  fk_var1    = rc1$fk_var,
  fk_var2    = rc2$fk_var,
  t_stats    = rc1$likelihood + (rc2$likelihood * scale_factor) - h0_likelihood
)

# Refit genes with negative T stats
idx   <- fk_tbl$t_stats < 0
h0_fk <- map(em_res, "fk")

em_hc <- pmap(
  list(h0_fk[idx], Xk1[idx], beta_int1[idx], chi_hat1[idx]),
  function(x, y, z, k) {
    tryCatch(
      main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
              beta_int = z, chi_hat = k, max_itr = max_itr, tor = tor),
      error = function(err)
        list("beta" = NA, "Yk" = NA, "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
    )
  }
)

em_ht <- pmap(
  list(h0_fk[idx], Xk2[idx], beta_int2[idx], chi_hat2[idx]),
  function(x, y, z, k) {
    tryCatch(
      main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
              beta_int = z, chi_hat = k, max_itr = max_itr, tor = tor),
      error = function(err)
        list("beta" = NA, "Yk" = NA, "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
    )
  }
)

h1_likelihood1 <- map_dbl(em_hc, ~ .x$likelihoods[[length(.x$likelihoods)]])
h1_likelihood2 <- map_dbl(em_ht, ~ .x$likelihoods[[length(.x$likelihoods)]])

fk_tbl_idx <- tibble(
  gene_id    = names(em_hc),
  beta1      = map_dbl(em_hc, "beta"),
  beta2      = map_dbl(em_ht, "beta"),
  lfc        = log2(beta2 / beta1),
  fk_mean1   = map_dbl(em_hc, "fk_mean"),
  fk_mean2   = map_dbl(em_ht, "fk_mean"),
  h0_fk      = map(em_res, "fk")[idx],
  h0_fk_mean = map_dbl(em_res, "fk_mean")[idx],
  h0_fk_var  = map_dbl(em_res, "fk_var")[idx],
  fk_var1    = map_dbl(em_hc, "fk_var"),
  fk_var2    = map_dbl(em_ht, "fk_var"),
  t_stats    = h1_likelihood1 + (h1_likelihood2 * scale_factor) - h0_likelihood[idx]
)

fk_tbl <- bind_rows(fk_tbl[!idx, ], fk_tbl_idx)

fk_tbl <- fk_tbl %>%
  mutate(
    p          = pchisq(2 * t_stats, df = 2, ncp = 0, lower.tail = FALSE),
    padj       = p.adjust(p, method = "BH"),
    delta_mean = fk_mean2 - fk_mean1,
    delta_sd   = (sqrt(fk_var2) - sqrt(fk_var1)) / (sqrt(fk_var2) + sqrt(fk_var1)),
    significant = case_when(
      (padj < sig_p) & ((abs(delta_sd) > 0.2) | (abs(delta_mean) > 20)) ~ "Yes",
      TRUE ~ "No"
    ),
    category = factor(significant, levels = c("Yes", "No"))
  )


# ============================================================
# VISUALIZATION
# ============================================================

# Beta LRT: volcano plot (lfc vs -log10 adjusted p-value)
p_beta_volcano <- beta_tbl %>%
  ggplot(aes(x = lfc, y = -log10(padj_plot), color = category)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#E41A1C", "gray", "#377EB8")) +
  labs(
    y     = expression(-log[10] * "(adjusted p-value)"),
    x     = expression(log[2] * "FC(Treated/Control)"),
    color = expression(beta * " change")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  cowplot::theme_cowplot()

# fk LRT: delta-mean vs delta-SD scatter with marginal histograms
p_fk_scatter <- fk_tbl %>%
  ggplot(aes(x = delta_mean, y = delta_sd, color = significant)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("grey", "#377EB8")) +
  geom_vline(xintercept = c(-20, 20), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey60") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    x     = expression(mu[Treat] - mu[Control]),
    y     = expression((sigma[Treat] - sigma[Control]) / (sigma[Treat] + sigma[Control])),
    color = "Significant"
  ) +
  cowplot::theme_cowplot() +
  theme(legend.position = c(0.75, 0.1))

p_fk_scatter <- ggMarginal(p_fk_scatter, type = "histogram",
                            bins = 50, fill = "grey60", color = "black")

# Metaplot: sum of observed read density across all genes
make_sum_distribution_df <- function(rc, condition_name) {
  total_counts <- reduce(rc$Xk, `+`)
  tibble(k = seq_along(total_counts) - 1, count = total_counts) %>%
    mutate(density = count / sum(count), condition = condition_name)
}

meta_sum_df <- bind_rows(
  make_sum_distribution_df(rc1, "Control"),
  make_sum_distribution_df(rc2, "Treated")
)

p_meta <- ggplot(meta_sum_df,
                 aes(x = k, y = density, color = condition, fill = condition)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("Control" = "grey35", "Treated" = "#377EB8")) +
  scale_fill_manual(values  = c("Control" = "grey35", "Treated" = "#377EB8")) +
  labs(
    x     = "Position relative to TSS (bp)",
    y     = "Normalized observed density",
    color = NULL,
    fill  = NULL
  ) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 12), legend.position = c(0.75, 0.8))

# fk LRT: joint beta-lfc vs delta-SD 2D density
p_fk_beta <- ggplot(fk_tbl, aes(x = lfc, y = delta_sd)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom  = "polygon",
    bins  = 6,
    alpha = 0.8
  ) +
  coord_cartesian(xlim = c(-10, 1), ylim = c(-0.4, 0.4)) +
  scale_fill_viridis_c(name = "Density") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = expression(log[2] ~ FC(beta[Treated] / beta[Control])),
    y = expression((sigma[Treat] - sigma[Control]) / (sigma[Treat] + sigma[Control]))
  ) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 14))


# ============================================================
# OUTPUT
# ============================================================

ggsave(file.path(result_dir, "beta_lfc_volcano.pdf"),
       plot = p_beta_volcano, width = 4.5, height = 2.5)

ggsave(file.path(result_dir, "lrt_mean_vs_std_of_k_pointdensity.pdf"),
       plot = p_fk_scatter, width = 5.5, height = 5)

ggsave(file.path(result_dir, "lrt_mean_vs_std_metaplot_raw.pdf"),
       plot = p_meta, width = 4.5, height = 2.5)

ggsave(file.path(result_dir, "lrt_beta_fk.pdf"),
       plot = p_fk_beta, width = 5.5, height = 4)

write_csv(fk_tbl %>% select(-h0_fk),
          file = file.path(result_dir, "lrt_fk.csv"))
