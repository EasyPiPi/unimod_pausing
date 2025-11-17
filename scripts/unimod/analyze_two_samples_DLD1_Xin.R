#### note ####
# this script currently takes care of samples with or without spike-ins
# note that it uses sum of reads across samples

#### load packages ####
library(tidyverse)
library(ggpubr)
library(ggpointdensity)
library(viridis)
library(ggExtra)

#### testing files ####
root_dir <- "/Users/zeng/Desktop/project/YiXin_Likelihood"

helper_tc_in <- file.path(root_dir, "codes/result1/helper_function_em_two_condition.R")
helper_pr_in <- file.path(root_dir, "codes/result1/helper_function_em_pause_release.R")
source(helper_tc_in)
source(helper_pr_in)

# rc1_in <- file.path(root_dir, "data/result1/PROseq-DLD1-aoi-NELFC_NVP2_Ctrl-SE/pause_release/rate.RDS")
# rc2_in <- file.path(root_dir, "data/result1/PROseq-DLD1-aoi-NELFC_NVP2-SE/pause_release/rate.RDS")


# rc1_in <- file.path(root_dir, "data/result1/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE/pause_release/rate.RDS")
# rc2_in <- file.path(root_dir, "data/result1/PROseq-DLD1-aoi-NELFC_Fp-SE/pause_release/rate.RDS")

rc1_in <- file.path(root_dir, "data/result1/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE/pause_release/rate.RDS")
rc2_in <- file.path(root_dir, "data/result1/PROseq-DLD1-aoi-NELFC_Auxin-SE/pause_release/rate.RDS")

spike_in <- file.path(root_dir, "data/result1/scaling_factor.csv")

quantile_normalization <- "identity"

# Null is for summing reads across samples, bioreplicate is for comparison between
# biological replicates to ensure not too many false positive are seen when there is
# no real changes

replicates <- "all"
# replicates <- "-bioreplicate"

# result_dir <-
#   file.path(root_dir, "outputs/result1/",
#             paste0("PROseq-DLD1-aoi-NELFC_NVP2-SE", "-", quantile_normalization, "-", replicates))
# result_dir <-
#   file.path(root_dir, "outputs/result1/",
#             paste0("PROseq-DLD1-aoi-NELFC_Fp-SE", "-", quantile_normalization, "-", replicates))

result_dir <-
  file.path(root_dir, "outputs/result1/",
            paste0("PROseq-DLD1-aoi-NELFC_Auxin-SE", "-", quantile_normalization, "-", replicates))

# read in summary statistics from the varied pause site model 
rc1 <- readRDS(rc1_in)
rc2 <- readRDS(rc2_in)

# Number of genes in each sample
NROW(rc1)
NROW(rc2)

rc1$max_count <- map_dbl(rc1$Xk, max)
rc1 <- rc1[rc1$max_count > 8,]

rc2$max_count <- map_dbl(rc2$Xk, max)
rc2 <- rc2[rc2$max_count > 8,]
# union set of genes being analyzed
gn_union <- intersect(rc1$gene_id, rc2$gene_id)
rc1 <- rc1[match(gn_union, rc1$gene_id), ]
rc2 <- rc2[match(gn_union, rc2$gene_id), ]


scale_tbl <- read_csv(spike_in, show_col_types = FALSE)
# subset the right table
scale_tbl <- scale_tbl[str_detect(rc2_in, scale_tbl$sample), ]


lambda1 <- scale_tbl$control_1 + ifelse(is.na(scale_tbl$control_2), 0, scale_tbl$control_2)
lambda2 <- scale_tbl$treated_1 + ifelse(is.na(scale_tbl$treated_2), 0, scale_tbl$treated_2)


# set up parameters
k <- 50
kmin <- 1
kmax <- 200 # also used as k on the poisson case

rnap_size <- 50
zeta <- 2000

# criteria for significance, p.adj < 0.05 and 2 fold differences
sig_p <- 0.05
beta_lfc1 <- log2(1.2/1)
beta_lfc2 <- log2(0.8/1)



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

scale_factor  <- lambda1/lambda2
# chi_hat for control and test sets
chi_hat1 <- rc1$chi
chi_hat2 <- rc2$chi
# max iterations and tolerance for EM
max_itr = 500
tor = 1e-6
# run EM for multiple combinations of parameters
em_res <- pmap(list(Xk1, Xk2, beta_int, chi_hat, chi_hat1, chi_hat2),
               function(x, y, z, k, m, n) {
                 tryCatch(main_EM_h0(fk_int, Xk1 = x, Xk2 = y, kmin, kmax, beta_int = z,
                                     chi_hat = k,  chi_hat1 = m, chi_hat2 = n,
                                     max_itr, tor),
                          error = function(err) {
                            # handling the error, one of the cases is when 
                            # there is no read counts in the pause region
                            list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
                          })
               }
)
# calculate t1 and t2 for H0
# t1_h0 <- map_dbl(em_res, ~ sum(.x[["Yk1"]]), .default = NA)
# t2_h0 <- map_dbl(em_res, ~ sum(.x[["Yk2"]]), .default = NA)

h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])

beta_tbl <-
  tibble(gene_id = rc1$gene_id,
         beta1 = rc1$beta,
         beta2 = rc2$beta,
         lfc = log2( beta2 / beta1),
         fk_mean1 = rc1$fk_mean,
         fk_mean2 = rc2$fk_mean,
         fk_var1 = rc1$fk_var,
         fk_var2 = rc2$fk_var,
         # use eq (25) instead of (31) to compute T stats
         t_stats = rc1$likelihood + rc2$likelihood - h0_likelihood)

# some genes with negative T stats, fix them
idx <- beta_tbl$t_stats < 0

# use parameter estimates from h0 as initial values for EM in h1
h0_beta <- map_dbl(em_res, "beta")
h0_fk1 <- map(em_res, "fk1") 
h0_fk2 <- map(em_res, "fk2")

em_hc <- pmap(list(h0_fk1[idx], Xk1[idx], h0_beta[idx], chi_hat1[idx]),
              function(x, y, z, k) {
                tryCatch(main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
                                 beta_int = z, chi_hat = k,
                                 max_itr = max_itr, tor = tor),
                         error = function(err) {
                           list("beta" = NA, "Yk" = NA,
                                "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
                         })
              }
)

em_ht <- pmap(list(h0_fk2[idx], Xk2[idx], h0_beta[idx], chi_hat2[idx]),
              function(x, y, z, k) {
                tryCatch(main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
                                 beta_int = z, chi_hat = k,
                                 max_itr = max_itr, tor = tor),
                         error = function(err) {
                           list("beta" = NA, "Yk" = NA,
                                "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
                         })
              }
)

h1_likelihood1 <- map_dbl(em_hc, ~ .x$likelihoods[[length(.x$likelihoods)]])
h1_likelihood2 <- map_dbl(em_ht, ~ .x$likelihoods[[length(.x$likelihoods)]])

beta_tbl_idx <- 
  tibble(gene_id = names(em_hc),
         beta1 = map_dbl(em_hc, "beta"),
         beta2 = map_dbl(em_ht, "beta"),
         lfc = log2( beta2 / beta1),
         fk_mean1 = map_dbl(em_hc, "fk_mean"),
         fk_mean2 = map_dbl(em_ht, "fk_mean"),
         fk_var1 = map_dbl(em_hc, "fk_var"),
         fk_var2 = map_dbl(em_ht, "fk_var"),
         t_stats = h1_likelihood1 + h1_likelihood2 - h0_likelihood[idx])

# beta_tbl <- beta_tbl %>%
#   bind_cols(
#     bind_rows(
#       pmap(list(s1, s2, t1_h0, t2_h0, t1_h1, t2_h1),
#            function(x, y, z, k, m, n) {
#              beta_lrt(s1 = x, s2 = y, t1_h0 = z, t2_h0 = k, t1_h1 = m, t2_h1 = n)
#            })
#     ))

beta_tbl <- bind_rows(beta_tbl[!idx, ], beta_tbl_idx)

beta_tbl <- beta_tbl %>%
  mutate(p = pchisq(2 * t_stats, df = 1, ncp = 0, lower.tail = F, log.p = FALSE))

beta_tbl <- beta_tbl %>% mutate(padj = p.adjust(p, method = "BH"))

beta_tbl <- beta_tbl %>%
  mutate(beta = (beta1 + beta2) / 2,
         logbeta_zeta = log2(beta * zeta),
         category =
           case_when(
             (padj < 0.05) & (lfc > beta_lfc1) ~ "Up",
             (padj < 0.05) & (lfc < beta_lfc2) ~ "Down",
             TRUE ~ "Others"
           ),
         category = factor(category, levels = c("Up", "Down", "Others")))

p <- beta_tbl %>%
  ggplot(aes(x = logbeta_zeta, y = lfc, color = category)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values=c("#E41A1C", "#377EB8", "gray")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ylim(-6, 6) +
  xlim(-10, 10) +
  labs(y = expression(log[2]*"FC(Treated/Control)"),
       x = expression(log[2]*"mean("*beta*zeta*")"),
       color = "Category") +
  cowplot::theme_cowplot()

ggsave(file.path(result_dir, "beta_mean_vs_lfc.pdf"), plot = p,
       width = 5, height = 3)





## LRT for beta distribution##
# need to jointly do EM one more time for H0, which assume betas are the same
# between conditions
# initialize fk with some reasonable values based on heuristic


fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
fk_int <- fk_int / sum(fk_int)
# try uniform distribution which gives similar results
# fk_int <- rep(1 / kmax, kmax)

# collect and construct stats
s1 <- rc1$s
#s2 <- rc2$s * lambda1/lambda2
s2 <- rc2$s

Xk1 <- rc1$Xk
Xk2 <- rc2$Xk
#Xk2 <- lapply(rc2$Xk, function(x) x * lambda1 / lambda2)
# gene body length, assumed to be same between conditions
M <- rc1$N

# some values inherit from EM for H1, could be further integrated into EM here
chi_hat <- (s1 + s2) / M
# beta_int <- chi_hat / (t1_h1 + t2_h1)
beta_int1 <- rc1$beta
beta_int2 <- rc2$beta

# chi_hat for control and test sets
chi_hat1 <- rc1$chi
chi_hat2 <- rc2$chi
# max iterations and tolerance for EM
max_itr = 500
tor = 1e-6
# run EM for multiple combinations of parameters
em_res <- pmap(list(Xk1, Xk2, beta_int1, beta_int2, chi_hat, chi_hat1, chi_hat2),
               function(x, y, z1,z2, k, m, n) {
                 tryCatch(main_EM_fk_h0(fk_int, Xk1 = x, Xk2 = y, kmin, kmax, beta_int1 = z1, beta_int2=z2,
                                     chi_hat = k,  chi_hat1 = m, chi_hat2 = n,scale_factor=scale_factor,
                                     max_itr, tor),
                          error = function(err) {
                            # handling the error, one of the cases is when 
                            # there is no read counts in the pause region
                            list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
                          })
               }
)
# calculate t1 and t2 for H0
# t1_h0 <- map_dbl(em_res, ~ sum(.x[["Yk1"]]), .default = NA)
# t2_h0 <- map_dbl(em_res, ~ sum(.x[["Yk2"]]), .default = NA)

h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])
h0_fk <- map(em_res, ~ .x$fk)
h0_fk_mean <- map_dbl(em_res, ~ .x$fk_mean[[length(.x$fk_mean)]])
h0_fk_var <- map_dbl(em_res, ~ .x$fk_var[[length(.x$fk_var)]])

beta_tbl <-
  tibble(gene_id = rc1$gene_id,
         beta1 = rc1$beta,
         beta2 = rc2$beta,
         lfc = log2( beta2 / beta1),
         fk_mean1 = rc1$fk_mean,
         fk_mean2 = rc2$fk_mean,
         h0_fk=h0_fk,
         h0_fk_mean=h0_fk_mean,
         h0_fk_var=h0_fk_var,
         fk_var1 = rc1$fk_var,
         fk_var2 = rc2$fk_var,
         # use eq (25) instead of (31) to compute T stats
         t_stats = rc1$likelihood + (rc2$likelihood * scale_factor) - h0_likelihood)

# some genes with negative T stats, fix them
idx <- beta_tbl$t_stats < 0

# use parameter estimates from h0 as initial values for EM in h1
h0_fk <- map(em_res, "fk") 

em_hc <- pmap(list(h0_fk[idx], Xk1[idx], beta_int1[idx], chi_hat1[idx]),
              function(x, y, z, k) {
                tryCatch(main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
                                 beta_int = z, chi_hat = k,
                                 max_itr = max_itr, tor = tor),
                         error = function(err) {
                           list("beta" = NA, "Yk" = NA,
                                "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
                         })
              }
)

em_ht <- pmap(list(h0_fk[idx], Xk2[idx], beta_int2[idx], chi_hat2[idx]),
              function(x, y, z, k) {
                tryCatch(main_EM(fk_int = x, Xk = y, kmin = kmin, kmax = kmax,
                                 beta_int = z, chi_hat = k,
                                 max_itr = max_itr, tor = tor),
                         error = function(err) {
                           list("beta" = NA, "Yk" = NA,
                                "fk_mean" = NA, "fk_var" = NA, "likelihoods" = NA)
                         })
              }
)

h1_likelihood1 <- map_dbl(em_hc, ~ .x$likelihoods[[length(.x$likelihoods)]])
h1_likelihood2 <- map_dbl(em_ht, ~ .x$likelihoods[[length(.x$likelihoods)]])

beta_tbl_idx <- 
  tibble(gene_id = names(em_hc),
         beta1 = map_dbl(em_hc, "beta"),
         beta2 = map_dbl(em_ht, "beta"),
         lfc = log2( beta2 / beta1),
         fk_mean1 = map_dbl(em_hc, "fk_mean"),
         fk_mean2 = map_dbl(em_ht, "fk_mean"),
         h0_fk = map(em_res, "fk")[idx], 
         h0_fk_mean = map_dbl(em_res, "fk_mean")[idx], 
         h0_fk_var= map_dbl(em_res, "fk_var")[idx],
         fk_var1 = map_dbl(em_hc, "fk_var"),
         fk_var2 = map_dbl(em_ht, "fk_var"),
         t_stats = h1_likelihood1 + h1_likelihood2 - h0_likelihood[idx])

# beta_tbl <- beta_tbl %>%
#   bind_cols(
#     bind_rows(
#       pmap(list(s1, s2, t1_h0, t2_h0, t1_h1, t2_h1),
#            function(x, y, z, k, m, n) {
#              beta_lrt(s1 = x, s2 = y, t1_h0 = z, t2_h0 = k, t1_h1 = m, t2_h1 = n)
#            })
#     ))

beta_tbl <- bind_rows(beta_tbl[!idx, ], beta_tbl_idx)

beta_tbl <- beta_tbl %>%
  mutate(p = pchisq(2 * t_stats, df = 2, ncp = 0, lower.tail = F, log.p = FALSE))

beta_tbl <- beta_tbl %>% mutate(padj = p.adjust(p, method = "BH"))

beta_tbl$delta_mean <- beta_tbl$fk_mean2 - beta_tbl$fk_mean1
beta_tbl$delta_sd <- (sqrt(beta_tbl$fk_var2) - sqrt(beta_tbl$fk_var1))/(sqrt(beta_tbl$fk_var2) + sqrt(beta_tbl$fk_var1))


beta_tbl <- beta_tbl %>%
  mutate(
         significant =
           case_when(
             (padj < 0.05) & ((abs(delta_sd) > 0.2) | (abs(delta_mean) > 20))   ~ "Yes",
             TRUE ~ "No"
           ),
         category = factor(significant, levels = c('Yes', 'No')))



p <- beta_tbl %>%
  ggplot(aes(x = delta_mean, y = delta_sd, color = significant)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values=c("#E41A1C", "#377EB8")) +
  geom_vline(xintercept = -20, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey60") +
  labs(
         x = "Differences in mean of k",
         y = "Normalized difference in sd of k",
       color = "Significant") +
  cowplot::theme_cowplot()+
  theme(
    legend.position = c(0.8, 0.05))

p <- ggMarginal(p, type = "histogram", bins = 50, fill = "grey60", color = "black")

ggsave(filename = file.path(result_dir, "lrt_mean_vs_std_of_k_pointdensity.pdf"), plot = p,
       width = 6, height = 5)

beta_tbl <- beta_tbl %>% select(-h0_fk)
write_csv(beta_tbl, file =  file.path(result_dir, "lrt_fk.csv"))










