#### load packages ####
library(GenomicRanges)
library(tidyverse)

#### snakemake files ####

#### testing files ####
root_dir <- '~/Desktop/project/YiXin_Likelihood'
# root_dir <- "~/projects/Snakemake_projects/unimod_human"

table_dir <- file.path(root_dir, "pausing_distribution/simulation/data/subsampling/high")
figure_dir <- file.path(root_dir, "outputs/simulation/figures")

#meta_in <- file.path(root_dir, "metadata/simulation_params_lrt.csv")
helper_in <- file.path(root_dir, "codes/result1/helper_function_em_two_condition.R")

compute_in_range_proportion <- function(mean, sd, lower = 0, upper = 200) {
  p_upper <- pnorm((upper - mean) / sd)
  p_lower <- pnorm((lower - mean) / sd)
  proportion <- p_upper - p_lower
  return(proportion)
}

# min and max position of pause site
kmin = 1
kmax = 200
matched_gb_len <- 2e4 - kmax

#### end of parsing arguments ####

source(helper_in)

file_names <- list.files(table_dir)

rate_tbls <- data.frame(
  id = file_names,
  k = as.numeric(str_extract(file_names, "(?<=k)\\d+(?=ksd)")),
  ksd = as.numeric(str_extract(file_names, "(?<=ksd)\\d+(?=kmin)")),
  a = rep(1, length(file_names)),
  b = rep(1, length(file_names)),
  z = rep(2000, length(file_names)),
  t = rep(40, length(file_names)),
  n = rep(20000, length(file_names)),
  s = rep(33, length(file_names)),
  h = rep(17, length(file_names)),
  l = rep(1950, length(file_names))
  )

rate_tbls$tbl <-
  map(file.path(table_dir, file_names), readRDS)

k_range <- c(30, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90)
ksd_range <- c(10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70)
 
k_combinations <- expand.grid(k = k_range, ksd = ksd_range)
k_combinations <- k_combinations %>%
    rowwise() %>%
    mutate(proportion = compute_in_range_proportion(k, ksd),
           above_0.9 = proportion > 0.9) %>%
    ungroup()

valid_combinations <- k_combinations %>%
  filter(above_0.9) %>%
  select(k, ksd)

rate_tbls_filtered <- rate_tbls %>%
    semi_join(valid_combinations, by = c("k", "ksd"))

rate_tbls_filtered <- rate_tbls_filtered %>%
  mutate(k = factor(k, levels = k_range),
         ksd = factor(ksd, levels = ksd_range))

rate_tbls_filtered <- rate_tbls_filtered %>% unnest(cols = tbl)
rate_tbls_filtered$t_h1 <- map_dbl(rate_tbls_filtered$Yk, sum)

#### LRT for beta ####
# do multiple LRTs between a control setting and varied parameters
get_beta_lrt_results <- function(ctrl_tbl, test_tbl) {
  # number of combinations
  n_rep <- NROW(test_tbl) / NROW(ctrl_tbl)
  # note how to check integer in R
  if (n_rep == round(n_rep)) {
    ctrl_tbl <- bind_rows(rep(list(ctrl_tbl), n_rep))
    } else stop("check combination numbers for comparison")
  # initialize fk with some reasonable values based on heuristic
  fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
  fk_int <- fk_int / sum(fk_int)
  # collect and construct stats
  s1 <- ctrl_tbl %>% pull(rc_gb)
  s2 <- test_tbl %>% pull(rc_gb)
  t1_h1 <- ctrl_tbl %>% pull(t_h1)
  t2_h1 <- test_tbl %>% pull(t_h1)
  Xk1 <- ctrl_tbl %>% pull(Xk)
  Xk2 <- test_tbl %>% pull(Xk)
  # gene body length
  if (!is.null(matched_gb_len)) {
    M <- matched_gb_len
  } else {
    M <- 1800
  }
  # some values inherit from EM for H1, could be further integrated into EM here
  chi_hat <- (s1 + s2) / M
  beta_int <- chi_hat / (t1_h1 + t2_h1)

  beta_int1 <- ctrl_tbl %>% pull(beta_adp)
  beta_int2 <- test_tbl %>% pull(beta_adp)
  # chi_hat for control and test sets
  chi_hat1 <- ctrl_tbl$chi
  chi_hat2 <- test_tbl$chi
  # max iterations and tolerance for EM
  max_itr = 500
  tor = 1e-6
  # run EM for multiple combinations of parameters
  em_res <- pmap(list(Xk1, Xk2, beta_int1, beta_int2, chi_hat, chi_hat1, chi_hat2),
                 function(x, y, z1, z2, k, m, n) {
                   tryCatch(main_EM_fk_h0(fk_int, Xk1 = x, Xk2 = y, kmin, kmax, beta_int1 = z1, beta_int2=z2,
                              chi_hat = k,  chi_hat1 = m, chi_hat2 = n, 1,
                              max_itr, tor),
                            error = function(err) {
                              # handling the error, one of the cases is when 
                              # there is no read counts in the pause region
                              list("beta" = NA, "Yk1" = NA, "Yk2" = NA)
                              })
                   }
                 )
  # calculate t1 and t2 for H0
  t1_h0 <- map_dbl(em_res, ~ sum(.x[["Yk1"]]), .default = NA)
  t2_h0 <- map_dbl(em_res, ~ sum(.x[["Yk2"]]), .default = NA)
  
  h0_likelihood <- map_dbl(em_res, ~ .x$likelihoods[[length(.x$likelihoods)]])

  # calculate T stats
  # # use eq (31) to compute T stats
  # test_res <- pmap(list(s1, s2, t1_h0, t2_h0, t1_h1, t2_h1),
  #                  function(x, y, z, k, m, n) {
  #                    beta_lrt(s1 = x, s2 = y, t1_h0 = z, t2_h0 = k, t1_h1 = m, t2_h1 = n)
  #                  }) %>%
  #   bind_rows()
  # # summarize results
  # test_res <-
  #   test_res %>% bind_cols(test_tbl) %>%
  #   select(a, b, chi, t_stats, p)
  
  # use eq (25) to compute T stats
  test_res <- test_tbl
  test_res <- 
    test_res %>%
    mutate(t_stats = ctrl_tbl$likelihood + test_tbl$likelihood - h0_likelihood,
           p = pchisq(2 * t_stats, df = 2, ncp = 0, lower.tail = F, log.p = FALSE)) %>%
    select(k, ksd, fk_mean, fk_var, chi, t_stats, p)
    
  test_res <- test_res %>%
  mutate(fk_mean_ctrl = ctrl_tbl$fk_mean,
         fk_var_ctrl = ctrl_tbl$fk_var)
  test_res <-
    test_res %>% mutate(
      # adjust p values after multiple comparisons
      sig = (p.adjust(p, method = "fdr") < .05),
      logp = log(p),
      t1_h0 = t1_h0,
      t2_h0 = t2_h0,
      t1_h1 = t1_h1,
      t2_h1 = t2_h1
    )
  
  return(test_res)
}


#nrow(k_combinations)
for (i in 1:nrow(valid_combinations)) {
  k_val <- valid_combinations$k[i]
  ksd_val <- valid_combinations$ksd[i]

  ctrl_tbl <- rate_tbls_filtered %>% filter(k == k_val, ksd == ksd_val)
  
  if (i == 1){
    lrt_res_all <- get_beta_lrt_results(ctrl_tbl, rate_tbls_filtered)
    lrt_res_all$k_ctrl <-  k_val
    lrt_res_all$ksd_ctrl <-  ksd_val
  }
  else{
    lrt_res <- get_beta_lrt_results(ctrl_tbl, rate_tbls_filtered)
    lrt_res$k_ctrl <-  k_val
    lrt_res$ksd_ctrl <-  ksd_val
 

    lrt_res_all <- rbind(lrt_res_all, lrt_res)
  }
  if (i %% 10 == 0){print(i)}
}

for (i in 1:nrow(valid_combinations)) {
  k_val <- valid_combinations$k[i]
  ksd_val <- valid_combinations$ksd[i]
  lrt_res_sbuset <- lrt_res_all %>% filter(k_ctrl == k_val, ksd_ctrl==ksd_val)
  power_by_k_ksd <- lrt_res_sbuset %>%
    group_by(k,ksd) %>%
    summarise(
        n = n(),
        power = mean(sig, na.rm = TRUE),.groups = "drop" 
    )

  p1 <- ggplot(power_by_k_ksd, aes(x = factor(ksd), y = factor(k), fill = power)) +
    geom_tile() +
    geom_text(aes(label = round(power, 2)), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(x = "ksd", y = "k", fill = "Power") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
    output_name <- paste0('k',k_val,'ksd', ksd_val, '_filtered.pdf')
    ggsave(paste0(figure_dir, '/',output_name),p1, width = 6, height = 4)
  }

lrt_res_all <- lrt_res_all %>%
  mutate(
    nd_k = as.numeric(as.character(k)) - as.numeric(as.character(k_ctrl)),
    nd_ksd = (as.numeric(as.character(ksd)) - as.numeric(as.character(ksd_ctrl))) / 
             (as.numeric(as.character(ksd)) + as.numeric(as.character(ksd_ctrl)))
  )

# lrt_res_all <- lrt_res_all %>%
#     mutate(
#         nd_k_pred = fk_mean - fk_mean_ctrl,
#         nd_ksd_pred = (sqrt(fk_var) - sqrt(fk_var_ctrl)) / 
#             (sqrt(fk_var) + sqrt(fk_var_ctrl)))


power_by_k_ksd_normalized <- lrt_res_all %>%
    group_by(nd_k,nd_ksd) %>%
    summarise(
        n = n(),
        power = mean(sig, na.rm = TRUE),.groups = "drop" 
    )

p1 <- ggplot(power_by_k_ksd_normalized, aes(x = nd_k, y = nd_ksd, color = power)) +
      geom_point(size = 1, alpha = 0.8) +
      scale_color_viridis_c(option = "plasma") +
      theme_minimal() +
       theme(             
          axis.ticks = element_line(color = "black"), 
          axis.line = element_line(color = "black")
         ) +
      labs(x = "Differences in mean of k", y = "Normalized difference in sd of k", color = "Statistical\npower")

ggsave(paste0(figure_dir, '/', 'statistical_power_filtered.pdf'),p1, width = 5, height = 4)

