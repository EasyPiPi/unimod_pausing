#### log file ####
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
helper_in <- snakemake@params[["helper"]]
table_dir <- snakemake@params[["table_dir"]]
figure_dir <- snakemake@params[["figure_dir"]]

#### load packages ####
library(GenomicRanges)
library(tidyverse)

#### testing files ####
root_dir <- "~/Desktop/github/unimod_pausing"

# table_dir <- file.path(root_dir, "outputs/simulation/tables/lrt_pause_escape")
# figure_dir <- file.path(root_dir, "outputs/simulation/figures/lrt_pause_escape")

# helper_in <- file.path(root_dir, "scripts/unimod/helper_function_em_two_condition.R")

#### end of parsing arguments ####
meta_in <- file.path(root_dir, "metadata/simulation_params_lrt.csv")

suffix <- "_subsample_cells"

# min and max position of pause site
kmin <- 1
kmax <- 200
matched_gb_len <- 2e4 - kmax

walk(c(table_dir, figure_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

source(helper_in)

# read in and clean up rate tibbles using meta dataframe
rate_tbls <- read_csv(meta_in, show_col_types = FALSE)

colnames(rate_tbls) <- c("id", "k", "ksd", "a", "b", "z", "t", "n", "s", "h", "l")

read_density <- c("lrt_high", "lrt_median", "lrt_low")

rate_tbls <- tibble(
  read_density = read_density,
  tbl = list(rate_tbls)
) %>%
  unnest(cols = tbl)

rate_tbls$tbl <-
  map(file.path(table_dir, paste0(rate_tbls$id, "_", rate_tbls$read_density, ".RDS")), readRDS)

rate_ranges <- c(0.1, 0.2, 0.5, 0.8, 1, 1.2, 2, 5, 10)
rate_tbls <- rate_tbls %>%
  mutate(
    a = factor(a, levels = rate_ranges),
    b = factor(b, levels = rate_ranges),
    r = ifelse(ksd == 25, "peak", "spike"),
    r = factor(r, levels = c("peak", "spike")),
    read_density = str_remove(read_density, "lrt_"),
    read_density = factor(read_density, levels = unique(read_density))
  )

rc_tbl <- rate_tbls %>% unnest(cols = tbl)

#### LRT for omega ####
# Formulas below are all based on MS version 5
# <t>, formula is located between (17) and (18), and below (20)
rc_tbl$t_h1 <- map_dbl(rc_tbl$Yk, sum)

# do multiple LRTs between a control setting and varied parameters
get_omega_lrt_results <- function(ctrl_tbl, test_tbl) {
  s1 <- ctrl_tbl %>% pull(rc_gb)
  s2 <- test_tbl %>% pull(rc_gb)

  # number of combinations
  n_rep <- length(s2) / length(s1)
  # note how to check integer in R
  if (n_rep == round(n_rep)) s1 <- rep(s1, times = n_rep) else stop("check combination numbers for comparison")

  test_res <- map2(s1, s2, ~ omega_lrt(.x, .y, tao1 = 1 / 2, tao2 = 1 / 2)) %>% bind_rows()
  test_res <-
    test_res %>%
    bind_cols(test_tbl) %>%
    select(a, b, chi, t_stats, p) %>%
    # adjust p values after multiple comparisons
    mutate(
      sig = (p.adjust(p, method = "fdr") < .05),
      logp = log(p)
    )

  return(test_res)
}

#### LRT for beta ####
# do multiple LRTs between a control setting and varied parameters
get_beta_lrt_results <- function(ctrl_tbl, test_tbl) {
  # number of combinations
  n_rep <- NROW(test_tbl) / NROW(ctrl_tbl)
  # note how to check integer in R
  if (n_rep == round(n_rep)) {
    ctrl_tbl <- bind_rows(rep(list(ctrl_tbl), n_rep))
  } else {
    stop("check combination numbers for comparison")
  }
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
  # chi_hat for control and test sets
  chi_hat1 <- ctrl_tbl$chi
  chi_hat2 <- test_tbl$chi
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
    mutate(
      t_stats = ctrl_tbl$likelihood + test_tbl$likelihood - h0_likelihood,
      p = pchisq(2 * t_stats, df = 1, ncp = 0, lower.tail = F, log.p = FALSE)
    ) %>%
    select(a, b, chi, t_stats, p)

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

# perform LRTs
do_lrt <- function(sub_tbl, a, b, fix = NULL, lrt = NULL) {
  # fix both alpha and beta as control
  ctrl_tbl <- sub_tbl %>% filter(a == {{ a }}, b == {{ b }})
  # fix either alpha or beta as test
  if (fix == "a") {
    # alpha is fixed, varied beta
    test_tbl <- sub_tbl %>% filter(a == {{ a }})
  } else if (fix == "b") {
    # beta is fixed, varied alpha
    test_tbl <- sub_tbl %>% filter(b == {{ b }})
  } else {
    stop("Either fix a (initation) or b (pause release) for a LRT")
  }
  # do the test
  if (lrt == "a") {
    lrt_res <- get_omega_lrt_results(ctrl_tbl, test_tbl)
  } else if (lrt == "b") {
    lrt_res <- get_beta_lrt_results(ctrl_tbl, test_tbl)
  } else {
    stop("Either do LRT for a (initation) or b (pause release)")
  }
  return(lrt_res)
}

rc_tbl_gb <- rc_tbl %>%
  group_by(r, read_density) %>%
  nest()

rc_tbl_gb <- rc_tbl_gb %>%
  mutate(
    lrt_a_fix_a = map(data, ~ do_lrt(.x, a = 1, b = 1, fix = "a", lrt = "a")),
    lrt_a_fix_b = map(data, ~ do_lrt(.x, a = 1, b = 1, fix = "b", lrt = "a")),
    lrt_b_fix_a = map(data, ~ do_lrt(.x, a = 1, b = 1, fix = "a", lrt = "b")),
    lrt_b_fix_b = map(data, ~ do_lrt(.x, a = 1, b = 1, fix = "b", lrt = "b"))
  )

# use chi square test for comparison
chi_tbl_gb <- rc_tbl_gb %>% select(read_density, r, data)

get_chi_results <- function(sub_tbl) {
  # In the chi square test, we only focus on differential pausing
  ctrl_tbl <- sub_tbl %>% filter(a == 1, b == 1)
  test_tbl <- sub_tbl %>% filter(a == 1)

  n_rep <- NROW(test_tbl) / NROW(ctrl_tbl)
  # note how to check integer in R
  if (n_rep == round(n_rep)) {
    ctrl_tbl <- bind_rows(rep(list(ctrl_tbl), n_rep))
  } else {
    stop("check combination numbers for comparison")
  }

  tbl <- bind_cols(
    test_tbl,
    ctrl_tbl %>%
      select(rc_tss, rc_gb) %>%
      rename_with(~ paste0(.x, "_ctrl"), everything())
  )

  tbl <- tbl %>% mutate(
    pval = pmap_dbl(
      list(rc_tss, rc_gb, rc_tss_ctrl, rc_gb_ctrl),
      ~ chisq.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2))$p.value
    ),
    padj = p.adjust(pval, method = "fdr"),
    sig = (padj < .05)
  )

  return(tbl)
}

chi_tbl_gb <- chi_tbl_gb %>%
  mutate(chisq = map(data, get_chi_results))

# # debug
# data <- rc_tbl_gb$data[[1]]
# ctrl_tbl <- data %>% filter(a == 1, b == 1)
# test_tbl <- data %>% filter(b == 1)
#
# n <- 106
# main_EM_h0(fk_int, Xk1 = Xk1[[n]], Xk2 = y[[n]], kmin, kmax, beta_int = beta_int[[n]],
#            chi_hat = k,  chi_hat1 = chi_hat1[[n]], chi_hat2 = chi_hat2[[n]],
#            max_itr, tor, percent = 0.01)
# debug(main_EM_h0)

#### visualization ####
prep_plot_tbl <- function(rc_tbl_gb, lrt_col, gb_col) {
  plot_tbl <- rc_tbl_gb %>%
    select(r, read_density, {{ lrt_col }}) %>%
    unnest(cols = {{ lrt_col }}) %>%
    group_by(r, read_density, {{ gb_col }}) %>%
    summarise(sig = mean(sig, na.rm = TRUE), .groups = "drop")
  return(plot_tbl)
}

theme_set(cowplot::theme_cowplot())

# power curves
p <- rc_tbl_gb %>%
  prep_plot_tbl(lrt_a_fix_b, a) %>%
  ggplot(aes(x = a, y = sig, group = interaction(r, read_density))) +
  # geom_line(aes(color = read_density, linetype = r)) +
  geom_line(aes(color = read_density)) +
  labs(
    y = "Statistical Power", x = expression("True " * alpha * zeta),
    color = "Exp. Level", linetype = "Category",
    title = expression("LRT on " * omega)
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("lrt_omega_fix_beta", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

p <- rc_tbl_gb %>%
  prep_plot_tbl(lrt_a_fix_a, b) %>%
  ggplot(aes(x = b, y = sig, group = interaction(r, read_density))) +
  # geom_line(aes(color = read_density, linetype = r)) +
  geom_line(aes(color = read_density)) +
  labs(
    y = "Statistical Power", x = expression("True " * beta * zeta),
    color = "Exp. Level", linetype = "Category",
    title = expression("LRT on " * omega)
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("lrt_omega_fix_alpha", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

pause_lrt_tbl <-
  rc_tbl_gb %>%
  prep_plot_tbl(lrt_b_fix_a, b)

p <- pause_lrt_tbl %>%
  ggplot(aes(x = b, y = sig, group = interaction(r, read_density))) +
  # geom_line(aes(color = read_density, linetype = r)) +
  geom_line(aes(color = read_density)) +
  labs(
    y = "Statistical Power", x = expression("True " * beta * zeta),
    color = "Exp. Level", linetype = "Category",
    title = expression("LRT on " * beta)
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("lrt_beta_fix_alpha", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

# figure for chi square test
chi_tbl <-
  chi_tbl_gb %>%
  prep_plot_tbl(chisq, b)

p <- chi_tbl %>%
  ggplot(aes(x = b, y = sig, group = interaction(r, read_density))) +
  # geom_line(aes(color = read_density, linetype = r)) +
  geom_line(aes(color = read_density)) +
  labs(
    y = "Statistical Power", x = expression("True " * beta * zeta),
    color = "Exp. Level", linetype = "Category",
    title = expression("LRT on " * beta)
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("chi_sq_test", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

# For beta estimates, compare LRT test and chi square test
pause_lrt_tbl <- pause_lrt_tbl %>% mutate(method = "LRT")
chi_tbl <- chi_tbl %>% mutate(method = "Chi-sq")

test_cpr <- bind_rows(pause_lrt_tbl, chi_tbl)
test_cpr_w <- test_cpr %>%
  filter(b == 0.8 | b == 1.2) %>%
  pivot_wider(
    id_cols = c(r, read_density, b),
    names_from = method,
    values_from = sig,
    names_glue = "{method}_{.value}"
  ) %>%
  mutate(
    diff = LRT_sig - `Chi-sq_sig`
  )

p <- test_cpr %>%
  filter(b != 1) %>%
  ggplot(aes(x = b, y = sig, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~read_density, nrow = 3) +
  labs(
    x = expression("True " * beta * zeta),
    y = "Statistical Power"
  )

ggsave(file.path(figure_dir, paste0("test_comparison", suffix, ".png")),
  plot = p,
  width = 8, height = 8
)

p <- test_cpr %>%
  filter(b != 1) %>%
  ggplot(aes(x = b, y = sig, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = sig),
    position = position_dodge(width = 0.8),
    vjust = 0.5,
    size = 3,
    color = "black"
  ) +
  facet_wrap(~read_density, nrow = 3) +
  labs(
    x = expression("True " * beta * zeta),
    y = "Statistical Power"
  )

ggsave(file.path(figure_dir, paste0("test_comparison_w_num", suffix, ".png")),
  plot = p,
  width = 8, height = 8
)

p <- rc_tbl_gb %>%
  prep_plot_tbl(lrt_b_fix_a, b) %>%
  filter(!b %in% c(0.1, 10)) %>%
  ggplot(aes(x = b, y = sig, group = interaction(r, read_density))) +
  # geom_line(aes(color = read_density, linetype = r)) +
  geom_line(aes(color = read_density)) +
  labs(
    y = "Statistical Power", x = expression("True " * beta * zeta),
    color = "", linetype = "Category"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.72, 0.3)
  )

ggsave(file.path(figure_dir, paste0("lrt_beta_fix_alpha", suffix, "_grant.png")),
  plot = p,
  width = 4, height = 2.5
)

p <- rc_tbl_gb %>%
  prep_plot_tbl(lrt_b_fix_b, a) %>%
  ggplot(aes(x = a, y = sig, group = interaction(r, read_density))) +
  # geom_line(aes(color = read_density, linetype = r)) +
  geom_line(aes(color = read_density)) +
  labs(
    y = "False Positive Rate", x = expression("True " * alpha * zeta),
    color = "Exp. Level", linetype = "Category",
    title = expression("LRT on " * beta)
  ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, paste0("lrt_beta_fix_beta", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

# Chi distribution for sanity check
y_int <-
  tibble(
    "read_density" = str_remove(read_density, "lrt_"),
    "y" = c(0.0902, 0.0489, 0.0233)
  ) %>%
  mutate(read_density = factor(read_density, levels = read_density))

p <- rc_tbl %>%
  filter(b == 1) %>%
  ggplot() +
  geom_boxplot(aes(x = a, y = chi, color = read_density)) +
  geom_hline(aes(yintercept = y, color = read_density),
    data = y_int,
    linetype = "dashed"
  ) +
  labs(x = expression("True " * alpha * zeta), y = expression(chi), color = "Exp. Level")

ggsave(file.path(figure_dir, paste0("chi_distribution_with_targeted_read_density", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

# t stats
p <- rc_tbl_gb %>%
  select(r, read_density, lrt_a_fix_b) %>%
  unnest(lrt_a_fix_b) %>%
  ggplot(aes(x = a, y = t_stats, color = read_density)) +
  geom_boxplot() +
  labs(y = "T statistic", x = expression("True " * alpha * zeta), color = "Exp. Level") +
  # facet_grid(r ~ .) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("lrt_omega_fix_beta_t_stats", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

p <- rc_tbl_gb %>%
  select(r, read_density, lrt_a_fix_a) %>%
  unnest(lrt_a_fix_a) %>%
  ggplot(aes(x = b, y = t_stats, color = read_density)) +
  geom_boxplot() +
  labs(y = "T statistic", x = expression("True " * beta * zeta), color = "Exp. Level") +
  # facet_grid(r ~ .) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("lrt_omega_fix_alpha_t_stats", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

p <- rc_tbl_gb %>%
  select(r, read_density, lrt_b_fix_a) %>%
  unnest(lrt_b_fix_a) %>%
  mutate(t_stats = ifelse(t_stats < 0, 0, t_stats)) %>%
  ggplot(aes(x = b, y = t_stats, color = read_density)) +
  geom_boxplot() +
  labs(y = "T statistic", x = expression("True " * beta * zeta), color = "Exp. Level") +
  # facet_grid(r ~ .) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("lrt_beta_fix_alpha_t_stats", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)

p <- rc_tbl_gb %>%
  select(r, read_density, lrt_b_fix_b) %>%
  unnest(lrt_b_fix_b) %>%
  mutate(t_stats = ifelse(t_stats < 0, 0, t_stats)) %>%
  ggplot(aes(x = a, y = t_stats, color = read_density)) +
  geom_boxplot() +
  labs(y = "T statistic", x = expression("True " * alpha * zeta), color = "Exp. Level") +
  # facet_grid(r ~ .) +
  cowplot::theme_cowplot()

ggsave(file.path(figure_dir, paste0("lrt_beta_fix_beta_t_stats", suffix, ".png")),
  plot = p,
  width = 8, height = 3
)