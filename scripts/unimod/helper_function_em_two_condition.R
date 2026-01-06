#### function for adapted model ####
# two conditions, e.g., comparisons between treatments, cell types, or species
# functions for EM based on Gaussian distributed k, and beta is tied in two conditions
get_expectation <- function(fk, Xk, beta) {
  Yk <- Xk / (1 - beta + beta / fk)

  return(Yk)
}

get_fk <- function(Xk, Yk, fk, kmin, kmax) {
  t <- sum(Yk)
  u <- sum(Yk * seq(kmin, kmax))
  v <- sum(Yk * seq(kmin, kmax)^2)

  w <- sum(fk / (1 - fk) * (Xk - Yk), na.rm = TRUE)
  z <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax), na.rm = TRUE)
  r <- sum(fk / (1 - fk) * (Xk - Yk) * seq(kmin, kmax)^2, na.rm = TRUE)

  fk_mean <- (u - z) / (t - w)
  fk_var <- (v - r) / (t - w) - fk_mean^2

  # avoid negative values
  if (fk_var < 1e-10) {
    fk[1:length(fk)] <- 0
    # sometimes it looks like an integer but actually it's not
    fk[round(fk_mean)] <- 1
    return(list("t" = t, "fk" = fk))
  }

  fk <- dnorm(kmin:kmax, mean = fk_mean, sd = fk_var^0.5)
  fk <- fk / sum(fk)

  return(list("t" = t, "fk" = fk, "fk_mean" = fk_mean, "fk_var" = fk_var))
}

get_maximization_h0 <- function(chi_hat, Xk1, Xk2, Yk1, Yk2, fk1, fk2, kmin, kmax, scale_factor) {
  fk1_ls <- get_fk(Xk1, Yk1, fk1, kmin, kmax)
  fk2_ls <- get_fk(Xk2, Yk2, fk2, kmin, kmax)

  beta <- chi_hat / (fk1_ls$t + fk2_ls$t * scale_factor)

  return(list("beta" = beta, "fk1" = fk1_ls$fk, "fk2" = fk2_ls$fk))
}



get_likelihood <- function(beta, chi, Xk, Yk, fk) {
  # part of the original likelihood function associated with beta, Xk and Yk
  # used as criteria to terminate EM
  t <- sum(Yk)
  # take care of the 0s
  idx_1 <- fk != 0
  idx_2 <- 1 - fk != 0
  likelihood <- -t * log(beta) - chi / beta +
    sum(Yk[idx_1] * log(fk[idx_1])) + sum((Xk - Yk)[idx_2] * log(1 - fk[idx_2]))
  return(likelihood)
}

# EM function to estimate parameters when assuming beta is the same between conditions
main_EM_h0 <- function(fk_int, Xk1, Xk2, kmin, kmax, beta_int,
                       chi_hat, chi_hat1, chi_hat2, scale_factor,
                       max_itr = 100, tor = 1e-3) {

  # lists to record changes of likelihood and betas in iterations
  betas <- list()
  likelihoods <- list()
  yks <- list()
  # default flag is normal
  flag <- "normal"

  for (i in 1:max_itr) {
    if (i == 1) {
      Yk1 <- get_expectation(fk_int, Xk1, beta_int)
      Yk2 <- get_expectation(fk_int, Xk2, beta_int)
      hats <- get_maximization_h0(
        chi_hat, Xk1, Xk2, Yk1, Yk2,
        fk_int, fk_int, kmin, kmax, scale_factor
      )
    }
    if (i != 1) {
      Yk1 <- get_expectation(hats$fk1, Xk1, hats$beta)
      Yk2 <- get_expectation(hats$fk2, Xk2, hats$beta)
      hats <- get_maximization_h0(
        chi_hat, Xk1, Xk2, Yk1, Yk2,
        hats$fk1, hats$fk2, kmin, kmax, scale_factor
      )
    }

    likelihoods[[i]] <-
      get_likelihood(beta = hats$beta, chi = chi_hat1, Xk = Xk1, Yk = Yk1, fk = hats$fk1) +
      get_likelihood(beta = hats$beta, chi = chi_hat2, Xk = Xk2, Yk = Yk2, fk = hats$fk2) * scale_factor

    betas[[i]] <- hats$beta
    # for troubleshoot
    # yks[[i]] <- list("Yk1" = sum(Yk1), "Yk2" = sum(Yk2))

    if (any(hats$fk1 == 1) & any(hats$fk2 == 1)) {
      hats$beta <- chi_hat / (Xk1[which(hats$fk1 == 1)] + Xk2[which(hats$fk2 == 1)])
      flag <- "single_site"
      break
    }

    if (i > 1) {
      diff <- likelihoods[[i]] - likelihoods[[i - 1]]
      if (diff <= tor) break
    }
  }

  if (i == max_itr) flag <- "max_iteration"

  return(list(
    "beta" = hats$beta, "Yk1" = Yk1, "Yk2" = Yk2,
    "fk1" = hats$fk1, "fk2" = hats$fk2, "betas" = betas,
    "likelihoods" = likelihoods, "yks" = yks, "flag" = flag
  ))
}

#### Functions for LRTs ####
# formulas are based on the unified model preprint v5
# formula (27), calculate t stats for omega
omega_lrt <- function(s1, s2, tao1, tao2) {
  # compute T statistic and p values
  t_stats <- s1 * log(s1 / (tao1 * (s1 + s2))) + s2 * log(s2 / (tao2 * (s1 + s2)))
  p <- pchisq(2 * t_stats, df = 1, ncp = 0, lower.tail = F, log.p = FALSE)
  return(c("t_stats" = t_stats, "p" = p))
}

# formula (31), calculate t stats for beta
# this equation is problematic, switch to eq (25) for computing T stats
# beta_lrt <- function(s1, s2, t1_h0, t2_h0, t1_h1, t2_h1) {
#
#   t_stats <-
#     s1 * log(s1) + s2 * log(s2) + t1_h1 * log(t1_h1) + t2_h1 * log(t2_h1) -
#     (s1 + s2) * log(s1 + s2) - (s1 + t1_h0) * log(s1 + t1_h0) -
#     (s2 + t2_h0) * log(s2 + t2_h0) - (t1_h0 + t2_h0) * log(t1_h0 + t2_h0) +
#     (s1 + s2 + t1_h0 + t2_h0) * log(s1 + s2 + t1_h0 + t2_h0)
#
#   p <- pchisq(2 * t_stats, df = 1, ncp = 0, lower.tail = F, log.p = FALSE)
#
#   return(c("t_stats" = t_stats, "p" = p))
# }


## adjust for doing LRT in fk
main_EM_fk_h0 <- function(fk_init, Xk1, Xk2, kmin, kmax, beta_int1, beta_int2,
                          chi_hat, chi_hat1, chi_hat2, scale_factor,
                          max_itr = 200, tor = 1e-3) {
  likelihoods <- list()
  flag <- "normal"

  for (i in 1:max_itr) {
    if (i == 1) {
      Yk1 <- get_expectation(fk_init, Xk1, beta_int1)
      Yk2 <- get_expectation(fk_init, Xk2, beta_int2)
      fk_ls <- get_fk(Xk1 + Xk2 * scale_factor, Yk1 + Yk2 * scale_factor, fk_init, kmin, kmax)
    } else {
      Yk1 <- get_expectation(fk, Xk1, beta1)
      Yk2 <- get_expectation(fk, Xk2, beta2)
      fk_ls <- get_fk(Xk1 + Xk2 * scale_factor, Yk1 + Yk2 * scale_factor, fk, kmin, kmax)
    }

    fk <- fk_ls$fk
    t_total <- fk_ls$t

    beta1 <- chi_hat1 / sum(Yk1)
    beta2 <- chi_hat2 / sum(Yk2)

    if (any(fk == 1)) {
      beta1 <- chi_hat1 / (Xk1[which(fk == 1)])
      beta2 <- chi_hat2 / (Xk2[which(fk == 1)])
      flag <- "single_site"
      fk_ls$fk_mean <- which(fk == 1)
      fk_ls$fk_var <- 0
    }
    l1 <- get_likelihood(beta1, chi_hat1, Xk1, Yk1, fk)
    l2 <- get_likelihood(beta2, chi_hat2, Xk2, Yk2, fk)

    likelihoods[[i]] <- l1 + l2 * scale_factor
    if (i > 1 && abs(likelihoods[[i]] - likelihoods[[i - 1]]) < tor) break
  }

  if (i == max_itr) flag <- "max_iteration"
  return(list(
    beta1 = beta1, beta2 = beta2, fk = fk, fk_mean = fk_ls$fk_mean, fk_var = fk_ls$fk_var, likelihoods = likelihoods,
    "flag" = flag
  ))
}