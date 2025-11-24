#### log file ####
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
rds_in <- snakemake@input[["rds"]]

helper_in <- snakemake@params[["helper"]]
sample_cell <- as.integer(snakemake@params[["sample_cell"]])
sample_n <- snakemake@params[["sample_n"]] # number of times to sample
lambda <- snakemake@params[["lambda_exp"]] # scaling factor to match simulation to coverage in experimental data
matched_len <- snakemake@params[["matched_len"]]
sel_sample <- snakemake@params[["sel_sample"]]
count_rnap <- snakemake@params[["count_rnap"]]

rate_out <- snakemake@output[["rate_tbl"]]
bw_out <- snakemake@output[["bw"]]

# args <- commandArgs(trailingOnly = TRUE)
# sel_sample <- args[1]
# lambda <- args[2]

#### load packages ####
library(GenomicRanges)
library(tidyverse)
library(rtracklayer)

#### testing files ####
# root_dir <- "~/Desktop/github/unimod_pausing"
# helper_in <-
#   "/home/yixin/Desktop/github/unimod_pausing/scripts/unimod/helper_function_em_pause_release.R"

# sel_sample <-
#   "k30ksd55kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t40n20000s33h17"
# lambda <- "lrt_high"

# rate_dir <-
#   file.path(root_dir, "outputs/simulation/tables/lrt_pause_distribution")
# rds_in <-
#   file.path(
#     root_dir, "outputs/simulation/data_fk/simPol",
#     paste0(sel_sample, "/positions/position_matrix_400000.csv")
#   )
# rate_out <-
#   file.path(rate_dir, paste0(sel_sample, "_", lambda, ".RDS"))
# bw_out <- file.path(rate_dir, paste0(sel_sample, "_", lambda, ".bw"))

# dir.create(rate_dir, showWarnings = FALSE, recursive = TRUE)

# sample_n <- 100
# sample_cell <- 5000
# matched_len <- 2e4
# count_rnap <- FALSE
#### end of parsing arguments ####

source(helper_in)
#### set parameters ####
kmin <- 1
kmax <- 200
matched_gb_len <- matched_len - kmax

# calculate total spacing from sample name
spacing <-
  sum(as.integer(
    str_split(
      str_split(sel_sample, "[:digit:]s", simplify = TRUE)[, 2],
      "h"
    )[[1]]
  ))

k <- as.integer(str_remove(str_split(sel_sample, "ksd", simplify = TRUE)[, 1], "k"))

start_point <- 0.99 * 1e6 # set a start coordinate for the simulated gene

# set lambda according to the read coverage of PRO-seq in the control samples
# from Dukler et al. 2017
if (is.null(lambda)) {
  lambda <- NULL
} else if (lambda == "extreme_high") {
  lambda <- 100000
} else if (lambda == "high") {
  lambda <- 188.3
} else if (lambda == "median") {
  lambda <- 102.1
} else if (lambda == "low") {
  lambda <- 48.7
} else if (str_detect(lambda, "lrt")) {
  lambda_lrt <- 1800 * 5000 / 2670 * c(0.0233, 0.0489, 0.0902)
  names(lambda_lrt) <- c("lrt_low", "lrt_median", "lrt_high")
  lambda <- lambda_lrt[[lambda]]
} else {
  stop("lambda should be either high, median or low.")
}

# retrieve Pol II positions (the last step) for re-sampling
polII_pos <- read.csv(rds_in, header = FALSE)
polII_pos <- t(polII_pos[, -ncol(polII_pos)])
# total number of cells
total_cell <- NCOL(polII_pos)
# gene length of the simulated gene
gene_len <- NROW(polII_pos) - 1

#### initiation and pause release rate estimates ####
# generate regions for read counting
gn_rng <-
  GRanges(
    seqnames = rep("chr1", 3),
    IRanges(
      start = c(1, kmax + 1, 1),
      end = c(kmax, gene_len, spacing)
    )
  )

gn_rng <- shift(gn_rng, shift = start_point)

region_names <- c("tss", "gb", "landing")
names(gn_rng) <- region_names

len <- as.list(width(gn_rng))
names(len) <- region_names

# set seeds for random sampling
seeds <- seq(from = 2013, by = 1, length.out = sample_n)
# a list to Granges for PolII positions
polII_grng <- list()
# a list recording number of RNAPs at or before pause site
if (count_rnap) rnap_n_ls <- list()

for (i in 1:sample_n) {
  set.seed(seeds[[i]])
  sel_cells <- sample(1:total_cell, size = sample_cell, replace = TRUE)
  res_pos <- polII_pos[, sel_cells]
  # get rid of position 1, which is always 1
  res_pos <- res_pos[-1, ]
  res_pos <- as.matrix(res_pos)
  if (count_rnap) {
    # get pause sites
    pause_site <- idx[sel_cells, 1] - 1
    # generate data mask
    # inspired by https://stackoverflow.com/questions/47732085/sum-of-some-positions-in-a-row-r
    res_shape <- dim(res_pos)
    after_pause_len <- res_shape[1] - pause_site
    mask_mx <- map2(
      pause_site, after_pause_len,
      function(x, y) c(rep(TRUE, x), rep(FALSE, y))
    )
    mask_mx <- Matrix::Matrix(unlist(mask_mx), nrow = res_shape[1], ncol = res_shape[2])
    # calculate rnap number before pause site for every cell
    rnap_n_ls[[i]] <- colSums(res_pos * mask_mx)
  }
  # combine PolII positions across all cells
  res_all <- rowSums(res_pos)
  # generate bigwig for positive strand
  polII_grng[[i]] <-
    GRanges(
      seqnames = "chr1",
      IRanges(
        start = (1 + start_point):(gene_len + start_point),
        width = 1
      ),
      score = res_all,
      strand = "+",
      seqlengths = c("chr1" = gene_len * 10) + start_point
    )

  rm(res_pos, res_all)
}

# read counting
summarise_bw <-
  function(bw, grng) {
    rc <- grng %>%
      plyranges::group_by_overlaps(bw) %>%
      plyranges::group_by(query) %>%
      plyranges::summarise(score = sum(score))
    if (!1 %in% rc$query) {
      rc <- rbind(DataFrame(list(query = 1, score = 0)), rc)
    }
    rc <- as.list(rc$score)
    names(rc) <- region_names
    return(rc)
  }

bw_dfs <- tibble(trial = 1:sample_n)
bw_dfs$rc_region <- map(polII_grng, ~ summarise_bw(.x, gn_rng))

bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
bw_dfs$rc_gb <- map_dbl(bw_dfs$rc_region, "gb")
bw_dfs$rc_landing <- map_dbl(bw_dfs$rc_region, "landing")

#### empirical way to calculate steric hindrance at pause site ####
bw_dfs <- bw_dfs %>%
  mutate( # number of RNAPs per cell per gene
    R = (rc_tss + rc_gb) / sample_cell,
    # number of RNAPs in the pause peak per cell per gene
    R_pause = rc_tss / sample_cell,
    # proportion of landing pad being occupied by RNAP, i.e., empirical phi
    polII_prop = rc_landing / sample_cell
  )

# whether to match the simulated number of RNAPs to read coverage in experimental data or not
# here match RNAP number within kmin to kmax, RNAP in gene body will be taken care
# of afterwards
if (!is.null(lambda)) {
  polII_grng <- map(polII_grng, function(grng) {
    grng$score[kmin:kmax] <-
      rpois(length(kmin:kmax), grng$score[kmin:kmax] / sample_cell * lambda)
    # first 20bp get removed because they are usually not seen in sequencing
    # for testing!!!!!
    # grng$score[1:20] <- 0
    return(grng)
  })
  bw_dfs$rc_region <- map(polII_grng, ~ summarise_bw(.x, gn_rng))
  bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
  # bw_dfs$rc_gb <-map_dbl(bw_dfs$rc_region, "gb")
  # bw_dfs$rc_landing <-map_dbl(bw_dfs$rc_region, "landing")
}

# match RNAP number within gene bodies to desired read coverage
if (!is.null(lambda)) {
  pois_mean <- (lambda * bw_dfs$rc_gb / sample_cell) * (matched_gb_len / len$gb)
  bw_dfs$rc_gb <- rpois(length(pois_mean), pois_mean)
  # assign matched gene body length as gene body length
  len$gb <- matched_gb_len
}

# get read counts on each position within pause peak (from kmin to kmax)
bw_dfs$Xk <- map(
  polII_grng,
  ~ .x[(start(.x) >= 990000 + kmin) & (start(.x) <= 990000 + kmax), ]$score
)

# Original model: Poisson based MLEs
# use read count within gene body to pre-estimate chi hat
bw_dfs$chi <- bw_dfs$rc_gb / len$gb

# take care of single pause site or variable pause sites (pause peak)
if (str_detect(sel_sample, "ksd0")) {
  bw_dfs$beta_org <- bw_dfs$chi / map_dbl(bw_dfs$Xk, k)
} else {
  bw_dfs$beta_org <- bw_dfs$chi / (bw_dfs$rc_tss / len$tss)
}

bw_dfs$beta_max_rc <- bw_dfs$chi / map_dbl(bw_dfs$Xk, max)

# Adapted model allows uncertainty in the pause site and steric hindrance
bw_dfs <- bw_dfs %>%
  mutate(alpha_empirical = chi / (1 - polII_prop))

# initialize beta using sum of read counts within pause peak
bw_dfs$Xk_sum <- sapply(bw_dfs$Xk, sum)
bw_dfs$beta_int <- bw_dfs$chi / bw_dfs$Xk_sum

# initialize fk with some reasonable values based on heuristic
fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
fk_int <- fk_int / sum(fk_int)

# estimate rates using EM
em_ls <- list()
main_EM <- possibly(main_EM, otherwise = NA)

for (i in 1:NROW(bw_dfs)) {
  # message("Dealing with the ", i, " row...")
  rc <- bw_dfs[i, ]
  # message("This is gene ", rc$gene_id)
  em_ls[[i]] <- main_EM(
    fk_int = fk_int, Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
    beta_int = rc$beta_int[[1]], chi_hat = rc$chi,
    max_itr = 500, tor = 1e-3
  )
}

# get rate estimates and posterior distribution of pause sites
bw_dfs$beta_adp <- map_dbl(em_ls, "beta", .default = NA)
bw_dfs$Yk <- map(em_ls, "Yk", .default = NA)
bw_dfs$fk <- map(em_ls, "fk", .default = NA)
bw_dfs$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
bw_dfs$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
# calculate Yk / Xk
bw_dfs$proportion_Yk <- sapply(bw_dfs$Yk, sum) / sapply(bw_dfs$Xk, sum)
# whether to count rnap or not also distinguish it's rate estimations or LRT
if (!count_rnap) bw_dfs$likelihood <- map_dbl(em_ls, ~ .x$likelihoods[[length(.x$likelihoods)]])
bw_dfs$flag <- map_chr(em_ls, "flag", .default = NA)

# Combine all GRanges and average the scores, used for bigwig output
combined_grng <- polII_grng[[1]]
combined_grng$score <- rowMeans(sapply(polII_grng, function(x) x$score))

#### save outputs for further analyses ####
saveRDS(bw_dfs, rate_out)
export.bw(object = combined_grng, con = bw_out)