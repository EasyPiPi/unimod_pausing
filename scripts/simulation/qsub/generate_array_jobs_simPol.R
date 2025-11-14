library(tidyverse)

# Remove scientific notation in the entire R session
options(scipen = 100)

# set up parameters for the model
threads <- 1

k_range <- c(50, 70)
ksd_range <- 25
kmin <- 17
kmax <- 200

geneLen <- 2000
m <- 250
l <- geneLen - k_range - m

a_range <- c(0.1, 1, 10, 50, 100, 1000)
b_range <- c(0.1, 1, 10, 50, 100, 1000)
g <- 1

z <- 2000
zeta_sd <- 1000
# zeta_sd <- 0
zeta_min <- 1500
zeta_max <- 2500

s <- 33
add_space_range <- c(0, 17, 37)

n <- 20000
t <- 20

# dirs in the github repo
root_dir <- "~/Desktop/github/unimod_human"
script_dir <- file.path(root_dir, "scripts/simulation/qsub")
metadata_dir <- file.path(root_dir, "metadata")
# output dir on the cluster
output_dir <- "~/projects/Snakemake_projects/unimod_human/outputs/simulation/data"
# a table to handle some failed jobs
partial_st_tbl <- read_csv(file.path(root_dir, "metadata/simulation_params_steric_hindrance_partial.csv"))
partial_pr_tbl <- read_csv(file.path(root_dir, "metadata/simulation_params_pause_release_partial.csv"))

# generate meta dataframes for snakemake
generate_metadataframe <-
  function(k_range, ksd_range, l, m, a_range, b_range, g, z, t, n, s, add_space_range) {
    param_tbl <-
      tidyr::expand_grid(k_range, ksd_range, m, a_range, b_range, g, z, t, n, s, add_space_range)
    param_tbl <- dplyr::mutate(param_tbl,
                               l = geneLen - k_range - m,
                               param_id = paste0("k", k_range, "ksd", ksd_range, "kmin", kmin, "kmax", kmax,
                                                 "l", l, "m", m, "a", a_range, "b", b_range,
                                                 "g", g, "z", z, "zsd", zeta_sd, "zmin", zeta_min,
                                                 "zmax", zeta_max, 
                                                 "t", t, "n", n, "s", s, "h", add_space_range))
    
    # param_tbl$param_id <- stringr::str_remove_all(param_tbl$param_id, "\\.")
    param_tbl <- dplyr::select(param_tbl, param_id, tidyr::everything())
    return(param_tbl)
  } 

## meta-dataframe for steric hindrance ##
param_st_tbl <-
  generate_metadataframe(k_range = k_range, ksd_range = ksd_range, l = l, m = m,
                         a_range = a_range, b_range = b_range,
                         g = g, z = z, t = t, n = n, s = s,
                         add_space_range = add_space_range)

readr::write_csv(x = param_st_tbl, file = file.path(metadata_dir, "simulation_params_steric_hindrance.csv"))

param_st_tbl <- param_st_tbl %>%
  mutate(cmd = paste("./simPol.R", "-p", threads, "-n", n,
                     "-a", a_range, "-b", b_range, "-g", g, "-z", z,
                     "-t", t, "-s", s, "-k", k_range,
                     paste0("--kSd=", ksd_range), paste0("--addSpace=", add_space_range),
                     paste0("--geneLen=", geneLen), paste0("--zetaSd=", zeta_sd),
                     paste0("--zetaMax=", zeta_max), paste0("--zetaMin=", zeta_min),
                     "-d", output_dir))

write(param_st_tbl$cmd, file.path(script_dir, paste0("simulation_rate_steric_hindrance.sh")))

# filter out files have been generated, overwrite simulation_rate_steric_hindrance.sh
if (NROW(partial_st_tbl)) {
  todo_st_tbl <- param_st_tbl %>%
    filter(!param_id %in% partial_st_tbl$param_id)
  # write a single sh file for the array job
  write(todo_st_tbl$cmd, file.path(script_dir, paste0("simulation_rate_steric_hindrance.sh")))
}

## meta-dataframe for pause release ##
k_range <- 50
ksd_range <- c(0, 5, 15)
add_space_range <- 17

param_pr_tbl <-
  generate_metadataframe(k_range = k_range, ksd_range = ksd_range, l = l, m = m,
                         a_range = a_range, b_range = b_range,
                         g = g, z = z, t = t, n = n, s = s,
                         add_space_range = add_space_range)

readr::write_csv(x = param_pr_tbl, file = file.path(metadata_dir, "simulation_params_pause_release.csv"))

param_pr_tbl <- param_pr_tbl %>%
  mutate(cmd = paste("./simPol.R", "-p", threads, "-n", n,
                     "-a", a_range, "-b", b_range, "-g", g, "-z", z,
                     "-t", t, "-s", s, "-k", k_range,
                     paste0("--kSd=", ksd_range), paste0("--addSpace=", add_space_range),
                     paste0("--geneLen=", geneLen), paste0("--zetaSd=", zeta_sd),
                     paste0("--zetaMax=", zeta_max), paste0("--zetaMin=", zeta_min),
                     "-d", output_dir))

write(param_pr_tbl$cmd, file.path(script_dir, paste0("simulation_rate_pause_release.sh")))

# filter out files have been generated
if (NROW(partial_pr_tbl)) {
  todo_pr_tbl <- param_pr_tbl %>%
    filter(!param_id %in% partial_pr_tbl$param_id)
  # write a single sh file for the array job
  write(todo_pr_tbl$cmd, file.path(script_dir, paste0("simulation_rate_pause_release.sh")))
}

#### continue to run some jobs with small alpha and beta ####
param_st_tbl <- param_st_tbl %>%
  mutate(cmd = paste("./simPol.R", "-p", threads, "-n", n,
                     "-a", a_range, "-b", b_range, "-g", g, "-z", z,
                     "-t", t, "-s", s, "-k", k_range,
                     paste0("--kSd=", ksd_range), paste0("--addSpace=", add_space_range),
                     paste0("--geneLen=", geneLen), paste0("--zetaSd=", zeta_sd),
                     paste0("--zetaMax=", zeta_max), paste0("--zetaMin=", zeta_min),
                     "--continue=TRUE", "-d", output_dir))

# param_st_tbl <- param_st_tbl %>% filter(a_range == 0.1, b_range == 0.1)

write(param_st_tbl$cmd, file.path(script_dir, paste0("simulation_rate_steric_hindrance.sh")))

param_pr_tbl <- param_pr_tbl %>%
  mutate(cmd = paste("./simPol.R", "-p", threads, "-n", n,
                     "-a", a_range, "-b", b_range, "-g", g, "-z", z,
                     "-t", t, "-s", s, "-k", k_range,
                     paste0("--kSd=", ksd_range), paste0("--addSpace=", add_space_range),
                     paste0("--geneLen=", geneLen), paste0("--zetaSd=", zeta_sd),
                     paste0("--zetaMax=", zeta_max), paste0("--zetaMin=", zeta_min),
                     "--continue=TRUE", "-d", output_dir))

# param_pr_tbl <- param_pr_tbl %>% filter(a_range == 0.1, b_range == 0.1)

write(param_pr_tbl$cmd, file.path(script_dir, paste0("simulation_rate_pause_release.sh")))

#### generate table for LRT ####
# Note SimPol have been updated 

# generate meta dataframes for snakemake
generate_metadataframe_lrt <-
  function(k_range, ksd_range, l, a_range, b_range, z, t, n, s, add_space_range) {
    param_tbl <-
      tidyr::expand_grid(k_range, ksd_range, a_range, b_range, z, t, n, s, add_space_range)
    param_tbl <- dplyr::mutate(param_tbl,
                               l = geneLen - k_range,
                               param_id = paste0("k", k_range, "ksd", ksd_range, "kmin", kmin, "kmax", kmax,
                                                 "l", l, "a", a_range, "b", b_range,
                                                 "z", z, "zsd", zeta_sd, "zmin", zeta_min,
                                                 "zmax", zeta_max, 
                                                 "t", t, "n", n, "s", s, "h", add_space_range))
    
    # param_tbl$param_id <- stringr::str_remove_all(param_tbl$param_id, "\\.")
    param_tbl <- dplyr::select(param_tbl, param_id, tidyr::everything())
    return(param_tbl)
  } 

k_range <- 50
ksd_range <- 25
a_range <- 1
b_range <- c(0.1, 0.2, 0.5, 0.8, 1, 1.2, 2, 5, 10)
t <- 40

param_beta_tbl <-
  generate_metadataframe_lrt(k_range = k_range, ksd_range = ksd_range, l = l,
                         a_range = a_range, b_range = b_range,
                         z = z, t = t, n = n, s = s,
                         add_space_range = add_space_range)

a_range <- c(0.1, 0.2, 0.5, 0.8, 1, 1.2, 2, 5, 10)
b_range <- 1

param_alpha_tbl <-
  generate_metadataframe_lrt(k_range = k_range, ksd_range = ksd_range, l = l,
                         a_range = a_range, b_range = b_range,
                         z = z, t = t, n = n, s = s,
                         add_space_range = add_space_range)

readr::write_csv(x = bind_rows(param_alpha_tbl, param_beta_tbl),
                 file = file.path(metadata_dir, "simulation_params_lrt.csv"))

# x <- param_alpha_tbl %>% filter(!a_range %in% c(0.1, 1, 10))
# y <- param_beta_tbl %>% filter(!b_range %in% c(0.1, 1, 10))

# zeta_sd <- 0 

param_lrt_tbl <-
  # x %>% bind_rows(y) %>% 
  param_alpha_tbl %>% bind_rows(param_beta_tbl) %>% 
  mutate(cmd = paste("./simPol.R", "-n", n,
                     "-a", a_range, "-b", b_range, "-z", z,
                     "-t", t, "-s", s, "-k", k_range,
                     paste0("--kSd=", ksd_range), paste0("--addSpace=", add_space_range),
                     paste0("--geneLen=", geneLen), paste0("--zetaSd=", zeta_sd),
                     paste0("--zetaMax=", zeta_max), paste0("--zetaMin=", zeta_min),
                     "-d", output_dir)) %>%
  unique()

write(param_lrt_tbl$cmd, file.path(script_dir, paste0("simulation_rate_lrt.sh")))
