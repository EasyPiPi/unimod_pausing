#### load packages ####
library(tidyverse)

#### parse paths and parameters ####
root_dir <- Sys.getenv("PROJECT_ROOT", ".")
simpol_bin <- Sys.getenv("SIMPOL_BIN", "simPol_Release")

meta_in <- file.path(root_dir, "metadata/simulation_params_lrt.csv")

n <- 20000
t <- 40
geneLen <- 2000
zeta_sd <- 1000
zeta_min <- 1500
zeta_max <- 2500

meta <- read_csv(meta_in, show_col_types = FALSE)
meta <- meta %>%
  mutate(cmd = paste(
    simpol_bin, "-n", n,
    "-a", a_range, "-b", b_range, "-z", z,
    "-t", t, "-s", s, "-k", k_range,
    paste0("--kSd=", ksd_range), paste0("--addSpace=", add_space_range),
    paste0("--geneLen=", geneLen), paste0("--zetaSd=", zeta_sd),
    paste0("--zetaMax=", zeta_max), paste0("--zetaMin=", zeta_min),
    "-d", paste0("./", param_id)
  ))

output_sh <- file.path(root_dir, "scripts/simulation/pause_escape", "simulation_rate_lrt.sh")
write(meta$cmd, output_sh)