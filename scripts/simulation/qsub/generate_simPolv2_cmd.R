#### load packages ####
library(tidyverse)

#### testing files ####
root_dir <- "~/Desktop/github/unimod_pausing"

meta_in <- file.path(root_dir, "metadata/simulation_params_lrt.csv")

#### end of parsing arguments ####

n <- 20000
t <- 40
geneLen <- 2000
zeta_sd <- 1000
zeta_min <- 1500
zeta_max <- 2500

meta <- read_csv(meta_in, show_col_types = FALSE)
meta <- meta %>%
  mutate(cmd = paste(
    "~/Desktop/github/SimPolv2/bin/simPol_Release", "-n", n,
    "-a", a_range, "-b", b_range, "-z", z,
    "-t", t, "-s", s, "-k", k_range,
    paste0("--kSd=", ksd_range), paste0("--addSpace=", add_space_range),
    paste0("--geneLen=", geneLen), paste0("--zetaSd=", zeta_sd),
    paste0("--zetaMax=", zeta_max), paste0("--zetaMin=", zeta_min),
    "-d", paste0("./", param_id)
  ))

write(meta$cmd, file.path(root_dir, "scripts/simulation", "simulation_rate_lrt_simPolv2.sh"))