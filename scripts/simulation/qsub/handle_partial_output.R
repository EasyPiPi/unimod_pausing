#### load packages ####
library(tidyverse)

#### testing files ####
root_dir <- "~/Desktop/github/unimod_human"

# tbl_in <- file.path(root_dir, "metadata/simulation_params_steric_hindrance_partial.csv")
tbl_in <- file.path(root_dir, "metadata/simulation_params_pause_release_partial.csv")

par_tbl <- read_csv(tbl_in)
par_tbl <- par_tbl %>% mutate(param_id = str_remove(param_id, ".RDS"))
pars <- str_split(str_remove(par_tbl$param_id, "k"),
          "ksd|kmin|kmax|l|m|a|b|g|zsd|zmin|zmax|t|n|h|z|s", simplify = TRUE) 
colnames(pars) <- c("k", "ksd", "kmin", "kmax", "l", "m", "a", "b",
                    "g", "z", "zsd", "zmin", "zmax", "t", "n", "s", "h")

par_tbl[, 2:NCOL(par_tbl)] <- pars[, c("k", "ksd", "m", "a", "b", "g", "z", "t", "n", "s", "h", "l")]

par_tbl %>% write_csv(tbl_in)
