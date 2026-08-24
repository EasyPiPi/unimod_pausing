library(STADyUM)
library(GenomicRanges)
library(tidyverse)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/betweenCellType/LRT/")

if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
}

# ==============================================================
# HELPER FUNCTIONS  (inlined from helper_function.R)
# ==============================================================

get_beta_lrt_info <- function(lrt_result, lfc_threshold = 0.5,
                               prefix1 = 'control', prefix2 = 'treatment') {
  beta_lrt <- betaTbl(lrt_result)

  sample1_beta <- rates(transcriptionRates1(lrt_result))
  sample1_beta <- sample1_beta[, c('geneId', 'chi', 'fkSD', 'betaGroup', 'chiGroup', 'sdGroup')]
  colnames(sample1_beta)[2:6] <- paste0(prefix1, '_', colnames(sample1_beta)[2:6])

  sample2_beta <- rates(transcriptionRates2(lrt_result))
  sample2_beta <- sample2_beta[, c('geneId', 'chi', 'fkSD', 'betaGroup', 'chiGroup', 'sdGroup')]
  colnames(sample2_beta)[2:6] <- paste0(prefix2, '_', colnames(sample2_beta)[2:6])

  beta_lrt <- left_join(beta_lrt, sample1_beta, by = 'geneId')
  beta_lrt <- left_join(beta_lrt, sample2_beta, by = 'geneId')

  # tss_type <- as.data.frame(pauseRegions(transcriptionRates1(lrt_result)))
  # tss_type <- tss_type[, c('gene_id', 'tss_type')]
  # colnames(tss_type) <- c('geneId', 'tssType')
  # beta_lrt <- left_join(beta_lrt, tss_type, by = 'geneId')

  beta_lrt <- beta_lrt %>%
    mutate(
      betaCategory = case_when(
        (padj < 0.05) & (lfc >  lfc_threshold) ~ "Up",
        (padj < 0.05) & (lfc < -lfc_threshold) ~ "Down",
        TRUE ~ "Others"
      ),
      betaCategory = factor(betaCategory, levels = c("Down", "Others", "Up"))
    )

  return(beta_lrt)
}


get_lrt_chi_stats <- function(lrt_df, control, treat, scale_factor, lfc_threshold) {
  control <- rlang::ensym(control)
  treat   <- rlang::ensym(treat)

  lrt_df <- lrt_df %>%
    mutate(
      chi_mean = ((!!treat * scale_factor) + !!control) / 2,
      chi_mean_group = cut(
        chi_mean,
        breaks = quantile(chi_mean, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
        labels = c("Low", "Medium", "High"),
        include.lowest = TRUE
      ),
      chi_lfc = log2((!!treat * scale_factor) / !!control),
      chi_lfc_group = case_when(
        chi_lfc <= -lfc_threshold ~ "Decrease",
        chi_lfc >=  lfc_threshold ~ "Increase",
        TRUE ~ "Unchanged"
      ),
      chi_lfc_group = factor(chi_lfc_group, levels = c("Increase", "Unchanged", "Decrease"))
    )

  return(lrt_df)
}

# ==============================================================
# SPECIES CONFIGURATION
# ==============================================================

species_config <- list(

  human = list(
    cd4_rate_in  = file.path(.paths$outputs, 'publish/singleCellType/1_transcriptionRate/human/cd4_rate.RDS'),
    cd14_rate_in = file.path(.paths$outputs, 'publish/singleCellType/1_transcriptionRate/human/cd14_rate.RDS'),
    scale_factor = 38359083 / 33294862,
    result_subdir = 'human'
  ),

  rhesus = list(
    cd4_rate_in  = file.path(.paths$outputs, 'publish/singleCellType/1_transcriptionRate/rhesus/cd4_rate.RDS'),
    cd14_rate_in = file.path(.paths$outputs, 'publish/singleCellType/1_transcriptionRate/rhesus/cd14_rate.RDS'),
    scale_factor = 58191605 / 32664636,
    result_subdir = 'rhesus'
  ),

  baboon = list(
    cd4_rate_in  = file.path(.paths$outputs, 'publish/singleCellType/1_transcriptionRate/baboon/cd4_rate.RDS'),
    cd14_rate_in = file.path(.paths$outputs, 'publish/singleCellType/1_transcriptionRate/baboon/cd14_rate.RDS'),
    scale_factor = 89558782 / 63574639,
    result_subdir = 'baboon'
  )

)

# ==============================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================

run_lrt_analysis <- function(species, cfg) {
  message("Processing: ", species)

  cd4_rate     <- readRDS(cfg$cd4_rate_in)
  cd14_rate    <- readRDS(cfg$cd14_rate_in)
  scale_factor <- cfg$scale_factor

  # Core LRT
  lrt      <- likelihoodRatioTest(cd4_rate, cd14_rate, scale_factor)
  beta_lrt <- get_beta_lrt_info(lrt, prefix1 = 'CD4', prefix2 = 'CD14')
  beta_lrt <- get_lrt_chi_stats(beta_lrt, CD4_chi, CD14_chi, scale_factor, 0.8)
  beta_lrt$pause_change    <- paste0(beta_lrt$CD4_sdGroup, '_to_', beta_lrt$CD14_sdGroup)

  # Species-specific post-processing
    # Rhesus and baboon: compute deltaSD (normalized) and deltaMean directly
  beta_lrt$deltaSD   <- (beta_lrt$CD14_fkSD - beta_lrt$CD4_fkSD) /
                        (beta_lrt$CD14_fkSD + beta_lrt$CD4_fkSD)
  beta_lrt$deltaMean <- beta_lrt$fkMean2 - beta_lrt$fkMean1

  out_dir <- paste0(result_dir, cfg$result_subdir, "/")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  out_path <- paste0(out_dir, '/lrt_result.txt')

  write.table(beta_lrt, out_path, sep = '\t', quote = FALSE, row.names = FALSE)
  message("Written: ", out_path)
}

for (sp in names(species_config)) {
  run_lrt_analysis(sp, species_config[[sp]])
}
