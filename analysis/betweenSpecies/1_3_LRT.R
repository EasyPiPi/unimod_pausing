# ============================================================
# LRT_between_species.R
#
# Between-species LRT analysis for all six pairwise comparisons.
# For each comparison, outputs one _lrt.txt to
#
# Comparisons:
#   CD4:  human-rhesus, human-baboon, rhesus-baboon
#   CD14: human-rhesus, human-baboon, rhesus-baboon
#
#
# Self-contained: no external R scripts are sourced.
# Dependencies: STADyUM, GenomicRanges, tidyverse (R packages only).
# ============================================================

library(STADyUM)
library(GenomicRanges)
library(tidyverse)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))


# ── Ortholog table (loaded once, shared across all comparisons) ───────────────

ortho_genes <- read.table(
  file.path(.paths$outputs, 'publish/betweenSpecies/primate_conserved_genes.txt'),
  sep = '\t', header = TRUE
)


# ── Comparison configurations ─────────────────────────────────────────────────
#
# scale_factor = reads_species1 / reads_species2
#   (scales species2 counts to species1 sequencing depth in likelihoodRatioTest)
#
# Library sizes used:
#   human  cd4:  38,359,083   baboon cd4:  89,558,782   rhesus cd4:  58,191,605
#   human cd14:  33,294,862   baboon cd14: 63,574,639   rhesus cd14: 32,664,636

comparisons <- list(

  list(
    name         = "human_rhesus_cd4",
    species1     = "rhesus",
    species2     = "human",
    cell_type    = "cd4",
    scale_factor = 58191605 / 38359083,
    tss_file     = "primate_cd4_conserved_tss.txt",
    tss_col      = "hr_tss_category"
  ),

  list(
    name         = "human_baboon_cd4",
    species1     = "baboon",
    species2     = "human",
    cell_type    = "cd4",
    scale_factor = 89558782 / 38359083,
    tss_file     = "primate_cd4_conserved_tss.txt",
    tss_col      = "hb_tss_category"
  ),

  list(
    name         = "rhesus_baboon_cd4",
    species1     = "baboon",
    species2     = "rhesus",
    cell_type    = "cd4",
    scale_factor = 89558782 / 58191605,
    tss_file     = "primate_cd4_conserved_tss.txt",
    tss_col      = "br_tss_category"
  ),

  list(
    name         = "human_rhesus_cd14",
    species1     = "rhesus",
    species2     = "human",
    cell_type    = "cd14",
    scale_factor = 32664636 / 33294862,
    tss_file     = "primate_cd14_conserved_tss.txt",
    tss_col      = "hr_tss_category"
  ),

  list(
    name         = "human_baboon_cd14",
    species1     = "baboon",
    species2     = "human",
    cell_type    = "cd14",
    scale_factor = 63574639 / 33294862,
    tss_file     = "primate_cd14_conserved_tss.txt",
    tss_col      = "hb_tss_category"
  ),

  list(
    name         = "rhesus_baboon_cd14",
    species1     = "baboon",
    species2     = "rhesus",
    cell_type    = "cd14",
    scale_factor = 63574639 / 32664636,
    tss_file     = "primate_cd14_conserved_tss.txt",
    tss_col      = "br_tss_category"
  )

)


# ── Helper functions ──────────────────────────────────────────────────────────
#
# Section A: STADyUM result extraction (copied from helper_function.R verbatim).
# Section B: Pipeline steps (load rates, run LRT, annotate, assemble output).

# --- Section A: STADyUM result extraction ------------------------------------

get_beta_lrt_info <- function(lrt_result, lfc_threshold=0.5, prefix1='control', prefix2='treatment'){
	# Summarize per-gene beta LRT results and rate groups for two samples.
	# Add beta/chi/sd group annotations for `prefix1` and `prefix2`.
	# Attach TSS type information and return a combined data.frame.
	beta_lrt <- betaTbl(lrt_result)

	sample1_beta <- rates(transcriptionRates1(lrt_result))
	sample1_beta <- sample1_beta[,c('geneId', 'chi', 'fkSD', 'betaGroup', 'chiGroup', 'sdGroup')]
	colnames(sample1_beta)[2:6] <- paste0(prefix1, '_', colnames(sample1_beta)[2:6])

	sample2_beta <- rates(transcriptionRates2(lrt_result))
	sample2_beta <- sample2_beta[,c('geneId', 'chi', 'fkSD', 'betaGroup', 'chiGroup', 'sdGroup')]
	colnames(sample2_beta)[2:6] <- paste0(prefix2, '_', colnames(sample2_beta)[2:6])

	beta_lrt <- left_join(beta_lrt, sample1_beta, by = 'geneId')
	beta_lrt <- left_join(beta_lrt, sample2_beta, by = 'geneId')

	beta_lrt <- beta_lrt %>%
  		mutate(
         betaCategory =
           case_when(
             (padj < 0.05) & (lfc > lfc_threshold) ~ "Up",
             (padj < 0.05) & (lfc < -lfc_threshold) ~ "Down",
             TRUE ~ "Others"
           ),
         betaCategory = factor(betaCategory, levels = c("Down", "Others", "Up")))

	return(beta_lrt)
}


get_lrt_chi_stats <- function(lrt_df, control, treat, scale_factor, lfc_threshold) {

  control <- rlang::ensym(control)
  treat <- rlang::ensym(treat)

  lrt_df <- lrt_df %>%
    mutate(
      chi_mean = ((!!treat * scale_factor) + !!control) / 2,
      chi_mean_group = cut(
        chi_mean,
        breaks = quantile(chi_mean,
                          probs = c(0, 1/3, 2/3, 1),
                          na.rm = TRUE),
        labels = c("Low", "Medium", "High"),
        include.lowest = TRUE
      ),
      chi_lfc = log2((!!treat * scale_factor) / !!control),
      chi_lfc_group = case_when(
        chi_lfc <= -lfc_threshold ~ "Decrease",
        chi_lfc >=  lfc_threshold ~ "Increase",
        TRUE          ~ "Unchanged"
      ),
      chi_lfc_group = factor(
        chi_lfc_group,
        levels = c("Decrease", "Unchanged", "Increase")
      )
    )

  return(lrt_df)
}


process_fk_lrt <- function(lrt_df){
  lrt_df$deltaMean <-  lrt_df$fkMean2 - lrt_df$fkMean1
  lrt_df$deltaSD <-  (sqrt(lrt_df$fkVar2) - sqrt(lrt_df$fkVar1))/(sqrt(lrt_df$fkVar2) + sqrt(lrt_df$fkVar1))

  return(lrt_df)
}


# --- Section B: Pipeline steps -----------------------------------------------

# Load rate RDS and remap species-specific gene IDs to ortholog group IDs.
load_and_map_rates <- function(species, cell_type, ortho_genes, outputs_dir) {
  path <- file.path(outputs_dir,
                    'publish/singleCellType/1_transcriptionRate',
                    species, paste0(cell_type, '_rate.RDS'))
  rate_obj <- readRDS(path)

  rate_df <- rates(rate_obj)
  # Rename geneId -> species name using base R to avoid plyr::rename masking.
  names(rate_df)[names(rate_df) == "geneId"] <- species

  rate_df <- rate_df %>%
    inner_join(
      ortho_genes %>% dplyr::select(all_of(c(species, "group"))),
      by = species
    ) %>%
    dplyr::rename(geneId = group)

  rate_obj@rates <- rate_df
  return(rate_obj)
}


# Run LRT, extract beta and chi stats, and compute deltaSD/deltaMean.
run_species_lrt <- function(rate1, rate2, scale_factor, prefix1, prefix2, ortho_genes) {
  lrt <- likelihoodRatioTest(rate1, rate2, scale_factor)

  # STADyUM bug: when baboon is species1, betaTbl may retain raw ENS IDs
  # instead of ortholog group IDs. Remap them here.
  if (startsWith(prefix1, "baboon")) {
    gene_id  <- lrt@betaTbl$geneId
    idx_ens  <- startsWith(gene_id, "ENS")
    if (any(idx_ens)) {
      baboon2group        <- setNames(ortho_genes$group, ortho_genes$baboon)
      gene_id[idx_ens]   <- baboon2group[gene_id[idx_ens]]
      lrt@betaTbl$geneId <- gene_id
    }
  }

  beta_lrt <- get_beta_lrt_info(lrt, prefix1 = prefix1, prefix2 = prefix2)

  # Add deltaSD and deltaMean (species2 minus species1 in the LRT call order).
  beta_lrt <- process_fk_lrt(beta_lrt)

  # Add chi-based stats using species1 as control and species2 as treatment.
  # rlang::inject() splices symbols into the call so ensym() receives bare names.
  beta_lrt <- rlang::inject(get_lrt_chi_stats(
    beta_lrt,
    !!rlang::sym(paste0(prefix1, "_chi")),
    !!rlang::sym(paste0(prefix2, "_chi")),
    scale_factor,
    2
  ))

  return(beta_lrt)
}


# Join conserved TSS category, fill missing entries, and rename to tssType.
add_tss_category <- function(beta_lrt, tss_file, tss_col, outputs_dir) {
  tss_path      <- file.path(outputs_dir, 'publish/betweenSpecies', tss_file)
  conserved_tss <- read.table(tss_path, sep = '\t', header = TRUE,
                               na.strings = c("", "NA"))

  beta_lrt <- left_join(
    beta_lrt,
    conserved_tss[, c('group', tss_col)],
    by = c('geneId' = 'group')
  )

  # Standardize to tssType; fill genes without a conserved TSS.
  beta_lrt$tssType <- beta_lrt[[tss_col]]
  beta_lrt$tssType[is.na(beta_lrt$tssType)] <- 'non-Conserved TSS'
  beta_lrt[[tss_col]] <- NULL

  return(beta_lrt)
}


# Add pause_change column (species1_sdGroup -> species2_sdGroup)
# and join species-specific gene IDs from the ortholog table.
process_lrt_result <- function(beta_lrt, prefix1, prefix2, ortho_genes) {
  beta_lrt$pause_change <- paste0(
    beta_lrt[[paste0(prefix1, "_sdGroup")]],
    '_to_',
    beta_lrt[[paste0(prefix2, "_sdGroup")]]
  )

  beta_lrt <- left_join(beta_lrt, ortho_genes, by = c('geneId' = 'group'))

  return(beta_lrt)
}


# ── Main loop: run all six comparisons ────────────────────────────────────────

for (cfg in comparisons) {
  message("Processing: ", cfg$name)

  cell_upper <- toupper(cfg$cell_type)            # "CD4" or "CD14"
  prefix1    <- paste0(cfg$species1, cell_upper)  # e.g. "rhesusCD4"
  prefix2    <- paste0(cfg$species2, cell_upper)  # e.g. "humanCD4"
  result_dir <- file.path(.paths$outputs, 'publish/betweenSpecies/LRT', cfg$name)

  # Load rates and remap gene IDs to ortholog group IDs.
  rate1 <- load_and_map_rates(cfg$species1, cfg$cell_type, ortho_genes, .paths$outputs)
  rate2 <- load_and_map_rates(cfg$species2, cfg$cell_type, ortho_genes, .paths$outputs)

  # Run LRT and compute beta/chi/fk stats.
  beta_lrt <- run_species_lrt(rate1, rate2, cfg$scale_factor, prefix1, prefix2, ortho_genes)

  # Annotate TSS conservation category.
  beta_lrt <- add_tss_category(beta_lrt, cfg$tss_file, cfg$tss_col, .paths$outputs)

  # Add pause_change and attach species-specific gene IDs.
  beta_lrt <- process_lrt_result(beta_lrt, prefix1, prefix2, ortho_genes)

  write.table(beta_lrt, paste0(result_dir, '_lrt.txt'),
              sep = '\t', quote = FALSE, row.names = FALSE)

  message("Done: ", cfg$name)
}
