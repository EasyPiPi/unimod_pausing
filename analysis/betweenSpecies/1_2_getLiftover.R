# ============================================================
# Generate sorted TSS BED files from STADyUM rate objects
# for human, rhesus, and baboon.
#
# Reads the RDS files produced by 1_EstimateRates and writes
# one BED file per sample.
# ============================================================


# ============================================================
# LIBRARIES
# ============================================================

library(STADyUM)   # rates(), @pauseRegions, @counts
library(plyranges) # anchor_5p(), mutate() on GRanges


# ============================================================
# CONFIGURATION
# ============================================================

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/")

# Each entry maps one input RDS file to one output BED file.
sample_configs <- list(
  list(rds = file.path(result_dir, "rhesus/cd4_rate.RDS"),  bed = file.path(result_dir, "rhesus/cd4_pause_region.bed")),
  list(rds = file.path(result_dir, "rhesus/cd14_rate.RDS"), bed = file.path(result_dir, "rhesus/cd14_pause_region.bed")),
  list(rds = file.path(result_dir, "baboon/cd4_rate.RDS"),  bed = file.path(result_dir, "baboon/cd4_pause_region.bed")),
  list(rds = file.path(result_dir, "baboon/cd14_rate.RDS"), bed = file.path(result_dir, "baboon/cd14_pause_region.bed"))
)


library(STADyUM)
library(plyranges)
library(rtracklayer)
library(GenomeInfoDb)
library(S4Vectors)

export_pause_region_bed <- function(rds_file, bed_file) {
  message("Reading: ", rds_file)

  rate_obj <- readRDS(rds_file)

  pause_region <- rate_obj@pauseRegions

  mcols(pause_region) <- DataFrame(
    name = pause_region$gene_id
  )

  seqlevelsStyle(pause_region) <- "UCSC"

  pause_region <- sort(pause_region)

  dir.create(dirname(bed_file), showWarnings = FALSE, recursive = TRUE)

  message("Writing: ", bed_file)
  export.bed(pause_region, bed_file)
}

invisible(lapply(sample_configs, function(x) {
  export_pause_region_bed(
    rds_file = x$rds,
    bed_file = x$bed
  )
}))
