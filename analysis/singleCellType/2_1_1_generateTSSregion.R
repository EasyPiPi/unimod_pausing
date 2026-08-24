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
  Sys.getenv("PROJECT_ROOT", "."),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/")

# Each entry maps one input RDS file to one output BED file.
sample_configs <- list(
  list(rds = file.path(result_dir, "human/cd4_rate.RDS"),   bed = file.path(result_dir, "human/cd4_tss_sort.bed")),
  list(rds = file.path(result_dir, "human/cd14_rate.RDS"),  bed = file.path(result_dir, "human/cd14_tss_sort.bed")),
  list(rds = file.path(result_dir, "rhesus/cd4_rate.RDS"),  bed = file.path(result_dir, "rhesus/cd4_tss_sort.bed")),
  list(rds = file.path(result_dir, "rhesus/cd14_rate.RDS"), bed = file.path(result_dir, "rhesus/cd14_tss_sort.bed")),
  list(rds = file.path(result_dir, "baboon/cd4_rate.RDS"),  bed = file.path(result_dir, "baboon/cd4_tss_sort.bed")),
  list(rds = file.path(result_dir, "baboon/cd14_rate.RDS"), bed = file.path(result_dir, "baboon/cd14_tss_sort.bed"))
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Extract a sorted, BED-like data frame of TSS positions from a STADyUM object.
# TSS is defined as the 5' end of each pause region (single base, width = 1).
# Pause counts and gene body counts are attached from the counts slot.
# Returned columns (in order): seqnames, start, end, gene_id, gCount, strand, pCount
make_tss_bed_df_from_stadyum <- function(rate_obj) {
  stopifnot(!is.null(rate_obj))

  rate_df <- rates(rate_obj)

  # Anchor pause regions at 5' end to get single-base TSS positions
  tss_gr <- rate_obj@pauseRegions %>%
    plyranges::anchor_5p() %>%
    mutate(width = 1)

  # Keep only genes present in the rates table
  tss_gr <- tss_gr[tss_gr$gene_id %in% rate_df$geneId]

  # Attach pause and gene body counts
  idx <- match(tss_gr$gene_id, rate_obj@counts$gene_id)
  tss_gr$pCount <- rate_obj@counts[idx, ]$summarizedPauseCounts
  tss_gr$gCount <- rate_obj@counts[idx, ]$summarizedGbCounts

  # Sort and select BED columns in required order
  bed_df <- as.data.frame(sort(tss_gr))
  bed_df[, c("seqnames", "start", "end", "gene_id", "gCount", "strand", "pCount")]
}


# ============================================================
# MAIN PROCESSING LOOP
# ============================================================

for (s in sample_configs) {
  rate_obj <- readRDS(s$rds)
  bed_df   <- make_tss_bed_df_from_stadyum(rate_obj)
  write.table(bed_df, file = s$bed,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
