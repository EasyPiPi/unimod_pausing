# ============================================================
# Estimate transcription rates for human, rhesus, and baboon.
#
# Set SPECIES in the CONFIGURATION section, then source this
# script. All visualization code has been removed; the script
# stops after saving the STADyUM rate objects as RDS files.
# ============================================================


# ============================================================
# LIBRARIES
# ============================================================

library(STADyUM)
library(GenomicRanges)
library(plyranges)
library(dplyr)
library(pracma)   # for findpeaks() inside find_valley_threshold()


# ============================================================
# CONFIGURATION
# Modify this section to select a species and update paths.
# ============================================================

# Which species to process: "human", "rhesus", or "baboon"
SPECIES <- "human"

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/")

species_configs <- list(

  human = list(
    sample1_name = "CD4",
    sample2_name = "CD14",
    tsn1_file    = file.path(.paths$data, "proseq/tidgrng/PROseq-HUMAN-CD4_tsn.RDS"),
    tsn2_file    = file.path(.paths$data, "proseq/tidgrng/PROseq-HUMAN-CD14_tsn.RDS"),
    tq_file      = file.path(.paths$data, "proseq/tq/HUMAN.RDS"),
    bw1_plus     = file.path(.paths$data, "proseq/bigwig/p3/PROseq-HUMAN-CD4_plus.bw"),
    bw1_minus    = file.path(.paths$data, "proseq/bigwig/p3/PROseq-HUMAN-CD4_minus.bw"),
    bw2_plus     = file.path(.paths$data, "proseq/bigwig/p3/PROseq-HUMAN-CD14_plus.bw"),
    bw2_minus    = file.path(.paths$data, "proseq/bigwig/p3/PROseq-HUMAN-CD14_minus.bw"),
    out_subdir   = "human"
  ),

  rhesus = list(
    sample1_name = "CD4",
    sample2_name = "CD14",
    tsn1_file    = file.path(.paths$data, "proseq/tidgrng/PROseq-RHESUS-CD4_tsn.RDS"),
    tsn2_file    = file.path(.paths$data, "proseq/tidgrng/PROseq-RHESUS-CD14_tsn.RDS"),
    tq_file      = file.path(.paths$data, "proseq/tq/RHESUS.RDS"),
    bw1_plus     = file.path(.paths$data, "proseq/bigwig/p3/PROseq-RHESUS-CD4_plus.bw"),
    bw1_minus    = file.path(.paths$data, "proseq/bigwig/p3/PROseq-RHESUS-CD4_minus.bw"),
    bw2_plus     = file.path(.paths$data, "proseq/bigwig/p3/PROseq-RHESUS-CD14_plus.bw"),
    bw2_minus    = file.path(.paths$data, "proseq/bigwig/p3/PROseq-RHESUS-CD14_minus.bw"),
    out_subdir   = "rhesus"
  ),

  baboon = list(
    sample1_name = "CD4",
    sample2_name = "CD14",
    tsn1_file    = file.path(.paths$data, "proseq/tidgrng/PROseq-BABOON-CD4_tsn.RDS"),
    tsn2_file    = file.path(.paths$data, "proseq/tidgrng/PROseq-BABOON-CD14_tsn.RDS"),
    tq_file      = file.path(.paths$data, "proseq/tq/BABOON.RDS"),
    bw1_plus     = file.path(.paths$data, "proseq/bigwig/p3/PROseq-BABOON-CD4_plus.bw"),
    bw1_minus    = file.path(.paths$data, "proseq/bigwig/p3/PROseq-BABOON-CD4_minus.bw"),
    bw2_plus     = file.path(.paths$data, "proseq/bigwig/p3/PROseq-BABOON-CD14_plus.bw"),
    bw2_minus    = file.path(.paths$data, "proseq/bigwig/p3/PROseq-BABOON-CD14_minus.bw"),
    out_subdir   = "baboon"
  )
)


# ============================================================
# HELPER FUNCTIONS
# Copied from helper_function.R — only the functions needed
# by this script are included so it runs independently.
# ============================================================

# Select the single most upstream TSS per gene (strand-aware).
# On "+" strand the most upstream TSS has the smallest coordinate;
# on "-" strand it has the largest coordinate.
keep_upstream_tss <- function(tsn) {
  as.data.frame(tsn) %>%
    group_by(ensembl_gene_id, strand) %>%
    slice_min(order_by = ifelse(strand == "+", start, -start),
              n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    as_granges()
}

# Compute TSS distances between two samples for QC.
# Returns a named numeric vector of absolute distances for genes present in both samples.
# TSS coordinates are not modified.
compute_tss_distances <- function(tsn1, tsn2) {
  common   <- intersect(tsn1$ensembl_gene_id, tsn2$ensembl_gene_id)
  tsn1i    <- tsn1[match(common, tsn1$ensembl_gene_id)]
  tsn2i    <- tsn2[match(common, tsn2$ensembl_gene_id)]
  d        <- abs(start(tsn1i) - start(tsn2i))
  names(d) <- tsn1i$ensembl_gene_id
  d
}

# Build promoter-proximal pause regions and trimmed gene body regions
# for read counting. Filters TSNs by score, constructs pause windows of
# width `kmax`, and retains gene bodies longer than `gb_min_length`
# after trimming `trim_len` from each end.
build_readcount_regions <- function(bw_tsn,
                                    transcripts,
                                    tsn_cutoff    = 5,
                                    gb_min_length = 6000,
                                    trim_len      = 2000,
                                    kmax          = 200) {
  if (!inherits(bw_tsn, "GRanges"))
    stop("bw_tsn must be a GRanges.")
  if (!("ensembl_gene_id" %in% colnames(mcols(bw_tsn))))
    stop("bw_tsn must contain 'ensembl_gene_id'.")
  if (!("score" %in% colnames(mcols(bw_tsn))))
    stop("bw_tsn must contain 'score'.")
  if (!inherits(transcripts, "GRanges"))
    stop("transcripts must be a GRanges.")
  if (!("ensembl_gene_id" %in% colnames(mcols(transcripts))))
    stop("transcripts must contain 'ensembl_gene_id'.")

  # Filter low-signal TSNs and build pause windows
  bw_tsn   <- bw_tsn[bw_tsn$score >= tsn_cutoff]
  bw_pause <- promoters(bw_tsn, upstream = 0, downstream = kmax)
  rm_cols  <- intersect(c("score", "type"), colnames(mcols(bw_pause)))
  if (length(rm_cols) > 0) mcols(bw_pause)[rm_cols] <- NULL

  # Derive one TTS position (3' end, width = 1) per gene from the transcript annotation
  gngrng <- transcripts %>%
    plyranges::group_by(ensembl_gene_id) %>%
    plyranges::reduce_ranges_directed() %>%
    sort()
  bw_tts <- gngrng %>%
    plyranges::anchor_3p() %>%
    mutate(width = 1)

  # Build gene bodies spanning from TSN to TTS, apply length filter and trim
  generate_gene_body <- function(bw_tsn, bw_tts, bw_pause, gb_min_length, trim_len) {
    idx_tts <- match(bw_pause$ensembl_gene_id, bw_tts$ensembl_gene_id)
    idx_tsn <- match(bw_pause$ensembl_gene_id, bw_tsn$ensembl_gene_id)
    keep    <- !is.na(idx_tts) & !is.na(idx_tsn)
    if (!any(keep)) stop("No overlapping genes between TSN/TTS and pause set.")

    bw_tts_filtered <- bw_tts[idx_tts[keep]]
    bw_tsn_filtered <- bw_tsn[idx_tsn[keep]]

    gb         <- punion(bw_tsn_filtered, bw_tts_filtered, fill.gap = TRUE)
    gb$gene_id <- bw_tsn_filtered$ensembl_gene_id
    gb_filt    <- gb[width(gb) > gb_min_length]
    gb_filt    <- gb_filt - trim_len
    gb_filt
  }

  bw_gb_filtered <- generate_gene_body(bw_tsn, bw_tts, bw_pause, gb_min_length, trim_len)

  # Match pause regions to the filtered gene bodies
  match_idx         <- match(bw_gb_filtered$gene_id, bw_pause$ensembl_gene_id)
  keep2             <- !is.na(match_idx)
  bw_gb_filtered    <- bw_gb_filtered[keep2]
  bw_pause_filtered <- bw_pause[match_idx[keep2]]
  names(mcols(bw_pause_filtered))[names(mcols(bw_pause_filtered)) == "ensembl_gene_id"] <- "gene_id"

  stats <- list(
    tsn_after_cutoff = length(bw_tsn),
    gene_body_kept   = length(bw_gb_filtered)
  )
  message("Genes kept after length filter: ", stats$gene_body_kept)

  list(gene_body = bw_gb_filtered, pause = bw_pause_filtered, stats = stats)
}

# Find the valley between the Sharp and Broad peaks in the fkSD histogram.
# Searches for the lowest local minimum between `from` and `to`.
find_valley_threshold <- function(x, from = 20, to = 50, bins = 100, smooth_k = 5) {
  hist_data     <- hist(x, breaks = bins, plot = FALSE)
  mids          <- hist_data$mids
  counts        <- hist_data$counts

  # smooth histogram counts
  smooth_counts <- as.numeric(stats::filter(
    counts,
    rep(1 / smooth_k, smooth_k),
    sides = 2
  ))
  smooth_counts[is.na(smooth_counts)] <- counts[is.na(smooth_counts)]

  candidate_idx <- which(mids > from & mids < to)
  minima_idx    <- findpeaks(-smooth_counts[candidate_idx])[, 2]
  mids[candidate_idx[minima_idx[which.min(smooth_counts[candidate_idx[minima_idx]])]]]
}

# Add quantile group labels (betaGroup, fkMeanGroup, chiGroup) and a
# Sharp/Broad SD classification to the rates data frame produced by STADyUM.
process_sample_rate <- function(df, from = 20, to = 50, bins = 100) {
  df$betaGroup <- cut(
    df$betaAdp,
    breaks = quantile(df$betaAdp, probs = seq(0, 1, 0.2), na.rm = TRUE),
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    include.lowest = TRUE
  )
  df$fkMeanGroup <- cut(
    df$fkMean,
    breaks = quantile(df$fkMean, probs = seq(0, 1, 0.2), na.rm = TRUE),
    labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
    include.lowest = TRUE
  )
  df$chiGroup <- cut(
    df$chi,
    breaks = quantile(df$chi, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE
  )

  df$fkSD   <- sqrt(df$fkVar)
  threshold <- find_valley_threshold(df$fkSD, from, to, bins)
  print(threshold)
  df$sdGroup <- ifelse(df$fkSD <= threshold, "Sharp", "Broad")
  attr(df, "threshold") <- threshold
  return(df)
}


# ============================================================
# MAIN PROCESSING FUNCTION
# ============================================================

estimate_rates_for_species <- function(cfg, result_dir) {

  out_dir <- file.path(result_dir, cfg$out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Read TSN files ---
  tsn1 <- readRDS(cfg$tsn1_file)
  tsn2 <- readRDS(cfg$tsn2_file)

  # --- Keep one upstream TSS per gene ---
  tsn1 <- keep_upstream_tss(tsn1)
  tsn2 <- keep_upstream_tss(tsn2)

  # --- Save TSS distance table (QC only; coordinates are not modified) ---
  dist_df   <- data.frame(Distance = compute_tss_distances(tsn1, tsn2))
  dist_file <- file.path(out_dir, paste0(tolower(cfg$sample1_name), "_",
                                         tolower(cfg$sample2_name), "_tss_dist.txt"))
  write.table(dist_df, dist_file, row.names = TRUE, sep = ",", quote = FALSE)

  # --- Build read-count regions using sample-specific TSSs ---
  tq      <- readRDS(cfg$tq_file)
  region1 <- build_readcount_regions(tsn1, tq@transcripts)
  region2 <- build_readcount_regions(tsn2, tq@transcripts)

  # --- Run STADyUM ---
  rate1 <- estimateTranscriptionRates(
    cfg$bw1_plus, cfg$bw1_minus,
    region1$pause, region1$gene_body, cfg$sample1_name
  )
  rate2 <- estimateTranscriptionRates(
    cfg$bw2_plus, cfg$bw2_minus,
    region2$pause, region2$gene_body, cfg$sample2_name
  )

  # --- Add quantile groups and Sharp/Broad classification ---
  rate1@rates <- process_sample_rate(rates(rate1))
  rate2@rates <- process_sample_rate(rates(rate2))

  # --- Save rate objects ---
  saveRDS(rate1, file.path(out_dir, paste0(tolower(cfg$sample1_name), "_rate.RDS")))
  saveRDS(rate2, file.path(out_dir, paste0(tolower(cfg$sample2_name), "_rate.RDS")))

  message("Done: ", cfg$out_subdir,
          " (", cfg$sample1_name, " and ", cfg$sample2_name, ")")
}


# ============================================================
# RUN
# ============================================================

cfg <- species_configs[[SPECIES]]
estimate_rates_for_species(cfg, result_dir)
