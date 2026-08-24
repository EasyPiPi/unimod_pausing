library(tidyverse)
library(plyranges)
library(genomation)
library(STADyUM)
library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# ── Paths ──────────────────────────────────────────────────────────────────────
root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning")

# ── Helper functions ───────────────────────────────────────────────────────────

get_1d_signal_from_contact_df <- function(contacts_df_path, seq_info) {
  contacts_df <- read.table(
    contacts_df_path, sep = "\t", header = FALSE,
    col.names = c("chrom", "start", "end", "name", "score", "strand", "tss"),
    stringsAsFactors = FALSE
  )
  gr <- makeGRangesFromDataFrame(
    contacts_df,
    seqnames.field = "chrom", start.field = "start", end.field = "end",
    strand.field = "strand", keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )
  seqinfo(gr) <- seq_info[seqnames(seqinfo(gr))]
  gr
}

loess_smooth <- function(y) {
  x_pos <- seq(from = -1000, to = 1000, length.out = 2000)
  fit <- loess(y ~ x_pos, span = 0.05)
  predict(fit, x_pos)
}

get_tid_from_stadyum_obj <- function(stadyum_obj) {
  tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1)
  tid <- resize(tsn, width = 2000, fix = "center")
  rate_df <- rates(stadyum_obj)
  matched_idx <- match(rate_df$geneId, tid$gene_id)
  tid[matched_idx[!is.na(matched_idx)], ]
}

get_matrix_from_nucleosome_signal <- function(rate_object, nucleosome_signal_in,
                                               seq_info, smooth_func = loess_smooth) {
  gr <- get_1d_signal_from_contact_df(nucleosome_signal_in, seq_info)
  tid <- get_tid_from_stadyum_obj(rate_object)
  score_matrix <- ScoreMatrix(gr, windows = tid,
    weight.col = "score", strand.aware = TRUE)
  score_matrix <- as.matrix(score_matrix@.Data)
  score_matrix_smooth <- t(apply(score_matrix, 1, smooth_func))
  list(raw = score_matrix, smooth = score_matrix_smooth, tid)
}

# Process one cell type through the full nucleosome positioning workflow.
process_cell_type <- function(rate_in, signal_in, seqinfo, cell_type,
                               out_dir, fix_na_sdgroup = FALSE) {
  rate      <- readRDS(rate_in)
  ns_matrix <- get_matrix_from_nucleosome_signal(rate, signal_in, seqinfo)
  rate_df   <- rate@rates

  # Baboon-specific: sdGroup has NAs that must be set to "Sharp" before saving
  if (fix_na_sdgroup) {
    rate_df[is.na(rate_df$sdGroup), "sdGroup"] <- "Sharp"
  }

  plus1_mean_200 <- rowMeans(ns_matrix$raw[, 1001:1200])
  rate_df[[paste0(cell_type, "_nuc_mean_200")]] <- plus1_mean_200

  ns_matrix$rate_df <- rate_df
  saveRDS(ns_matrix, file.path(out_dir, paste0(cell_type, "_ns_matrix.RDS")))
  invisible(ns_matrix)
}

# ── Seqinfo ────────────────────────────────────────────────────────────────────

# Human: from TxDb annotation package
seqlevelsStyle(TxDb.Hsapiens.UCSC.hg38.knownGene) <- "Ensembl"
seqinfo_human <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Rhesus: from genome FASTA index
seqinfo_rhesus <- seqinfo(scanFaIndex(
  file.path(.paths$genome, "rheMac10.fa")
))
seqlevelsStyle(seqinfo_rhesus) <- "Ensembl"

# Baboon: from genome FASTA index
seqinfo_baboon <- seqinfo(scanFaIndex(
  file.path(.paths$genome, "papAnu4.fa")
))
seqlevelsStyle(seqinfo_baboon) <- "Ensembl"

# ── Species configuration ──────────────────────────────────────────────────────
# fix_na_sdgroup: baboon has NA values in sdGroup that must be filled with "Sharp"

species_configs <- list(
  human = list(
    seqinfo        = seqinfo_human,
    fix_na_sdgroup = FALSE,
    cell_types = list(
      cd4  = list(
        rate_in   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd4_rate.RDS"),
        signal_in = file.path(.paths$microc_1d, "human_cd4_1d_signal.bed")
      ),
      cd14 = list(
        rate_in   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd14_rate.RDS"),
        signal_in = file.path(.paths$microc_1d, "human_cd14_1d_signal.bed")
      )
    )
  ),
  rhesus = list(
    seqinfo        = seqinfo_rhesus,
    fix_na_sdgroup = FALSE,
    cell_types = list(
      cd4  = list(
        rate_in   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/rhesus/cd4_rate.RDS"),
        signal_in = file.path(.paths$microc_1d, "rhesus_cd4_1d_signal.bed")
      ),
      cd14 = list(
        rate_in   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/rhesus/cd14_rate.RDS"),
        signal_in = file.path(.paths$microc_1d, "rhesus_cd14_1d_signal.bed")
      )
    )
  ),
  baboon = list(
    seqinfo        = seqinfo_baboon,
    fix_na_sdgroup = TRUE,
    cell_types = list(
      cd4  = list(
        rate_in   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/baboon/cd4_rate.RDS"),
        signal_in = file.path(.paths$microc_1d, "baboon_cd4_1d_signal.bed")
      ),
      cd14 = list(
        rate_in   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/baboon/cd14_rate.RDS"),
        signal_in = file.path(.paths$microc_1d, "baboon_cd14_1d_signal.bed")
      )
    )
  )
)

# ── Run ────────────────────────────────────────────────────────────────────────

for (species in names(species_configs)) {
  cfg     <- species_configs[[species]]
  out_dir <- file.path(result_dir, species)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (cell_type in names(cfg$cell_types)) {
    ct <- cfg$cell_types[[cell_type]]
    message("Processing ", species, " / ", cell_type, " ...")
    process_cell_type(
      rate_in        = ct$rate_in,
      signal_in      = ct$signal_in,
      seqinfo        = cfg$seqinfo,
      cell_type      = cell_type,
      out_dir        = out_dir,
      fix_na_sdgroup = cfg$fix_na_sdgroup
    )
  }
}

message("Done.")
