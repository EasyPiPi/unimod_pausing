library(tidyverse)
library(STADyUM)

# ── Paths ─────────────────────────────────────────────────────────────────────

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

lrt_dir    <- file.path(.paths$outputs, "publish/betweenCellType/LRT/human")
rate_dir   <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human")
tss_dir    <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human")
result_dir <- file.path(.paths$outputs, "publish/betweenCellType/EP")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parameters ────────────────────────────────────────────────────────────────

# Sequencing depth normalisation: CD14 / CD4 (matches 1_1_LRT_unified.R human config)
scale_factor  <- 38359083 / 33294862
lfc_threshold <- 2

input_paths <- list(
  lrt      = file.path(lrt_dir,  "lrt_result.txt"),
  cd4_rate = file.path(rate_dir, "cd4_rate.RDS"),
  cd14_rate = file.path(rate_dir, "cd14_rate.RDS"),
  cd4_tss  = file.path(tss_dir,  "cd4_tss_sort_1000bp_extend.bed"),
  cd14_tss = file.path(tss_dir,  "cd14_tss_sort_1000bp_extend.bed")
)

# ── Input file checks ─────────────────────────────────────────────────────────

for (nm in names(input_paths)) {
  if (!file.exists(input_paths[[nm]])) {
    stop("Required input file not found: ", input_paths[[nm]])
  }
}

# ── Helper functions ──────────────────────────────────────────────────────────

# Reads a TSS BED file and returns a two-column table: geneId → file_name.
# file_name is built as <chromosome_no_chr_prefix>_<center>.csv.
get_tss_filename_table <- function(tss_bed_path) {
  tss <- read.table(tss_bed_path, sep = "\t")
  tss$chromosome <- sub("^chr", "", tss$V1)
  tss$center     <- (tss$V2 + tss$V3) / 2
  tss$file_name  <- paste0(tss$chromosome, "_", tss$center, ".csv")
  data.frame(geneId = tss$V4, file_name = tss$file_name, stringsAsFactors = FALSE)
}

# Reads a STADyUM rate RDS and returns a two-column table: geneId → <count_col>.
# pauseCount = sum(actualPauseSiteCounts) per gene.
read_rate_pause_count <- function(rate_rds_path, count_col) {
  rate_df <- rates(readRDS(rate_rds_path))
  rate_df$pauseCount <- vapply(rate_df$actualPauseSiteCounts, sum, numeric(1))
  result <- rate_df[, c("geneId", "pauseCount")]
  colnames(result)[2] <- count_col
  result
}

# Reports how many shared genes have identical file_name derived from CD4 and CD14 TSS BEDs.
# Only used for sanity checking; the final output always uses the CD4 file_name.
check_tss_filename_consistency <- function(cd4_tss_table, cd14_tss_table) {
  merged <- inner_join(
    cd4_tss_table,
    cd14_tss_table,
    by     = "geneId",
    suffix = c("_cd4", "_cd14")
  )
  n_match <- sum(merged$file_name_cd4 == merged$file_name_cd14, na.rm = TRUE)
  n_total <- nrow(merged)
  message(sprintf(
    "TSS file_name consistency: %d / %d shared genes have identical CD4 and CD14 file_name.",
    n_match, n_total
  ))
  if (n_match < n_total) {
    warning(n_total - n_match, " gene(s) have mismatched CD4 / CD14 file_name.")
  }
  invisible(merged)
}

# Adds pauseCount_lfc = log2(CD14 * scale_factor / CD4). Warns if any value is zero.
add_pause_count_lfc <- function(df, cd4_col, cd14_col, scale_factor) {
  zeros <- sum(df[[cd4_col]] == 0 | df[[cd14_col]] == 0, na.rm = TRUE)
  if (zeros > 0) {
    warning(zeros, " gene(s) have zero pauseCount in CD4 or CD14; log2 fold-change will be ±Inf or NaN.")
  }
  df$pauseCount_lfc <- log2(df[[cd14_col]] * scale_factor / df[[cd4_col]])
  df
}

# Assigns pauseCount_lfc_group: Down / Others / Up based on lfc_threshold.
add_pause_count_lfc_group <- function(df, lfc_threshold) {
  df$pauseCount_lfc_group <- factor(
    case_when(
      df$pauseCount_lfc <= -lfc_threshold ~ "Down",
      df$pauseCount_lfc >=  lfc_threshold ~ "Up",
      TRUE                                ~ "Others"
    ),
    levels = c("Down", "Others", "Up")
  )
  df
}

# ── Load inputs ───────────────────────────────────────────────────────────────

human_lrt      <- read.table(input_paths$lrt, sep = "\t", header = TRUE)
cd4_tss_table  <- get_tss_filename_table(input_paths$cd4_tss)
cd14_tss_table <- get_tss_filename_table(input_paths$cd14_tss)

message(sprintf("LRT result loaded: %d genes.", nrow(human_lrt)))

# ── Sanity: TSS file_name consistency across cell types ───────────────────────

check_tss_filename_consistency(cd4_tss_table, cd14_tss_table)

# ── Join CD4 file_name into LRT result ────────────────────────────────────────

human_lrt <- left_join(human_lrt, cd4_tss_table, by = "geneId")

n_matched    <- sum(!is.na(human_lrt$file_name))
n_missing_fn <- nrow(human_lrt) - n_matched
message(sprintf("%d / %d LRT genes have a matched file_name.", n_matched, nrow(human_lrt)))
if (n_missing_fn > 0) {
  warning(n_missing_fn, " LRT gene(s) have missing file_name after join.")
}

# ── Compute and join pauseCount for each cell type ────────────────────────────

cd4_pause  <- read_rate_pause_count(input_paths$cd4_rate,  "CD4_pauseCount")
cd14_pause <- read_rate_pause_count(input_paths$cd14_rate, "CD14_pauseCount")

human_lrt <- left_join(human_lrt, cd4_pause,  by = "geneId")
human_lrt <- left_join(human_lrt, cd14_pause, by = "geneId")

# ── Compute pauseCount LFC and group label ────────────────────────────────────

human_lrt <- add_pause_count_lfc(human_lrt, "CD4_pauseCount", "CD14_pauseCount", scale_factor)
human_lrt <- add_pause_count_lfc_group(human_lrt, lfc_threshold)

# ── Write output ──────────────────────────────────────────────────────────────

out_path <- file.path(result_dir, "human_lrt_withFilename.txt")
write.table(human_lrt, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: ", out_path)
