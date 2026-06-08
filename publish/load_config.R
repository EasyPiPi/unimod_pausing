# load_config.R — defines .paths for all R scripts in the publish pipeline.
#
# Set root_dir before sourcing, then .paths is available in the calling script.
# Override large-data paths via environment variables if they live elsewhere:
#   ENCODE_BAM_ROOT, ENCODE_BW_ROOT, GENOME_ROOT, MICROC_1D_ROOT

stopifnot(exists("root_dir"))

.e <- function(var, default) {
  v <- Sys.getenv(var, "")
  if (nzchar(v)) normalizePath(v, mustWork = FALSE) else default
}

.paths <- list(
  root       = root_dir,
  data       = file.path(root_dir, "data", "publish"),
  outputs    = file.path(root_dir, "outputs"),
  encode_bam = .e("ENCODE_BAM_ROOT", file.path(root_dir, "data", "result2", "encode", "bam")),
  encode_bw  = .e("ENCODE_BW_ROOT",  file.path(root_dir, "data", "result2", "encode", "bw_clean")),
  genome     = .e("GENOME_ROOT",     file.path(root_dir, "codes", "micro-c")),
  microc_1d  = .e("MICROC_1D_ROOT",  file.path(root_dir, "data", "micro-c"))
)

rm(.e)
