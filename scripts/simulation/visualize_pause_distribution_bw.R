#### log file ####
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
table_dir <- snakemake@params[["table_dir"]]
figure_dir <- snakemake@params[["figure_dir"]]
lambda <- snakemake@params[["lambda_exp"]]

#### load packages ####
library(tidyverse)
library(rtracklayer)
library(Gviz)

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_pausing"

# lambda <- "lrt_high"

# table_dir <-
#   file.path(root_dir, "outputs/simulation/tables/lrt_pause_distribution")
# figure_dir <-
#   file.path(root_dir, "outputs/simulation/figures/lrt_pause_distribution", lambda)

#### end of parsing arguments ####
bar_colors <- c("#1868B2", "#F3A332")

walk(
  c(file.path(figure_dir, "k"), file.path(figure_dir, "ksd")),
  dir.create,
  showWarnings = FALSE, recursive = TRUE
)

kmin <- 1
kmax <- 200

file_names <- list.files(table_dir)
# only use high coverage setting for LRT
file_names <- str_subset(file_names, paste0(lambda, ".bw"))

rate_tbls <- tibble(
  id = file_names,
  k = as.numeric(str_extract(file_names, "(?<=k)\\d+(?=ksd)")),
  ksd = as.numeric(str_extract(file_names, "(?<=ksd)\\d+(?=kmin)")),
  a = rep(1, length(file_names)),
  b = rep(1, length(file_names)),
  z = rep(2000, length(file_names)),
  t = rep(40, length(file_names)),
  n = rep(20000, length(file_names)),
  s = rep(33, length(file_names)),
  h = rep(17, length(file_names)),
  l = rep(1950, length(file_names))
)

rate_tbls$bws <-
  map(
    file.path(table_dir, rate_tbls$id),
    rtracklayer::import
  )

make_gviz_tracks <- function(rate_tbl) {
  # Build DataTrack objects from GRanges and plot
  bw_list <- rate_tbl$bws
  names(bw_list) <- rate_tbl$ksd

  # choose a chromosome to plot (common across all, or fallback to first)
  all_chrs <- lapply(bw_list, function(gr) unique(as.character(seqnames(gr))))
  common_chr <- Reduce(intersect, all_chrs)
  chr <- if (length(common_chr) > 0) common_chr[1] else all_chrs[[1]][1]

  # restrict to chosen chromosome
  bw_list_chr <- lapply(bw_list, function(gr) gr[as.character(seqnames(gr)) == chr])

  # compute window
  from <-
    min(unlist(lapply(bw_list_chr, function(gr) if (length(gr)) min(start(gr)) else NA)), na.rm = TRUE)
  to <-
    max(unlist(lapply(bw_list_chr, function(gr) if (length(gr)) max(end(gr)) else NA)), na.rm = TRUE)

  # compute global y-axis limits
  all_scores <-
    unlist(lapply(bw_list_chr, function(gr) if (length(gr)) gr$score else numeric(0)))
  ylim <- c(-0.05 * max(all_scores, na.rm = TRUE), max(all_scores, na.rm = TRUE))

  # Create gradient colors for "#F3A332"
  n_tracks <- length(bw_list_chr)
  gradient_colors <- colorRampPalette(c("#F3A332", "#FDE3C1"))(n_tracks)

  tracks <- Map(function(gr, nm, col) {
    DataTrack(
      range = gr, data = gr$score, type = "hist", name = nm, ylim = ylim,
      fill = col, col = col, col.histogram = col, fill.histogram = col,
      background.title = "white", col.title = "grey85", col.axis = "grey85",
      cex.axis = 0.7
    )
  }, bw_list_chr, names(bw_list_chr), gradient_colors)

  return(list("tracks" = tracks, "chr" = chr, "from" = from, "to" = to))
}

save_tracks_plot <-
  function(tracks, chr, from, to, lambda, figure_dir, width, height) {
    out_file <- file.path(
      figure_dir,
      paste0("tracks_", chr, "_", lambda, "_", from, "-", to, ".pdf")
    )

    pdf(out_file, width = width, height = height) # 1800/200=9in, 1200/200=6in
    plotTracks(
      c(list(GenomeAxisTrack()), tracks),
      chromosome = chr,
      from = from,
      to = to
    )
    dev.off()
  }

# subet to k = 80 for visualization
rate_subset_k <- rate_tbls %>%
  filter(k == 80, ksd %in% c(10, 20, 30, 35, 40, 45, 50, 60, 70)) %>%
  arrange(ksd)

tracks_k <- make_gviz_tracks(rate_subset_k)

# Call the function to save the plot
# Whole Gene
save_tracks_plot(
  tracks_k$tracks, tracks_k$chr, tracks_k$from, tracks_k$to, lambda,
  file.path(figure_dir, "k"),
  width = 6, height = 6
)
# Pause Region
save_tracks_plot(
  tracks_k$tracks, tracks_k$chr, tracks_k$from + kmin - 1,
  tracks_k$from + kmax - 1, lambda,
  file.path(figure_dir, "k"),
  width = 2.5, height = 6
)

rate_subset_ksd <- rate_tbls %>%
  filter(ksd == 20, k %in% c(30, 40, 50, 55, 60, 65, 70, 80, 90)) %>%
  arrange(k)

tracks_ksd <- make_gviz_tracks(rate_subset_ksd)

# Call the function to save the plot
# Whole Gene
save_tracks_plot(
  tracks_ksd$tracks, tracks_ksd$chr, tracks_ksd$from, tracks_ksd$to, lambda,
  file.path(figure_dir, "ksd"),
  width = 6, height = 6
)
# Pause Region
save_tracks_plot(
  tracks_ksd$tracks, tracks_ksd$chr, tracks_ksd$from + kmin - 1,
  tracks_ksd$from + kmax - 1, lambda,
  file.path(figure_dir, "ksd"),
  width = 2.5, height = 6
)