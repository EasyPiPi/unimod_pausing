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

#### testing files ####
# root_dir <- "."
# table_dir <- file.path(root_dir, "outputs/simulation/tables/lrt_pause_escape")
# figure_dir <- file.path(root_dir, "outputs/simulation/figures/lrt_pause_escape")
# lambda <- "lrt_high"

#### end of parsing arguments ####
bar_colors <- c("#1868B2", "#F3A332")

kmin <- 1
kmax <- 200

meta_in <- "metadata/simulation_params_lrt.csv"

rate_tbls <- read_csv(meta_in)
colnames(rate_tbls) <-
  c("id", "k", "ksd", "a", "b", "z", "t", "n", "s", "h", "l")

rate_tbls$bws <-
  map(
    file.path(table_dir, paste0(rate_tbls$id, "_", lambda, ".bw")),
    rtracklayer::import
  )

rate_subset <- rate_tbls %>%
  filter(a == 1) %>%
  arrange(b)

# Build DataTrack objects from GRanges and plot
bw_list <- rate_subset$bws
names(bw_list) <- rate_subset$b

# choose a chromosome to plot (common across all, or fallback to first)
all_chrs <- lapply(bw_list, function(gr) unique(as.character(seqnames(gr))))
common_chr <- Reduce(intersect, all_chrs)
chr <- if (length(common_chr) > 0) common_chr[1] else all_chrs[[1]][1]

# restrict to chosen chromosome
bw_list_chr <- lapply(bw_list, function(gr) gr[as.character(seqnames(gr)) == chr])

# compute window
from <- min(unlist(lapply(bw_list_chr, function(gr) if (length(gr)) min(start(gr)) else NA)), na.rm = TRUE)
to <- max(unlist(lapply(bw_list_chr, function(gr) if (length(gr)) max(end(gr)) else NA)), na.rm = TRUE)

# compute global y-axis limits
all_scores <- unlist(lapply(bw_list_chr, function(gr) if (length(gr)) gr$score else numeric(0)))
ylim <- c(-0.05 * max(all_scores, na.rm = TRUE), max(all_scores, na.rm = TRUE))

# Create gradient colors from #1868B2
n_tracks <- length(bw_list_chr)
gradient_colors <- colorRampPalette(c("#1868B2", "#87CEEB"))(n_tracks)

tracks <- Map(function(gr, nm, col) {
  DataTrack(
    range = gr, data = gr$score, type = "hist", name = nm, ylim = ylim,
    fill = col, col = col, col.histogram = col, fill.histogram = col,
    background.title = "white", col.title = "grey85", col.axis = "grey85",
    cex.axis = 0.7
  )
}, bw_list_chr, names(bw_list_chr), gradient_colors)

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

# Call the function to save the plot
# Whole Gene
save_tracks_plot(
  tracks, chr, from, to, lambda, figure_dir,
  width = 6, height = 6
)
# TSS Region
save_tracks_plot(
  tracks, chr, from + kmin - 81,
  from + kmax - 21, lambda, figure_dir,
  width = 1.5, height = 6
)

tracks2 <- Map(function(gr, nm, col) {
  DataTrack(
    range = gr, data = gr$score, type = "hist", name = nm, ylim = c(-0.025, 0.5),
    fill = col, col = col, col.histogram = col, fill.histogram = col,
    background.title = "white", col.title = "grey85", col.axis = "grey85",
    cex.axis = 0.7
  )
}, bw_list_chr, names(bw_list_chr), gradient_colors)
# Gene Body
save_tracks_plot(
  tracks2, chr, from + 500 - 1,
  to, lambda, figure_dir,
  width = 6, height = 6
)