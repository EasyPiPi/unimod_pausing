#### log file ####
log <- file(snakemake@log[[1]], open = "wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)

#### snakemake files ####
bwp_rep1_in <- snakemake@input[["bwp_rep1"]]
bwm_rep1_in <- snakemake@input[["bwm_rep1"]]

bwp_rep2_in <- snakemake@input[["bwp_rep2"]]
bwm_rep2_in <- snakemake@input[["bwm_rep2"]]

gtf_in <- snakemake@input[["gtf"]]

rc_cutoff <- snakemake@params[["rc_cutoff"]]
threads <- snakemake@threads[[1]]

max_tsn_gn_out <- snakemake@output[["max_tsn_gn"]]

root_dir <- "."
result_dir <- file.path(root_dir, "outputs/read_dt")
proseq_in <- list.files(file.path(root_dir, "outputs/bigwig/p5"))
message("Bigwig files being analyzed:")
print(proseq_in)

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_pausing"
# result_dir <- file.path(root_dir, "outputs/read_dt")
#
# gtf_in <- file.path(result_dir, "human_transcript_granges.rds")
#
# bwp_rep1_in <- file.path(root_dir, "ext_data/copro/hg38/GSM4296337_PRO-cap-01-NELFC-AID-untreated-rep1.plus.bw")
# bwm_rep1_in <- file.path(root_dir, "ext_data/copro/hg38/GSM4296337_PRO-cap-01-NELFC-AID-untreated-rep1.minus.bw")
#
# bwp_rep2_in <- file.path(root_dir, "ext_data/copro/hg38/GSM4296339_PRO-cap-03-NELFC-AID-untreated-rep2.plus.bw")
# bwm_rep2_in <- file.path(root_dir, "ext_data/copro/hg38/GSM4296339_PRO-cap-03-NELFC-AID-untreated-rep2.minus.bw")
#
# proseq_in <- list.files(file.path(root_dir, "outputs/bigwig/p5"))
#
# max_tsn_gn_out <- file.path(result_dir, "max_tsn_per_gene_DLD1.rds")
# # read count cutoff for coPRO-cap 5' end read count
# rc_cutoff <- 10
# threads <- 12

#### end of parsing arguments ####
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

gtf <- readRDS(gtf_in)
promoter_grng <- promoters(gtf, upstream = 250, downstream = 250)

# function for read in bigwig files
read_bw <- function(bw_in, strand) {
  bw <- rtracklayer::import.bw(bw_in)
  seqlevelsStyle(bw) <- "Ensembl"
  bw <- bw[bw$score != 0]
  bw <- keepStandardChromosomes(bw, pruning.mode = "coarse")
  strand(bw) <- strand
  if (strand == "-") bw$score <- abs(bw$score)
  bw <- BRGenomics::makeGRangesBRG(bw)
  return(bw)
}

# read in PRO-cap 5p end read counts
bwp_rep1 <- read_bw(bwp_rep1_in, strand = "+")
bwm_rep1 <- read_bw(bwm_rep1_in, strand = "-")

bwp_rep2 <- read_bw(bwp_rep2_in, strand = "+")
bwm_rep2 <- read_bw(bwm_rep2_in, strand = "-")

bw_procap <- BRGenomics::mergeGRangesData(bwp_rep1, bwm_rep1, bwp_rep2, bwm_rep2,
  field = "score",
  multiplex = FALSE,
  makeBRG = TRUE,
  exact_overlaps = FALSE,
  ncores = getOption("mc.cores", threads)
)

# reduce TSS regions for searching PRO-cap signals
promoter_rd <-
  GenomicRanges::reduce(promoter_grng,
    with.revmap = TRUE,
    drop.empty.ranges = TRUE, ignore.strand = FALSE
  )

# A much faster implementation
# ideas from https://stackoverflow.com/questions/60311568/get-the-longest-ranges-per-seqnames
get_max_tsn <- function(copro_5p, promoter) {
  # find overlaps between the promoter region and read counts
  copro_ovp <- findOverlaps(promoter, copro_5p, ignore.strand = FALSE)
  # split read counts into groups based on overlaps
  copro_5p_ovp <- copro_5p[subjectHits(copro_ovp)]
  copro_prom_rc <- split(copro_5p_ovp, queryHits(copro_ovp))
  # split read counts into groups for getting positions with max values
  copro_5p_rc <- splitAsList(copro_5p_ovp$score, queryHits(copro_ovp))
  # get max values with ties, when ties happen, use the most upstream position
  copro_max_tsn_all <- copro_prom_rc[which(copro_5p_rc == max(copro_5p_rc))]
  # get most upstream position according to strand
  if (unique(strand(copro_5p)) == "+") {
    # note the global argument gives the global position across all positions
    copro_end_position <- which.min(start(copro_max_tsn_all), global = TRUE)
  } else if (unique(strand(copro_5p)) == "-") {
    copro_end_position <- which.max(end(copro_max_tsn_all), global = TRUE)
  } else {
    stop("strand is not + or -")
  }
  # concat all positions together. note using "c" with do.call() is too slow
  # while using unlist the group name is not in the metadata frame
  copro_max_tsn_all <- stack(copro_max_tsn_all, index.var = "group")
  # subset and return single position per group
  copro_max_tsn <- copro_max_tsn_all[copro_end_position]

  return(copro_max_tsn)
}

map_max_tsn2gene <- function(bw, promoter_rd) {
  bwp_5p <- bw[strand(bw) == "+"]
  bwm_5p <- bw[strand(bw) == "-"]

  max_tsn_plus <- get_max_tsn(bwp_5p, promoter = promoter_rd)
  max_tsn_minus <- get_max_tsn(bwm_5p, promoter = promoter_rd)

  # This is the max tsn for every promoter
  max_tsn <- sort(c(max_tsn_plus, max_tsn_minus))

  # Map genes names back to promoter tsn
  max_tsn$gene_id <-
    promoter_grng[min(promoter_rd[as.integer(as.vector(max_tsn$group)), ]$revmap)]$gene_id
  # get the maximum tsn for every gene
  max_tsn_gn <- splitAsList(max_tsn$score, max_tsn$gene_id)
  max_tsn_gn <-
    unlist(split(
      max_tsn,
      max_tsn$gene_id
    )[max_tsn_gn == max(max_tsn_gn)])
  # filter out promoters with counts lower than a certain threshold
  max_tsn_gn <- max_tsn_gn[max_tsn_gn$score > rc_cutoff]
  # deduplicate genes with ties in the max tsn, only very few genes
  max_tsn_gn <- max_tsn_gn[!duplicated(max_tsn_gn$gene_id)]
  return(max_tsn_gn)
}

max_tsn_gn <- map_max_tsn2gene(bw_procap, promoter_rd)

# read in other PRO-seq data
proseq_names <- unique(stringr::str_remove(proseq_in, "_plus.bw|_minus.bw"))
# keep selected cell types
proseq_names <- stringr::str_subset(proseq_names, "DLD1")
# select a subst for further analyses
proseq_names <- stringr::str_subset(proseq_names, "NELFC_Auxin|NELFC_Fp|NELFC_NVP2")

proseq_5p_list <- list()

for (proseq_name in proseq_names) {
  bwp <- read_bw(file.path(
    root_dir, "outputs/bigwig/p5",
    paste0(proseq_name, "_plus.bw")
  ), strand = "+")
  bwm <- read_bw(file.path(
    root_dir, "outputs/bigwig/p5",
    paste0(proseq_name, "_minus.bw")
  ), strand = "-")
  # bwm$score <- abs(bwm$score)
  bw <- c(bwp, bwm)
  proseq_5p_list[[proseq_name]] <- bw
}

# rename sample for better visualization
names(proseq_5p_list) <-
  sapply(
    proseq_names,
    function(x) paste0(stringr::str_split(x, "-")[[1]][c(4)], collapse = "-")
  )

# clean up
rm(bwp, bwm, bw)

pro_max_tsn_gn_ls <-
  purrr::map(proseq_5p_list, map_max_tsn2gene, promoter_rd = promoter_rd)

compare_copro_and_proseq <- function(pro_max_tsn_gn) {
  uniset_gn <- intersect(max_tsn_gn$gene_id, pro_max_tsn_gn$gene_id)
  uniset_grng <- punion(max_tsn_gn[uniset_gn],
    pro_max_tsn_gn[uniset_gn],
    fill.gap = TRUE
  )
  return(uniset_grng)
}

compare_max_tsn_gn_ls <- purrr::map(pro_max_tsn_gn_ls, compare_copro_and_proseq)

message("Number of genes are compared between PRO-cap and PRO-seq for max TSN:")
print(lengths(compare_max_tsn_gn_ls))

p <- tibble::tibble(
  sample = names(compare_max_tsn_gn_ls),
  width = purrr::map(compare_max_tsn_gn_ls, width)
) %>%
  tidyr::unnest(cols = c(width)) %>%
  ggplot2::ggplot(ggplot2::aes(x = log2(width))) +
  ggplot2::geom_histogram(color = "black", fill = "white") +
  ggplot2::facet_wrap(. ~ sample, nrow = 4) +
  ggplot2::geom_vline(ggplot2::aes(xintercept = 7, color = "red"),
    linetype = "dashed", show.legend = FALSE
  ) +
  ggplot2::labs(x = "log2(Distance to PRO-cap max TSN)", y = "Number of genes") +
  cowplot::theme_cowplot()

ggplot2::ggsave(
  file.path(
    result_dir,
    "proseq_max_tsn_distance_to_procap_max_tsn_DLD1.png"
  ),
  plot = p, width = 16, height = 10
)

saveRDS(max_tsn_gn, max_tsn_gn_out)

# go 100 bp upstream and downstream to check 5p end signal
max_tsn_rn <- resize(max_tsn_gn, width = 101, fix = "center")
max_tsn_rn <- sort(max_tsn_rn, by = ~score, decreasing = TRUE)

sm_5p <- ScoreMatrixList(c(list("PRO-cap" = bw_procap), proseq_5p_list),
  windows = max_tsn_rn,
  strand.aware = TRUE, weight.col = "score", cores = 2
)

save_png <- function(file_name, plot_fun, width = 2000, height = 600) {
  png(filename = file_name, width = width, height = height, units = "px")
  plot_fun
  invisible(dev.off())
}

# heatmap
# http://zvfak.blogspot.com/2015/10/summary-of-new-features-of-genomation_15.html
# exclude top 2% for better visualization
# skip plotting "PROseq-K562-chivu-treated-PE" cause too few reads in it
save_png(
  file_name = paste0(file.path(result_dir, "procap_5p_read_count_heatmap"), ".png"),
  plot_fun = multiHeatMatrix(sm_5p,
    xcoords = c(-50, 50),
    winsorize = c(0, 98), order = TRUE,
    xlab = "bases around TSS"
  )
)

# metaplot
paired_color <- RColorBrewer::brewer.pal(name = "Paired", n = length(proseq_names))

meta_plot <- function(sm, smlcolors, meta.rescale = TRUE, xcoords = c(-50, 50),
                      dispersion = NULL, ylab = "average score",
                      kmin = NULL, kmax = NULL, legend = "topright") {
  plotMeta(sm,
    xcoords = xcoords, meta.rescale = meta.rescale,
    line.col = smlcolors, dispersion = dispersion,
    xlab = "bases around TSS", ylab = ylab,
    dispersion.col = rainbow(length(sm), alpha = 0.5)
  )
  legend(legend, names(sm), lty = c(1, 1), lwd = c(2.5, 2.5), col = smlcolors, cex = 1.5)
  if (!(is.null(kmin) & is.null(kmax))) {
    abline(v = c(kmin, kmax), col = c("red", "red"), lty = c(2, 2), lwd = c(2, 2))
  }
}

save_png(
  file_name = paste0(file.path(result_dir, "procap_5p_read_count_scaled_metaplot"), ".png"),
  plot_fun = meta_plot(sm_5p, c("purple", paired_color)),
  width = 1000, height = 600
)

save_png(
  file_name = paste0(file.path(result_dir, "procap_5p_read_count_noscale_metaplot"), ".png"),
  plot_fun = meta_plot(sm_5p, c("purple", paired_color), meta.rescale = FALSE),
  width = 1000, height = 600
)

# calculate 3p density
proseq_3p_list <- list()

for (proseq_name in proseq_names) {
  bwp <- read_bw(file.path(
    root_dir, "outputs/bigwig/p3",
    paste0(proseq_name, "_plus.bw")
  ), strand = "+")
  bwm <- read_bw(file.path(
    root_dir, "outputs/bigwig/p3",
    paste0(proseq_name, "_minus.bw")
  ), strand = "-")
  # bwm$score <- abs(bwm$score)
  bw <- c(bwp, bwm)
  proseq_3p_list[[proseq_name]] <- bw
}

rm(bwp, bwm, bw)

# rename sample for better visualization
names(proseq_3p_list) <- sapply(
  proseq_names,
  function(x) paste0(str_split(x, "-")[[1]][c(4)], collapse = "-")
)

max_tsn_rn_3p <- resize(max_tsn_rn, width = 501, fix = "center")

sm_3p <- ScoreMatrixList(proseq_3p_list,
  windows = max_tsn_rn_3p,
  strand.aware = TRUE, weight.col = "score", cores = 2
)

# heatmap
# exclude top 1% for better visualization
save_png(
  file_name = paste0(file.path(
    result_dir,
    "procap_3p_read_count_heatmap"
  ), ".png"),
  plot_fun = multiHeatMatrix(sm_3p,
    xcoords = c(-250, 250),
    winsorize = c(0, 99), order = TRUE, xlab = "bases around TSS"
  )
)

# metaplot
save_png(
  file_name = paste0(file.path(
    result_dir,
    "procap_3p_read_count_scaled_metaplot"
  ), ".png"),
  plot_fun = meta_plot(sm_3p, paired_color,
    xcoords = c(-250, 250), legend = "topleft"
  ),
  width = 1000, height = 600
)

save_png(
  file_name = paste0(file.path(result_dir, "procap_3p_read_count_noscale_metaplot"), ".png"),
  plot_fun = meta_plot(sm_3p, paired_color,
    xcoords = c(-250, 250), meta.rescale = FALSE, legend = "topleft"
  ),
  width = 1000, height = 600
)

# summarize distances for pausing peaks
pause_sites <-
  sapply(lapply(
    sm_3p@.Data,
    function(x) apply(x[, 251:501], 1, which.max)
  ), summary)
colnames(pause_sites) <- names(sm_3p)

message("Some statistics for the distance of maximum 3' end signals from the TSN")
print(pause_sites)

# calculate k0, which will be used in Poisson based model adaptation
pause_regions <- lapply(sm_3p@.Data, function(x) x[, 251:501])
names(pause_regions) <- names(sm_3p)

kmin <- 21
kmax <- 200

save_png(
  file_name = paste0(file.path(result_dir, "procap_3p_read_count_noscale_metaplot_wt_kmin_kmax"), ".png"),
  plot_fun = meta_plot(sm_3p, paired_color,
    xcoords = c(-250, 250),
    meta.rescale = FALSE, kmin = kmin, kmax = kmax,
    legend = "topleft"
  ),
  width = 1000, height = 600
)
