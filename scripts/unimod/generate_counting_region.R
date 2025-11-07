#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

#### snakemake files ####
gtf_in <- snakemake@input[["gtf"]]
tsn_in <- snakemake@input[["tsn"]]

tss_length <- snakemake@params[["tss_length"]]
tts_length <- snakemake@params[["tts_length"]] # parameter m
gb_min_length <- snakemake@params[["gb_min_length"]]
gb_max_length <- snakemake@params[["gb_max_length"]]
dist_to_tss <- snakemake@params[["dist_to_tss"]] 

grng_out <- snakemake@output[["grng"]]
tss_out <- snakemake@output[["tss"]]

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_human"

# gtf_in <- file.path(root_dir, "outputs/read_dt/human_transcript_granges.rds")
# tsn_in <- file.path(root_dir, "outputs/read_dt/copro_max_tsn_per_gene.rds")

# tss_length <- 250 # parameter k
# tts_length <- 250 # parameter m
# gb_min_length <- 6e3 # filter out genes shorter than this length
# # when gene is longer than this length, only use length up to this cutoff for read counting
# gb_max_length <- 9e4
# dist_to_tss <- 1000

# tss_out <- file.path(root_dir, "outputs/read_dt/tss.bed")
# grng_out <- file.path(root_dir, "outputs/read_dt/granges_for_read_counting.RData")
#### end of parsing arguments ####

#### generate regions for read counting ####
tsn <- readRDS(tsn_in)
tsn$group <- NULL

# get TSNs downstream regions for pause peak
bw_pause <- promoters(tsn, upstream = 0, downstream = tss_length)
bw_pause$score <- NULL

# # get gene body region by fixed length
# bw_gb <- bw_pause %>% plyranges::shift_downstream(gb_start)
# bw_gb <- bw_gb %>% plyranges::anchor_5p() %>% mutate(width = gb_length)

# get TTS by using annotated gene regions
txgrng <- readRDS(gtf_in)

# export tss coordinates for other use
tssgrng <- promoters(txgrng, upstream = 0, downstream = 1)
tssgrng$score <- 0
tssgrng$name <- tssgrng$transcript_id
export(tssgrng, tss_out, format = "BED")

txgrng <- txgrng[txgrng$gene_id %in% bw_pause$gene_id]
gngrng <- txgrng %>%
  plyranges::group_by(gene_id) %>%
  plyranges::reduce_ranges_directed() %>%
  sort()

seqlevels(gngrng) <- seqlevelsInUse(gngrng)

bw_tts <- gngrng %>% plyranges::anchor_3p() %>% mutate(width = tts_length)

# get gene body region by pause and termination sites
bw_pause_end <- bw_pause %>% plyranges::anchor_3p() %>% mutate(width = 1)
bw_tts_end  <- bw_tts %>% plyranges::anchor_5p() %>% mutate(width = 1)

seqlevels(bw_pause_end) <- seqlevels(bw_tts_end)

bw_pause_end <- bw_pause_end[order(bw_pause_end$gene_id)]
bw_tts_end <- bw_tts_end[order(bw_tts_end$gene_id)]

bw_gb <- punion(bw_pause_end, bw_tts_end, fill.gap = TRUE)
bw_gb$gene_id <- bw_pause_end$gene_id

bw_gb_filtered <- bw_gb[width(bw_gb) > gb_min_length]

# trim either end to avoid pausing and termination peaks
bw_gb_filtered <- bw_gb_filtered - dist_to_tss

# cap gene length up tp certain cutoff 
bw_gb_filtered_long <- bw_gb_filtered[width(bw_gb_filtered) >= gb_max_length]
bw_gb_filtered <- bw_gb_filtered[width(bw_gb_filtered) < gb_max_length]

bw_gb_filtered_long_pl <- bw_gb_filtered_long[strand(bw_gb_filtered_long) == "+"]
bw_gb_filtered_long_mn <- bw_gb_filtered_long[strand(bw_gb_filtered_long) == "-"]

end(bw_gb_filtered_long_pl) <- start(bw_gb_filtered_long_pl) + gb_max_length - 1
start(bw_gb_filtered_long_mn) <- end(bw_gb_filtered_long_mn) - gb_max_length + 1

bw_gb_filtered <- c(bw_gb_filtered, bw_gb_filtered_long_pl, bw_gb_filtered_long_mn)
rm(bw_gb_filtered_long, bw_gb_filtered_long_pl, bw_gb_filtered_long_mn)

bw_pause_filtered <-
  bw_pause[bw_pause$gene_id %in% bw_gb_filtered$gene_id]
bw_tts_filtered <-
  bw_tts[bw_tts$gene_id %in% bw_gb_filtered$gene_id]

# # make sure gene body doesn't exceed TTS
# bw_gb_end <- bw_gb %>% plyranges::anchor_3p() %>% mutate(width = 1)
# 
# bw_gb_end$indicator <- (start(bw_tts) > start(bw_gb_end))
# 
# gn_filter <-
#   ifelse((strand(bw_gb_end) == "+" & bw_gb_end$indicator) |
#            (strand(bw_gb_end) == "-" & !bw_gb_end$indicator),
#        TRUE, FALSE)
# 
# # get three regions for read counting
# bw_gb_filtered <- bw_gb[gn_filter]
# 
# bw_tts_filtered <-
#   bw_tts[bw_tts$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id] %>%
#   plyranges::anchor_3p() %>%
#   mutate(width = tts_length)
# 
# bw_pause_filtered <- bw_pause[bw_pause$ensembl_gene_id %in% bw_gb_filtered$ensembl_gene_id]

count_grng <- sort(c(bw_pause_filtered, bw_gb_filtered, bw_tts_filtered))
save(bw_pause_filtered, bw_gb_filtered, bw_tts_filtered, count_grng, file = grng_out)
