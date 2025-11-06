#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
non_olp_gn_in <- snakemake@input[["non_olp_gn"]]
gtf_in <- snakemake@input[["gtf"]]

gtf_out <- snakemake@output[["gtf"]]

#### load packages ####
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_human"
# result_dir <- file.path(root_dir, "outputs/read_dt")

# non_olp_gn_in <- file.path(root_dir, "ext_data/annotation/non_overlapping_coding_genes.csv")
# gtf_in <- file.path(root_dir, "ext_data/annotation/transcript.gtf")

# gtf_out <- file.path(result_dir, "human_transcript_granges.rds")

#### end of parsing arguments ####
# read non-overlapping genes
non_olp_gn <- read.csv(non_olp_gn_in)
# get human transcripts
gtf <- gffToGRanges(gtf_in)
gtf <-
  gtf[gtf$type == "transcript" & gtf$gene_id %in% non_olp_gn$ensembl_gene_id]
gtf <- keepStandardChromosomes(gtf, pruning.mode="coarse")

saveRDS(gtf, gtf_out)

