#### note ####
# Generate visualization for select gene with simulated data

#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(GenomicRanges)
library(rtracklayer)
library(Gviz)

#### snakemake files ####

#### testing files ####
root_dir <- "~/Desktop/github/unimod_human"

grng_in <- file.path(root_dir, "outputs/read_dt/human_transcript_granges.rds")

bwp_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_plus.bw")
bwm_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_minus.bw")

# bwp_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-vihervaara-control-SE_plus.bw")
# bwm_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-vihervaara-control-SE_minus.bw")

# bwp_sim_in <- file.path(root_dir, "outputs/simulation/data/k50l19700m250r20a1b5g1z1000t30n5000s50h20.plus.bw")
# bwm_sim_in <- file.path(root_dir, "outputs/simulation/data/k50l19700m250r20a1b5g1z1000t30n5000s50h20.minus.bw")

# simulation from SimPolv2
bwp_sim_in <- file.path(root_dir, "ext_data/simPolv2/results/750/combined_cell_data.csv")

result_dir <-
  file.path(root_dir, "outputs/simulation/figures", "PROseq-K562-dukler-control-SE")

# result_dir <-
#   file.path(root_dir, "outputs/simulation/figures", "PROseq-K562-vihervaara-control-SE")

#### end of parsing arguments ####
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

gn_grng <- readRDS(grng_in)
seqlevelsStyle(gn_grng) <- "UCSC"
tx_grng <- gn_grng[gn_grng$gene_name == "DNAJA1"]


bwp_p3 <- import(bwp_p3_in)
bwp_p3 <- keepStandardChromosomes(bwp_p3, pruning.mode = "coarse")
seqlevelsStyle(bwp_p3) <- "UCSC"

# bwp_sim <- import(bwp_sim_in)

bwp_sim <- read.csv(bwp_sim_in, col.names = "score")

axis_track <- GenomeAxisTrack(tx_grng)
tx_track <- GeneRegionTrack(tx_grng)
tx_track@dp@pars$shape <- "arrow"
data_track <- DataTrack(bwp_p3, type = "h", window = -1)

tx_start <- min(start(tx_grng))
bwp_sim_plt <- GRanges(seqnames = "chr9",
        ranges= IRanges(start = tx_start : (tx_start + 2e4 - 1),
                        width = 1),
        strand = "+")
bwp_sim_plt$score <- bwp_sim$score

sample_cell <- 5000
lambda <- 1000

bwp_sim_plt$score[1:200] <-
  rpois(length(bwp_sim_plt[1:200]),
        bwp_sim_plt$score[1:200] / sample_cell * lambda)

bwp_sim_plt$score[201:length(bwp_sim_plt)] <-
  rpois(length(bwp_sim_plt[201:length(bwp_sim_plt)]),
        bwp_sim_plt$score[201:length(bwp_sim_plt)] / sample_cell * lambda)

sim_track <- DataTrack(bwp_sim_plt, type = "h", window = -1)

plotTracks(list(axis_track, tx_track, data_track, sim_track),
           from = min(start(tx_grng)) - 5e3, to = max(end(tx_grng)) + 5e3)

genome(bwp_sim_plt) <- "GRCh38"
seqlengths(bwp_sim_plt) <- 1e9
  
export(bwp_sim_plt, file.path(result_dir, "sim_track_SimPolv2.bw"))

