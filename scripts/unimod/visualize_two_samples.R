#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####
grng_in <- snakemake@input[["grng"]]
gtf_in <- snakemake@input[["gtf"]]

spike_in <- snakemake@input[["spike_in"]]

bwp1_p3_in <- snakemake@input[["bwp1_p3"]]
bwm1_p3_in <- snakemake@input[["bwm1_p3"]]

bwp2_p3_in <- snakemake@input[["bwp2_p3"]]
bwm2_p3_in <- snakemake@input[["bwm2_p3"]]

beta_in <- snakemake@input[["beta"]]
omega_in <- snakemake@input[["omega"]]
alpha_in <- snakemake@input[["alpha"]]
tf_target_in <- snakemake@input[["tf_target"]]

result_dir <- snakemake@params[["result_dir"]]

#### load packages ####
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(ggbeeswarm)
library(rtracklayer)
library(genomation)
library(plyranges)
library(Gviz)
library(ggpubr)
library(ggpointdensity)
library(viridis)

options(ucscChromosomeNames=FALSE)

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_human"
# 
# # grng_in <- file.path(root_dir, "outputs/read_dt/granges_for_read_counting_K562.RData")
# grng_in <- file.path(root_dir, "outputs/read_dt/granges_for_read_counting_DLD1.RData")
# 
# gtf_in <- file.path(root_dir, "outputs/read_dt/human_transcript_granges.rds")
# 
# spike_in <- file.path(root_dir, "metadata/scaling_factor.csv")
# 
# bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_plus.bw")
# bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_Auxin_Ctrl-SE_minus.bw")
# 
# bwp2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_Auxin-SE_plus.bw")
# bwm2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_Auxin-SE_minus.bw")
# 
# result_dir <- file.path(root_dir, "outputs/between_samples", "PROseq-DLD1-aoi-NELFC_Auxin-SE-identity-all")
# 
# bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_NVP2_Ctrl-SE_plus.bw")
# bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_NVP2_Ctrl-SE_minus.bw")
# 
# bwp2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_NVP2-SE_plus.bw")
# bwm2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-DLD1-aoi-NELFC_NVP2-SE_minus.bw")
# 
# result_dir <- file.path(root_dir, "outputs/between_samples", "PROseq-DLD1-aoi-NELFC_NVP2-SE-identity-all")
# 
# # bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_plus.bw")
# # bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_minus.bw")
# # 
# # bwp2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-treated-SE_plus.bw")
# # bwm2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-treated-SE_minus.bw")
# # 
# # result_dir <- file.path(root_dir, "outputs/between_samples", "PROseq-K562-dukler-SE-identity-all")
# 
# # bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-control-PE_plus.bw")
# # bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-control-PE_minus.bw")
# #
# # bwp2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-treated-PE_plus.bw")
# # bwm2_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-treated-PE_minus.bw")
# #
# # result_dir <- file.path(root_dir, "outputs/between_samples", "PROseq-K562-chivu-PE-identity-all")
# 
# omega_in <- file.path(result_dir, "omega.csv")
# alpha_in <- file.path(result_dir, "alpha.csv")
# beta_in <- file.path(result_dir, "beta.csv")
# 
# tf_target_in <- file.path(root_dir, "outputs/between_samples/table/gsea_targets.csv")

#### end of parsing arguments ####
fig_dir <- file.path(result_dir, "gviz")

walk(c(result_dir, fig_dir, file.path(fig_dir, "omega"),
       file.path(fig_dir, "beta", c("up", "down")),
       file.path(fig_dir, "combine",
                 c("beta_up_omega_up", "beta_down_omega_down",
                   "beta_others_omega_up", "beta_others_omega_down",
                   "beta_down_omega_up", "select")),
       file.path(fig_dir, "tf_target",
                 c("beta_up_omega_up", "beta_down_omega_down",
                   "beta_others_omega_up", "beta_others_omega_down"))), 
     dir.create, showWarnings = FALSE, recursive = TRUE)
# dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

load(grng_in)
gtf <- readRDS(gtf_in)

# read in number of spike-in or total number of mappable reads
# use them as scaling factor
scale_tbl <- read_csv(spike_in, show_col_types = FALSE)
# subset the right table
scale_tbl <- scale_tbl[str_detect(bwp1_p3_in, scale_tbl$sample), ]

scale_factor <-
  (scale_tbl$control_1 + ifelse(is.na(scale_tbl$control_2), 0, scale_tbl$control_2)) /
  (scale_tbl$treated_1 + ifelse(is.na(scale_tbl$treated_2), 0, scale_tbl$treated_2))

# process data
process_bw <- function(bwp_in, bwm_in, scale = 1) {
    bwp <- import.bw(bwp_in)
    bwm <- import.bw(bwm_in)
    strand(bwp) <- "+"
    strand(bwm) <- "-"
    bwp <- bwp[abs(bwp$score) > 0]
    bwm <- bwm[abs(bwm$score) > 0]
    bwp <- BRGenomics::makeGRangesBRG(bwp)
    bwm <- BRGenomics::makeGRangesBRG(bwm)
    bwp$score <- bwp$score * scale
    bwm$score <- bwm$score * scale
    return(list("bwp" = bwp, "bwm" = bwm))
}

bw_ctrl <- process_bw(bwp1_p3_in, bwm1_p3_in) 
bw_trtd <- process_bw(bwp2_p3_in, bwm2_p3_in, scale = scale_factor)

# reduce to gene level. but not pretty
# gn_grng <- count_grng %>%
#     plyranges::group_by(gene_id) %>%
#     plyranges::summarise(start = min(start), end = max(end))
# 
# gn_grng <- DataFrame(gene_id = count_grng$gene_id,
#           strand = strand(count_grng),
#           seqnames = seqnames(count_grng)) %>%
#     unique() %>%
#     merge(gn_grng, by = "gene_id") %>%
#     plyranges::as_granges()

# extent pausing and termination region a bit for better visualization 
bw_pause_ext <- bw_pause_filtered %>% plyranges::anchor_3p() %>% plyranges::stretch(251)

# concat plus and minus strand for visualization
bw_ctrl_bs <- c(bw_ctrl$bwp, bw_ctrl$bwm) 
bw_ctrl_bs$score <- abs(bw_ctrl_bs$score)

bw_trtd_bs <- c(bw_trtd$bwp, bw_trtd$bwm)
bw_trtd_bs$score <- abs(bw_trtd_bs$score)

# function to save figure and do metaplot
save_png <- function(file_name, plot_fun, width = 600, height = 400) {
    png(filename = file_name, width = width, height = height, units = "px")
    plot_fun
    invisible(dev.off())
}

meta_plot <- function(sm, smlcolors, meta.rescale = FALSE, xcoords = c(-250, 250),
                      centralTend = "mean", dispersion = NULL, ylab = "Average Read Counts", xlab = "Distance from TSS",
                      ylim = c(0, 0.15)) {
    plotMeta(sm, xcoords = xcoords, meta.rescale = meta.rescale,
             line.col = smlcolors[c(2, 6)],
             centralTend = centralTend, dispersion = dispersion,
             xlab = xlab, ylab = ylab,
             dispersion.col = smlcolors[c(1, 5)], ylim = ylim, cex.lab = 1.5, cex.axis = 1.2)
    legend("topright", names(sm), lty=c(1,1), lwd=c(2.5,2.5), col=smlcolors[c(2, 6)], cex = 1.5)
}

## plot PRO-seq density around TSS and within gene body ##
# metaplot
beta_tbl <- read_csv(beta_in, show_col_types = FALSE)
beta_up <- beta_tbl %>% filter(category == "Up") %>% pull(gene_id)
beta_down <- beta_tbl %>% filter(category == "Down") %>% pull(gene_id)

paired_color <- RColorBrewer::brewer.pal(name = "Paired", n = 6)

sm_pause <- ScoreMatrixList(target = GRangesList("Control" = bw_ctrl_bs, "Treated" = bw_trtd_bs),
               windows = bw_pause_ext, strand.aware = TRUE, weight.col = "score")

save_png(file_name = file.path(result_dir, "proseq_signal_around_tss.png"),
         plot_fun = meta_plot(sm_pause, paired_color, ylim = NULL, dispersion = "se"))

sm_pause_up <- 
  ScoreMatrixList(target = GRangesList("Control" = bw_ctrl_bs, "Treated" = bw_trtd_bs),
                  windows = bw_pause_ext[bw_pause_ext$gene_id %in% beta_up, ],
                  strand.aware = TRUE, weight.col = "score")

save_png(file_name = file.path(result_dir, "proseq_signal_around_tss_beta_up.png"),
         plot_fun = meta_plot(sm_pause_up, paired_color, ylim = NULL, dispersion = "se"))

sm_pause_down <- 
  ScoreMatrixList(target = GRangesList("Control" = bw_ctrl_bs, "Treated" = bw_trtd_bs),
                  windows = bw_pause_ext[bw_pause_ext$gene_id %in% beta_down, ],
                  strand.aware = TRUE, weight.col = "score")

save_png(file_name = file.path(result_dir, "proseq_signal_around_tss_beta_down.png"),
         plot_fun = meta_plot(sm_pause_down, paired_color, ylim = NULL, dispersion = "se"))

sm_gb <- ScoreMatrixList(target = GRangesList("Control" = bw_ctrl_bs, "Treated" = bw_trtd_bs),
                          bin.num = 100, windows = bw_gb_filtered, strand.aware = TRUE,
                          weight.col = "score")

save_png(file_name = file.path(result_dir, "proseq_signal_within_genebody.png"),
         plot_fun = meta_plot(sm_gb, paired_color, xlab = "Bins within gene body",
                              ylim = NULL, xcoords = c(0, 100), dispersion = "se"),
         width = 800, height = 400)

sm_gb_up <- ScoreMatrixList(target = GRangesList("Control" = bw_ctrl_bs, "Treated" = bw_trtd_bs),
                         bin.num = 100, windows = bw_gb_filtered[bw_gb_filtered$gene_id %in% beta_up, ],
                         strand.aware = TRUE, weight.col = "score")

save_png(file_name = file.path(result_dir, "proseq_signal_within_genebody_beta_up.png"),
         plot_fun = meta_plot(sm_gb_up, paired_color, xlab = "Bins within gene body",
                              ylim = NULL, xcoords = c(0, 100), dispersion = "se"),
         width = 800, height = 400)

if (length(beta_down) > 5) {
  sm_gb_down <- ScoreMatrixList(target = GRangesList("Control" = bw_ctrl_bs, "Treated" = bw_trtd_bs),
                                bin.num = 100, windows = bw_gb_filtered[bw_gb_filtered$gene_id %in% beta_down, ],
                                strand.aware = TRUE, weight.col = "score")
  
  save_png(file_name = file.path(result_dir, "proseq_signal_within_genebody_beta_down.png"),
           plot_fun = meta_plot(sm_gb_down, paired_color, xlab = "Bins within gene body",
                                ylim = NULL, xcoords = c(0, 100), dispersion = "se"),
           width = 800, height = 400)
}

# heatmap
# sm_pause_scaled <- scaleScoreMatrixList(sm_pause)
sort_pause_site <- function(sml) {
  sort_name <- names(sort(apply(sml[[1]]@.Data[, 251:501], MARGIN = 1, FUN = which.max)))
  return(as(lapply(sml, function(x) x@.Data[sort_name, ]), "ScoreMatrixList"))
 }

save_png(file_name = file.path(result_dir, "proseq_heatmap_around_tss.png"),
         plot_fun = multiHeatMatrix(sort_pause_site(sm_pause), xcoords = c(-250, 250),
                                    col = RColorBrewer::brewer.pal("Reds", n = 9),
                                    winsorize = c(0, 99), common.scale = TRUE, xlab="Distance from TSS",
                                    cex.axis = 0.8, cex.lab = 1.2),
         width = 600, height = 600)

save_png(file_name = file.path(result_dir, "proseq_heatmap_around_tss_beta_up.png"),
         plot_fun = multiHeatMatrix(sort_pause_site(sm_pause_up), xcoords = c(-250, 250),
                                    col = RColorBrewer::brewer.pal("Reds", n = 9),
                                    winsorize = c(0, 99), common.scale = TRUE, xlab="Distance from TSS",
                                    cex.axis = 0.8, cex.lab = 1.2),
         width = 600, height = 400)

if (length(beta_down) > 5) {
  save_png(file_name = file.path(result_dir, "proseq_heatmap_around_tss_beta_down.png"),
           plot_fun = multiHeatMatrix(sort_pause_site(sm_pause_down), xcoords = c(-250, 250),
                                      col = RColorBrewer::brewer.pal("Reds", n = 9),
                                      winsorize = c(0, 99), common.scale = TRUE, xlab="Distance from TSS",
                                      cex.axis = 0.8, cex.lab = 1.2),
           width = 600, height = 600)
}

save.image(file = file.path(result_dir, "data.RData"))

# load(file.path(result_dir, "data.RData"))
if (stringr::str_detect(omega_in, "chivu")) quit(save = "no", status = 0, runLast = FALSE)

zeta <- 2000 

theme_set(cowplot::theme_cowplot())

# visualize some interesting tf targets
omega_tbl <- read_csv(omega_in, show_col_types = FALSE)
alpha_tbl <- read_csv(alpha_in, show_col_types = FALSE)
gsea_target <- read_csv(tf_target_in, show_col_types = FALSE)

target_rate <- gsea_target %>%
  left_join(beta_tbl[c("gene_id", "beta1", "beta2", "category", "lfc")],
            by = c("ensembl_gene_id" = "gene_id")) %>%
  na.omit()

target_rate <- target_rate %>%
  left_join(omega_tbl[c("gene_id", "chi1", "chi2", "category", "lfc")],
            by = c("ensembl_gene_id" = "gene_id"),
            suffix = c("_beta", "_omega")) %>%
  left_join(alpha_tbl[c("gene_id", "alpha_zeta1", "alpha_zeta2", "category", "lfc")],
            by = c("ensembl_gene_id" = "gene_id")) %>%
  dplyr::rename("lfc_alpha" = "lfc", "category_alpha" = "category")

target_rate <- target_rate %>%
  mutate(across(contains("category"), ~ factor(.x,  levels = c("Up", "Down", "Others"))))

median_rates <- target_rate %>% 
  summarise(lfc_beta = median(lfc_beta),
            lfc_omgega = median(lfc_omega),
            lfc_alpha = median(lfc_alpha))

p1 <- target_rate %>%
  select(tf, category_omega, lfc_omega) %>%
  ggplot() +
  geom_beeswarm(aes(x = tf, y = lfc_omega, color = category_omega), alpha = 1, size = 0.7) +
  scale_color_manual(values=c("#E41A1C", "#377EB8", "gray")) +
  stat_summary(aes(x = tf, y = lfc_omega), fun = "median", geom = "crossbar", width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ylim(c(-4.5, 5.5)) +
  labs(y = expression(log[2]*"("*chi["Treated"]/chi["Control"]*")"),
       color = "Category") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p2 <- target_rate %>%
  select(tf, category_beta, lfc_beta) %>%
  ggplot() +
  geom_beeswarm(aes(x = tf, y = lfc_beta, color = category_beta), alpha = 1, size = 0.7) +
  scale_color_manual(values=c("#E41A1C", "#377EB8", "gray")) +
  stat_summary(aes(x = tf, y = lfc_beta), fun = "median", geom = "crossbar", width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  ylim(c(-4.5, 5.5)) +
  labs(y = expression(log[2]*"("*beta["Treated"]/beta["Control"]*")"),
       color = "Category") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

p3 <- target_rate %>%
  select(tf, category_alpha, lfc_alpha) %>%
  ggplot() +
  geom_beeswarm(aes(x = tf, y = lfc_alpha, color = category_alpha), alpha = 1, size = 0.7) +
  scale_color_manual(values=c("#E41A1C", "#377EB8", "gray")) +
  stat_summary(aes(x = tf, y = lfc_alpha), fun = "median", geom = "crossbar", width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  #ylim(c(-4.5, 5.5)) +
  labs(y = expression(log[2]*"("*alpha["Treated"]/alpha["Control"]*")"),
       color = "Category") +
  theme(axis.title.x = element_blank())

p <- cowplot::plot_grid(p1, p2, p3, axis = "tblr", rel_widths = c(2.4, 2.4, 3), nrow = 1)
ggsave(filename = file.path(result_dir, "tf_target_rates.png"), plot = p,
       width = 12, height = 4)

# plot genes with changes
beta_down <- beta_tbl %>% arrange(padj) %>% filter(category == "Down") %>% pull(gene_id)
beta_up <- beta_tbl %>% arrange(padj) %>% filter(category == "Up") %>% pull(gene_id)

count_region <- c(bw_pause_filtered, bw_gb_filtered)

# functions for gviz plots
gviz_plot <- function(gene_sel, extend_region, category, dir = "beta") {

  make_data_track <- function(bw, bsize, s, chrom, condition) {
    data_track <- Gviz::DataTrack(
      range = bw,
      type = "h",
      window = -1,
      windowSize = bsize,
      name = paste0(condition, " (", s, ")"),
      col = strand_col[s],
      strand = s,
      chromosome = chrom)
    return(data_track)
  }

  set_datatrack_ylim <- function (data_track, ylim) {
    s <- Gviz::strand(data_track)
    s <- factor(as.character(s), levels = c("+", "-", "*"))
    scale <- do.call(c(I, rev, I)[as.numeric(s)][[1]], list(ylim))
    Gviz::displayPars(data_track)$ylim <- scale
    return(data_track)
  }

  gene_tx <- gtf[gtf$gene_id == gene_sel]
  gene_range <-
    GRanges(seqnames = seqnames(gene_tx[1]),
            ranges = IRanges(start = min(start(gene_tx)),
                             end = max(end(gene_tx))))
  gene_range <- plyranges::stretch(gene_range, extend_region)

  region_tx <- subsetByOverlaps(gtf, gene_range)
  region_gn <- region_tx %>%
    group_by(gene_id) %>%
    reduce_ranges_directed()

  bw_range <-
    GRanges(seqnames = seqnames(region_gn[1]),
            ranges = IRanges(start = min(start(region_gn)),
                             end = max(end(region_gn))))
  bw_range <- plyranges::stretch(bw_range, extend_region)

  axis_track <- GenomeAxisTrack()
  tx_track <-
    GeneRegionTrack(region_gn, name = "Gene", shape = "arrow",
                    chromosome = seqnames(gene_tx[1]))
  count_track <-
    AnnotationTrack(subsetByOverlaps(count_region, gene_range), shape = "box")

  strand_col <- c(`+` = "blue", `-` = "red")
  bsize <- 10

  ctrl_bwp <- subsetByOverlaps(bw_ctrl$bwp, bw_range)
  ctrl_bwm <- subsetByOverlaps(bw_ctrl$bwm, bw_range)
  ctrl_bwm$score <- abs(ctrl_bwm$score)

  trtd_bwp <- subsetByOverlaps(bw_trtd$bwp, bw_range)
  trtd_bwm <- subsetByOverlaps(bw_trtd$bwm, bw_range)
  trtd_bwm$score <- abs(trtd_bwm$score)

  ctrl_plus <- make_data_track(ctrl_bwp, bsize, "+", seqnames(gene_tx[1]), condition = "Control")
  ctrl_minus <- make_data_track(ctrl_bwm, bsize, "-", seqnames(gene_tx[1]), condition = "Control")
  trtd_plus <- make_data_track(trtd_bwp, bsize, "+", seqnames(gene_tx[1]), condition = "Treated")
  trtd_minus <- make_data_track(trtd_bwm, bsize, "-", seqnames(gene_tx[1]), condition = "Treated")

  # Use rolling mean to determine ylim
  ylim <-
    ceiling(max(map_dbl(list(ctrl_bwp, ctrl_bwm, trtd_bwp, trtd_bwm),
                        function(x) {
                          max(zoo::rollmean(x$score, bsize))
                        })))

  save_png(file.path(fig_dir, dir, category, paste0(gene_sel, ".png")),
           plotTracks(
             c(list(axis_track, tx_track, count_track),
               map(list(ctrl_plus, ctrl_minus, trtd_plus, trtd_minus), set_datatrack_ylim, c(0, ylim)))
           ),
           width = 600, height = 400)
}

# the most significant genes
for (gene_sel in beta_up[1:10]) {
  gviz_plot(gene_sel, 5000, "up")
}

for (gene_sel in beta_down[1:10]) {
  gviz_plot(gene_sel, 5000, "down")
}

# plots for some TF targets
beta_up_omega_up <- target_rate %>%
  filter(category_beta == "Up", category_omega == "Up") %>%
  pull(ensembl_gene_id)

if (length(beta_up_omega_up) > 0) {
  for (gene_sel in beta_up_omega_up) {
    gviz_plot(gene_sel, 5000, "beta_up_omega_up", dir = "tf_target")
  }
}

beta_down_omega_down <- target_rate %>%
  filter(category_beta == "Down", category_omega == "Down") %>%
  pull(ensembl_gene_id)

if (length(beta_down_omega_down) > 0) {
  for (gene_sel in beta_down_omega_down) {
    gviz_plot(gene_sel, 5000, "beta_down_omega_down", dir = "tf_target")
  }
}

beta_others_omega_up <- target_rate %>%
  filter(category_beta == "Others", category_omega == "Up") %>%
  pull(ensembl_gene_id)

if (length(beta_others_omega_up) > 0) {
  for (gene_sel in beta_others_omega_up) {
    gviz_plot(gene_sel, 5000, "beta_others_omega_up", dir = "tf_target")
  }
}

beta_others_omega_down <- target_rate %>%
  filter(category_beta == "Others", category_omega == "Down") %>%
  pull(ensembl_gene_id)

if (length(beta_others_omega_down) > 0) {
  for (gene_sel in beta_others_omega_down) {
    gviz_plot(gene_sel, 5000, "beta_others_omega_down", dir = "tf_target")
  }
}

# plots for gene may be driven by initiation or pause-escape
rate_tbl <- beta_tbl %>%
  select(gene_id, beta1,  beta2, category, t_stats, lfc, padj) %>%
  left_join(omega_tbl %>% select(gene_id, chi1, chi2, category, t_stats, lfc, padj),
            by = "gene_id", suffix = c("_beta", "_omega"))

beta_up_omega_up <- rate_tbl %>%
  filter(category_beta == "Up", category_omega == "Up") %>%
  arrange(desc(t_stats_omega)) %>% 
  slice_head(n = 20) %>%
  pull(gene_id)

if (length(beta_up_omega_up) > 0) {
  for (gene_sel in beta_up_omega_up) {
    gviz_plot(gene_sel, 5000, "beta_up_omega_up", dir = "combine")
  }
}

beta_down_omega_down <- rate_tbl %>%
  filter(category_beta == "Down", category_omega == "Down") %>%
  arrange(desc(t_stats_omega)) %>% 
  slice_head(n = 20) %>%
  pull(gene_id)

if (length(beta_down_omega_down) > 0) {
  for (gene_sel in beta_down_omega_down) {
    gviz_plot(gene_sel, 5000, "beta_down_omega_down", dir = "combine")
  }
}

beta_others_omega_up <- rate_tbl %>%
  filter(category_beta == "Others", category_omega == "Up") %>%
  arrange(desc(t_stats_omega)) %>% 
  slice_head(n = 20) %>%
  pull(gene_id)

if (length(beta_others_omega_up) > 0) {
  for (gene_sel in beta_others_omega_up) {
    gviz_plot(gene_sel, 5000, "beta_others_omega_up", dir = "combine")
  }
}

beta_others_omega_down <- rate_tbl %>%
  filter(category_beta == "Others", category_omega == "Down") %>%
  arrange(desc(t_stats_omega)) %>% 
  slice_head(n = 20) %>%
  pull(gene_id)

if (length(beta_others_omega_down) > 0) {
  for (gene_sel in beta_others_omega_down) {
    gviz_plot(gene_sel, 5000, "beta_others_omega_down", dir = "combine")
  }
}

beta_down_omega_up <- rate_tbl %>%
  filter(category_beta == "Down", category_omega == "Up") %>%
  arrange(desc(t_stats_omega)) %>% 
  slice_head(n = 20) %>%
  pull(gene_id)

if (length(beta_down_omega_up) > 0) {
  for (gene_sel in beta_down_omega_up) {
    gviz_plot(gene_sel, 5000, "beta_down_omega_up", dir = "combine")
  }
}

# visualize select genes
sel_genes <-
  c("ENSG00000099194", "ENSG00000132182",
    "ENSG00000133454", "ENSG00000139793",
    "ENSG00000110958", "ENSG00000170606",
    "ENSG00000144381", "ENSG00000151929")

for (gene_sel in sel_genes) {
  try(gviz_plot(gene_sel, 5000, "select", dir = "combine"))
}

# lfc scatter plot 
rate_tbl %>%
  count(category_beta, category_omega) %>%
  pivot_wider(id_cols = category_omega, names_from = category_beta, values_from = n)

theme_set(cowplot::theme_cowplot())

sel_symbol <-
  AnnotationDbi::select(org.Hs.eg.db, keys = sel_genes,
                        keytype = "ENSEMBL", columns = "SYMBOL")

gene_tbl <- rate_tbl %>%
  filter(gene_id %in% sel_genes) %>%
  left_join(sel_symbol, by = c("gene_id" = "ENSEMBL"))

gene_tbl %>% select(SYMBOL, lfc_omega, padj_omega, lfc_beta, padj_beta)

p <- rate_tbl %>%
  ggplot(aes(x = lfc_omega, y = lfc_beta)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  # geom_abline(slope = 1) +
  stat_cor() +
  geom_point(data = gene_tbl, aes(x = lfc_omega, y = lfc_beta), color = "red") +
  # geom_text(data = gene_tbl, aes(label=SYMBOL), size=3) +
  ggrepel::geom_text_repel(data = gene_tbl, aes(label=SYMBOL),
                           size=3, color = "#595959") +
  labs(x = expression(log[2]*"("*chi["Treated"]/chi["Control"]*")"),
       y = expression(log[2]*"("*beta["Treated"]/beta["Control"]*")")) +
  ylim(-6, 6) +
  xlim(-5, 7)

ggsave(file.path(result_dir, "lfc_alpha_vs_beta_density.png"), plot = p,
       width = 7, height = 5)

# rate_tbl %>%
#   # mutate(category_omega = ifelse(lfc_omega > 0, "Up", "Down")) %>%
#   ggplot(aes(x = lfc_omega, y = lfc_beta, color = category_omega)) +
#   geom_point() +
#   geom_smooth(aes(color = category_omega), method=lm, se=FALSE, fullrange=TRUE)

