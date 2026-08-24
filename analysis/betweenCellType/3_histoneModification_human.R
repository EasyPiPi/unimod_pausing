# ==============================================================
# LIBRARIES
# ==============================================================

library(tidyverse)
library(rtracklayer)
library(DiffBind)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(plyranges)
library(STADyUM)

# ==============================================================
# GLOBAL PATHS AND CONFIGURATION
# ==============================================================

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/betweenCellType/histoneModification")

if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
}
bam_dir   <- .paths$encode_bam

cd4_rate_in       <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd4_rate.RDS")
cd14_rate_in      <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd14_rate.RDS")
tss_dist_in       <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd4_cd14_tss_dist.txt")
lrt_results_in    <- file.path(.paths$outputs, "publish/betweenCellType/LRT/human/lrt_result.txt")
bam_metadata_path <- file.path(.paths$data, "encode/meta_bam.RDS")

proximal_region_len <- 500
histones <- c("H3K4me3", "H3K27ac")

# ==============================================================
# HELPER FUNCTIONS
# ==============================================================

run_diffbind_histone <- function(histone_mark, sample_sheet, result_dir) {
    sample_sheet_sub <- sample_sheet[sample_sheet$Factor == histone_mark, ]
    if (nrow(sample_sheet_sub) < 2) stop("Need at least two samples after filtering by histone mark.")

    chipseq <- dba(
        sampleSheet = sample_sheet_sub,
        config = data.frame(
            AnalysisMethod = DBA_DESEQ2, th = 1,
            DataType = DBA_DATA_GRANGES, RunParallel = TRUE,
            minQCth = 15, fragmentSize = 125,
            bCorPlot = FALSE, reportInit = "DBA",
            bUsePval = FALSE, design = TRUE,
            doBlacklist = TRUE, doGreylist = FALSE
        )
    )
    chipseq$config$cores <- 16
    chipseq <- dba.count(chipseq, summits = 0, filter = 0)
    chipseq <- dba.normalize(chipseq)
    chipseq <- dba.contrast(chipseq, reorderMeta = list(Tissue = "CD4_memory"))
    chipseq <- dba.analyze(chipseq)
    print(chipseq)

    pdf(file.path(result_dir, paste0(histone_mark, "_chipseq_sample_info.pdf")), width = 8, height = 6)
    plot(chipseq)
    dev.off()
    saveRDS(chipseq, file.path(result_dir, paste0(histone_mark, "_chipseq_sample_info.RDS")))

    chipseq_gr <- dba.report(chipseq, bAll = TRUE)
    seqlevelsStyle(chipseq_gr) <- "ucsc"
    return(chipseq_gr)
}

merge_diffbind_info <- function(region_of_interest, diffbind_result) {
    seqlevelsStyle(region_of_interest) <- seqlevelsStyle(diffbind_result)
    olp <- findOverlaps(region_of_interest, diffbind_result)
    bind_cols(
        as_tibble(mcols(region_of_interest[queryHits(olp), ])),
        as_tibble(mcols(diffbind_result[subjectHits(olp), ]))
    )
}

my_boxplot_style <- function() {
    list(
        geom_boxplot(outlier.shape = NA, linewidth = 0.3),
        theme_cowplot(),
        theme(
            legend.position = "none",
            axis.title.y = element_text(size = 14),
            axis.title.x = element_text(size = 18),
            axis.text    = element_text(size = 14),
            strip.text   = element_text(size = 14),
            axis.line    = element_line(size = 0.3),
            axis.ticks   = element_line(size = 0.3)
        )
    )
}

# Boxplot of histone LFC grouped by chi_lfc_group
plot_histone_boxplot <- function(data, y_var, y_label) {
    ggplot(data, aes(x = chi_lfc_group, y = .data[[y_var]], fill = chi_lfc_group)) +
        stat_compare_means(
            label.x = 1.2, label.y = 2.5,
            aes(label = gsub("<", "p < ", after_stat(p.format)))
        ) +
        labs(x = expression(chi ~ group), y = y_label) +
        coord_cartesian(ylim = c(-1, 3)) +
        my_boxplot_style() +
        scale_fill_brewer(palette = "Greens") +
        theme(
            axis.text.x  = element_text(angle = 30, hjust = 1),
            axis.title.x = element_blank()
        )
}

# Violin + boxplot of histone LFC by betaCategory, faceted by chi_lfc_group
plot_histone_violin <- function(data, y_var, y_label, ylim_range) {
    ggplot(data, aes(x = betaCategory, y = .data[[y_var]], fill = betaCategory)) +
        geom_violin(width = 0.9, alpha = 0.5, color = NA, trim = TRUE) +
        geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", linewidth = 0.3) +
        facet_wrap(~ chi_lfc_group, ncol = 3) +
        stat_compare_means(label.x = 1.2, label.y = 3, label = "p") +
        ylab(y_label) +
        xlab("") +
        coord_cartesian(ylim = ylim_range) +
        scale_fill_manual(values = c("Others" = "#BFBFBF", "Up" = "#E45756", "Down" = "#4C78A8")) +
        labs(fill = expression(beta * " change")) +
        cowplot::theme_cowplot() +
        theme(
            strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
            axis.text.x      = element_blank(),
            axis.ticks.x     = element_blank(),
            legend.position  = "bottom",
            legend.justification = "center"
        )
}

# ==============================================================
# LOAD INPUT DATA
# ==============================================================

human_lrt <- read.table(lrt_results_in, sep = "\t", header = TRUE)
cd4_rate  <- readRDS(cd4_rate_in)
cd14_rate <- readRDS(cd14_rate_in)

# ==============================================================
# GENERATE PROMOTER REGIONS
# ==============================================================

cd4_tsn <- cd4_rate@pauseRegions %>%
    anchor_5p() %>%
    mutate(width = 1, tss_type = "cd4_specific") %>%
    filter(gene_id %in% rates(cd4_rate)$geneId)

cd14_tsn <- cd14_rate@pauseRegions %>%
    anchor_5p() %>%
    mutate(width = 1, tss_type = "cd14_specific") %>%
    filter(gene_id %in% rates(cd14_rate)$geneId)

# All TSS proximal regions (500 bp centered on TSS)
tsn_all <- unique(c(cd4_tsn, cd14_tsn))
all_tss_proximal <- resize(tsn_all, proximal_region_len, fix = "center")
colnames(mcols(all_tss_proximal))[1] <- "name"
seqlevelsStyle(all_tss_proximal) <- "ucsc"
all_tss_proximal$score <- 10
export.bed(all_tss_proximal, file.path(result_dir, "all_promoter.bed"))

# Proximal regions restricted to genes with identical TSS across cell types
tss_dist <- read.table(tss_dist_in, header = TRUE, sep = ',')

id_tss_gene_names <- rownames(tss_dist)[tss_dist$Distance == 0]
id_tss_proximal <- cd4_rate@pauseRegions[cd4_rate@pauseRegions$gene_id %in% id_tss_gene_names] %>%
    anchor_5p() %>%
    mutate(width = 1)
id_tss_proximal <- resize(id_tss_proximal, proximal_region_len, fix = "center")
seqlevelsStyle(id_tss_proximal) <- "ucsc"

# ==============================================================
# DIFFBIND ANALYSIS
# ==============================================================

meta_bam <- readRDS(bam_metadata_path)
meta_bam <- meta_bam %>%
    filter(Target != "Control") %>%
    left_join(
        meta_bam %>% filter(Target == "Control"),
        by = c("Donor", "Biosample"),
        suffix = c(".chip", ".input")
    )

sample_sheet <- meta_bam %>%
    filter(Donor != "ENCDO157PKH") %>%  # excluded: low read fraction within peak
    mutate(
        bamReads   = file.path(bam_dir, paste0(Accession.chip, ".bam")),
        bamControl = file.path(bam_dir, paste0(Accession.input, ".bam")),
        Replicate  = case_when(
            Donor == "ENCDO836ICQ" ~ 1,
            Donor == "ENCDO922TOC" ~ 2,
            Donor == "ENCDO553CCA" ~ 3,
            TRUE                   ~ 5
        ),
        Peaks      = file.path(result_dir, "all_promoter.bed"),
        PeakCaller = "bed",
        Biosample  = str_replace(Biosample, " ", "_")
    ) %>%
    dplyr::select(Accession.chip, Biosample, Target.chip, Replicate,
                  bamReads, bamControl, PeakCaller, Peaks)

colnames(sample_sheet)[1:3] <- c("SampleID", "Tissue", "Factor")
sample_sheet <- as.data.frame(sample_sheet)

# Run DiffBind for each histone mark and merge with promoter regions
tss_diffbind_list <- list()
for (histone in histones) {
    diffbind_result              <- run_diffbind_histone(histone, sample_sheet, result_dir)
    tss_diffbind_list[[histone]] <- merge_diffbind_info(id_tss_proximal, diffbind_result)
}

# ==============================================================
# MERGE HISTONE SIGNALS WITH LRT RESULTS
# ==============================================================



merge_diffbind_df <- dplyr::left_join(
    tss_diffbind_list$H3K4me3 %>%
        dplyr::select(gene_id, H3K4me3_lfc = Fold),

    tss_diffbind_list$H3K27ac %>%
        dplyr::select(gene_id, H3K27ac_lfc = Fold),

    by = "gene_id"
)

human_lrt_histone <- human_lrt %>%
    filter(geneId %in% id_tss_gene_names) %>%
    left_join(merge_diffbind_df, by = c("geneId" = "gene_id"))

write.table(human_lrt_histone,
            file.path(result_dir, "lrt_histone.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

human_lrt_histone <- human_lrt_histone %>%
    mutate(chi_lfc_group = factor(chi_lfc_group, levels = c("Decrease", "Unchanged", "Increase"))) %>%
    droplevels()

# ==============================================================
# VISUALIZATION
# ==============================================================

# Boxplots: transcription rate change (chi_lfc_group) vs histone LFC
p_box_h3k27ac <- plot_histone_boxplot(
    human_lrt_histone,
    y_var   = "H3K27ac_lfc",
    y_label = expression(log[2] * FC ~ "(H3K27ac)")
)
ggsave(file.path(result_dir, "chi_lfc_h3K27ac.pdf"), p_box_h3k27ac, width = 2.5, height = 4)

p_box_h3k4me3 <- plot_histone_boxplot(
    human_lrt_histone,
    y_var   = "H3K4me3_lfc",
    y_label = expression(log[2] * FC ~ "(H3K4me3)")
)
ggsave(file.path(result_dir, "chi_lfc_h3K4me3.pdf"), p_box_h3k4me3, width = 2.5, height = 4)

# Recode group labels for violin plots
human_lrt_histone <- human_lrt_histone %>%
    mutate(chi_lfc_group = recode(chi_lfc_group,
        "Decrease"  = "Decreased expression",
        "Unchanged" = "Unchanged expression",
        "Increase"  = "Increased expression"
    ))

# Violin plots: beta change category vs histone LFC, faceted by chi_lfc_group
p_vln_h3k27ac <- plot_histone_violin(
    human_lrt_histone,
    y_var      = "H3K27ac_lfc",
    y_label    = expression(Log[2] * "FC (H3K27ac signal)"),
    ylim_range = c(-1.5, 3.5)
)
ggsave(file.path(result_dir, "pause_lfc_H3K27ac.pdf"), p_vln_h3k27ac, width = 8, height = 3)

p_vln_h3k4me3 <- plot_histone_violin(
    human_lrt_histone,
    y_var      = "H3K4me3_lfc",
    y_label    = expression(Log[2] * "FC (H3K4me3 signal)"),
    ylim_range = c(-1, 3.5)
)
ggsave(file.path(result_dir, "pause_lfc_H34me3.pdf"), p_vln_h3k4me3, width = 8, height = 4)
