# ============================================================
# 2_2_nucleosomePositioning_cd4.R
#
# Compare nucleosome positioning (entropy of MNase signal)
# across primate species for CD4 cells, stratified by TSS type.
#
# Output:
#   <result_dir>/all_cd4_comparison.pdf
#
# Self-contained: no external R scripts are sourced.
# Dependencies: ggpubr, tidyverse, cowplot, rstatix.
# ============================================================

library(ggpubr)
library(tidyverse)
library(cowplot)
library(rstatix)
library(patchwork)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", "."),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

# ── Path configuration ────────────────────────────────────────────────────────
lrt_root   <- file.path(.paths$outputs, "publish/betweenSpecies/LRT/")
nuc_root   <- file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/")
result_dir <- file.path(.paths$outputs, "publish/betweenSpecies/NucleosomePositioning/")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)


# ── Helper functions ──────────────────────────────────────────────────────────

# Normalized Shannon entropy of a signal vector; returns NA for zero-sum input.
compute_entropy_norm <- function(signal, eps = 1e-8) {
  signal <- pmax(signal, 0)
  if (sum(signal) == 0) return(NA_real_)
  p <- signal / sum(signal)
  H <- -sum(p * log(p + eps))
  H / log(length(p))
}

# Match rows of two matrices via gene-ID lookup in map_df; scale matrix2 by scale_factor.
get_match_matrix <- function(map_df, key1, mat1_key, mat1,
                              key2, mat2_key, mat2, scale_factor = 1) {
  i1 <- match(map_df[[key1]], mat1_key)
  i2 <- match(map_df[[key2]], mat2_key)
  list(
    matrix1 = mat1[i1, , drop = FALSE],
    matrix2 = mat2[i2, , drop = FALSE] * scale_factor
  )
}

# Compute log2FC of (1 - entropy) for sp2 vs sp1 using MNase matrix columns.
# Positive values indicate sp2 is more organized (sharper positioning) than sp1.
compute_nuc_entropy_lfc <- function(tss_df,
                                    sp1_name, sp1_geneId, sp1_raw,
                                    sp2_name, sp2_geneId, sp2_raw,
                                    scale_factor,
                                    cols = 1001:1501,
                                    eps = 1e-8) {
  matched <- get_match_matrix(
    tss_df,
    sp1_name, sp1_geneId, sp1_raw[, cols, drop = FALSE],
    sp2_name, sp2_geneId, sp2_raw[, cols, drop = FALSE],
    scale_factor
  )
  sp1_entropy <- apply(matched$matrix1, 1, compute_entropy_norm)
  sp2_entropy <- apply(matched$matrix2, 1, compute_entropy_norm)
  log2(pmax(1 - sp2_entropy, eps) / pmax(1 - sp1_entropy, eps))
}

my_boxplot_style <- function() {
  list(
    geom_boxplot(outlier.shape = NA, linewidth = 0.3),
    theme_cowplot(),
    theme(
      legend.position = "none",
      axis.title.y    = element_text(size = 14),
      axis.title.x    = element_text(size = 18),
      axis.text       = element_text(size = 14),
      strip.text      = element_text(size = 14),
      axis.line       = element_line(size = 0.3),
      axis.ticks      = element_line(size = 0.3)
    )
  )
}


# ── Constants ─────────────────────────────────────────────────────────────────

# Sequencing depth per species (reads); scale_factor = sp1_reads / sp2_reads
# normalises the sp2 matrix to sp1 sequencing depth before entropy computation.
HUMAN_READS  <- 577010542
RHESUS_READS <- 485755186
BABOON_READS <- 522837005

# pause_switch groups compared against "Unchanged" baseline
MY_COMP <- list(c("Broad_to_Sharp", "Unchanged"), c("Sharp_to_Broad", "Unchanged"))

FILL_COLORS <- c(
  "Broad_to_Sharp" = "#FC8D62",
  "Sharp_to_Broad" = "#66C2A5",
  "Unchanged"      = "grey80"
)

Y_BREAKS <- c(-1, 0, 1)
Y_LABELS <- c("-1.0", "0", "1.0")

# ── Load nucleosome data (one RDS per species) ────────────────────────────────
human_cd4_nuc  <- readRDS(file.path(nuc_root, "human/cd4_ns_matrix.RDS"))
rhesus_cd4_nuc <- readRDS(file.path(nuc_root, "rhesus/cd4_ns_matrix.RDS"))
baboon_cd4_nuc <- readRDS(file.path(nuc_root, "baboon/cd4_ns_matrix.RDS"))


# ── 1. Human vs Rhesus ───────────────────────────────────────────────────────
# sp1 = rhesus (denominator), sp2 = human (numerator):
#   positive nuc_positioning_lfc → human more organised than rhesus
#   pause_switch = rhesusCD4_sdGroup → humanCD4_sdGroup
human_rhesus_lrt <- read.table(file.path(lrt_root, "human_rhesus_cd4_lrt.txt"),
                                sep = "\t", header = TRUE)

human_rhesus_lrt$pause_switch <- paste0(human_rhesus_lrt$rhesusCD4_sdGroup, "_to_",
                                         human_rhesus_lrt$humanCD4_sdGroup)
human_rhesus_lrt$nuc_positioning_lfc <- compute_nuc_entropy_lfc(
  tss_df    = human_rhesus_lrt,
  sp1_name  = "rhesus", sp1_geneId = rhesus_cd4_nuc$rate_df$geneId, sp1_raw = rhesus_cd4_nuc$raw,
  sp2_name  = "human",  sp2_geneId = human_cd4_nuc$rate_df$geneId,  sp2_raw = human_cd4_nuc$raw,
  scale_factor = RHESUS_READS / HUMAN_READS
)

human_rhesus_lrt$pause_switch[
  human_rhesus_lrt$pause_switch %in% c("Broad_to_Broad", "Sharp_to_Sharp")] <- "Unchanged"
human_rhesus_lrt <- human_rhesus_lrt %>%
  mutate(pause_switch = factor(pause_switch, levels = c("Sharp_to_Broad", "Unchanged", "Broad_to_Sharp")))

human_rhesus_stat <- human_rhesus_lrt %>%
  group_by(tssType) %>%
  pairwise_wilcox_test(nuc_positioning_lfc ~ pause_switch, comparisons = MY_COMP) %>%
  ungroup()
human_rhesus_stat$y.position <- c(1, 1, 1, 1)

human_rhesus_hline_df <- human_rhesus_lrt %>%
  group_by(tssType) %>%
  summarise(yint = median(nuc_positioning_lfc[pause_switch == "Unchanged"], na.rm = TRUE),
            .groups = "drop")

human_rhesus_delta_p <- ggplot(human_rhesus_lrt,
    aes(x = pause_switch, y = nuc_positioning_lfc, fill = pause_switch)) +
  facet_wrap(~ tssType, ncol = 1) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.6) +
  stat_pvalue_manual(human_rhesus_stat, label = "p.adj.signif") +
  ylab(expression(Log[2] * "FC(1-Entropy)")) +
  xlab("") +
  scale_fill_manual(values = FILL_COLORS) +
  coord_cartesian(ylim = c(-1, 1.4)) +
  scale_y_continuous(breaks = Y_BREAKS, labels = Y_LABELS) +
  my_boxplot_style() +
  geom_hline(data = human_rhesus_hline_df, aes(yintercept = yint),
             linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
  ggtitle("Human/Rhesus") +
  theme(legend.position  = "none",
        strip.background = element_blank(),
        strip.text       = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5, face = "plain"))


# ── 2. Human vs Baboon ───────────────────────────────────────────────────────
# sp1 = baboon (denominator), sp2 = human (numerator):
#   positive nuc_positioning_lfc → human more organised than baboon
#   pause_switch = baboonCD4_sdGroup → humanCD4_sdGroup
human_baboon_lrt <- read.table(file.path(lrt_root, "human_baboon_cd4_lrt.txt"),
                                sep = "\t", header = TRUE)

human_baboon_lrt$baboonCD4_sdGroup[is.na(human_baboon_lrt$baboonCD4_sdGroup)] <- "Sharp"
human_baboon_lrt$pause_switch <- paste0(human_baboon_lrt$baboonCD4_sdGroup, "_to_",
                                         human_baboon_lrt$humanCD4_sdGroup)
human_baboon_lrt$nuc_positioning_lfc <- compute_nuc_entropy_lfc(
  tss_df    = human_baboon_lrt,
  sp1_name  = "baboon", sp1_geneId = baboon_cd4_nuc$rate_df$geneId, sp1_raw = baboon_cd4_nuc$raw,
  sp2_name  = "human",  sp2_geneId = human_cd4_nuc$rate_df$geneId,  sp2_raw = human_cd4_nuc$raw,
  scale_factor = BABOON_READS / HUMAN_READS
)

human_baboon_lrt$pause_switch[
  human_baboon_lrt$pause_switch %in% c("Broad_to_Broad", "Sharp_to_Sharp")] <- "Unchanged"
human_baboon_lrt <- human_baboon_lrt %>%
  mutate(pause_switch = factor(pause_switch, levels = c("Sharp_to_Broad", "Unchanged", "Broad_to_Sharp")))

human_baboon_stat <- human_baboon_lrt %>%
  group_by(tssType) %>%
  pairwise_wilcox_test(nuc_positioning_lfc ~ pause_switch, comparisons = MY_COMP) %>%
  ungroup()
human_baboon_stat$y.position <- c(1, 1, 1, 1)

human_baboon_hline_df <- human_baboon_lrt %>%
  group_by(tssType) %>%
  summarise(yint = median(nuc_positioning_lfc[pause_switch == "Unchanged"], na.rm = TRUE),
            .groups = "drop")

human_baboon_delta_p <- ggplot(human_baboon_lrt,
    aes(x = pause_switch, y = nuc_positioning_lfc, fill = pause_switch)) +
  facet_wrap(~ tssType, ncol = 1) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.6) +
  stat_pvalue_manual(human_baboon_stat, label = "p.adj.signif") +
  labs(x = "", y = "", fill = expression(sigma * " change")) +
  scale_fill_manual(values = FILL_COLORS) +
  coord_cartesian(ylim = c(-1, 1.4)) +
  my_boxplot_style() +
  ggtitle("Human/Baboon") +
  geom_hline(data = human_baboon_hline_df, aes(yintercept = yint),
             linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
  theme(legend.position  = "bottom",
        strip.background = element_blank(),
        strip.text       = element_blank(),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5, face = "plain"))


# ── 3. Rhesus vs Baboon ──────────────────────────────────────────────────────
# sp1 = baboon (denominator), sp2 = rhesus (numerator):
#   positive nuc_positioning_lfc → rhesus more organised than baboon
#   pause_switch = baboonCD4_sdGroup → rhesusCD4_sdGroup
rhesus_baboon_lrt <- read.table(file.path(lrt_root, "rhesus_baboon_cd4_lrt.txt"),
                                 sep = "\t", header = TRUE)

rhesus_baboon_lrt$baboonCD4_sdGroup[is.na(rhesus_baboon_lrt$baboonCD4_sdGroup)] <- "Sharp"
rhesus_baboon_lrt$pause_switch <- paste0(rhesus_baboon_lrt$baboonCD4_sdGroup, "_to_",
                                          rhesus_baboon_lrt$rhesusCD4_sdGroup)
rhesus_baboon_lrt$nuc_positioning_lfc <- compute_nuc_entropy_lfc(
  tss_df    = rhesus_baboon_lrt,
  sp1_name  = "baboon", sp1_geneId = baboon_cd4_nuc$rate_df$geneId, sp1_raw = baboon_cd4_nuc$raw,
  sp2_name  = "rhesus", sp2_geneId = rhesus_cd4_nuc$rate_df$geneId, sp2_raw = rhesus_cd4_nuc$raw,
  scale_factor = BABOON_READS / RHESUS_READS
)

rhesus_baboon_lrt$pause_switch[
  rhesus_baboon_lrt$pause_switch %in% c("Broad_to_Broad", "Sharp_to_Sharp")] <- "Unchanged"
rhesus_baboon_lrt <- rhesus_baboon_lrt %>%
  mutate(pause_switch = factor(pause_switch, levels = c("Sharp_to_Broad", "Unchanged", "Broad_to_Sharp")))

rhesus_baboon_stat <- rhesus_baboon_lrt %>%
  group_by(tssType) %>%
  pairwise_wilcox_test(nuc_positioning_lfc ~ pause_switch, comparisons = MY_COMP) %>%
  ungroup()
rhesus_baboon_stat$y.position <- c(1, 1, 1, 1)

rhesus_baboon_hline_df <- rhesus_baboon_lrt %>%
  group_by(tssType) %>%
  summarise(yint = median(nuc_positioning_lfc[pause_switch == "Unchanged"], na.rm = TRUE),
            .groups = "drop")

rhesus_baboon_delta_p <- ggplot(rhesus_baboon_lrt,
    aes(x = pause_switch, y = nuc_positioning_lfc, fill = pause_switch)) +
  facet_wrap(~ tssType, ncol = 1, strip.position = "right") +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.6) +
  stat_pvalue_manual(rhesus_baboon_stat, label = "p.adj.signif") +
  labs(x = "", y = "", fill = expression(sigma * " change")) +
  scale_fill_manual(values = FILL_COLORS) +
  coord_cartesian(ylim = c(-1, 1.4)) +
  my_boxplot_style() +
  geom_hline(data = rhesus_baboon_hline_df, aes(yintercept = yint),
             linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE) +
  ggtitle("Rhesus/Baboon") +
  theme(axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        legend.position  = "none",
        strip.background = element_blank(),
        strip.text.y.right = element_text(size = 10, angle = 270),
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5, face = "plain"))


# ── Save combined figure ──────────────────────────────────────────────────────
ggsave(filename = file.path(result_dir, "all_cd4_comparison.pdf"),
       plot = human_rhesus_delta_p | human_baboon_delta_p | rhesus_baboon_delta_p,
       width = 6, height = 4)
