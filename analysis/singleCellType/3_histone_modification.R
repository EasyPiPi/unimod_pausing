# ============================================================
# 4_histone_modification.R
# Analyzes promoter histone modification signals for
# CD4 regulatory cells and CD14 monocytes.
# Self-contained — no external helper files required.
# ============================================================

# ---- Libraries ----
library(STADyUM)
library(rtracklayer)
library(plyranges)
library(tidyverse)
library(cowplot)
library(ggpubr)

# ---- Paths ----
root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))

cd4_rate_in     <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd4_rate.RDS")
cd14_rate_in    <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd14_rate.RDS")
histone_meta_in <- file.path(.paths$data, "encode/meta_fc.RDS")
result_dir      <- file.path(.paths$outputs, "publish/singleCellType/4_histoneModification/")

# ---- Input checks ----
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)
for (f in c(cd4_rate_in, cd14_rate_in, histone_meta_in)) {
  if (!file.exists(f)) stop("Input file not found: ", f)
}

# ---- Constants ----
HISTONE_MARKS <- c("H3K4me3", "H3K27ac")
PROMOTER_LEN  <- 500
CHI_LABELS    <- c("Low expression", "Medium expression", "High expression")

# ---- Cell-type configuration ----
# To add a new cell type, append an entry here.
# strict_p: use "p < 2e-16" floor instead of ggpubr's default p.format.
cell_configs <- list(
  cd4 = list(
    rate_in   = cd4_rate_in,
    biosample = "CD4 regulatory",
    prefix    = "cd4",
    save_rds  = TRUE,
    ylims     = list(H3K4me3 = c(2, 9),    H3K27ac = c(-0.5, 6.5)),
    label_x   = list(H3K4me3 = 1.5,        H3K27ac = 1.5),
    strict_p  = list(H3K4me3 = FALSE,       H3K27ac = FALSE)
  ),
  cd14 = list(
    rate_in   = cd14_rate_in,
    biosample = "CD14 monocyte",
    prefix    = "cd14",
    save_rds  = FALSE,
    ylims     = list(H3K4me3 = c(2, 10),   H3K27ac = c(2, 8)),
    label_x   = list(H3K4me3 = 1.2,        H3K27ac = 1.5),
    strict_p  = list(H3K4me3 = TRUE,        H3K27ac = FALSE)
  )
)

# ============================================================
# Helper functions
# ============================================================

get_promoter_df_from_stadyum_obj <- function(stadyum_obj, promoter_len = 500) {
  tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1)
  tsn <- tsn[tsn$gene_id %in% rates(stadyum_obj)$geneId]
  promoter_grng <- resize(tsn, promoter_len, fix = "center")

  promoter_grng_df <- as.data.frame(promoter_grng) %>%
    left_join(rates(stadyum_obj), by = c("gene_id" = "geneId"))

  promoter_grng_df$pauseCount <- sapply(promoter_grng_df$actualPauseSiteCounts, sum)
  promoter_grng_df <- promoter_grng_df %>% dplyr::select(-actualPauseSiteCounts)

  list(df = promoter_grng_df, grang = promoter_grng)
}

get_histone_modification_sum <- function(promoter_df, promoter_grng, meta_sel) {
  for (i in seq_len(nrow(meta_sel))) {
    histone_bw <- import.bw(meta_sel$bw[i])
    hits       <- findOverlaps(promoter_grng, histone_bw)

    sum_scores <- data.frame(
      promoter_idx  = queryHits(hits),
      histone_score = mcols(histone_bw)$score[subjectHits(hits)]
    ) %>%
      group_by(promoter_idx) %>%
      summarize(histone_sum_scores = mean(histone_score), .groups = "drop")

    promoter_df[, meta_sel$Target[i]] <- sum_scores$histone_sum_scores
  }
  promoter_df
}

my_boxplot_style <- function() {
  list(
    geom_boxplot(outlier.shape = NA, linewidth = 0.3),
    scale_fill_brewer(palette = "Blues"),
    theme_cowplot(),
    theme(
      legend.position  = "none",
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      axis.title.y     = element_text(size = 12),
      axis.title.x     = element_text(size = 12),
      axis.text        = element_text(size = 11),
      strip.text       = element_text(size = 9),
      axis.line        = element_line(size = 0.3),
      axis.ticks       = element_line(size = 0.3)
    )
  )
}

# ============================================================
# Analysis functions
# ============================================================

process_cell_type <- function(cfg, meta_all) {
  meta_sel <- meta_all %>%
    filter(Biosample %in% cfg$biosample, Target %in% HISTONE_MARKS)

  rate_obj <- readRDS(cfg$rate_in)
  promoter <- get_promoter_df_from_stadyum_obj(rate_obj, promoter_len = PROMOTER_LEN)
  promoter$df <- get_histone_modification_sum(promoter$df, promoter$grang, meta_sel)

  levels(promoter$df$chiGroup) <- CHI_LABELS
  promoter
}

plot_histone_boxplot <- function(df, histone, cfg, result_dir) {
  if (!histone %in% colnames(df)) {
    stop("Column '", histone, "' not found in promoter data frame for ", cfg$prefix, ".")
  }

  ylim   <- cfg$ylims[[histone]]
  lx     <- cfg$label_x[[histone]]
  strict <- cfg$strict_p[[histone]]

  p_label_aes <- if (strict) {
    aes(label = ifelse(
      after_stat(p) < 2e-16,
      "p < 2e-16",
      paste0("p = ", after_stat(p.format))
    ))
  } else {
    aes(label = gsub("<", "p < ", after_stat(p.format)))
  }

  y_labels <- list(
    H3K4me3 = expression(log[2] ~ "(H3K4me3)"),
    H3K27ac = expression(log[2] ~ "(H3K27ac)")
  )

  p <- ggplot(df, aes(x = betaGroup, y = log2(.data[[histone]]), fill = betaGroup)) +
    facet_wrap(~chiGroup, nrow = 1, ncol = 3) +
    stat_compare_means(label.x = lx, mapping = p_label_aes) +
    labs(
      x    = expression(beta ~ "Groups"),
      y    = y_labels[[histone]],
      fill = expression(beta ~ group)
    ) +
    coord_cartesian(ylim = ylim) +
    my_boxplot_style() +
    theme(
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      legend.position = "right"
    )

  out_path <- file.path(result_dir, paste0(cfg$prefix, "_", histone, "_beta_chi.pdf"))
  ggsave(out_path, p, width = 8, height = 2.5)
  invisible(p)
}

# ============================================================
# Main
# ============================================================

meta_fc <- readRDS(histone_meta_in)
# Override absolute bw paths stored in the RDS with portable equivalents
meta_fc$bw <- file.path(.paths$encode_bw, basename(meta_fc$bw))

for (cfg in cell_configs) {
  message("Processing: ", cfg$prefix)
  promoter <- process_cell_type(cfg, meta_fc)

  if (cfg$save_rds) {
    saveRDS(promoter, file.path(result_dir, paste0(cfg$prefix, "_histone.rds")))
  }

  for (histone in HISTONE_MARKS) {
    plot_histone_boxplot(promoter$df, histone, cfg, result_dir)
  }
}

message("Done.")
