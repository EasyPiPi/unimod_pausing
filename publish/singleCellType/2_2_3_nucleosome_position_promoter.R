library(tidyverse)
library(plyranges)
library(cowplot)
library(RColorBrewer)
library(STADyUM)
library(UpSetR)

# ---------------------------------------------------------------------------- #
#  Paths
# ---------------------------------------------------------------------------- #

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/")

human_cd4_rate_in       <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd4_rate.RDS")
human_cd14_rate_in      <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd14_rate.RDS")
human_cd4_ns_matrix_in  <- file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/human/cd4_ns_matrix.RDS")
human_cd14_ns_matrix_in <- file.path(.paths$outputs, "publish/singleCellType/2_nucleosomePositioning/human/cd14_ns_matrix.RDS")

# ---------------------------------------------------------------------------- #
#  Helper functions
# ---------------------------------------------------------------------------- #

get_promoter_motif_from_stadyum_obj <- function(stadyum_obj, promoter_df) {
  # This function links pause regions from a STADyUM object to nearby gene promoters.
  # For each pause site, it identifies the closest TSS within 200 bp on the same strand.
  # Promoters are then classified into motif classes based on core promoter elements.

  tsn <- stadyum_obj@pauseRegions %>% anchor_5p() %>% mutate(width = 1)
  tsn_df <- as.data.frame(tsn) %>%
    as_tibble() %>%
    dplyr::rename(seqname = seqnames) %>%
    mutate(
      seqname = paste0("chr", seqname),
      tsn_pos = start,
      tsn_id  = dplyr::row_number()
    )

  promoter_df <- promoter_df %>%
    dplyr::rename(gene_id = gene_ensembl_merged)

  pairs <- tsn_df %>%
    inner_join(promoter_df, by = c("seqname", "gene_id", "strand")) %>%
    mutate(dist_bp = abs(tsn_pos - TSS))

  nearest <- pairs %>%
    group_by(tsn_id) %>%
    slice_min(order_by = dist_bp, with_ties = FALSE) %>%
    ungroup() %>%
    filter(dist_bp <= 200)

  nearest <- nearest %>%
    mutate(
      motif_class = case_when(
        TATA.box == 1 ~ "TATA-box",
        CCAAT.box == 1 & GC.box == 1 & Inr == 0 & TATA.box == 0 ~ "GC-CCAAT",
        CCAAT.box == 1 & GC.box == 0 & Inr == 0 & TATA.box == 0 ~ "CCAAT only",
        GC.box == 1 & Inr == 0 & TATA.box == 0 & CCAAT.box == 0 ~ "GC-box only",
        GC.box == 1 & Inr == 1 & TATA.box == 0 & CCAAT.box == 0 ~ "GC-Inr",
        Inr == 1 & GC.box == 0 & TATA.box == 0 & CCAAT.box == 0 ~ "Inr only",
        TRUE ~ "Other"
      )
    )
  return(nearest)
}

mycolor <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB",
             "#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00")

plot_beta_group_line_promoter <- function(count_matrix, rate_df, group_col, range = NULL) {

  df_long <- make_beta_group_df_long(
    count_matrix = count_matrix,
    rate_df      = rate_df,
    group_col    = group_col,
    range        = range,
    drop_groups  = c("Unknown", "Other")
  )

  plot_beta_group_line_core(
    df_long       = df_long,
    color_values  = mycolor,
    legend_inside = c(0.6, 0.8)
  )
}


# Plots nucleosome profiles faceted by motif class, with lines colored by beta group.
plot_beta_group_line_faceted_with_promoter <- function(count_matrix,
                                                        rate_df,
                                                        group_col    = "beta_group",
                                                        class_col    = "motif_class",
                                                        range        = NULL,
                                                        colors       = NULL,
                                                        x_label      = c(0, 200),
                                                        shift        = 1001,
                                                        legend_title = expression(beta ~ group),
                                                        y_label      = "Mean Micro-C signal") {
  stopifnot(nrow(count_matrix) == nrow(rate_df))

  df       <- as.data.frame(count_matrix)
  df$group <- rate_df[[group_col]]
  df$class <- factor(rate_df[[class_col]])

  group_means_df <- df %>%
    group_by(class, group) %>%
    summarise(across(starts_with("V"), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  df_long <- group_means_df %>%
    pivot_longer(cols = starts_with("V"), names_to = "position_str", values_to = "value") %>%
    mutate(position = as.numeric(str_remove(position_str, "V")) - shift)

  if (!is.null(range)) {
    df_long <- df_long %>% filter(position >= range[1], position <= range[2])
  }
  df_long <- df_long %>% filter(!class %in% c("Unknown", "Other") & !is.na(class))

  if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(6, "Blues")[2:6]
  }

  ggplot(df_long, aes(x = position, y = value, color = group)) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ class, ncol = 6) +
    labs(x = "Distance from TSS (bp)", y = y_label, color = legend_title) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = x_label) +
    scale_y_continuous(n.breaks = 4) +
    theme_cowplot() +
    theme(
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      legend.position  = "right",
      axis.ticks       = element_line(color = "black", size = 0.3),
      axis.line        = element_line(color = "black", size = 0.3),
      legend.box       = "horizontal",
      axis.text        = element_text(size = 11),
      strip.text       = element_text(size = 9)
    )
}

# ---------------------------------------------------------------------------- #
#  Per-cell-type analysis
# ---------------------------------------------------------------------------- #

run_promoter_nucleosome_analysis <- function(cell_type, rate_rds, ns_matrix_rds,
                                              promoter_motif, out_dir) {
  ns_matrix <- readRDS(ns_matrix_rds)
  rate_obj  <- readRDS(rate_rds)

  # Extract promoter motif annotations for this cell type
  promoter_info <- get_promoter_motif_from_stadyum_obj(rate_obj, promoter_motif)

  # Join motif info into the nucleosome rate data frame
  rate_df <- ns_matrix$rate_df %>%
    left_join(
      promoter_info %>%
        dplyr::select(
          gene_id, TSS_nearest = TSS, dist_bp,
          promoter_index, TATA.box, Inr, CCAAT.box, GC.box, motif_class
        ),
      by = c("geneId" = "gene_id")
    )

  # UpSet plot for promoter motif combinations
  promoter_info_plot <- promoter_info %>%
    dplyr::rename(
      "TATA box"  = TATA.box,
      "Initiator" = Inr,
      "CCAAT"     = CCAAT.box,
      "GC box"    = GC.box
    )
  pdf(
    paste0(file.path(out_dir, paste0("human/", cell_type, "_motif_info")), ".pdf"),
    width = 4, height = 3
  )
  print(upset(
    as.data.frame(promoter_info_plot),
    sets        = c("TATA box", "Initiator", "CCAAT", "GC box"),
    nintersects = NA,
    order.by    = "freq"
  ))
  dev.off()

  # Motif-grouped nucleosome plot (single panel, colored by motif class)
  rate_df$motif_class <- fct_relevel(
    rate_df$motif_class,
    "CCAAT only", "GC-CCAAT", "GC-box only", "GC-Inr", "Inr only", "TATA-box"
  )
  p_motif <- plot_beta_group_line_promoter(
    ns_matrix$smooth, rate_df, 'motif_class', range = c(0, 1000)
  )
  ggsave(
    paste0(file.path(out_dir, paste0("human/", cell_type, "_motif_nuc_merged")), ".pdf"),
    p_motif, width = 5, height = 4
  )

  # Beta-group faceted nucleosome plot (faceted by motif class, colored by beta group)
  p_beta <- plot_beta_group_line_faceted_with_promoter(
    ns_matrix$smooth, rate_df,
    range     = c(0, 1000),
    x_label   = c(0, 500),
    group_col = 'betaGroup'
  )
  ggsave(
    paste0(file.path(out_dir, paste0("human/", cell_type, "_Micro-C_beta_gb_motif")), ".pdf"),
    p_beta, width = 12, height = 3
  )
}

# ---------------------------------------------------------------------------- #
#  Shared annotation data
# ---------------------------------------------------------------------------- #

promoter_motif <- read.table(
  file.path(.paths$data, "promoter_annotation/merge_promoter_motif.txt"),
  sep = '\t', header = TRUE
)

# ---------------------------------------------------------------------------- #
#  Run analysis for each cell type
# ---------------------------------------------------------------------------- #

run_promoter_nucleosome_analysis(
  cell_type     = "cd4",
  rate_rds      = human_cd4_rate_in,
  ns_matrix_rds = human_cd4_ns_matrix_in,
  promoter_motif = promoter_motif,
  out_dir       = result_dir
)

run_promoter_nucleosome_analysis(
  cell_type     = "cd14",
  rate_rds      = human_cd14_rate_in,
  ns_matrix_rds = human_cd14_ns_matrix_in,
  promoter_motif = promoter_motif,
  out_dir       = result_dir
)
