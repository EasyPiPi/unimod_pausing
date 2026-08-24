# ===========================================================================
# 3_4_GCinfo_metaplot.R
# Metaplot of mean GC content and GC skew around TSS across species/cell types
# ===========================================================================


# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------

library(ggplot2)
library(cowplot)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", "."),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

seq_dir    <- file.path(.paths$outputs, "publish/betweenSpecies/sequence")
result_dir <- file.path(.paths$outputs, "publish/betweenSpecies")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- GC content input files (one per species / cell-type combination) ---

gc_content_files <- list(
  "Human CD4"   = file.path(seq_dir, "human_cd4_tss_sort_1000bp_extend_100_gc_content.csv"),
  "Human CD14"  = file.path(seq_dir, "human_cd14_tss_sort_1000bp_extend_100_gc_content.csv"),
  "Rhesus CD4"  = file.path(seq_dir, "rhesus_cd4_tss_sort_1000bp_extend_100_gc_content.csv"),
  "Rhesus CD14" = file.path(seq_dir, "rhesus_cd14_tss_sort_1000bp_extend_100_gc_content.csv"),
  "Baboon CD4"  = file.path(seq_dir, "baboon_cd4_tss_sort_1000bp_extend_100_gc_content.csv"),
  "Baboon CD14" = file.path(seq_dir, "baboon_cd14_tss_sort_1000bp_extend_100_gc_content.csv")
)

# --- GC skew input files ---

gc_skew_files <- list(
  "Human CD4"   = file.path(seq_dir, "human_cd4_tss_sort_1000bp_extend_100_gc_skew.csv"),
  "Human CD14"  = file.path(seq_dir, "human_cd14_tss_sort_1000bp_extend_100_gc_skew.csv"),
  "Rhesus CD4"  = file.path(seq_dir, "rhesus_cd4_tss_sort_1000bp_extend_100_gc_skew.csv"),
  "Rhesus CD14" = file.path(seq_dir, "rhesus_cd14_tss_sort_1000bp_extend_100_gc_skew.csv"),
  "Baboon CD4"  = file.path(seq_dir, "baboon_cd4_tss_sort_1000bp_extend_100_gc_skew.csv"),
  "Baboon CD14" = file.path(seq_dir, "baboon_cd14_tss_sort_1000bp_extend_100_gc_skew.csv")
)

# --- Shared plotting constants ---

# 2000 bp window centered on TSS, 100 bp sliding step -> 1901 windows
TSS_POSITIONS <- -950:950

# Colors assigned by name so order in the legend is always explicit
SAMPLE_COLORS <- c(
  "Human CD4"   = "#DC143C",
  "Human CD14"  = "#0000FF",
  "Rhesus CD4"  = "#20B2AA",
  "Rhesus CD14" = "#FFA500",
  "Baboon CD4"  = "#9370DB",
  "Baboon CD14" = "#98FB98"
)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Read a sliding-window GC CSV (gene names in col 1) and return column means
read_gc_means <- function(path) {
  df <- read.table(path, sep = ",")
  colMeans(df[, -1])
}

# Combine a named list of per-window means into a long data frame for ggplot
build_metaplot_df <- function(means_list, positions) {
  n <- length(means_list[[1]])
  data.frame(
    Sample   = rep(names(means_list), each = n),
    Value    = unlist(means_list, use.names = FALSE),
    Position = rep(positions, times = length(means_list))
  )
}

# Build a TSS metaplot from a long data frame produced by build_metaplot_df()
make_metaplot <- function(df, y_label, colors) {
  ggplot(df, aes(x = Position, y = Value, color = Sample)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = colors) +
    theme_cowplot() +
    labs(x = "Position relative to TSS (bp)", y = y_label) +
    theme(
      axis.title   = element_text(size = 12),
      axis.text    = element_text(size = 12),
      strip.text   = element_text(size = 12),
      axis.line    = element_line(linewidth = 0.3),
      axis.ticks   = element_line(linewidth = 0.3),
      legend.title = element_blank()
    )
}


# ---------------------------------------------------------------------------
# GC content: load, summarize, and build plot data frame
# ---------------------------------------------------------------------------

gc_content_means <- lapply(gc_content_files, read_gc_means)
gc_content_df    <- build_metaplot_df(gc_content_means, TSS_POSITIONS)


# ---------------------------------------------------------------------------
# GC skew: load, summarize, and build plot data frame
# ---------------------------------------------------------------------------

gc_skew_means <- lapply(gc_skew_files, read_gc_means)
gc_skew_df    <- build_metaplot_df(gc_skew_means, TSS_POSITIONS)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

p_gc_content <- make_metaplot(gc_content_df, y_label = "GC content", colors = SAMPLE_COLORS)
p_gc_skew    <- make_metaplot(gc_skew_df,    y_label = "GC skew",    colors = SAMPLE_COLORS)


# ---------------------------------------------------------------------------
# Save figures
# ---------------------------------------------------------------------------

ggsave(file.path(result_dir, "all_gc_content.pdf"), p_gc_content, width = 5, height = 3)
ggsave(file.path(result_dir, "all_gc_skew.pdf"),    p_gc_skew,    width = 5, height = 3)
