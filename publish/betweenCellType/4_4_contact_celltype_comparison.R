library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# ============================================================
# Enhancer-promoter contact comparison: CD4 vs CD14
#
# Normalization strategy follows Barshad et al. 2023 (Nature Genetics)
# K562 vs Jurkat cell-line comparison logic, NOT the treatment logic.
#
# Key distinction (from paper Methods and README):
#   Treatment comparison (DMSO/FLV/TRP, same cell line):
#     keep pairs with CPB > threshold in AT LEAST ONE condition (OR)
#   Cell-line comparison (K562 vs Jurkat, or here CD4 vs CD14):
#     keep pairs with CPB > threshold in BOTH cell types (AND)
#     "to avoid ascertainment bias" — Barshad et al. 2023
# ============================================================

# ---- Paths ------------------------------------------------------------------

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/betweenCellType/EP/")
human_dir  <- file.path(result_dir, "human")
dir.create(human_dir, recursive = TRUE, showWarnings = FALSE)

cd4_ep_in  <- file.path(.paths$data, "ep_contacts/cd4_EP_and_BG_contacts.txt")
cd14_ep_in <- file.path(.paths$data, "ep_contacts/cd14_EP_and_BG_contacts.txt")
lrt_in     <- file.path(result_dir, "human_lrt_withFilename.txt")

# ---- Helper functions -------------------------------------------------------

my_boxplot_style <- function() {
  list(
    geom_boxplot(outlier.shape = NA, linewidth = 0.3),
    theme_cowplot(),
    theme(legend.position = "none",
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.text  = element_text(size = 14),
      strip.text = element_text(size = 14),
      axis.line = element_line(size = 0.3),
      axis.ticks = element_line(size = 0.3)))
}

# ---- Input checks -----------------------------------------------------------

for (f in c(cd4_ep_in, cd14_ep_in, lrt_in)) {
  if (!file.exists(f)) stop("Required input file not found: ", f)
}

# ---- Parameters -------------------------------------------------------------

ep_contact_threshold <- 8   # contacts per billion (CPB) minimum

# Sequencing depths: total cis contacts with MAPQ >= 30
# (matches the denominator used in Barshad et al. Python script: / 10,000,000,000)
cd4_total_contacts  <- 485936844
cd14_total_contacts <- 392861612

cd4_scale_factor  <- cd4_total_contacts  / 1e10
cd14_scale_factor <- cd14_total_contacts / 1e10

# ---- Load and label EP contact files ----------------------------------------

ep_col_names <- c("seqname", "promoter", "enhancer", "ep_count", "bg_count")

cd4_ep <- read.table(cd4_ep_in, sep = "\t", col.names = ep_col_names, colClasses = c(seqname = "character")) %>%
  mutate(ep_bg_ratio = ep_count / bg_count)

cd14_ep <- read.table(cd14_ep_in, sep = "\t", col.names = ep_col_names, colClasses = c(seqname = "character")) %>%
  mutate(ep_bg_ratio = ep_count / bg_count)

# ---- Merge on shared anchor pairs -------------------------------------------
# inner join: only pairs present in both cell types

merged_ep <- inner_join(
  cd4_ep,
  cd14_ep,
  by     = c("seqname", "promoter", "enhancer"),
  suffix = c("_cd4", "_cd14")
)

message(sprintf("E-P pairs present in both cell types (before CPB filter): %d", nrow(merged_ep)))

# ---- CPB normalization ------------------------------------------------------
# CPB = contacts per billion  (EP count / scale factor, where scale = depth / 1e10)

merged_ep <- merged_ep %>%
  mutate(
    cd4_ep_cpb  = ep_count_cd4  / cd4_scale_factor,
    cd14_ep_cpb = ep_count_cd14 / cd14_scale_factor
  )

# ---- Cell-line CPB filter (AND logic) ---------------------------------------
#
# IMPORTANT: Barshad et al. use AND here, not OR.
# From the paper Methods (highlighted):
#   "Because TSS calling data (PRO-seq and coPRO-capped) were more abundant
#    for K562 than Jurkat, when comparing K562 and Jurkat libraries we
#    considered enhancer-promoter pairs with at least eight CPB in BOTH
#    cell lines, to avoid ascertainment bias."
#
# The treatment-comparison script (Compering_EP_contacts_between_treatments.py)
# uses OR because one condition may legitimately have depleted contacts.
# For a cell-type comparison there is no such asymmetry, so AND is correct.

filtered_ep <- merged_ep %>%
  dplyr::filter(
    cd4_ep_cpb  > ep_contact_threshold &   # AND — both cell types must pass
    cd14_ep_cpb > ep_contact_threshold
  ) %>%
  drop_na()

message(sprintf("E-P pairs after CPB > %d AND filter:               %d",
                ep_contact_threshold, nrow(filtered_ep)))
message(sprintf("Pairs removed by CPB filter:                        %d",
                nrow(merged_ep) - nrow(filtered_ep)))

# ---- Log2 fold change of EP/background ratio --------------------------------
# Positive values = stronger contact in CD14 relative to CD4

filtered_ep <- filtered_ep %>%
  mutate(
    ep_lfc   = log2(ep_bg_ratio_cd14 / ep_bg_ratio_cd4),
    file_name = paste0(seqname, "_", promoter, ".csv")
  )

# ---- Load LRT annotations and join ------------------------------------------

human_lrt <- read.table(lrt_in, sep = "\t", header = TRUE)

ep_lrt <- left_join(filtered_ep, human_lrt, by = "file_name")

n_missing_lrt <- sum(is.na(ep_lrt$chi_lfc_group))
message(sprintf("Rows missing LRT annotation after left_join:        %d / %d",
                n_missing_lrt, nrow(ep_lrt)))

# ---- Required columns check -------------------------------------------------

required_cols <- c("chi_lfc_group", "betaCategory", "pauseCount_lfc_group",
                   "CD4_sdGroup", "CD14_sdGroup")
missing_cols <- base::setdiff(required_cols, colnames(ep_lrt))
if (length(missing_cols) > 0) {
  stop("Required columns missing after merge: ", paste(missing_cols, collapse = ", "))
}

# ---- Factor ordering for chi_lfc_group --------------------------------------

ep_lrt <- ep_lrt %>%
  mutate(
    chi_lfc_group = factor(chi_lfc_group,
                           levels = c("Decrease", "Unchanged", "Increase"))
  ) %>%
  droplevels()

# ---- Drop NAs in columns used for plotting ----------------------------------

ep_lrt <- ep_lrt %>%
  tidyr::drop_na(ep_lfc, chi_lfc_group, betaCategory,
                 pauseCount_lfc_group, CD4_sdGroup, CD14_sdGroup)

message(sprintf("Rows available for plotting (after NA removal):     %d", nrow(ep_lrt)))

# ============================================================
# PLOTS
# ============================================================

# ---- Plot 0: chi_lfc_group vs ep_lfc ----------------------------------------

p0 <- ep_lrt %>%
  ggplot(aes(x = chi_lfc_group, y = ep_lfc, fill = chi_lfc_group)) +
  stat_compare_means(
    label.x = 1.2, label.y = 1.5,
    aes(label = gsub("<", "p < ", after_stat(p.format)))
  ) +
  labs(x = "", y = expression(log[2]*FC ~ "(E-P contact)")) +
  coord_cartesian(ylim = c(-1.5, 2)) +
  my_boxplot_style() +
  scale_fill_brewer(palette = "Greens") +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank()
  )

ggsave(file.path(human_dir, "chi_EP.pdf"), p0, width = 2.5, height = 4)

# ---- Recode chi_lfc_group labels for remaining plots ------------------------

ep_lrt <- ep_lrt %>%
  mutate(
    chi_lfc_group = recode(chi_lfc_group,
      "Decrease"  = "Decreased expression",
      "Unchanged" = "Unchanged expression",
      "Increase"  = "Increased expression"
    )
  )

# ---- Plot 1: betaCategory vs ep_lfc, faceted by expression group ------------

p1 <- ep_lrt %>%
    ggplot(aes(x = betaCategory, y = ep_lfc, fill = betaCategory)) +
    geom_violin(width = 0.9, alpha = 0.5, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", linewidth = 0.3) +
    facet_wrap(~ chi_lfc_group, ncol = 3) +
    ylab(expression(log[2]*FC ~ "(E-P contact)")) +
    scale_fill_manual(values = c("Others" = "#BFBFBF",
                                 "Up"     = "#E45756",
                                 "Down"   = "#4C78A8")) +
    labs(fill = expression(beta * " change")) +
    xlab("") +
    coord_cartesian(ylim = c(-3, 3)) +
    cowplot::theme_cowplot() +
    theme(
      strip.background     = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      axis.text.x          = element_blank(),
      axis.ticks.x         = element_blank(),
      legend.position      = "bottom",
      legend.justification = "center"
    )

ggsave(file.path(human_dir, "beta_groupby_chi_EP.pdf"), p1, width = 7, height = 3)

# ---- Plot 2: pauseCount_lfc_group vs ep_lfc, faceted by expression group ----

p2 <- ep_lrt %>%
    ggplot(aes(x = pauseCount_lfc_group, y = ep_lfc, fill = pauseCount_lfc_group)) +
    geom_violin(width = 0.9, alpha = 0.5, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, color = "black", linewidth = 0.3) +
    stat_compare_means(label.x = 1.2, label.y = 2.5, label = "p") +
    facet_wrap(~ chi_lfc_group, ncol = 3) +
    ylab(expression(log[2]*FC ~ "(E-P contact)")) +
    coord_cartesian(ylim = c(-3, 3)) +
    scale_fill_brewer(palette = "RdPu")  +
    labs(fill = "Paused count\nchange") +
    xlab("") +
    cowplot::theme_cowplot() +
    theme(
      strip.background     = element_rect(fill = "white", colour = "black", linewidth = 0.4),
      axis.text.x          = element_blank(),
      axis.ticks.x         = element_blank(),
      legend.position      = "bottom",
      legend.justification = "center"
    )

ggsave(file.path(human_dir, "pauseCount_groupby_chi_EP.pdf"), p2, width = 7, height = 3)


message("Done. Plots saved to: ", human_dir)
