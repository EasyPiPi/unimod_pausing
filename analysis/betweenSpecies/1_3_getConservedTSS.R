library(tidyverse)
library(STADyUM)
library(GenomicRanges)
library(cowplot)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

out_dir <- file.path(.paths$outputs, "publish/betweenSpecies")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper functions ──────────────────────────────────────────────────────────

# Base-pair overlap between two genomic intervals; 0 if on different chromosomes
overlap_len <- function(chr1, s1, e1, chr2, s2, e2) {
  ifelse(chr1 == chr2, pmax(0L, pmin(e1, e2) - pmax(s1, s2)), 0L)
}

# Histogram of overlap distances with a 150 bp classification threshold
plot_histogram <- function(overlap_vec, save_path) {
  p <- ggplot(data.frame(Overlap = overlap_vec), aes(x = Overlap)) +
    geom_histogram(bins = 100) +
    ylab("Count") +
    geom_vline(xintercept = 150, linetype = "dashed", linewidth = 0.6) +
    theme_cowplot()
  ggsave(save_path, p, width = 4, height = 2)
}

# Extract pause regions from a STADyUM rate RDS; join with ortholog groups
get_human_pause_regions <- function(rate_path, ortho_genes) {
  gr <- readRDS(rate_path)@pauseRegions
  mcols(gr)          <- gr$gene_id
  names(mcols(gr))   <- "human"
  seqlevelsStyle(gr) <- "UCSC"
  as.data.frame(gr) %>%
    inner_join(ortho_genes[, c("group", "human")], by = "human")
}

# Parse a liftover BED file: keep 1-to-1 genes and join with ortholog groups.
# The BED name field (V4) has format chr:range:gene_id.
parse_liftover_bed <- function(bed_path, ortho_genes, species_col) {
  read.table(bed_path) %>%
    separate(V4, into = c("chr_target", "range", "gene"), sep = ":", remove = FALSE) %>%
    add_count(gene) %>%
    dplyr::filter(n == 1) %>%
    dplyr::select(-n) %>%
    inner_join(ortho_genes[, c("group", species_col)], by = c("gene" = species_col))
}

# Join human pause regions with one species' liftover BED, compute overlap,
# save histogram, classify TSS, and return a tidy coordinate + category table.
#   species_col : "baboon" or "rhesus"  (becomes the gene-ID column name)
#   cat_col     : "hb_tss_category" or "hr_tss_category"
#   prefix      : "b" or "r"  (prefixes the lifted coordinate columns)
compute_species_overlap <- function(human_pause, liftover_df,
                                    species_col, cat_col, prefix,
                                    hist_path) {
  joined <- inner_join(human_pause, liftover_df, by = "group") %>%
    mutate(dist = overlap_len(seqnames, start, end, V1, V2, V3))

  plot_histogram(joined$dist, hist_path)

  classified <- joined %>%
    mutate(!!sym(cat_col) := case_when(
      dist > 150 ~ "Conserved TSS",
      TRUE       ~ "non-Conserved TSS"
    ))

  # Rename with base R to avoid S4Vectors masking dplyr::rename
  names(classified)[names(classified) == "gene"] <- species_col
  names(classified)[names(classified) == "V1"]   <- paste0(prefix, "_chr")
  names(classified)[names(classified) == "V2"]   <- paste0(prefix, "_start")
  names(classified)[names(classified) == "V3"]   <- paste0(prefix, "_end")

  classified %>%
    dplyr::select(group, human,
                  all_of(c(species_col, cat_col,
                           paste0(prefix, c("_chr", "_start", "_end")))))
}

# Full-join baboon and rhesus single-species tables, compute baboon-rhesus
# overlap, save histogram, and return the final 7-column conserved-TSS table.
combine_overlaps <- function(baboon_overlap, rhesus_overlap, hist_path) {
  combined <- full_join(baboon_overlap, rhesus_overlap, by = "group") %>%
    mutate(dist = overlap_len(b_chr, b_start, b_end, r_chr, r_start, r_end))

  plot_histogram(combined$dist, hist_path)

  combined %>%
    mutate(
      br_tss_category = case_when(
        is.na(dist) ~ NA_character_,
        dist > 150  ~ "Conserved TSS",
        TRUE        ~ "non-Conserved TSS"
      ),
      human = dplyr::coalesce(human.x, human.y)
    ) %>%
    dplyr::select(-human.x, -human.y) %>%
    dplyr::select(group, human, baboon, rhesus,
                  hb_tss_category, hr_tss_category, br_tss_category)
}

# End-to-end pipeline for one cell type: parse, overlap, classify, write output
process_cell_type <- function(rate_path, baboon_bed, rhesus_bed,
                              ortho_genes, cell_type, out_dir) {
  human_pause <- get_human_pause_regions(rate_path, ortho_genes)
  baboon_lift <- parse_liftover_bed(baboon_bed, ortho_genes, "baboon")
  rhesus_lift <- parse_liftover_bed(rhesus_bed, ortho_genes, "rhesus")

  hb_overlap <- compute_species_overlap(
    human_pause, baboon_lift,
    species_col = "baboon", cat_col = "hb_tss_category", prefix = "b",
    hist_path   = file.path(out_dir, paste0("human_baboon_", cell_type, "_tssOverlap.pdf"))
  )

  hr_overlap <- compute_species_overlap(
    human_pause, rhesus_lift,
    species_col = "rhesus", cat_col = "hr_tss_category", prefix = "r",
    hist_path   = file.path(out_dir, paste0("human_rhesus_", cell_type, "_tssOverlap.pdf"))
  )

  result <- combine_overlaps(
    hb_overlap, hr_overlap,
    hist_path = file.path(out_dir, paste0("baboon_rhesus_", cell_type, "_tssOverlap.pdf"))
  )

  write.table(
    result,
    file.path(out_dir, paste0("primate_", cell_type, "_conserved_tss.txt")),
    sep = "\t", na = "", quote = FALSE, row.names = FALSE
  )

  invisible(result)
}

# ── Main ──────────────────────────────────────────────────────────────────────

ortho_genes <- read.table(
  file.path(.paths$outputs, "publish/betweenSpecies/primate_conserved_genes.txt"),
  sep = "\t", header = TRUE
)

# CD4
process_cell_type(
  rate_path   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd4_rate.RDS"),
  baboon_bed  = file.path(.paths$data, "liftover/baboon_cd4_to_human.bed"),
  rhesus_bed  = file.path(.paths$data, "liftover/rhesus_cd4_to_human.bed"),
  ortho_genes = ortho_genes,
  cell_type   = "cd4",
  out_dir     = out_dir
)

# CD14
process_cell_type(
  rate_path   = file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human/cd14_rate.RDS"),
  baboon_bed  = file.path(.paths$data, "liftover/baboon_cd14_to_human.bed"),
  rhesus_bed  = file.path(.paths$data, "liftover/rhesus_cd14_to_human.bed"),
  ortho_genes = ortho_genes,
  cell_type   = "cd14",
  out_dir     = out_dir
)
