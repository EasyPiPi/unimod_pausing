library(org.Hs.eg.db)
library(clusterProfiler)
library(STADyUM)
library(cowplot)
library(ggplot2)
library(stringr)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", path.expand("~/Desktop/project/YiXin_Likelihood")),
  mustWork = FALSE)
source(file.path(root_dir, "codes", "publish", "load_config.R"))

result_dir <- file.path(.paths$outputs, "publish/singleCellType/5_GO/")

# --- Configuration -----------------------------------------------------------

rate_base <- file.path(.paths$outputs, "publish/singleCellType/1_transcriptionRate/human")
tsn_base  <- file.path(.paths$data, "proseq/tidgrng")

cell_configs <- list(
  cd4 = list(
    name       = "CD4",
    rate_rds   = file.path(rate_base, "cd4_rate.RDS"),
    tsn_rds    = file.path(tsn_base,  "PROseq-HUMAN-CD4_tsn.RDS"),
    output_pdf = file.path(result_dir, "human/cd4_beta_GO_chiHigh.pdf"),
    width      = 8,
    height     = 9
  ),
  cd14 = list(
    name       = "CD14",
    rate_rds   = file.path(rate_base, "cd14_rate.RDS"),
    tsn_rds    = file.path(tsn_base,  "PROseq-HUMAN-CD14_tsn.RDS"),
    output_pdf = file.path(result_dir, "human/cd14_beta_GO_chiHigh.pdf"),
    width      = 8,
    height     = 9
  )
)

# --- Helper functions --------------------------------------------------------

run_beta_group_go <- function(rate_rds, tsn_rds) {
  for (f in c(rate_rds, tsn_rds)) {
    if (!file.exists(f)) stop("Input file not found: ", f)
  }

  rate_obj <- readRDS(rate_rds)
  tsn      <- readRDS(tsn_rds)

  rate_df <- rates(rate_obj)
  required_cols <- c("chiGroup", "betaGroup", "geneId")
  missing <- setdiff(required_cols, colnames(rate_df))
  if (length(missing) > 0)
    stop("Missing columns in rate dataframe: ", paste(missing, collapse = ", "))
  if (!"ensembl_gene_id" %in% colnames(mcols(tsn)))
    stop("Column 'ensembl_gene_id' not found in TSN file")

  rate_df <- rate_df[rate_df$chiGroup == "High", ]
  gene_groups <- split(rate_df$geneId, rate_df$betaGroup)

  go_result <- compareCluster(
    gene_groups,
    fun            = enrichGO,
    OrgDb          = org.Hs.eg.db,
    keyType        = "ENSEMBL",
    ont            = "BP",
    pAdjustMethod  = "BH",
    qvalueCutoff   = 0.05,
    universe       = tsn$ensembl_gene_id
  )

  simplify(go_result, cutoff = 0.9, by = "p.adjust", select_fun = min)
}

plot_go_dotplot <- function(go_simple) {
  p <- dotplot(go_simple, showCategory = 5) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    theme_cowplot() +
    labs(x = NULL)
  p@theme$axis.text.y <- element_text(size = 13, hjust = 1)
  p@theme$axis.text.x <- element_text(size = 16)
  p
}

process_one_cell_type <- function(cfg) {
  message("Processing ", cfg$name, " ...")
  dir.create(dirname(cfg$output_pdf), recursive = TRUE, showWarnings = FALSE)

  go_simple <- run_beta_group_go(cfg$rate_rds, cfg$tsn_rds)
  p         <- plot_go_dotplot(go_simple)

  ggsave(cfg$output_pdf, p, dpi = 300, width = cfg$width, height = cfg$height)
  message("Saved: ", cfg$output_pdf)
}

# --- Run ---------------------------------------------------------------------

for (cfg in cell_configs) {
  process_one_cell_type(cfg)
}
