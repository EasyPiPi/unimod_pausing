library(tidyverse)

root_dir <- normalizePath(
  Sys.getenv("PROJECT_ROOT", "."),
  mustWork = FALSE)
source(file.path(root_dir, "analysis", "load_config.R"))

out_dir <- file.path(.paths$outputs, "publish/betweenSpecies")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read an OrthoFinder pairwise ortholog table, keep 1-to-1 orthologs,
# strip version suffixes from all gene ID columns, and drop Orthogroup.
read_one2one <- function(path, cols) {
  df <- read.table(path, sep = "\t", header = TRUE)
  # Discard rows where any species column contains a comma (multi-gene rows)
  is_one2one <- !apply(df[, cols], 1, function(x) any(grepl(",", x)))
  df <- df[is_one2one, ]
  # Remove version suffixes (e.g. ENSG00000142065.14 -> ENSG00000142065)
  df[, cols] <- lapply(df[, cols], function(x) sub("\\..*$", "", x))
  df[, cols]   # drop Orthogroup by returning only the species columns
}

baboon_pairs <- read_one2one(
  file.path(.paths$data, "orthology/human__v__baboon.tsv"),
  cols = c("human", "baboon")
)

rhesus_pairs <- read_one2one(
  file.path(.paths$data, "orthology/human__v__rhesus.tsv"),
  cols = c("human", "rhesus")
)

# Full join on human gene ID so genes present in only one comparison are kept
primate_conserved_genes <- full_join(rhesus_pairs, baboon_pairs, by = "human")
primate_conserved_genes$group <- paste0("Group", seq_len(nrow(primate_conserved_genes)))

write.table(
  primate_conserved_genes,
  file.path(out_dir, "primate_conserved_genes.txt"),
  sep       = "\t",
  na        = "",
  quote     = FALSE,
  row.names = FALSE
)
