# Publish Code

Reproducibility code for the YiXin Likelihood manuscript. Each sub-folder corresponds to one analysis module; scripts are numbered in execution order.

## Setup

1. Edit `paths.yaml` to set `PROJECT_ROOT` to where you placed the project.
2. Run `python check_paths.py` to verify required data files are present.

## Folder overview

| Folder | Contents |
|--------|----------|
| `singleCellType/` | Rate estimation (β, χ, σ) from PRO-seq; nucleosome positioning; histone modification; GO enrichment — all for a single cell type within one species |
| `betweenCellType/` | LRT comparing transcriptional parameters between CD4 and CD14 cells; nucleosome and histone comparisons; enhancer-promoter contact analysis |
| `betweenSpecies/` | Cross-species LRT (human vs. rhesus, human vs. baboon); nucleosome occupancy comparisons; GC content analysis at TSS |
| `perturbation/` | EM-based two-condition model for perturbation experiments (DLD1 cells) |

## Path helpers (at this level)

| File | Purpose |
|------|---------|
| `paths.yaml` | Set `PROJECT_ROOT` here; documents optional path overrides |
| `load_config.R` | Sourced by every R script; defines the `.paths` list |
| `paths.sh` | Sourced by every shell script; exports `DATA`, `OUTPUTS`, and related variables |
| `check_paths.py` | Checks that all required input files are present |
