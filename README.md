# Comparative Analysis of Promoter-Proximal Pausing

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-STADyUM-blue.svg)](https://bioconductor.org/packages/release/bioc/html/STADyUM.html)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20598895-blue.svg)](https://doi.org/10.5281/zenodo.20598895)

This repository contains the analysis pipelines, Snakemake workflows, and figure-generation scripts for the manuscript:

> **A comparative analysis of promoter-proximal pausing reveals kinetic and distributional dimensions of variation**  
> Xin Zeng, Gilad Barshad, Rebecca Hassett, Edward J. Rice, Charles G. Danko, Adam Siepel\*, Yixin Zhao\*

---

## Overview

Promoter-proximal pausing of RNA polymerase II (Pol II) has separable kinetic ($\beta$, pause-escape rate) and distributional ($\sigma$, pause-site dispersion) dimensions that vary differently across biological contexts. This repository provides:

1. **Preprocessing & Snakemake Pipeline**: Automated preprocessing of PRO-seq/PRO-cap data, active TSS identification, parameter estimation, SimPol synthetic benchmarks, and perturbation LRT tests (Figs. 1, 2, S1, S2).
2. **Downstream Analysis Workflows (`analysis/`)**: Scripts for single-cell-type kinetics/dispersion, cross-cell-type comparisons, cross-species evolutionary analyses, +1 nucleosome integration (Micro-C), histone modification ChIP-seq integration, and promoter sequence GC content / skew analyses (Figs. 3–7, S3–S11).
3. **Statistical Package Implementation**: The underlying likelihood-ratio testing and EM algorithms are formally packaged in the Bioconductor R package **[STADyUM](https://bioconductor.org/packages/release/bioc/html/STADyUM.html)**.

---

## Repository Structure

| Directory / File | Description |
| :--- | :--- |
| `Snakefile` | Top-level Snakemake entry point for preprocessing, simulation sweeps, and perturbation models. |
| `rules/` | Snakemake rule sets: `proseq.smk` (read processing), `unimod.smk` (model fitting), `simulation.smk` (synthetic data sweeps). |
| `scripts/` | Helper scripts called by Snakemake rules for PRO-seq 2.0 processing, EM optimization, and simulation analysis. |
| `analysis/` | End-to-end analysis and figure-generation scripts for within-cell-type, between-cell-type, and between-species comparisons (Figs. 3–7). |
| `metadata/` | Sample sheets, contrast definitions, spike-in scaling factors, and simulation parameters. |
| `config.yml` | Configuration file for paths, chromosome info, genome indices, and analysis cutoffs. |
| `environment.yml` | Conda environment specification (Python, R, Bioconductor, and bioinformatics tools). |
| `figures.md` | Comprehensive mapping from each manuscript figure panel to its generating code. |

---

## Quick Start

### 1. Environment Setup

Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate unimod
```

### 2. Configuration

1. Review and adjust `config.yml` with your local reference paths if needed (reference genome FASTA, BWA index, chromosome sizes, and gene annotations; defaults use standard relative paths).
2. For downstream figure scripts in `analysis/`, set `PROJECT_ROOT` in `analysis/paths.yaml` or export the environment variable:
   ```bash
   export PROJECT_ROOT="/path/to/unimod_pausing"
   ```

### 3. Running the Snakemake Workflow

* **Dry run:**
  ```bash
  snakemake -np
  ```
* **Execute workflow:**
  ```bash
  snakemake --cores 8
  ```

---

## Figure Reproduction

Detailed instructions and file-level mappings for all panels (Main Figures 1–7 and Supplementary Figures S1–S11) are cataloged in [`figures.md`](figures.md).

* **Simulations & LRT Benchmark (Fig. 1, Supp. Fig. S1)**: Run via `rules/simulation.smk` or standalone scripts in `scripts/simulation/`.
* **Perturbation Experiments (Fig. 2, Supp. Fig. S2)**: Run via `rules/unimod.smk` or standalone scripts in `analysis/perturbation/`.
* **Single Cell-Type & Chromatin Features (Fig. 3, Supp. Figs. S3–S6)**: Run scripts in `analysis/singleCellType/`.
* **Between Cell-Type Comparisons (Fig. 4, Supp. Figs. S7–S8)**: Run scripts in `analysis/betweenCellType/`.
* **Between Species & Evolutionary Divergence (Fig. 5, Supp. Fig. S9)**: Run scripts in `analysis/betweenSpecies/`.
* **Chromatin & Sequence Determinants (Figs. 6, 7, Supp. Figs. S10–S11)**: Run integrative scripts in `analysis/betweenSpecies/` and `analysis/betweenCellType/`.

---

## Data & Software Availability

* **STADyUM Package**: [Bioconductor release](https://bioconductor.org/packages/release/bioc/html/STADyUM.html)
* **Processed Data & Vignettes**: [Zenodo (DOI: 10.5281/zenodo.20598895)](https://doi.org/10.5281/zenodo.20598895)
* **Raw Sequencing Data**:
  * Newly generated primate PRO-seq and Micro-C: NCBI GEO (GSE342014) and NIH dbGaP (phs002146).
  * Perturbation PRO-seq & PRO-cap (DLD-1 cells): NCBI GEO (GSE144786, Aoi et al., *Mol. Cell* 2020).
  * Histone modifications: ENCODE Consortium (see Supplementary Tables S1–S2).

---

## Citation

If you use this codebase or the STADyUM package in your research, please cite:

```bibtex
@article{zeng2026comparative,
  title={A comparative analysis of promoter-proximal pausing reveals kinetic and distributional dimensions of variation},
  author={Zeng, Xin and Barshad, Gilad and Hassett, Rebecca and Rice, Edward J and Danko, Charles G and Siepel, Adam and Zhao, Yixin},
  year={2026}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
