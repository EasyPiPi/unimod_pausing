---
TITLE: "Unimod Pausing Analysis"
AUTHORS: "Yixin Zhao; Xin Zeng"
DATE: 2026-06-08
---

# Unimod pausing analysis

This repository contains the Snakemake workflow, metadata, and R scripts used to preprocess PRO-seq/PRO-cap data, fit the UniMod models, generate simulated data, perform LRT tests, and reproduce the manuscript figures.

## Repository layout

| Path | Purpose |
| --- | --- |
| `Snakefile` | Entry point that loads metadata/configuration, sets up helper functions, and includes the modular rule files. |
| `rules/` | Snakemake rule collections for PRO-seq preprocessing (`proseq.smk`), unimod pause-release modeling (`unimod.smk`), and simulation post-processing (`simulation.smk`). |
| `scripts/` | R and bash helpers invoked by the rules (PRO-seq wrapper, active TSS calling, statistical models, simulation visualizations). |
| `metadata/` | CSV tables describing experiments (`metadata_aoi.csv`), contrasts (`metadata_comparison.csv`), PRO-cap samples (`copro_sample.csv`), spike-in scaling (`scaling_factor.csv`), and simulation parameter grids. |
| `config.yml` | Paths to references (genomes, annotations) and analysis parameters (window sizes, gene filters, etc.). |
| `environment.yml` | Conda environment mixing Python tooling (Snakemake, numpy/pandas) with the R/Bioconductor stack and command-line genomics utilities (bwa, bedtools, CrossMap, PRO-seq 2.0 dependencies). |
| `analysis/` | Downstream analysis scripts used to reproduce the main and supplementary figures presented in the manuscript. |
| `ext_data/` | Expected location for genomes, chain files, and PRO-cap bigWigs (currently symlinks to external paths). |
| `outputs/`, `tmp/`, `indicator/` | Generated results, intermediates, and completion markers created by the workflow. |
| `ms/` | Manuscript notes such as `pausedRNAP_nucleosome_quotes.md`. |

## Pipeline overview

The `Snakefile` orchestrates the complete analysis by:

1. Importing helper libraries and reading shared metadata tables up front so rules can reuse their columns as wildcards.
2. Defining `rule all`, which enumerates the final deliverables: processed sequencing bigWigs, pause-escape parameter estimates, likelihood-ratio test (LRT) tables, simulation summaries, and visualization files for the manuscript.
3. Including three focused rule modules:
   * **`rules/proseq.smk`** – indexes genomes, links GEO fastqs, runs the PRO-seq 2.0 bash pipeline, extracts 5′ end PRO-seq signals, normalizes coverage, and writes strand-specific bigWigs.
   * **`rules/unimod.smk`** – converts transcripts to counting regions, fits the unimod models per sample, compares two conditions, and generates figure-ready summaries.
   * **`rules/simulation.smk`** – subsamples RDS simulation outputs, summarizes parameter sweeps, and creates LRT tables that match experimental coverage.

## Scripts directory

* `scripts/proseq/`: PRO-seq 2.0 wrapper (`proseq2.0.bsh`) plus helper R scripts for transcript processing and active TSS detection (`find_active_tss_DLD1.R`).
* `scripts/unimod/`: Core statistical routines, including `analyze_one_sample_poisson_pause_release.R`, `analyze_two_samples_DLD1.R`, EM helpers, counting-region generation, and visualization (`visualize_two_samples.R`).
* `scripts/simulation/`: Post-processing utilities for simulation RDS files, coverage matching, LRT summaries, and BW visualizations.

## Metadata and configuration tips

* Update `config.yml` with absolute paths to your reference genomes, chromosome sizes, annotation GTF, and analysis-specific filters before running Snakemake.
* Inspect the CSV schemas in `metadata/` to understand how assays, groups, references, spike-in scaling factors, and simulation grids are encoded. These tables drive wildcard resolution everywhere in the workflow.
* Create the conda environment via `conda env create -f environment.yml` (or `mamba env create ...`) so both Python and R dependencies match the expected versions.

## Getting started

1. **Dry-run the workflow:** after editing `config.yml`, run `snakemake -np` to inspect the dependency graph without executing heavy jobs.
2. **Run the workflow:** `snakemake --cores 4` (adjust cores as needed); targets default to `rule all` defined in `Snakefile`.
3. **Trace a rule end-to-end:** e.g., follow `rules/unimod.smk:analyze_one_sample_pause_release` into `scripts/unimod/analyze_one_sample_poisson_pause_release.R` to see how pause-release rates are estimated.
4. **Review visualization scripts:** use `scripts/unimod/visualize_two_samples.R` to understand how figure panels are composed from processed outputs.
5. **Check manuscript mappings:** `figures.md` documents which scripts generate each figure; keep those references handy when recreating panels.

## Additional notes

* Supplemental figures combine Snakemake outputs with manual layout steps in Inkscape and IGV.
* Temporary outputs under `tmp/` can be safely removed between runs; Snakemake will regenerate them as needed.
* Simulation rules expect pre-existing inputs under `outputs/simulation/data_lrt/` and `outputs/simulation/data_fk/`; acquisition/generation of those inputs is TODO to document.
* `ext_data/copro/hg19/*.bw` and `ext_data/chain/hg19ToHg38.over.chain.gz` are required for PRO-cap liftover; ensure these exist or symlink to their locations.

## Citation
If you use this repository, please cite:
Zeng X, Barshad G, Hassett R, Rice EJ, Danko CG, Siepel A, Zhao Y.
*A comparative analysis of promoter-proximal pausing reveals kinetic and distributional dimensions of variation.* bioRxiv (2026).
Preprint: https://www.biorxiv.org/content/10.64898/2026.06.01.729264v1
