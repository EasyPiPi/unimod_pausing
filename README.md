---
TITLE: "Unimod Pausing Analysis"
AUTHORS: "Yixin Zhao; Xin Zeng"
DATE: 2025-11-17
---

# Unimod pausing analysis

This repository contains the Snakemake workflow, metadata, and R scripts used to preprocess PRO-seq/PRO-cap data, fit the UniMod models, generate simulated data, perform LRT tests, and reproduce the manuscript figures.

## Repository layout

| Path | Purpose |
| --- | --- |
| `Snakefile` | Entry point that loads metadata/configuration, sets up helper functions, and includes the modular rule files. |
| `rules/` | Snakemake rule collections for PRO-seq preprocessing (`proseq.smk`), unimod pause-release modeling (`unimod.smk`), and simulation post-processing (`simulation.smk`). |
| `scripts/` | R and bash helpers invoked by the rules (pipeline wrappers, statistical models, visualization scripts). |
| `metadata/` | CSV tables describing experiments, spike-ins, parameter grids, and simulation settings that populate Snakemake wildcards. |
| `config.yml` | Paths to references (genomes, annotations) and analysis parameters (window sizes, gene filters, etc.). |
| `environment.yml` | Conda environment mixing Python tooling (Snakemake, numpy/pandas) with the R/Bioconductor stack and command-line genomics utilities. |
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

* `scripts/proseq/`: PRO-seq 2.0 wrapper plus helper R scripts for active TSS detection (`find_active_tss_DLD1.R`) and transcript processing.
* `scripts/unimod/`: Core statistical routines, including `analyze_one_sample_poisson_pause_release.R`, `analyze_two_samples_DLD1.R`, expectation-maximization helpers, and `visualize_two_samples.R` for figure panels.
* `scripts/simulation/`: Post-processing utilities for simulation RDS files, coverage matching, and generation of supplemental figure tables.

## Metadata and configuration tips

* Update `config.yml` with absolute paths to your reference genomes, chromosome sizes, annotation GTF, and analysis-specific filters before running Snakemake.
* Inspect the CSV schemas in `metadata/` to understand how assays, groups, references, and spike-in scaling factors are encoded. These tables drive wildcard resolution everywhere in the workflow.
* Create the conda environment via `conda env create -f environment.yml` (or `mamba env create ...`) so both Python and R dependencies match the expected versions.

## Getting started

1. **Dry-run the workflow:** after editing `config.yml`, run `snakemake -np` to inspect the dependency graph without executing heavy jobs.
2. **Trace a rule end-to-end:** e.g., follow `rules/unimod.smk:analyze_one_sample_pause_release` into its corresponding R script to see how pause-release rates are estimated.
3. **Review visualization scripts:** use `scripts/unimod/visualize_two_samples.R` to understand how figure panels (e.g., Figure XXXX) are composed from the processed outputs.
4. **Check manuscript mappings:** figures.md documented which scripts generate each figure; keep those references handy when recreating panels in Inkscape or IGV.

## Additional notes

* Supplemental figures combine Snakemake outputs with manual layout steps in Inkscape and IGV.
* Temporary outputs under `tmp/` can be safely removed between runs; Snakemake will regenerate them as needed.
