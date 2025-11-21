# Unimod pausing analysis

This repository hosts the Snakemake workflow, metadata, and R helpers used to preprocess PRO-seq/PRO-cap data, fit pause-release models, explore synthetic datasets, and reproduce manuscript figures.

## Repository layout

| Path | Purpose |
| --- | --- |
| `Snakefile` | Entry point that loads configuration/metadata, defines `rule all`, and includes the rule modules. |
| `rules/` | Snakemake rule collections for PRO-seq preprocessing (`proseq.smk`), unimodal pause-release modeling (`unimod.smk`), and simulation post-processing (`simulation.smk`). |
| `scripts/` | R and bash helpers invoked by the rules: pipeline wrappers, statistical models, visualization scripts, and simulation utilities. |
| `metadata/` | CSV tables describing experiments, spike-ins, group comparisons, and simulation parameter grids that populate Snakemake wildcards. |
| `config.yml` | Paths to references (genome files, annotations) and analysis parameters (windows, filters). |
| `environment.yml` | Conda environment combining Python tooling (Snakemake, numpy/pandas) with the R/Bioconductor stack and command-line genomics utilities. |
| `ms/` | Manuscript notes such as `pausedRNAP_nucleosome_quotes.md`. |
| `tmp/` | Temporary outputs written by intermediate rules. |

## Pipeline overview

The workflow is orchestrated by the root `Snakefile`:

1. **Metadata loading.** Shared tables (e.g., `metadata/metadata_aoi.csv`, `metadata/metadata_comparison.csv`, simulation parameter grids) are read up front so their columns can be reused as wildcards throughout the rules.
2. **Targets.** `rule all` enumerates end-to-end outputs: processed sequencing bigWigs, within-sample pause-release rate estimates, between-sample LRTs plus visualizations, and simulation summaries/LRT tables matched to experimental coverage.
3. **Modules.**
   * **`rules/proseq.smk`** – indexes genomes, links GEO fastqs, runs the PRO-seq 2.0 bash pipeline, extracts 5′ ends, normalizes coverage, and writes strand-specific bigWigs.
   * **`rules/unimod.smk`** – converts transcripts to counting regions, fits pause-release models per sample, performs two-condition comparisons, and triggers visualization-ready summaries for figure panels.
   * **`rules/simulation.smk`** – subsamples RDS simulation outputs, summarizes pause-release and steric hindrance sweeps, and creates coverage-matched LRT tables using the helper scripts in `scripts/simulation/` (e.g., `lrt_matched_coverage.R`, `aux/lrt_power.R`).

## Scripts directory

* `scripts/proseq/`: PRO-seq 2.0 wrapper plus helper R scripts for active TSS detection (`find_active_tss_DLD1.R`) and transcript processing.
* `scripts/unimod/`: Core statistical routines, including `analyze_one_sample_poisson_pause_release.R`, `analyze_two_samples_DLD1.R`, expectation-maximization helpers, and `visualize_two_samples.R` for figure panels.
* `scripts/simulation/`: Post-processing utilities for simulation RDS files, coverage matching (`lrt_matched_coverage.R`), and power analyses under different pause/escape scenarios (`aux/lrt_power.R`).

## Metadata and configuration tips

* Update `config.yml` with absolute paths to reference genomes, chromosome sizes, annotations, and analysis-specific filters before running Snakemake.
* Inspect the CSV schemas in `metadata/` to understand how assays, groups, references, replicates, and spike-in scaling factors are encoded; these tables drive wildcard resolution everywhere in the workflow.
* Create the conda environment via `conda env create -f environment.yml` (or `mamba env create ...`) so both Python and R dependencies match the expected versions.

## Getting started

1. **Dry-run the workflow:** after editing `config.yml`, run `snakemake -np` to inspect the dependency graph without executing heavy jobs.
2. **Trace a rule end-to-end:** e.g., follow `rules/unimod.smk:analyze_one_sample_pause_release` into its corresponding R script to see how pause-release rates are estimated.
3. **Review visualization scripts:** use `scripts/unimod/visualize_two_samples.R` to understand how figure panels (e.g., Figure 1G–I) are composed from the processed outputs.
4. **Check manuscript mappings:** the original manuscript references which scripts generate each figure; keep those pointers handy when recreating panels in Inkscape or IGV.

## Additional notes

* Supplemental figures (e.g., Figure S1 panels) combine Snakemake outputs with manual layout steps in Inkscape and IGV.
* Temporary outputs under `tmp/` can be safely removed between runs; Snakemake will regenerate them as needed.
