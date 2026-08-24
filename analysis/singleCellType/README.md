# Single Cell Type Analysis Pipeline

## Project Summary

This directory contains the analysis scripts used to generate **Main Figure 3** and **Supplementary Figures S3–S6, S7B**.

- **Processed input data** used in this analysis are stored in `data/publish/` (downloadable from Zenodo).
- **Intermediate files** (rate objects, nucleosome matrices, signal BED files) are stored in `outputs/publish/singleCellType/`.
- **Final figure outputs** (PDF files) are stored in `outputs/publish/singleCellType/`.

The pipeline characterizes transcriptional dynamics in human, rhesus macaque, and baboon immune cells (CD4+ T cells and CD14+ monocytes) by integrating PRO-seq transcription rates, Micro-C nucleosome positioning signals, histone modification data, and promoter sequence elements.

---

## Step 1: Transcription Rate Estimation

**Biological/computational purpose:** Estimate per-gene transcriptional kinetic parameters — pause release rate (β), elongation rate (χ), and TSS fuzziness (μ/σ) — from PRO-seq data for each species and cell type using the STADyUM framework. Genes are grouped into quantile bins and classified into Sharp/Broad TSS categories.

---

### `1_1_EstimateRates_all_species.R`

**Purpose:**
Runs the STADyUM rate estimation pipeline for a user-specified species (human, rhesus, or baboon). For each species, it processes two cell types (CD4+ T cells and CD14+ monocytes). The script selects one upstream TSS per gene, constructs pause and gene body read-count windows, runs `estimateTranscriptionRates()`, annotates genes with quantile groups (betaGroup, chiGroup, fkMeanGroup) and a Sharp/Broad TSS classification based on a valley-threshold in the fkSD distribution, and saves the resulting STADyUM rate objects.

**Required input data:**
- PRO-seq TSN RDS files per species and cell type (`PROseq-{SPECIES}-{CELLTYPE}_tsn.RDS`)
- STADyUM transcript annotation RDS files (`{SPECIES}.RDS`)
- PRO-seq BigWig files (plus and minus strand) per sample

**Generated outputs:**
- `{species}/{celltype}_rate.RDS` — STADyUM rate object with annotated kinetic parameters (one per sample, six total across human/rhesus/baboon × CD4/CD14)
- `{species}/cd4_cd14_tss_dist.txt` — TSS distance QC table comparing TSS positions between cell types

**Related figures:**
- None (rate estimation only; figures are produced by the downstream visualization script)

---

### `1_2_EstimateRates_all_species_viz.R`

**Purpose:**
Reads the STADyUM rate RDS objects produced by `1_1_EstimateRates_all_species.R` and generates scatter and density plots of key kinetic parameter pairs for all six species × cell type combinations. Human CD4+ is shown as the main figure panel; the remaining five samples are arranged into supplementary multi-panel figures.

**Required input data:**
- `{species}/{celltype}_rate.RDS` — STADyUM rate objects from `1_1_EstimateRates_all_species.R`

**Generated outputs:**
- `human/cd4_beta_chi_scatter_new.pdf` — scatter plot of log₁₀(β) vs. log₁₀(χ) for human CD4+
- `human/beta_fk_new.pdf` — 2D density plot of log₁₀(β) vs. σ with Sharp/Broad boundary for human CD4+
- `beta_chi_scatter_all.pdf` — supplementary multi-panel β vs. χ scatter for all other samples
- `mu_sigma_scatter_all.pdf` — supplementary multi-panel μ vs. σ scatter for all six samples
- `beta_sigma_2d_all.pdf` — supplementary multi-panel β vs. σ density for all other samples

**Related figures:**
- Main Figure 3A, 3D (β vs. χ scatter; β vs. σ density — Human CD4+)
- Supplementary Figure S3A, S5A–B (β vs. χ, μ vs. σ, β vs. σ — all species/cell types)

---

## Step 2: Nucleosome Positioning and Micro-C Signal

**Biological/computational purpose:** Extract and visualize chromatin architecture signals (1D Micro-C) around transcription start sites to characterize nucleosome positioning as a function of transcriptional kinetic parameters and core promoter sequence elements.

---

### `2_1_1_generateTSSregion.R`

**Purpose:**
Extracts TSS genomic coordinates from STADyUM pause regions for all six samples (human, rhesus, baboon × CD4, CD14). The TSS is defined as the 5′ end of each pause region (single-base resolution). Pause and gene body read counts are appended, and the result is written as a sorted BED file for use as input to the 1D signal extraction scripts.

**Required input data:**
- `{species}/{celltype}_rate.RDS` — STADyUM rate objects from `1_1_EstimateRates_all_species.R`

**Generated outputs:**
- `{species}/{celltype}_tss_sort.bed` — sorted, single-base TSS BED file with pause and gene body counts (six files total)

**Related figures:**
- None (intermediate data preparation)

---

### `2_1_2_extract_1d_signal.py`

**Purpose:**
Python script that converts Micro-C paired-end contact data into a 1D signal at TSS regions. Reads a `.pairs.gz` file, retains reads with MAPQ ≥ 30, collapses both ends into a 1D table, shifts each end 75 bp downstream according to its strand, then intersects the resulting positions with ±1000 bp windows around each TSS. The output records per-position read counts for each gene's TSS window.

**Required input data:**
- A `.pairs.gz` Micro-C contact file for a single species and cell type
- A TSS BED file (from `2_1_1_generateTSSregion.R`)

**Generated outputs:**
- A BED file of strand-shifted 1D Micro-C signal counts intersected with TSS windows (one output per species × cell type)

**Related figures:**
- None (intermediate signal extraction)

---

### `2_1_3_batch_extract_1d_signal.sh`

**Purpose:**
Bash wrapper that calls `2_1_2_extract_1d_signal.py` for all six species × cell type combinations (human, rhesus, baboon × CD4, CD14), skipping any sample for which input files are missing.

**Required input data:**
- `{species}_{celltype}_merged.pairs.gz` — merged Micro-C pairs files per sample
- `{species}/{celltype}_tss_sort.bed` — TSS BED files from `2_1_1_generateTSSregion.R`
- `2_1_2_extract_1d_signal.py` — must be in the working directory

**Generated outputs:**
- `{species}_{celltype}_1d_signal.bed` — 1D Micro-C signal BED file for each sample (six files total)

**Related figures:**
- None (batch runner)

---

### `2_2_1_nucleosome_position.R`

**Purpose:**
Constructs nucleosome score matrices from the 1D Micro-C signal for each of the six samples. For each gene, a 2000 bp window centered on the TSS is extracted, and scores are aggregated into a matrix using `genomation::ScoreMatrix`. Each row is LOESS-smoothed. The mean signal over the +1 nucleosome region (positions 1001–1200 relative to the TSS) is added to the rate data frame as a summary statistic. Results are saved as RDS objects containing both raw and smoothed matrices alongside the annotated rate data frame.

**Required input data:**
- `{species}/{celltype}_rate.RDS` — STADyUM rate objects from `1_1_EstimateRates_all_species.R`
- `{species}_{celltype}_1d_signal.bed` — 1D Micro-C signal BED files from `2_1_3_batch_extract_1d_signal.sh`
- `TxDb.Hsapiens.UCSC.hg38.knownGene` R package (for human seqinfo)
- FASTA index files `rheMac10.fa` and `papAnu4.fa` (for rhesus and baboon seqinfo)

**Generated outputs:**
- `{species}/{celltype}_ns_matrix.RDS` — list containing raw score matrix, LOESS-smoothed matrix, TSS GRanges, and annotated rate data frame (six files total)

**Related figures:**
- None (intermediate matrix computation)

---

### `2_2_2_nucleosome_position_vis.R`

**Purpose:**
Reads nucleosome score matrices from `2_2_1_nucleosome_position.R` and generates line-profile plots of mean Micro-C signal over the 0–1000 bp downstream of the TSS, faceted by expression group (χ group) and colored by β or σ group. Human CD4+ produces the main figure panels; the remaining samples are assembled into supplementary multi-panel figures.

**Required input data:**
- `{species}/{celltype}_ns_matrix.RDS` — nucleosome matrix RDS files from `2_2_1_nucleosome_position.R`

**Generated outputs:**
- `human/cd4_Micro-C_beta_gb_chi.pdf` — nucleosome profile by β group, faceted by expression level (main figure)
- `human/cd4_Micro-C_fkSD_gb_chi.pdf` — nucleosome profile by σ group, faceted by expression level (main figure)
- `beta_chi_Micro_all.pdf` — supplementary 5-panel layout: β groups by expression, all samples except Human CD4+
- `sigma_chi_Micro_all.pdf` — supplementary 5-panel layout: σ groups by expression, all samples except Human CD4+
- `sigma_beta_Micro_all.pdf` — supplementary 6-panel layout: σ groups by β group, all six samples

**Related figures:**
- Main Figure 3B, 3E (nucleosome profiles by β and σ group — Human CD4+)
- Supplementary Figure S3B, S6A–C (nucleosome profiles — all other species/cell types; σ by β group — all samples)

---

### `2_2_3_nucleosome_position_promoter.R`

**Purpose:**
Links TSS pause sites to annotated core promoter elements (TATA-box, Initiator, CCAAT-box, GC-box from the EPD database) and examines how promoter motif class relates to nucleosome positioning. For each gene, the nearest annotated TSS within 200 bp on the same strand is identified, and a motif class is assigned. The script generates UpSet plots of motif co-occurrence patterns and nucleosome line-profile plots stratified by motif class and β group, for human CD4+ and CD14+ cells.

**Required input data:**
- `human/cd4_rate.RDS` and `human/cd14_rate.RDS` — STADyUM rate objects from `1_1_EstimateRates_all_species.R`
- `human/cd4_ns_matrix.RDS` and `human/cd14_ns_matrix.RDS` — nucleosome matrix RDS files from `2_2_1_nucleosome_position.R`
- `data/promoter_anno/epd/merge_promoter_motif.txt` — EPD promoter motif annotation table

**Generated outputs:**
- `human/cd4_motif_info.pdf` and `human/cd14_motif_info.pdf` — UpSet plots of promoter motif combinations
- `human/cd4_motif_nuc_merged.pdf` and `human/cd14_motif_nuc_merged.pdf` — nucleosome profiles colored by motif class
- `human/cd4_Micro-C_beta_gb_motif.pdf` and `human/cd14_Micro-C_beta_gb_motif.pdf` — nucleosome profiles faceted by motif class, colored by β group

**Related figures:**
- Supplementary Figure S4 (promoter motif UpSet plots and motif-stratified nucleosome profiles)

---

## Step 3: Histone Modification Analysis

**Biological/computational purpose:** Quantify promoter-proximal H3K4me3 and H3K27ac signals from ENCODE ChIP-seq data and test whether these active chromatin marks differ systematically across transcriptional kinetic groups (β quantile bins), stratified by expression level (χ group).

---

### `3_histone_modification.R`

**Purpose:**
Defines 500 bp promoter windows centered on each TSS using STADyUM pause regions. For each histone mark (H3K4me3 and H3K27ac), reads the corresponding BigWig fold-change files from ENCODE, computes the mean signal within each promoter window, and generates boxplots of log₂(signal) by β group, faceted by expression (χ) group. Statistical comparisons are annotated using Kruskal–Wallis tests. Processes human CD4+ and CD14+ cells independently.

**Required input data:**
- `human/cd4_rate.RDS` and `human/cd14_rate.RDS` — STADyUM rate objects from `1_1_EstimateRates_all_species.R`
- `data/result2/encode/meta_fc.RDS` — ENCODE metadata table linking biosample, histone mark, and BigWig file paths

**Generated outputs:**
- `cd4_H3K4me3_beta_chi.pdf` — boxplot of H3K4me3 by β group, faceted by expression level (CD4+)
- `cd4_H3K27ac_beta_chi.pdf` — boxplot of H3K27ac by β group, faceted by expression level (CD4+)
- `cd14_H3K4me3_beta_chi.pdf` — boxplot of H3K4me3 by β group, faceted by expression level (CD14+)
- `cd14_H3K27ac_beta_chi.pdf` — boxplot of H3K27ac by β group, faceted by expression level (CD14+)
- `cd4_histone.rds` — RDS object containing promoter data frame with histone signal for CD4+

**Related figures:**
- Main Figure 3C, Supplementary Figure S4 (histone modification signals by β and expression group)

---

## Step 4: Gene Ontology Enrichment Analysis

**Biological/computational purpose:** Identify biological processes enriched in genes from different β-group quantiles (Q1–Q5) among highly expressed genes, to characterize the functional consequences of transcriptional pause-release kinetics.

---

### `4_GO_human.R`

**Purpose:**
Filters human CD4+ and CD14+ genes to those with high expression (chiGroup == "High"), splits them by β group (Q1–Q5), and runs Gene Ontology Biological Process enrichment analysis using `clusterProfiler::compareCluster` with all detected TSN genes as the background universe. Redundant GO terms are removed by semantic similarity simplification. Results are visualized as dot plots showing the top 5 enriched terms per β group.

**Required input data:**
- `human/cd4_rate.RDS` and `human/cd14_rate.RDS` — STADyUM rate objects from `1_1_EstimateRates_all_species.R`
- `PROseq-HUMAN-CD4_tsn.RDS` and `PROseq-HUMAN-CD14_tsn.RDS` — TSN RDS files (used to define the gene universe)

**Generated outputs:**
- `human/cd4_beta_GO_chiHigh.pdf` — GO enrichment dot plot by β group for highly expressed CD4+ genes
- `human/cd14_beta_GO_chiHigh.pdf` — GO enrichment dot plot by β group for highly expressed CD14+ genes

**Related figures:**
- Supplementary Figure S7B (GO enrichment by β group in CD4+ and CD14+ cells)

---

## Pipeline Overview

The scripts form a linear dependency chain with a branching visualization step at each stage:

```
1_1  Estimate transcription rates (β, χ, σ) from PRO-seq data
      → 1_2  Visualize rate parameter distributions (main + supplementary figures)
      ↓
2_1_1  Extract TSS positions as BED files
      ↓
2_1_2 / 2_1_3  Extract 1D Micro-C nucleosome signal around TSS
      ↓
2_2_1  Compute per-gene nucleosome score matrices (LOESS-smoothed)
      → 2_2_2  Visualize nucleosome profiles by β / σ / χ group (main + supplementary)
      → 2_2_3  Visualize nucleosome profiles by core promoter motif class
      ↓
3_histone_modification.R  Quantify H3K4me3 / H3K27ac at promoters by β group
      ↓
4_GO_human.R              GO enrichment by β group in highly expressed genes
```

| Step | Key output | Purpose in paper |
|------|-----------|-----------------|
| 1 | `*_rate.RDS` | Kinetic parameter estimation and classification |
| 2 | `*_ns_matrix.RDS`, PDFs | Nucleosome positioning relative to TSS kinetics |
| 3 | Histone signal PDFs | Active chromatin marks vs. kinetic groups |
| 4 | GO dot plot PDFs | Functional annotation of pause-release kinetics |
