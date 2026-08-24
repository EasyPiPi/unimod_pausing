# Figure-to-Script Mapping (Manuscript 2026)

This document maps all figures and supplementary figures in the manuscript to their corresponding code in this repository and the companion packages.

---

## Main Figures

### Figure 1. Likelihood-ratio testing framework for quantifying pause-escape kinetics and pausing distribution
* **Panel A**: Methodological framework schematic (Adobe Illustrator / Inkscape).
* **Panel B**: Simulated PRO-seq profiles across pause-escape rates ($\beta$) $\rightarrow$ `scripts/simulation/visualize_pause_escape_bw.R`
* **Panel C**: Statistical power comparison between LRT and chi-squared test $\rightarrow$ `scripts/simulation/pause_escape/lrt_rates.R`
* **Panel D**: Representative simulated pause-site distributions ($\mu$, $\sigma$) $\rightarrow$ `scripts/simulation/visualize_pause_distribution_bw.R`
* **Panel E**: Statistical power surface of the distributional LRT $\rightarrow$ `scripts/simulation/pause_distribution/lrt_fk.R`

### Figure 2. Perturbation experiments reveal distinct regulatory effects on pause-escape kinetics and pausing distribution
* **Panels A, C**: Genome-wide LRT and $\Delta\beta$ vs $\Delta\sigma$ under P-TEFb inhibition (NVP-2) $\rightarrow$ `scripts/unimod/analyze_two_samples_DLD1.R` (or `analysis/perturbation/analyze_two_samples_DLD1.R`)
* **Panel B**: PRO-seq metaplots downstream of TSS under NVP-2 treatment $\rightarrow$ `scripts/unimod/visualize_two_samples.R`
* **Panels D, F**: Genome-wide LRT and $\Delta\beta$ vs $\Delta\sigma$ under NELF acute depletion $\rightarrow$ `scripts/unimod/analyze_two_samples_DLD1.R`
* **Panel E**: PRO-seq metaplots downstream of TSS under NELF depletion $\rightarrow$ `scripts/unimod/visualize_two_samples.R`

### Figure 3. Promoter-proximal pausing can be described using separable kinetic and distributional aspects
* **Panels A, D**: $\beta$ vs transcriptional activity $\chi$, and $\beta$ vs $\sigma$ orthogonality in human $\text{CD4}^+$ T cells $\rightarrow$ `analysis/singleCellType/1_1_EstimateRates_all_species.R` & `1_2_EstimateRates_all_species_viz.R`
* **Panels B, E**: +1 nucleosome occupancy stratified by $\beta$ quintiles and $\sigma$ categories $\rightarrow$ `analysis/singleCellType/2_2_1_nucleosome_position.R` & `2_2_2_nucleosome_position_vis.R`
* **Panel C**: Promoter $\text{H3K27ac}$ signal distributions across $\beta$ quintiles $\rightarrow$ `analysis/singleCellType/3_histone_modification.R`

### Figure 4. Cell-type differences in pause-escape rate and pausing distribution
* **Panel A**: Comparison of $\beta$ between human $\text{CD4}^+$ T cells and $\text{CD14}^+$ monocytes $\rightarrow$ `analysis/betweenCellType/1_1_LRT_unified.R` & `1_2_LRT_viz.R`
* **Panel B**: Comparison of pausing dispersion ($\sigma$) transitions between cell types $\rightarrow$ `analysis/betweenCellType/1_2_LRT_viz.R`
* **Panel C**: Joint distribution of changes $\Delta\beta$ vs $\Delta\sigma$ across cell types $\rightarrow$ `analysis/betweenCellType/1_2_LRT_viz.R`

### Figure 5. Cross-species changes in pause-escape rate and pausing distribution
* **Panel A**: Conserved vs non-conserved TSS classification across primate species $\rightarrow$ `analysis/betweenSpecies/1_1_getOrthoGenes.R`, `1_2_getLiftover.R`, `1_3_getConservedTSS.R`
* **Panel B**: Cross-species $\beta$ comparison $\rightarrow$ `analysis/betweenSpecies/1_4_LRT.R` & `1_5_LRT_vis.R`
* **Panel C**: Cross-species pausing dispersion ($\sigma$) comparison $\rightarrow$ `analysis/betweenSpecies/1_5_LRT_vis.R`
* **Panel D**: Joint cross-species changes $\Delta\beta$ vs $\Delta\sigma$ $\rightarrow$ `analysis/betweenSpecies/1_5_LRT_vis.R`
* **Panel E**: Representative orthologous gene loci tracks $\rightarrow$ IGV / `analysis/betweenSpecies/1_5_LRT_vis.R`

### Figure 6. Chromatin features are associated with variation in pause-escape kinetics and pausing distributions
* **Panel A**: +1 nucleosome occupancy changes vs $\beta$ changes across cell types $\rightarrow$ `analysis/betweenCellType/2_nucleosome.R`
* **Panel B**: $\text{H3K27ac}$ changes vs $\beta$ changes across cell types $\rightarrow$ `analysis/betweenCellType/3_histoneModification_human.R`
* **Panel C**: Cross-species +1 nucleosome occupancy changes vs $\beta$ changes $\rightarrow$ `analysis/betweenSpecies/2_1_nucleosomeOccupancy_cd4.R`
* **Panel D**: Cross-species nucleosome positioning entropy changes vs $\sigma$ changes $\rightarrow$ `analysis/betweenSpecies/2_2_nucleosomePositioning_cd4.R`

### Figure 7. Promoter sequence features are associated with variation in pause-escape kinetics and pausing distributions
* **Panels A, C**: Promoter GC content and GC skew metaplots $\rightarrow$ `analysis/betweenSpecies/3_4_GCinfo_metaplot.R`
* **Panel B**: GC content changes vs $\beta$ changes across species $\rightarrow$ `analysis/betweenSpecies/3_5_GCinfo_beta_vis.R`
* **Panel D**: GC skew changes vs $\sigma$ changes across species $\rightarrow$ `analysis/betweenSpecies/3_6_GCinfo_sigma_vis.R`

---

## Supplementary Figures

* **Fig. S1 (A–D)**: Simulation scheme and LRT power curves across expression levels $\rightarrow$ `scripts/simulation/`
* **Fig. S2 (A–C)**: Flavopiridol perturbation response $\rightarrow$ `scripts/unimod/analyze_two_samples_DLD1.R`
* **Fig. S3 (A–B)**: $\beta$ vs $\chi$ and +1 nucleosome occupancy across all three primate species $\rightarrow$ `analysis/singleCellType/1_2_EstimateRates_all_species_viz.R`, `2_2_2_nucleosome_position_vis.R`
* **Fig. S4 (A–I)**: Promoter motif classes, nucleosome profiles, and histone marks in $\text{CD4}^+$ and $\text{CD14}^+$ cells $\rightarrow$ `analysis/singleCellType/2_2_3_nucleosome_position_promoter.R`, `3_histone_modification.R`
* **Fig. S5 (A–B)**: $\beta$ vs $\sigma$ and $\mu$ vs $\sigma$ within cell types $\rightarrow$ `analysis/singleCellType/1_2_EstimateRates_all_species_viz.R`
* **Fig. S6 (A–C)**: Sharp vs broad pausing examples and nucleosome organization $\rightarrow$ `analysis/singleCellType/2_2_2_nucleosome_position_vis.R`
* **Fig. S7 (A–G)**: Cell-type kinetic ordering, GO enrichment, and rhesus/baboon LRT $\rightarrow$ `analysis/betweenCellType/1_2_LRT_viz.R`, `analysis/singleCellType/4_GO_human.R`
* **Fig. S8 (A–D)**: Cell-type $\sigma$ transitions and $\Delta\beta$ vs $\Delta\sigma$ in non-human primates $\rightarrow$ `analysis/betweenCellType/1_2_LRT_viz.R`
* **Fig. S9 (A–C)**: Cross-species $\beta$, $\sigma$, and $\Delta\beta$ vs $\Delta\sigma$ distributions $\rightarrow$ `analysis/betweenSpecies/1_5_LRT_vis.R`
* **Fig. S10 (A–H)**: Chromatin features stratified by $\chi$ and positioning entropy across cell types and species $\rightarrow$ `analysis/betweenCellType/2_nucleosome.R`, `3_histoneModification_human.R`, `analysis/betweenSpecies/2_1_nucleosomeOccupancy_cd4.R`, `2_2_nucleosomePositioning_cd4.R`
* **Fig. S11 (A–B)**: Cross-species GC content vs $\sigma$ and GC skew vs $\beta$ $\rightarrow$ `analysis/betweenSpecies/3_5_GCinfo_beta_vis.R`, `3_6_GCinfo_sigma_vis.R`
