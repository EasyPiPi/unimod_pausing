#!/usr/bin/env bash
# Prepares human CD4/CD14 potential regulatory elements and common promoters
# for EP-contact analysis.
# mamba activate unimodel

set -euo pipefail

# shellcheck source=../paths.sh
source "$(dirname "$0")/../paths.sh"

ENCODE="${DATA}/encode/peaks"
DREG="${DATA}/dreg"
OUTDIR="${OUTPUTS}/publish/betweenCellType/EP/activeElement"

TSS_DIR="${OUTPUTS}/publish/singleCellType/1_transcriptionRate/human"
CD4_TSS="${TSS_DIR}/cd4_tss_sort_1000bp_extend.bed"
CD14_TSS="${TSS_DIR}/cd14_tss_sort_1000bp_extend.bed"

mkdir -p "${OUTDIR}"

########################################################################
# 1. Sort H3K27ac ENCODE peak files
########################################################################

# CD4+ replicates
for acc in ENCFF880RNP ENCFF993CWV ENCFF353RHA ENCFF780UIQ; do
    sort -k1,1 -k2,2n "${ENCODE}/${acc}.bed" > "${ENCODE}/${acc}.sort.bed"
done

# CD14+ replicates
for acc in ENCFF472QFJ ENCFF634HCJ ENCFF465ZYQ ENCFF506FIB; do
    sort -k1,1 -k2,2n "${ENCODE}/${acc}.bed" > "${ENCODE}/${acc}.sort.bed"
done

########################################################################
# 2. Compute multi-replicate overlap (bedtools multiinter)
########################################################################

bedtools multiinter \
  -i "${ENCODE}"/ENCFF880RNP.sort.bed \
     "${ENCODE}"/ENCFF993CWV.sort.bed \
     "${ENCODE}"/ENCFF353RHA.sort.bed \
     "${ENCODE}"/ENCFF780UIQ.sort.bed \
  > "${OUTDIR}/human_cd4_H3K27ac_multiinter.bed"

bedtools multiinter \
  -i "${ENCODE}"/ENCFF472QFJ.sort.bed \
     "${ENCODE}"/ENCFF634HCJ.sort.bed \
     "${ENCODE}"/ENCFF465ZYQ.sort.bed \
     "${ENCODE}"/ENCFF506FIB.sort.bed \
  > "${OUTDIR}/human_cd14_H3K27ac_multiinter.bed"

########################################################################
# 3. Filter by replicate-overlap threshold; strip 'chr' prefix
#    Column 4 from multiinter = number of replicates containing the region.
#    Thresholds: 2 = at least 2/4, 3 = at least 3/4, 4 = all 4 replicates.
########################################################################

for ct in cd4 cd14; do
    for thresh in 2 3 4; do
        awk -v t="$thresh" 'BEGIN{OFS="\t"} $4 >= t {
            print substr($1, 4), $2, $3
        }' "${OUTDIR}/human_${ct}_H3K27ac_multiinter.bed" \
        > "${OUTDIR}/human_${ct}_H3K27ac_multiinter_filtered${thresh}.bed"
    done
done

########################################################################
# 4. Sort dREG peak files
########################################################################

sort -k1,1 -k2,2n "${DREG}/PROseq-HUMAN-CD4.dREG.peak.prob.bed" \
  > "${DREG}/PROseq-HUMAN-CD4.dREG.peak.prob.sort.bed"

sort -k1,1 -k2,2n "${DREG}/PROseq-HUMAN-CD14.dREG.peak.prob.bed" \
  > "${DREG}/PROseq-HUMAN-CD14.dREG.peak.prob.sort.bed"

########################################################################
# 5. For each threshold: dREG ∩ H3K27ac → add '+' strand → union preys
#    Preys = dREG peaks overlapping H3K27ac (includes promoter-proximal).
#    Threshold-4 outputs are also copied to the original file names so
#    downstream scripts run without changes.
########################################################################

for thresh in 2 3 4; do

    # CD4: dREG peaks overlapping filtered H3K27ac regions
    bedtools intersect \
      -a "${DREG}/PROseq-HUMAN-CD4.dREG.peak.prob.sort.bed" \
      -b "${OUTDIR}/human_cd4_H3K27ac_multiinter_filtered${thresh}.bed" \
      -u \
    | awk 'BEGIN{OFS="\t"} {$4="+"; print}' \
    > "${OUTDIR}/human_cd4_dREG_in_H3K27ac_threshold${thresh}_plus.bed"

    # CD14: dREG peaks overlapping filtered H3K27ac regions
    bedtools intersect \
      -a "${DREG}/PROseq-HUMAN-CD14.dREG.peak.prob.sort.bed" \
      -b "${OUTDIR}/human_cd14_H3K27ac_multiinter_filtered${thresh}.bed" \
      -u \
    | awk 'BEGIN{OFS="\t"} {$4="+"; print}' \
    > "${OUTDIR}/human_cd14_dREG_in_H3K27ac_threshold${thresh}_plus.bed"

    # Combine CD4 + CD14, sort, merge overlapping regions
    cat "${OUTDIR}/human_cd4_dREG_in_H3K27ac_threshold${thresh}_plus.bed" \
        "${OUTDIR}/human_cd14_dREG_in_H3K27ac_threshold${thresh}_plus.bed" \
      | sort -k1,1 -k2,2n -k3,3n \
      | bedtools merge -i stdin -c 4 -o distinct \
      > "${OUTDIR}/human_cd4_cd14_union_preys_plus_threshold${thresh}.bed"

done

# Preserve original output names (threshold-4) for downstream compatibility
cp "${OUTDIR}/human_cd4_dREG_in_H3K27ac_threshold4_plus.bed"  "${OUTDIR}/human_cd4_dREG_in_H3K27ac_plus.bed"
cp "${OUTDIR}/human_cd14_dREG_in_H3K27ac_threshold4_plus.bed" "${OUTDIR}/human_cd14_dREG_in_H3K27ac_plus.bed"
cp "${OUTDIR}/human_cd4_cd14_union_preys_plus_threshold4.bed"  "${OUTDIR}/human_cd4_cd14_union_preys_plus.bed"

########################################################################
# 6. Prepare TSS as baits
########################################################################

awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $1); print $1,$2,$3,$6}' \
  "${CD4_TSS}"  > "${OUTDIR}/cd4_tss_as_baits.bed"

awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $1); print $1,$2,$3,$6}' \
  "${CD14_TSS}" > "${OUTDIR}/cd14_tss_as_baits.bed"

########################################################################
# 7. Common promoters
########################################################################

# Sorted copies needed for comm and center matching (temporary)
sort "${OUTDIR}/cd4_tss_as_baits.bed"  > "${OUTDIR}/cd4_tss_as_baits.sorted.tmp.bed"
sort "${OUTDIR}/cd14_tss_as_baits.bed" > "${OUTDIR}/cd14_tss_as_baits.sorted.tmp.bed"

# 8a. Exact match: identical coordinates and strand in both cell types
comm -12 "${OUTDIR}/cd4_tss_as_baits.sorted.tmp.bed" \
         "${OUTDIR}/cd14_tss_as_baits.sorted.tmp.bed" \
  > "${OUTDIR}/human_common_promoters.bed"

# 8b. Relaxed match: promoters are considered common if their TSS centers
#     (midpoints) are within 200 bp of each other.
#
#     The bait regions are already TSS ±1000 bp windows (~2000 bp wide), so
#     applying bedtools window directly on the intervals would be overly
#     permissive — two promoters with TSS positions 1800 bp apart could still
#     be flagged as "common" due to interval overlap alone. Using midpoints
#     keeps the criterion anchored to the underlying TSS locations.

# Collapse each promoter region to its midpoint; tag rows with line number
# so we can recover the original coordinates after the window query.
awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); print $1, c, c+1, NR}' \
  "${OUTDIR}/cd4_tss_as_baits.sorted.tmp.bed"  > "${OUTDIR}/cd4_centers.tmp.bed"
awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); print $1, c, c+1, NR}' \
  "${OUTDIR}/cd14_tss_as_baits.sorted.tmp.bed" > "${OUTDIR}/cd14_centers.tmp.bed"

# Find which CD4 rows have a CD14 center within 200 bp, and vice versa
bedtools window -w 200 -a "${OUTDIR}/cd4_centers.tmp.bed" -b "${OUTDIR}/cd14_centers.tmp.bed" -u \
  | awk '{print $4}' | sort -un > "${OUTDIR}/cd4_center_rows.tmp.txt"
bedtools window -w 200 -a "${OUTDIR}/cd14_centers.tmp.bed" -b "${OUTDIR}/cd4_centers.tmp.bed" -u \
  | awk '{print $4}' | sort -un > "${OUTDIR}/cd14_center_rows.tmp.txt"

# Recover original promoter coordinates for matching rows
awk 'NR==FNR{keep[$1]=1; next} keep[FNR]' \
  "${OUTDIR}/cd4_center_rows.tmp.txt"  "${OUTDIR}/cd4_tss_as_baits.sorted.tmp.bed"  > "${OUTDIR}/cd4_relaxed.tmp.bed"
awk 'NR==FNR{keep[$1]=1; next} keep[FNR]' \
  "${OUTDIR}/cd14_center_rows.tmp.txt" "${OUTDIR}/cd14_tss_as_baits.sorted.tmp.bed" > "${OUTDIR}/cd14_relaxed.tmp.bed"

# Combine CD4 + CD14 relaxed matches, merge overlapping regions
cat "${OUTDIR}/cd4_relaxed.tmp.bed" "${OUTDIR}/cd14_relaxed.tmp.bed" \
  | sort -k1,1 -k2,2n -k3,3n \
  | bedtools merge -i stdin -c 4 -o distinct \
  > "${OUTDIR}/human_common_promoters_relaxed.bed"

########################################################################
# 9. Cleanup
########################################################################

# --- Temporary files: removed after use ---

# Sorted ENCODE peak files (originals retained; re-sortable if needed)
rm -f "${ENCODE}"/ENCFF{880RNP,993CWV,353RHA,780UIQ}.sort.bed \
      "${ENCODE}"/ENCFF{472QFJ,634HCJ,465ZYQ,506FIB}.sort.bed

# Raw multiinter outputs (superseded by per-threshold filtered files in OUTDIR)
rm -f "${OUTDIR}/human_cd4_H3K27ac_multiinter.bed" \
      "${OUTDIR}/human_cd14_H3K27ac_multiinter.bed"

# Sorted bait copies and all center/row index files used only in section 8
rm -f "${OUTDIR}/cd4_tss_as_baits.sorted.tmp.bed" \
      "${OUTDIR}/cd14_tss_as_baits.sorted.tmp.bed" \
      "${OUTDIR}/cd4_centers.tmp.bed"     "${OUTDIR}/cd14_centers.tmp.bed" \
      "${OUTDIR}/cd4_center_rows.tmp.txt" "${OUTDIR}/cd14_center_rows.tmp.txt" \
      "${OUTDIR}/cd4_relaxed.tmp.bed"     "${OUTDIR}/cd14_relaxed.tmp.bed"

# --- Final outputs retained in ${OUTDIR} ---
#
# Threshold-filtered H3K27ac regions (CD4 and CD14, thresholds 2/3/4):
#   human_{cd4,cd14}_H3K27ac_multiinter_filtered{2,3,4}.bed
#
# Per-threshold dREG/H3K27ac overlaps:
#   human_{cd4,cd14}_dREG_in_H3K27ac_threshold{2,3,4}_plus.bed
#
# Per-threshold union prey BEDs:
#   human_cd4_cd14_union_preys_plus_threshold{2,3,4}.bed
#
# Backward-compatible threshold-4 names (for downstream scripts):
#   human_cd4_dREG_in_H3K27ac_plus.bed
#   human_cd14_dREG_in_H3K27ac_plus.bed
#   human_cd4_cd14_union_preys_plus.bed
#
# TSS bait files:
#   cd4_tss_as_baits.bed
#   cd14_tss_as_baits.bed
#
# Common promoters:
#   human_common_promoters.bed         (exact match)
#   human_common_promoters_relaxed.bed (centers within 200 bp)
