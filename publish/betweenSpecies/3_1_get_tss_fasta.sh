#!/usr/bin/env bash
set -euo pipefail

# shellcheck source=../paths.sh
source "$(dirname "$0")/../paths.sh"

EXT=1000

extend_bed() {
  local in_bed="$1"
  local out_bed="$2"
  local ext="${3:-1000}"

  awk -v OFS="\t" -v EXT="$ext" '{
      start=$2-EXT; if(start<0) start=0;
      end=$3+EXT;
      chr=$1;
      if($1 !~ /^chr/) chr="chr"$1;
      print chr, start, end, $4, $5, $6
  }' "$in_bed" > "$out_bed"
}

get_fasta() {
  local fasta="$1"
  local bed="$2"
  local out_fa="$3"
  bedtools getfasta -fi "$fasta" -bed "$bed" -fo "$out_fa" -name -s
}

run_one() {
  local in_bed="$1"      # absolute path to input bed
  local fasta="$2"       # absolute path to fasta
  local out_dir="$3"     # absolute path to output directory
  local ext="${4:-1000}"

  mkdir -p "$out_dir"

  local base
  base=$(basename "$in_bed" .bed)

  local ext_bed="${out_dir}_${base}_${ext}bp_extend.bed"
  local out_fa="${out_dir}_${base}_${ext}bp_extend.fa"

  extend_bed "$in_bed" "$ext_bed" "$ext"
  get_fasta "$fasta" "$ext_bed" "$out_fa"
}

########################################
INPUT_DIR="${OUTPUTS}/publish/singleCellType/1_transcriptionRate"
RESULT_DIR="${OUTPUTS}/publish/betweenSpecies/sequence"
HG38="${GENOME_ROOT}/hg38.fa"
PAPANU4="${GENOME_ROOT}/papAnu4.fa"
RHEMAC10="${GENOME_ROOT}/rheMac10.fa"

########################################
##The input of this file is from 2_1_1_generateTSSregion.r
# human
run_one "${INPUT_DIR}/human/cd4_tss_sort.bed"  "$HG38"    "${RESULT_DIR}/human" "$EXT"
run_one "${INPUT_DIR}/human/cd14_tss_sort.bed" "$HG38"    "${RESULT_DIR}/human" "$EXT"

# baboon
run_one "${INPUT_DIR}/baboon/cd4_tss_sort.bed"  "$PAPANU4" "${RESULT_DIR}/baboon" "$EXT"
run_one "${INPUT_DIR}/baboon/cd14_tss_sort.bed" "$PAPANU4" "${RESULT_DIR}/baboon" "$EXT"

# rhesus
run_one "${INPUT_DIR}/rhesus/cd4_tss_sort.bed"  "$RHEMAC10" "${RESULT_DIR}/rhesus"  "$EXT"
run_one "${INPUT_DIR}/rhesus/cd14_tss_sort.bed" "$RHEMAC10" "${RESULT_DIR}/rhesus" "$EXT"


###