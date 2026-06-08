#!/usr/bin/env bash
# NOTE: This script runs on an HPC cluster where Micro-C .pairs.gz files live.
# Set MICROC_PAIRS_ROOT (location of pairs files) and PROJECT_ROOT before running,
# or configure them in codes/publish/paths.yaml.
set -euo pipefail

# shellcheck source=../paths.sh
source "$(dirname "$0")/../paths.sh"

# MICROC_PAIRS_ROOT must be set — no default exists because the pairs files are
# HPC-specific and are not distributed with the publish package.
if [[ -z "${MICROC_PAIRS_ROOT}" ]]; then
  echo "ERROR: MICROC_PAIRS_ROOT is not set. Set it in the environment or in" >&2
  echo "       codes/publish/paths.yaml before running this script." >&2
  exit 1
fi

species_list=("human" "rhesus" "baboon")
celltype_list=("cd4" "cd14")

for species in "${species_list[@]}"
do
    for celltype in "${celltype_list[@]}"
    do

        pair_path="${MICROC_PAIRS_ROOT}/${species}_${celltype}_merged.pairs.gz"

        tss_file="${OUTPUTS}/publish/singleCellType/1_transcriptionRate/${species}/${celltype}_tss_sort.bed"

        output_file="${MICROC_1D_ROOT}/${species}_${celltype}_1d_signal.bed"

        echo "Processing ${species} ${celltype}"

        if [[ ! -f "${pair_path}" ]]; then
            echo "Missing pair file: ${pair_path}"
            continue
        fi

        if [[ ! -f "${tss_file}" ]]; then
            echo "Missing tss file: ${tss_file}"
            continue
        fi

        python extract_1d_signal.py \
            --pair_path "${pair_path}" \
            --tss_file "${tss_file}" \
            --output_file "${output_file}"

    done
done