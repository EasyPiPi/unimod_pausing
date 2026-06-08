#!/usr/bin/env bash
# paths.sh — shell path helper for the publish pipeline.
#
# Source from any script in codes/publish/:
#   source "$(dirname "$0")/../paths.sh"
#
# Set PROJECT_ROOT in the environment before sourcing to override auto-detection.

# This file lives at codes/publish/paths.sh; two levels up is the project root.
_PATHS_SH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
: "${PROJECT_ROOT:="$(cd "${_PATHS_SH_DIR}/../.." && pwd)"}"
unset _PATHS_SH_DIR

DATA="${PROJECT_ROOT}/data/publish"
OUTPUTS="${PROJECT_ROOT}/outputs"

ENCODE_PEAKS="${DATA}/encode/peaks"
DREG_BED="${DATA}/dreg"
LIFTOVER_DIR="${DATA}/liftover"

: "${ENCODE_BAM_ROOT:="${PROJECT_ROOT}/data/result2/encode/bam"}"
: "${ENCODE_BW_ROOT:="${PROJECT_ROOT}/data/result2/encode/bw_clean"}"
: "${GENOME_ROOT:="${PROJECT_ROOT}/codes/micro-c"}"
: "${MICROC_1D_ROOT:="${DATA}/micro-c"}"
: "${MICROC_PAIRS_ROOT:=""}"
