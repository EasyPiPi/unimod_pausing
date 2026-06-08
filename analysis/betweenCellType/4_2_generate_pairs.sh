#!/usr/bin/env bash
# for EP-contact analysis.
# mamba activate unimodel


########################################################################
# Filter Micro-C pairs (CD4)
########################################################################

# Keep intra-chromosomal, deduplicated pairs with MAPQ >= 30
zcat human_cd4_merged.pairs.gz \
| awk 'BEGIN{OFS="\t"}
       $1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30 {
           print $2, $3, $4, $5, $6, $7, $8, $9, $10
       }' \
> human_cd4_merged.nodups_30_intra.pairs

########################################################################
# Filter Micro-C pairs (CD14)
########################################################################

# Keep intra-chromosomal, deduplicated pairs with MAPQ >= 30
zcat human_cd14_merged.pairs.gz \
| awk 'BEGIN{OFS="\t"}
       $1 == "." && $2 == $4 && $9 >= 30 && $10 >= 30 {
           print $2, $3, $4, $5, $6, $7, $8, $9, $10
       }' \
> human_cd14_merged.nodups_30_intra.pairs
