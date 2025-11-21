#!/bin/bash
#$ -l m_mem_free=50G
#$ -pe threads 8
#$ -o /grid/siepel/home/zeng/projects/Zhao_2025/code/Xin_code/simulation/batch_subsample_v2.log
#$ -j y

#$ -cwd

source ~/.bashrc
mamba activate compara_reg_r

DIR="/grid/siepel/home/zeng/projects/Zhao_2025/simulation"

for ITEM in "$DIR"/*; do
	NAME=$(basename "$ITEM")
	EXPR=(low median high)
	for L in "${EXPR[@]}"; do
		echo "Running sample $NAME with lambda = $L"
		Rscript subsample_simulation_pause_release_Xin.R "$NAME" "$L"
	done
done

#Rscript subsample_simulation_pause_release_Xin.R k75ksd20kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t40n20000s33h17 'high'