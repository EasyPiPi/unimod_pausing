output_file="simulation_pause_distribution.sh"

> "$output_file"

OUTPUT_DIR="${SIMULATION_OUTPUT_DIR:-outputs/simulation/data_fk}"
SIMPOL_SIF="${SIMPOL_SIF:-tools/simpol.sif}"
SIMPOL_BIN="${SIMPOL_BIN:-simPol_Release}"

for mean in 30 40 45 50 55 60 65 70 75 80 90
do
  for sd in 10 20 25 30 35 40 45 50 55 60 70  
  do
    if [[ -n "${USE_SINGULARITY:-}" ]]; then
      echo "singularity exec $SIMPOL_SIF $SIMPOL_BIN -k $mean --kSd=${sd} -n 20000 -a 1 -b 1 -z 2000 -t 40 -s 33 --addSpace=17 --geneLen=2000 --zetaSd=1000 --zetaMax=2500 --zetaMin=1500 -d $OUTPUT_DIR/k${mean}ksd${sd}kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t40n20000s33h17" >> "$output_file"
    else
      echo "$SIMPOL_BIN -k $mean --kSd=${sd} -n 20000 -a 1 -b 1 -z 2000 -t 40 -s 33 --addSpace=17 --geneLen=2000 --zetaSd=1000 --zetaMax=2500 --zetaMin=1500 -d $OUTPUT_DIR/k${mean}ksd${sd}kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t40n20000s33h17" >> "$output_file"
    fi
  done
done