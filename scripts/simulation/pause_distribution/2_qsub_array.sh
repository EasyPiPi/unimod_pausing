#$ -S /usr/bash
#$ -cwd
#$ -pe OpenMP 5
#$ -N simulation_array_v1
#$ -o simulation_array_v1.out
#$ -e simulation_array_v1.err
#$ -j y
#$ -l m_mem_free=2G
#$ -t 1-121
## Run


# Simulation with a wide range of parameters for rate estimates
# a string to record when the script was run
echo "Script is run on $(date +"%Y%m%d%I%M%S%p")"
# a single sh file contains all the CMDs
SHFILE=simulation_pause_distribution.sh
# each time pull out a single line to run
SH=$(cat $SHFILE | head -n $SGE_TASK_ID | tail -n 1)
# run the chosen CMD
eval "$SH"
