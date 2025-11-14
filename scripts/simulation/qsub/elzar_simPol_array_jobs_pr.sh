#$ -S /usr/bash
#$ -cwd
#$ -pe threads 1
#$ -N simulation_array
#$ -o simulation_array.out
#$ -e simulation_array.err
#$ -j y
#$ -l m_mem_free=8G
#$ -t 1-108
## Run

# Simulation with a wide range of parameters for rate estimates
# a string to record when the script was run
echo "Script is run on $(date +"%Y%m%d%I%M%S%p")"
# a single sh file contains all the CMDs
SHFILE=simulation_rate_pause_release.sh
# each time pull out a single line to run
SH=$(cat $SHFILE | head -n $SGE_TASK_ID | tail -n 1)
# run the chosen CMD
eval "$SH"
