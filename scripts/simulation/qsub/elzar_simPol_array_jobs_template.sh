#$ -S /usr/bash
#$ -cwd
#$ -pe threads 48
#$ -N simulation_array
#$ -o simulation_array.out
#$ -e simulation_array.err
#$ -j y
#$ -l m_mem_free=1G
#$ -t 1-9
## Run

# each sh contains some cmds to run
bash simulation_${SGE_TASK_ID}.sh
