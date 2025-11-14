#$ -S /usr/bash
#$ -cwd
#$ -pe threads 1
#$ -N simulation
#$ -o simulation.out
#$ -e simulation.err
#$ -j y
#$ -l m_mem_free=8G
## Run

# a string to record when the script was run
echo "Script is run on $(date +"%Y%m%d%I%M%S%p")"

./simPol.R -p 1 -n 20000 -a 50 -b 1 -g 1 -z 2000 -t 10 -s 33 -k 50 --kSd=0 --addSpace=37 --geneLen=2000 -d ~/projects/Snakemake_projects/unimod_human/outputs/simulation/data
