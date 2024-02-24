#!/bin/bash

#SBATCH --mail-user=isaac.yellan95@gmail.com
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --cpus-per-task=8   # number of MPI processes
#SBATCH --mem=90G    
#SBATCH --time=8:00:00         # time (DD-HH:MM)
#SBATCH --array=0-136

readarray -t TFs < /home/iyellan/scratch/tfbs_ages/TF_list.txt

module load StdEnv/2023 gcc/12.3 r/4.3.1 r-bundle-bioconductor/3.18
Rscript /home/iyellan/scripts/tfbs_liftover_parse.R \
/home/iyellan/scratch/tfbs_ages/liftover_results "${TFs[$SLURM_ARRAY_TASK_ID]}"

module load StdEnv/2020 gcc/9.3.0 kentutils/453 
find ~/scratch/tfbs_ages/liftover_beds -name *"${TFs[$SLURM_ARRAY_TASK_ID]}"_*.bed \
| parallel -j 11 ~/scripts/run_twoBitToFa.sh {}

module load StdEnv/2023 gcc/12.3 r/4.3.1 r-bundle-bioconductor/3.18
Rscript /home/iyellan/scripts/summarize_tfbs_ages.R \
/home/iyellan/scratch/tfbs_ages/liftover_results "${TFs[$SLURM_ARRAY_TASK_ID]}"

