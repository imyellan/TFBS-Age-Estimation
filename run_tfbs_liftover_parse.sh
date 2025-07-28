#!/bin/bash

#SBATCH --mail-user=isaac.yellan95@gmail.com
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --cpus-per-task=8   # number of MPI processes
#SBATCH --mem=15G    
#SBATCH --time=24:00:00         # time (DD-HH:MM)
#SBATCH --array=0-135
#SBATCH --requeue

# Handle memory scaling on restart
RESTART_COUNT=${SLURM_RESTART_COUNT:-0}
if [ $RESTART_COUNT -gt 0 ]; then
    # Calculate new memory: 15G * (1.5^restart_count)
    NEW_MEM=$(echo "15 * (1.5^$RESTART_COUNT)" | bc -l | cut -d. -f1)
    echo "Restart #$RESTART_COUNT: Using ${NEW_MEM}G memory"
    # Note: This requires resubmitting the job manually with new memory
fi
readarray -t TFs < /home/iyellan/scratch/tfbs_ages/TF_list_v2.txt

module load StdEnv/2023 gcc/12.3 r/4.3.1 r-bundle-bioconductor/3.18
Rscript /home/iyellan/tfbs_age_estimation/tfbs_liftover_parse.R \
/home/iyellan/scratch/tfbs_ages/liftover_results "${TFs[$SLURM_ARRAY_TASK_ID]}" \
nonancestral
Rscript /home/iyellan/tfbs_age_estimation/tfbs_liftover_parse.R \
/home/iyellan/scratch/tfbs_ages/liftover_results_ancestral "${TFs[$SLURM_ARRAY_TASK_ID]}" \
ancestral

module load StdEnv/2020 gcc/9.3.0 hal/2.2 kentutils/453 
# find ~/scratch/tfbs_ages/liftover_beds -name *"${TFs[$SLURM_ARRAY_TASK_ID]}"_*.bed \
# | parallel -j 11 ~/tfbs_age_estimation/run_twoBitToFa.sh {}
find ~/scratch/tfbs_ages/liftover_beds \
-name *"${TFs[$SLURM_ARRAY_TASK_ID]}"_*.bed \
| parallel -j 8 ~/tfbs_age_estimation/run_twoBitToFa.sh {}
find ~/scratch/tfbs_ages/liftover_beds_ancestral \
-name *"${TFs[$SLURM_ARRAY_TASK_ID]}"_*.bed \
| parallel -j 8 ~/tfbs_age_estimation/run_twoBitToFa.sh {}

module load StdEnv/2023 gcc/12.3 r/4.3.1 r-bundle-bioconductor/3.18
Rscript /home/iyellan/tfbs_age_estimation/summarize_tfbs_ages.R \
/home/iyellan/scratch/tfbs_ages/liftover_results "${TFs[$SLURM_ARRAY_TASK_ID]}"
Rscript /home/iyellan/tfbs_age_estimation/summarize_tfbs_ages_ancestral.R \
/home/iyellan/scratch/tfbs_ages/liftover_results_ancestral "${TFs[$SLURM_ARRAY_TASK_ID]}"

