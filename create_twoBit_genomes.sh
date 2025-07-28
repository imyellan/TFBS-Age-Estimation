#!/bin/bash

#SBATCH --mail-user=isaac.yellan95@gmail.com
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --cpus-per-task=1   # number of MPI processes
#SBATCH --mem=5G
#SBATCH --time=24:00:00         # time (DD-HH:MM)
#SBATCH --array=0-239

module load StdEnv/2020 gcc/9.3.0 hal/2.2 kentutils/453 

readarray -t genomes < /home/iyellan/scratch/tfbs_ages/genome_list.txt

## check if fasta file both exists and has a size greater than 0, if not, create it
if [[ ! -s ~/scratch/tfbs_ages/genomes/"${genomes[$SLURM_ARRAY_TASK_ID]}".fa ]]; then
    hal2fasta ~/scratch/tfbs_ages/241-mammalian-2020v2.hal "${genomes[$SLURM_ARRAY_TASK_ID]}" \
    > ~/scratch/tfbs_ages/genomes/"${genomes[$SLURM_ARRAY_TASK_ID]}".fa
    echo "${genomes[$SLURM_ARRAY_TASK_ID]} extracted"
fi

## check if twoBit file both exists and has a size greater than 16, if not, create it
if [[ ! -s ~/scratch/tfbs_ages/genomes/"${genomes[$SLURM_ARRAY_TASK_ID]}".2bit ]]; then
    faToTwoBit ~/scratch/tfbs_ages/genomes/"${genomes[$SLURM_ARRAY_TASK_ID]}".fa \
    ~/scratch/tfbs_ages/genomes/"${genomes[$SLURM_ARRAY_TASK_ID]}".2bit
    echo "${genomes[$SLURM_ARRAY_TASK_ID]} converted to 2bit"
fi