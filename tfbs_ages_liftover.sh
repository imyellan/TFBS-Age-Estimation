#!/bin/bash

#SBATCH --mail-user=isaac.yellan95@gmail.com
#SBATCH --mail-type=ARRAY_TASKS,FAIL
#SBATCH --cpus-per-task=3   # number of MPI processes
#SBATCH --mem-per-cpu=1G    # memory; default unit is megabytes
#SBATCH --time=72:00:00         # time (DD-HH:MM)
#SBATCH --array=0-11

module load StdEnv/2020 gcc/9.3.0 hal/2.2 kentutils/453

export proj_dir=/home/iyellan/scratch/tfbs_ages
genome_list="${proj_dir}"/genome_list.txt
genome_list_arr=($(cat "${genome_list}"))

hum_ancestral_genome_list="${proj_dir}"/hum_ancestral_genomes.txt
hum_ancestral_genome_list_arr=($(cut -d" " -f 1 "${hum_ancestral_genome_list}"))

mkdir -p "${proj_dir}"/tmp_beds
mkdir -p "${proj_dir}"/liftover_results

# run_liftover(){
#     in_bed=$1
#     target_spec=$2
#     TF_nm=$(basename "${in_bed}" | cut -d "_" -f 2)
#     proj_dir=/home/iyellan/scratch/tfbs_ages
    
#     halLiftover --outPSLWithName "${proj_dir}"/241-mammalian-2020v2.hal \
#     Homo_sapiens "${in_bed}" "${target_spec}" \
#     "${proj_dir}"/liftover_results/"${TF_nm}"."${target_spec}".txt
# }
# export -f run_liftover

# find "${proj_dir}"/tmp_beds/ -name "*_renamed.bed" \
# | parallel -j ${SLURM_CPUS_PER_TASK} --compress --tmpdir "${SLURM_TMPDIR}" \
# bash ~/scripts/run_liftover.sh {} "${genome_list_arr[$SLURM_ARRAY_TASK_ID]}"

mkdir -p "${proj_dir}"/liftover_results_ancestral
find "${proj_dir}"/tmp_beds/ -name "*_renamed.bed" \
| parallel -j ${SLURM_CPUS_PER_TASK} --compress --tmpdir "${SLURM_TMPDIR}" \
bash ~/scripts/run_liftover.sh {} \
"${hum_ancestral_genome_list_arr[$SLURM_ARRAY_TASK_ID]}" \
"${proj_dir}"/241-mammalian-2020v2.hal \
"${proj_dir}"/liftover_results_ancestral
