#!/bin/bash

in_bed=$1
target_spec=$2
hal_fil=$3
outdir=$4
# TF_nm=$(basename "${in_bed}" | cut -d "_" -f 2)
TF_nm=$(basename "${in_bed%%.bed}" | sed 's/_.*//g')
# proj_dir=/home/iyellan/scratch/tfbs_ages
    
halLiftover --bedType 4 --outPSLWithName "${hal_fil}" Homo_sapiens "${in_bed}" \
"${target_spec}" "${outdir}"/"${TF_nm}"."${target_spec}".txt