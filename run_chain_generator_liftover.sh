#!/bin/bash

in_bed=$1
in_bed_base=$(basename "${in_bed}")
target_spec=$2
proj_dir=/home/iyellan/scratch/tfbs_ages

halLiftover --outPSL "${proj_dir}"/241-mammalian-2020v2.hal Homo_sapiens \
"$in_bed" "${target_spec}" /dev/stdout \
| pslPosTarget stdin "${proj_dir}"/chains/"${in_bed_base}"_Homo_sapiens-to-"${target_spec}".psl