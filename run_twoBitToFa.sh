in_bed=$1
genome=~/scratch/tfbs_ages/genomes/$(basename $in_bed | sed -E 's/.*_([A-Z][A-Za-z]+_[A-Za-z]+(_[A-Za-z]+)?).bed/\1/g').2bit
twoBitToFa -noMask -bed="${in_bed}" "${genome}" "${in_bed%%.bed}".fa
