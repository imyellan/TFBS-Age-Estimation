in_bed=$1

if [[ "${in_bed}" =~ ancestral ]]; then
  genome=~/scratch/tfbs_ages/genomes/$(basename $in_bed | sed -E 's/.*_([A-Za-z0-9]+)_.*/\1/g').2bit
else
    genome=~/scratch/tfbs_ages/genomes/$(basename $in_bed | sed -E 's/.*_([A-Z][A-Za-z]+_[A-Za-z]+(_[A-Za-z]+)?).bed/\1/g').2bit
fi
# if [[ ! -f "${genome}" ]]; then
#     species=$(basename "${in_bed}" | sed -E 's/.*_([A-Za-z0-9]+)_.*/\1/g')
#     hal2fasta ~/scratch/tfbs_ages/241-mammalian-2020v2.hal "${species}" \
#     > ~/scratch/tfbs_ages/genomes/"${species}".fa
#     faToTwoBit ~/scratch/tfbs_ages/genomes/"${species}".fa \
#     ~/scratch/tfbs_ages/genomes/"${species}".2bit
# fi
twoBitToFa -noMask -bed="${in_bed}" "${genome}" "${in_bed%%.bed}".fa
