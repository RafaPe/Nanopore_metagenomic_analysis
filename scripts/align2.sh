no_chimera_dir="chimera_removed/filtered/"
out_dir="aligned2/"

mkdir -p ${out_dir}

for file in $(ls ${no_chimera_dir})
do
    echo ${file}
    minimap2 -cx map-ont db/operons.100.fa ${no_chimera_dir}${file} -z 70 > ${out_dir}${file}_align.paf
done