path="trimmed/"
out_dir="aligned/"

mkdir -p ${out_dir}
for file in $(ls ${path})
do
    echo ${file}
    minimap2 -cx map-ont db/operons.100.fa ${path}${file} -z 70 > ${out_dir}${file}_align.paf
done