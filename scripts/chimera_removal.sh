path_aligned="aligned/"
path_trimmed="trimmed/"
out_dir="chimera_removed/"

mkdir -p "${out_dir}report/"
mkdir -p "${out_dir}filtered/"


for file in $(ls ${path_aligned})
do
    echo ${file}
    yacrd -i ${path_aligned}${file} -o ${out_dir}report/$(basename "${file}" .fastq_align.paf).yacrd filter -i ${path_trimmed}$(basename "${file}" _align.paf) -o ${out_dir}filtered/$(basename "${file}" .fastq_align.paf).fasta
done

# echo ${files1:^files2}

# for file ${files1:^files2}
# do
#     echo ${f1}
#     # echo ${f2}
# done
    # OCI.py $f1 $f2 3/$f1:t:r+$f2:t
