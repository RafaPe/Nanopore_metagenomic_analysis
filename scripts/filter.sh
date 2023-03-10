path="porechop/"
out_dir="trimmed/"
for file in $(ls ${path})
do
    echo ${file}
    awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 1200 && length(seq) <= 1800) {print header, seq, qheader, qseq}}' < ${path}${file} > ${out_dir}$(basename "${file}" .fastq)_filt.fastq
done