path_to_porechopexe="../Porechop/"
path_to_raw_data="Whisky_metagenomics/3 Days & 5 Days Bacteria Raw Fastq"

# for file in ${path_to_raw_data}; do mv "${file}" `echo ${file} | tr ' ' '_'` ; done

for file in $(ls "${path_to_raw_data}")
do
    echo "${file}"
done