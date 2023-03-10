RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

path="porechop/"
dir_trimmed="trimmed/"
dir_aligned="aligned/"
dir_chimera="chimera_removed/"
dir_align2="aligned2/"
dir_no_chimera="chimera_removed/filtered/"

echo -e "${RED}Starting filtering process...${NC}"
mkdir -p ${dir_trimmed}

for file in $(ls ${path})
do
    awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 1200 && length(seq) <= 1800) {print header, seq, qheader, qseq}}' < ${path}${file} > ${dir_trimmed}$(basename "${file}" .fastq)_filt.fastq
done


echo -e "${RED}Starting alignment process...${NC}"
mkdir -p ${dir_aligned}
for file in $(ls ${dir_trimmed})
do
    minimap2 -cx map-ont db/operons.100.fa ${dir_trimmed}${file} -z 70 > ${dir_aligned}${file}_align.paf
done

echo -e "${RED}Starting chimera removal...${NC}"
mkdir -p "${dir_chimera}report/"
mkdir -p "${dir_chimera}filtered/"

for file in $(ls ${dir_aligned})
do
    yacrd -i ${dir_aligned}${file} -o ${dir_chimera}report/$(basename "${file}" .fastq_align.paf).yacrd filter -i ${dir_trimmed}$(basename "${file}" _align.paf) -o ${dir_chimera}filtered/$(basename "${file}" .fastq_align.paf).fasta
done


echo -e "${RED}Starting taxonomic asignment...${NC}"
mkdir -p ${dir_align2}

for file in $(ls ${dir_no_chimera})
do
    minimap2 -cx map-ont db/operons.100.fa ${dir_no_chimera}${file} -z 70 > ${dir_align2}${file}_align.paf
done