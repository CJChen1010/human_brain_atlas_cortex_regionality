REGION_BED="/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/all/BICCN_5K_regions.bed"
GENOME_FASTA="/tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa"
CHROMSIZES="/tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/hg38.nochrM.chrom.sizes"
DATABASE_PREFIX="BICCN_subclass"
SCRIPT_DIR="/tscc/nfs/home/biy022/softwares/create_cisTarget_databases"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
    ${GENOME_FASTA} \
    ${CHROMSIZES} \
    ${REGION_BED} \
    hg38_BICCN_subclass_with_1kb_bg_padding.fa \
    1000 \
    yes