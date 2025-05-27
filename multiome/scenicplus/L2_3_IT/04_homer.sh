ROOT=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L2_3_IT
INPUT_DIR=${ROOT}/corr_regions_ST_5e-1_5e-2_expr_filtered

cd ${INPUT_DIR}
mkdir -p sbatch_scripts
mkdir -p sbatch_logs

awk -F '[-]' '{print $1 "\t" $2 "\t" $3 "\tpositive_" NR}' L2_3_IT_positive.tsv >L2_3_IT_positive.bed
awk -F '[-]' '{print $1 "\t" $2 "\t" $3 "\tnegative_" NR}' L2_3_IT_negative.tsv >L2_3_IT_negative.bed

for pathname in *tive.bed; do
    subname=$(basename ${pathname##*_} .bed)
    mkdir -p ${subname}_homer
    cat >sbatch_scripts/${subname}.sh <<EOF
#! /bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH -J L2_3_IT_${subname}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem 100G
#SBATCH -t 4:00:00
#SBATCH -o $INPUT_DIR/sbatch_logs/${subname}.out
#SBATCH -e $INPUT_DIR/sbatch_logs/${subname}.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
conda activate pycistarget
cd ${INPUT_DIR}/

if [ ! -f "${subname}_homer/knownResults.txt" ]
then
    findMotifsGenome.pl \\
        ${INPUT_DIR}/${pathname} \\
        /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa \\
        ${subname}_homer \\
        -size given -mask -p 12
fi
EOF
done


## Use a different background file
root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L2_3_IT/corr_regions_ST_5e-1_5e-2_expr_filtered
cd $root_dir

awk 'BEGIN{FS="-";OFS="\t"}{print $1,$2,$3}' L2_3_IT_positive.tsv >L2_3_IT_positive_test.bed
awk 'BEGIN{FS="-";OFS="\t"}{print $1,$2,$3}' L2_3_IT_negative.tsv >L2_3_IT_negative_test.bed
awk 'BEGIN{FS="-";OFS="\t"}{print $1,$2,$3}' L2_3_IT_full_list.tsv >L2_3_IT_full_list.bed

bedtools intersect \
    -a L2_3_IT_full_list.bed \
    -b L2_3_IT_positive_test.bed L2_3_IT_negative_test.bed \
    -v >L2_3_IT_background.bed


for file in *test.bed; do
    subname=$(basename $file _test.bed)
    subname=${subname##L2_3_IT_}
    mkdir -p ${subname}_homer_bg
    cat >sbatch_scripts/${subname}_bg.sh <<EOF
#! /bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH -J L2_3_IT_${subname}_bg
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem 100G
#SBATCH -t 4:00:00
#SBATCH -o $root_dir/sbatch_logs/${subname}_bg.out
#SBATCH -e $root_dir/sbatch_logs/${subname}_bg.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
conda activate pycistarget
cd $root_dir

if [ ! -f "${subname}_homer_bg/knownResults.txt" ]
then
    findMotifsGenome.pl \\
        $file \\
        /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa \\
        ${subname}_homer_bg \\
        -size given -mask -p 12 \\
        -bg L2_3_IT_background.bed
fi
EOF
done


for file in *test.bed; do
    subname=$(basename $file _test.bed)
    subname=${subname##L2_3_IT_}
    mkdir -p ${subname}_homer_bg_selected_tfs
    cat >sbatch_scripts/${subname}_bg_selected_tfs.sh <<EOF
#! /bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH -J L2_3_IT_${subname}_bg_st
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem 100G
#SBATCH -t 4:00:00
#SBATCH -o $root_dir/sbatch_logs/${subname}_bg_st.out
#SBATCH -e $root_dir/sbatch_logs/${subname}_bg_st.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
conda activate pycistarget
cd $root_dir

if [ ! -f "${subname}_homer_bg_selected_tfs/knownResults.txt" ]
then
    findMotifsGenome.pl \\
        $file \\
        /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa \\
        ${subname}_homer_bg_selected_tfs \\
        -size given -mask \\
        -bg L2_3_IT_background.bed \\
        -p 20 -nomotif \\
        -mknown selected_tfs.motif
fi
EOF
done




