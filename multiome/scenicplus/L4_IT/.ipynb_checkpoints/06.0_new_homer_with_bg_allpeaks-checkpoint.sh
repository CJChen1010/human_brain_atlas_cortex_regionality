root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/corr_regions_5e-1_5e-2_expr_filtered
cd $root_dir

# convert tsv file to bed file
awk 'BEGIN{FS="-";OFS="\t"}{print $1,$2,$3}' L4_IT_positive.tsv >L4_IT_positive_test.bed
awk 'BEGIN{FS="-";OFS="\t"}{print $1,$2,$3}' L4_IT_negative.tsv >L4_IT_negative_test.bed
awk 'BEGIN{FS="-";OFS="\t"}{print $1,$2,$3}' L4_IT_full_list.tsv >L4_IT_full_list.bed

bedtools intersect \
    -a L4_IT_full_list.bed \
    -b L4_IT_positive_test.bed L4_IT_negative_test.bed \
    -v >L4_IT_background.bed

mkdir -p sbatch_scripts
for file in *test.bed; do
    subname=$(basename $file _test.bed)
    subname=${subname##L4_IT_}
    mkdir -p ${subname}_homer
    cat >sbatch_scripts/${subname}.sh <<EOF
#! /bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH -J L4_IT_${subname}
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 8
#SBATCH --mem 100G
#SBATCH -t 4:00:00
#SBATCH -o $root_dir/sbatch_logs/${subname}.out
#SBATCH -e $root_dir/sbatch_logs/${subname}.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
conda activate pycistarget
cd $root_dir

if [ ! -f "${subname}_homer/knownResults.txt" ]
then
    findMotifsGenome.pl \\
        $file \\
        /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa \\
        ${subname}_homer \\
        -size given -mask -p 12 \\
        -bg L4_IT_background.bed
fi
EOF
done

mkdir -p sbatch_logs

for file in *test.bed; do
    subname=$(basename $file _test.bed)
    subname=${subname##L4_IT_}
    mkdir -p ${subname}_homer_selected_tfs
    cat >sbatch_scripts/${subname}_selected_tfs.sh <<EOF
#! /bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH -J L4_IT_${subname}_st
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 8
#SBATCH --mem 100G
#SBATCH -t 4:00:00
#SBATCH -o $root_dir/sbatch_logs/${subname}_st.out
#SBATCH -e $root_dir/sbatch_logs/${subname}_st.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
conda activate pycistarget
cd $root_dir

if [ ! -f "${subname}_homer_selected_tfs/knownResults.txt" ]
then
    findMotifsGenome.pl \\
        $file \\
        /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa \\
        ${subname}_homer_selected_tfs \\
        -size given -mask -p 12 \\
        -bg L4_IT_background.bed \\
        -p 20 -nomotif \\
        -mknown selected_tfs.motif
fi
EOF
done


