RESDIR=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/rostral_caudal_spearman/peaks_bh_5e-1_5e-2_donor_replicate
ROOT=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/rostral_caudal_spearman
cd $ROOT
mkdir -p ${RESDIR}/pbs_scripts
mkdir -p ${RESDIR}/pbs_logs

for file in *_peaks_bh_5e-1_5e-2_donor_replicate.bed
do
    subname=$(basename $file _peaks_bh_5e-1_5e-2_donor_replicate.bed)
    mkdir -p ${RESDIR}/${subname}_correlation
    cat >${RESDIR}/pbs_scripts/${subname}.sh <<EOF
#! /bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH -J L4_IT_${subname}
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 8
#SBATCH --mem 50G
#SBATCH -t 4:00:00
#SBATCH -o ${RESDIR}/pbs_logs/${subname}.out
#SBATCH -e ${RESDIR}/pbs_logs/${subname}.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
conda activate pycistarget
cd ${RESDIR}/

if [ ! -f "${subname}_correlation/knownResults.txt" ]
then
    findMotifsGenome.pl \\
        ${ROOT}/${file} \\
        /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa \\
        ${subname}_correlation \\
        -size given -mask -p 12
fi
EOF
done


## Some scripts manually run for a set of TFs
findMotifsGenome.pl positive_peaks_bh_5e-1_5e-2_donor_replicate.bed /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/positive_correlation/ -size given -mask -p 20 -nomotif -mknown peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/selected_tfs.motif -find peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/selected_tfs.motif >peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/positive_correlation/selected_tf_occurences.tsv

findMotifsGenome.pl negative_peaks_bh_5e-1_5e-2_donor_replicate.bed /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/negative_correlation/ -size given -mask -p 20 -nomotif -mknown peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/selected_tfs.motif -find peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/selected_tfs.motif >peaks_bh_5e-1_5e-2_donor_replicate_selected_tfs/negative_correlation/selected_tf_occurences.tsv





## I rerun the correlation pipeline, and kept genes with padjusted < 0.01 without cutoffs on spearman scores
## The above script is for genes with padjusted < 0.05 and |score| > 0.5
RESDIR=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/rostral_caudal_spearman/peaks_bh_1e-2_donor_replicate
ROOT=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/rostral_caudal_spearman
cd $ROOT
mkdir -p ${RESDIR}/pbs_scripts
mkdir -p ${RESDIR}/pbs_logs

for file in *_peaks_bh_1e-2_donor_replicate.bed
do
    subname=$(basename $file _peaks_bh_1e-2_donor_replicate.bed)
    mkdir -p ${RESDIR}/${subname}_correlation
    cat >${RESDIR}/pbs_scripts/${subname}.sh <<EOF
#! /bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A csd772
#SBATCH -J L4_IT_${subname}
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 8
#SBATCH --mem 50G
#SBATCH -t 4:00:00
#SBATCH -o ${RESDIR}/pbs_logs/${subname}.out
#SBATCH -e ${RESDIR}/pbs_logs/${subname}.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
conda activate pycistarget
cd ${RESDIR}/

if [ ! -f "${subname}_correlation/knownResults.txt" ]
then
    findMotifsGenome.pl \\
        ${ROOT}/${file} \\
        /tscc/projects/ps-epigen/users/biy022/scmethylhic/genome/genome_hg38/hg38.fa \\
        ${subname}_correlation \\
        -size given -mask -p 12
fi
EOF
done



# Use expr filtered corr atac
# The cutoff is |corr| > 0.5 and padjusted < 0.05
# RC axis
RESDIR=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/rostral_caudal_spearman_expr_filtered/homer
ROOT=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/rostral_caudal_spearman_expr_filtered




