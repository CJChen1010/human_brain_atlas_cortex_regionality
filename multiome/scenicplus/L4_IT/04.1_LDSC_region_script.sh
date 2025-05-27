### Codes for LDSC analysis regional
srun -t 8:00:00 -N 1 -c 16 -A csd772 -p condo -q condo --mem 100G --pty bash
conda activate macs3

root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT
cd ${root_dir}
cd regional_atac_peaks

for file in ../regional_atac_fragments/*bed
do
    region_name=$(basename $file .bed)
    printf "Processing %s\n" ${region_name}
    macs3 callpeak --treatment $file \
        --extsize 150 --shift -75 --nomodel -g hs \
        --name ${region_name} -q 0.1 --call-summits -f BED
done

for file in *_summits.bed
do
    region_name=$(basename $file _summits.bed)
    printf "Processing %s\n" ${region_name}
    bedtools intersect \
        -a ../L4_IT_atac_regions.bed \
        -b $file \
        -wa -u >${region_name}_intersect_baseset.bed
done

for file in *_intersect_baseset.bed
do
    prefix=$(basename $file _intersect_baseset.bed)
    echo $prefix
    /tscc/projects/ps-renlab/nzemke/software/liftOver \
        -minMatch=0.99 \
        $file \
        /tscc/projects/ps-renlab/nzemke/software/liftover_annotations/hg38Tohg19.over.chain.gz \
        ./${prefix}_intersect_baseset_hg19.bed \
        ./${prefix}_intersect_baseset_unmapped.bed
done


cd ${root_dir}
plink=/tscc/projects/ps-renlab/yangli/resource/1000G_EUR_Phase3_plink
ldsc_home=~/softwares/ldsc/
conda activate ldsc

mkdir -p LDSC_region/sbatch_logs LDSC_region/sbatch_scripts LDSC_region/ld_score

function loadavg {
    while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt 50 ]
    do
        sleep 120
        date
    done
}

for file in regional_atac_peaks/*_intersect_baseset_hg19.bed
do
    region=$(basename ${file} _intersect_baseset_hg19.bed)
    echo $region
    mkdir -p LDSC_region/ld_score/${region}
    for j in {1..22}
    do
        python ${ldsc_home}/make_annot.py \
            --bed-file $file \
            --bimfile ${plink}/1000G.EUR.QC.$j.bim \
            --annot-file LDSC_region/ld_score/${region}/${region}.$j.annot.gz &
        sleep 3
        loadavg
    done
done

/tscc/projects/ps-renlab/nzemke/software/liftOver \
    -minMatch=0.99 \
    L4_IT_atac_regions.bed \
    /tscc/projects/ps-renlab/nzemke/software/liftover_annotations/hg38Tohg19.over.chain.gz \
    ./L4_IT_atac_regions_hg19.bed \
    ./L4_IT_atac_regions_unmapped.bed

mkdir -p LDSC_region/ld_score/Baseset
for j in {1..22}
do
    python ${ldsc_home}/make_annot.py \
        --bed-file L4_IT_atac_regions_hg19.bed \
        --bimfile ${plink}/1000G.EUR.QC.$j.bim \
        --annot-file LDSC_region/ld_score/Baseset/Baseset.$j.annot.gz &
    sleep 3
    loadavg
done


root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT
ldsc_home=~/softwares/ldsc/
DIR=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/LDSC_region
resources=/tscc/projects/ps-renlab/yangli/resource
union_peaks=L4_IT_atac_regions_hg19.bed

cd ${root_dir}

for file in regional_atac_peaks/*_intersect_baseset_hg19.bed
do
    region=$(basename ${file} _intersect_baseset_hg19.bed)
    echo ${region}
    cat >${DIR}/sbatch_scripts/${region}.sbatch <<EOF
#! /bin/bash
#SBATCH -J ${region}.ldscore
#SBATCH -A csd772
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 2
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH --mem 100G
#SBATCH -o ${DIR}/sbatch_logs/${region}.ldscore.out
#SBATCH -e ${DIR}/sbatch_logs/${region}.ldscore.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
cd $DIR
conda activate ldsc

for i in {1..22}
do
    python ${ldsc_home}/ldsc.py \\
        --l2 \\
        --bfile ${resources}/1000G_EUR_Phase3_plink/1000G.EUR.QC.\$i \\
        --ld-wind-cm 1 \\
        --annot ld_score/${region}/${region}.\$i.annot.gz \\
        --thin-annot \\
        --out ld_score/${region}/${region}.\$i \\
        --print-snps ${resources}/1000G_Phase3_hapmap3/1000G_Phase3_hapmap3_print_snps.\$i.snp
done

conda deactivate
EOF
done

cat >$DIR/sbatch_scripts/Baseset.sbatch <<EOF
#! /bin/bash
#SBATCH -J Baseset.ldscore
#SBATCH -A csd772
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 2
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH --mem 100G
#SBATCH -o ${DIR}/sbatch_logs/Baseset.ldscore.out
#SBATCH -e ${DIR}/sbatch_logs/Baseset.ldscore.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
cd $DIR
conda activate ldsc

for i in {1..22}
do
    python ${ldsc_home}/ldsc.py \\
        --l2 \\
        --bfile $resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.\$i \\
        --ld-wind-cm 1 \\
        --annot ld_score/Baseset/Baseset.\$i.annot.gz \\
        --thin-annot \\
        --out ld_score/Baseset/Baseset.\$i \\
        --print-snps $resources/1000G_Phase3_hapmap3/1000G_Phase3_hapmap3_print_snps.\$i.snp
done

conda deactivate
EOF

for file in LDSC_region/sbatch_scripts/*sbatch; do sbatch $file; done


root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT
cd ${root_dir}

mkdir -p LDSC_region/results
mkdir -p LDSC_region/sbatch_test_scripts
mkdir -p LDSC_region/sbatch_test_logs

for file in regional_atac_peaks/*_intersect_baseset_hg19.bed
do
    region=$(basename ${file} _intersect_baseset_hg19.bed)
    printf "${region}\tld_score/${region}/${region}.,ld_score/Baseset/Baseset.\n"
done >LDSC_region/regions.ldcts

cts_name=regions
DIR=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/LDSC_region
sumstats=/tscc/projects/ps-renlab/yangli/resource/GWAStraits
resources=/tscc/projects/ps-renlab/yangli/resource
ldsc_home=~/softwares/ldsc

for file in ${sumstats}/*sumstats.gz
do
    i=$(basename $file .sumstats.gz)
    echo $i
    cat >LDSC_region/sbatch_test_scripts/$i.$cts_name.sbatch <<EOF
#! /bin/bash
#SBATCH -J ${i}.${cts_name}.ldsc
#SBATCH -A csd772
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 2
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH --mem 100G
#SBATCH -o $DIR/sbatch_test_logs/$i.$cts_name.ldsc.out
#SBATCH -e $DIR/sbatch_test_logs/$i.$cts_name.ldsc.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
cd $DIR
conda activate ldsc

python ${ldsc_home}/ldsc.py \\
    --h2-cts $sumstats/$i.sumstats.gz \\
    --ref-ld-chr $resources/1000G_EUR_Phase3_baseline/baseline. \\
    --out results/${i}.$cts_name \\
    --ref-ld-chr-cts $cts_name.ldcts \\
    --w-ld-chr $resources/weights_hm3_no_hla/weights.
conda deactivate
EOF

done

for file in LDSC_region/sbatch_test_scripts/*sbatch; do sbatch $file; done








### Codes for LDSC analysis correlation
root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT

cd ${root_dir}
plink=/tscc/projects/ps-renlab/yangli/resource/1000G_EUR_Phase3_plink
ldsc_home=~/softwares/ldsc/
conda activate ldsc

mkdir -p LDSC_spearman_corr/sbatch_logs LDSC_spearman_corr/sbatch_scripts LDSC_spearman_corr/ld_score

cd rostral_caudal_spearman
for file in *_peaks_bh_5e-1_5e-2_donor_replicate.bed
do
    prefix=$(basename $file _peaks_bh_5e-1_5e-2_donor_replicate.bed)
    echo $prefix
    /tscc/projects/ps-renlab/nzemke/software/liftOver \
        -minMatch=0.99 \
        $file \
        /tscc/projects/ps-renlab/nzemke/software/liftover_annotations/hg38Tohg19.over.chain.gz \
        ./${prefix}_peaks_bh_5e-1_5e-2_donor_replicate_hg19.bed \
        ./${prefix}_peaks_bh_5e-1_5e-2_donor_replicate_unmapped.bed
done
cd ${root_dir}

for file in rostral_caudal_spearman/*_peaks_bh_5e-1_5e-2_donor_replicate_hg19.bed
do
    region=$(basename ${file} _peaks_bh_5e-1_5e-2_donor_replicate_hg19.bed)
    echo $region
    mkdir -p LDSC_spearman_corr/ld_score/${region}
    for j in {1..22}
    do
        python ${ldsc_home}/make_annot.py \
            --bed-file $file \
            --bimfile ${plink}/1000G.EUR.QC.$j.bim \
            --annot-file LDSC_spearman_corr/ld_score/${region}/${region}.$j.annot.gz &
        sleep 3
        loadavg
    done
done

root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT
ldsc_home=~/softwares/ldsc/
DIR=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/LDSC_spearman_corr
resources=/tscc/projects/ps-renlab/yangli/resource
union_peaks=L4_IT_atac_regions_hg19.bed

cd ${root_dir}

for file in rostral_caudal_spearman/*_peaks_bh_5e-1_5e-2_donor_replicate_hg19.bed
do
    region=$(basename ${file} _peaks_bh_5e-1_5e-2_donor_replicate_hg19.bed)
    echo ${region}
    cat >${DIR}/sbatch_scripts/${region}.sbatch <<EOF
#! /bin/bash
#SBATCH -J ${region}.ldscore
#SBATCH -A csd772
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 2
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH --mem 100G
#SBATCH -o ${DIR}/sbatch_logs/${region}.ldscore.out
#SBATCH -e ${DIR}/sbatch_logs/${region}.ldscore.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
cd $DIR
conda activate ldsc

for i in {1..22}
do
    python ${ldsc_home}/ldsc.py \\
        --l2 \\
        --bfile ${resources}/1000G_EUR_Phase3_plink/1000G.EUR.QC.\$i \\
        --ld-wind-cm 1 \\
        --annot ld_score/${region}/${region}.\$i.annot.gz \\
        --thin-annot \\
        --out ld_score/${region}/${region}.\$i \\
        --print-snps ${resources}/1000G_Phase3_hapmap3/1000G_Phase3_hapmap3_print_snps.\$i.snp
done

conda deactivate
EOF
done

root_dir=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT
cd ${root_dir}

mkdir -p LDSC_spearman_corr/results
mkdir -p LDSC_spearman_corr/sbatch_test_scripts
mkdir -p LDSC_spearman_corr/sbatch_test_logs

for file in rostral_caudal_spearman/*_peaks_bh_5e-1_5e-2_donor_replicate_hg19.bed
do
    region=$(basename ${file} _peaks_bh_5e-1_5e-2_donor_replicate_hg19.bed)
    printf "${region}\tld_score/${region}/${region}.,ld_score/Baseset/Baseset.\n"
done >LDSC_spearman_corr/regions.ldcts

cts_name=regions
DIR=/tscc/projects/ps-epigen/users/biy022/biccn/data/SNAREdata/scenicplus/L4_IT/LDSC_spearman_corr
sumstats=/tscc/projects/ps-renlab/yangli/resource/GWAStraits
resources=/tscc/projects/ps-renlab/yangli/resource
ldsc_home=~/softwares/ldsc

for file in ${sumstats}/*sumstats.gz
do
    i=$(basename $file .sumstats.gz)
    echo $i
    cat >LDSC_spearman_corr/sbatch_test_scripts/$i.$cts_name.sbatch <<EOF
#! /bin/bash
#SBATCH -J ${i}.${cts_name}.ldsc
#SBATCH -A csd772
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -N 2
#SBATCH -c 8
#SBATCH -t 8:00:00
#SBATCH --mem 100G
#SBATCH -o $DIR/sbatch_test_logs/$i.$cts_name.ldsc.out
#SBATCH -e $DIR/sbatch_test_logs/$i.$cts_name.ldsc.err
#SBATCH --mail-user biy022@health.ucsd.edu
#SBATCH --mail-type FAIL

source ~/.bashrc
cd $DIR
conda activate ldsc

python ${ldsc_home}/ldsc.py \\
    --h2-cts $sumstats/$i.sumstats.gz \\
    --ref-ld-chr $resources/1000G_EUR_Phase3_baseline/baseline. \\
    --out results/${i}.$cts_name \\
    --ref-ld-chr-cts $cts_name.ldcts \\
    --w-ld-chr $resources/weights_hm3_no_hla/weights.
conda deactivate
EOF

done

for file in LDSC_spearman_corr/sbatch_test_scripts/*sbatch; do sbatch $file; done