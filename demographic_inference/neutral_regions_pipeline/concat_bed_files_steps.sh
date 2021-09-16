#! /bin/bash
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/concat_bed_mask_files
#$ -l h_rt=04:00:00,h_data=5G,h_vmem=10G
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/concat_files.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/concat_files.err.txt
#$ -m abe

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools


# 1. concatenates bed files with CpG islands regions, repetitive regions and high density SNPs regions that will be latter removed from neutral regions.

cat /u/project/rwayne/jarobins/vaquita/reference_genomes/vaquita_genome/mPhoSin1.pri.cur.20190723_rename/mPhoSin1.pri.cur.20190723_rename_repeats_TRF_RM.bed /u/project/rwayne/jarobins/vaquita/reference_genomes/vaquita_genome/mPhoSin1.pri.cur.20190723_rename/mPhoSin1.pri.cur.20190723_rename_cpgIslands.bed /u/project/rwayne/jarobins/vaquita/reference_genomes/vaquita_genome/mPhoSin1.pri.cur.20190723_rename/vaquita_20_autos_PASS_10kbwins_1kbstep_highSNPdensity.bed > CpG_repeats_highSNPdensity_all_unsorted.bed

# 2. sorts the concatenated file containing all the regions explained above.

sort -k1,1 -k2,2n CpG_repeats_highSNPdensity_all_unsorted.bed > CpG_repeats_highSNPdensity_all_sorted.bed

# 3. Mergew the regions using bedtools

bedtools merge -i CpG_repeats_highSNPdensity_all_sorted.bed > CpG_repeats_highSNPdensity_all_merged.bed

conda deactivate
