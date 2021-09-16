#! /bin/bash
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/concat_bed_mask_files
#$ -l h_rt=04:00:00,h_data=5G,h_vmem=10G
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/exon_dist.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/exon_dist.err.txt
#$ -m abe

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

HOMEDIR=/u/project/rwayne/jarobins/vaquita
WORKDIR=${HOMEDIR}/Phsi_neutral_sfs/Neutral_Regions
REFDIR=${HOMEDIR}/reference_genomes/vaquita_genome/mPhoSin1.pri.cur.20190723_rename
GTF=${HOMEDIR}/Phsi_neutral_sfs/concat_bed_mask_files/GCF_008692025.1_mPhoSin1.pri_genomic_rename.gtf
EXONS=${HOMEDIR}/Phsi_neutral_sfs/concat_bed_mask_files/GCF_008692025.1_mPhoSin1.pri_genomic_rename_exons.bed
CDRS=${HOMEDIR}/Phsi_neutral_sfs/concat_bed_mask_files/GCF_008692025.1_mPhoSin1.pri_genomic_rename_allCDR.bed
hqSites=${WORKDIR}/HiQualCoords/all_HQCoords_sorted_merged.bed

mkdir -p ${WORKDIR}/DistanceFromExons

### Get exonic regions or coding regions (do only once) and make sure is sorted

awk '$3=="exon"' ${GTF} | awk '{OFS="\t";print $1,$4-1,$5,$9}' | sort -k1,1 -k2,2n > ${HOMEDIR}/Phsi_neutral_sfs/concat_bed_mask_files/GCF_008692025.1_mPhoSin1.pri_genomic_rename_exons.bed

## For vaquita we selected coding regions

# bedtools closest -d -a ${hqSites} -b ${EXONS} > ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt
bedtools closest -d -a ${hqSites} -b ${CDRS} > ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt


### Exploring and choosing different distances

## The last column (8) is the distance; we wanted the distance to be at least 10,000, and want to keep track of the distance. Collect all that are >10,000 away.

# awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min10kb_DistFromExons.0based.bed

awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' ${WORKDIR}/DistanceFromExons/all_HQCoords_DistFromCDR.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > ${WORKDIR}/DistanceFromExons/all_HQCoords_min10kb_DistFromCDR.0based.bed

## Note: 1,2,3 columns are the HQ SITES position, NOT the position of the exon. (If you mess up what is a and b in bedtools closest this would be messed up)


> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
for i in `ls ${WORKDIR}/DistanceFromExons/*DistFromExons.0based.bed`
do
echo $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromExons.txt
done

> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
for i in `ls ${WORKDIR}/DistanceFromExons/*DistFromCDR.0based.bed`
do
echo $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> ${WORKDIR}/DistanceFromExons/totalSequenceByDistanceFromCDRs.txt
done

conda deactivate

