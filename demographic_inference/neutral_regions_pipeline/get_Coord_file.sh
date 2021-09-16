#! /bin/bash
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs
#$ -l h_rt=23:00:00,h_data=4G,h_vmem=10G
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/getNeutralCoords.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/getNeutralCoords.err.txt
#$ -m abe
#$ -t 1-21

# Gets neutral coordinates to bed file
# Adapted from Annabel Beichman's script (to analyze exomes) by Sergio Nigenda to analyze whole genome data.
# Usage: qsub get_Coord_file.sh Vaquita


########## Setting environment

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -oe pipefail

########## Set variables, files and directories

REF=$1

HOMEDIR=/u/project/rwayne/jarobins/vaquita
WORKDIR=${HOMEDIR}/Phsi_neutral_sfs/Neutral_Regions
VCFDIR=${WORKDIR}/neutralVCFs/noCpGRepeats
OUTDIR=${WORKDIR}/neutralVCFs
SCRIPTDIR=${HOMEDIR}/Phsi_neutral_sfs/neutral_regions_pipeline
NeutralCoord_SCRIPT=${SCRIPTDIR}/obtain_noCpG_noRepetitive_coordinates.py
IDX=$(printf %01d ${SGE_TASK_ID})


##### make directories were this information will be stored

mkdir -p ${WORKDIR}/repeatRegions
mkdir -p ${WORKDIR}/get_fasta


##### Logging

cd ${OUTDIR}

mkdir -p ./logs
mkdir -p ./temp

echo "[$(date "+%Y-%m-%d %T")] Start creating bed files for ${REF} ${SGE_TASK_ID} JOB_ID: ${JOB_ID}"
echo "The qsub input"
echo "${REF} ${SGE_TASK_ID}"

PROGRESSLOG=./logs/create_beds_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


########## Creates a bed file for each vcf file

echo -e "[$(date "+%Y-%m-%d %T")]  creating bed files ... " >> ${PROGRESSLOG}
LOG=./logs/Step2_Creating_bed_files_${REF}_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

python ${NeutralCoord_SCRIPT} --VCF ${VCFDIR}/nocpg_repeats_SFS_${IDX}.vcf.gz --outfile ${WORKDIR}/repeatRegions/min10kb_DistFromCDRs_noCpGIsland_noRepeat_${IDX}.bed

bedtools merge -d 10 -i ${WORKDIR}/repeatRegions/min10kb_DistFromCDRs_noCpGIsland_noRepeat_${IDX}.bed > ${WORKDIR}/repeatRegions/min10kb_DistFromCDRs_noCpGIsland_noRepeat_mergedMaxDistance10_${IDX}.bed

cat ${WORKDIR}/repeatRegions/min10kb_DistFromCDRs_noCpGIsland_noRepeat_mergedMaxDistance10_*.bed | sort -k1,1 -k2,2n > ${WORKDIR}/get_fasta/min10kb_DistFromCDRs_noCpGIsland_noRepeat_mergedMaxDistance10_sorted.bed 


exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


conda deactivate

