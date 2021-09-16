#! /bin/bash
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs
#$ -l h_rt=23:00:00,h_data=20G,h_vmem=30G
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/NoCpGRepeats_sites.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/NoCpGRepeats_sites.err.txt
#$ -m abe
#$ -t 1-21

# Extract genomic regions that are at least 10 Kb apart from exons and are not within repetitive regions or cpg islands, to new vcf files per chromosome or scaffold
# Adapted from Annabel Beichman's script (for analyzing exomes) by Sergio Nigenda to analyze whole genome data.
# Usage: qsub Extract_noCpG_noRepetitive.sh [reference species]


########## Setting environement

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -oe pipefail


########## Set variables, files and directories

REF=${1}

HOMEDIR=/u/project/rwayne/jarobins/vaquita
WORKDIR=${HOMEDIR}/Phsi_neutral_sfs/Neutral_Regions
VCFDIR=${HOMEDIR}/vaquita_20/vcfs	
OUTDIR=${WORKDIR}/neutralVCFs/noCpGRepeats
REFDIR=${HOMEDIR}/reference_genomes/vaquita_genome/mPhoSin1.pri.cur.20190723_rename
REFERENCE=${REFDIR}/mPhoSin1.pri.cur.20190723_rename.fasta
CPG_REPEATS_HDENSNPS=${REFDIR}/CpG_repeats_highSNPdensity_all_merged.bed
IDX=$(printf %01d ${SGE_TASK_ID})
TenKb=${WORKDIR}/DistanceFromExons/all_HQCoords_min10kb_DistFromCDR.0based.bed

mkdir -p ${WORKDIR}
mkdir -p ${OUTDIR}

##### Logging

cd ${WORKDIR}/neutralVCFs

mkdir -p ./logs
mkdir -p ./temp

# echo the input
echo "[$(date "+%Y-%m-%d %T")] Start extracting not repetitive regions neither cpg islands for ${REF} ${SGE_TASK_ID} JOB_ID: ${JOB_ID}"
echo "The qsub input"
echo "${REF} ${SGE_TASK_ID}"

PROGRESSLOG=./logs/Extract_No_repeats_cpg_${REF}_${IDX}_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


########## Obtain vcf files per chromosome or scaffold that do not contain repeat regions or cpg islands and are at least 10 Kb apart from exons

echo -e "[$(date "+%Y-%m-%d %T")]  Extracting neutral regions with GATK SelecVariants... " >> ${PROGRESSLOG}
LOG=./logs/Step1_ExtractNeutralSites_${REF}_SelectVariants_${IDX}.log
date "+%Y-%m-%d %T" > ${LOG}

gatk3 -Xmx15g -Djava.io.tmpdir=./temp -T SelectVariants \
-R ${REFERENCE} \
--variant ${VCFDIR}/vaquita_20_chr${IDX}_simple_PASS.vcf.gz \
-XL ${CPG_REPEATS_HDENSNPS} \
-L ${TenKb} \
-o ${OUTDIR}/nocpg_repeats_SFS_${IDX}.vcf.gz &>> ${LOG}

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

echo "[$(date "+%Y-%m-%d %T")] Done extracting no repeats regions or cpg islands for ${REF} ${SGE_TASK_ID} Job ID: ${JOB_ID}"

conda deactivate

