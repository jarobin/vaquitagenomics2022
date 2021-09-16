#! /bin/bash
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral
#$ -l h_rt=20:00:00,h_data=15G,h_vmem=20G
#$ -N PreviewProjection
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/Preview_Neutral_Projection.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/Preview_Neutral_Projection.err.txt
#$ -t 1-21
#$ -m abe

# Usage: qsub Neutral_PreviewProjection_EasySFS_1D.sh
# This script will get and parse project preview for a given vcf file 
# @modification Mon Jul 13 21:44:56 2020
# @modification Final settings. Add per population maxHet filters at 0.75.

set -eo pipefail

###########################################################
## import packages 
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

###########################################################
## define functions
 
# select variants 
# $1 = input VCF name
# $2 = output VCF name
gatk_select_variants_onlyPASS(){
echo -e "[$(date "+%Y-%m-%d %T")] Start gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, select SNPs for ${1}" >> ${PROGRESSLOG}
gatk3 -Xmx4G -R $REFERENCE -T SelectVariants \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude_sample_name ${ss} "; done) \
--excludeFiltered \
--removeUnusedAlternates \
--excludeNonVariants \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
--variant ${1} -o ${2} &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Finish gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, select SNPs for ${1}" >> ${PROGRESSLOG}
}

# select invariant sites
# $1 = input VCF name
# $2 = output VCF name
gatk_select_invariants_onlyPASS(){
echo -e "[$(date "+%Y-%m-%d %T")] Start gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, Select Invariant sites for ${1}" >> ${PROGRESSLOG}
gatk3 -Xmx4G -R $REFERENCE -T SelectVariants \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude_sample_name ${ss} "; done) \
--excludeFiltered \
--removeUnusedAlternates
--selectTypeToInclude NO_VARIATION \
--selectTypeToExclude INDEL \
--variant ${1} -o ${2} &>> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Finishh gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]}, Select Invariant sites for ${1}" >> ${PROGRESSLOG}
}

# perform the projection; parse the projection and plot the projection
# $1 = input VCF name
# $2 = output file suffix (OUTDIRSNPS previously defined)
easySFS_projection_popHet75() {
# generate easySFS preview 
# -a keep all snps 
# -v verbose 
echo -e "[$(date "+%Y-%m-%d %T")] Start easySFS, perPop maxHet 0.75 for ${1}" >> ${PROGRESSLOG}
python ${WORKSCRIPT} -i ${1} -maxHetFilter ${maxHetFilter} -p ${POPFILE} --preview -a -v > ${OUTDIRSNPS}/PreviewProjection_${2}.txt 2>> ${LOG}

# parse the preview (the population had to match presupplied variables)
INFILE=${OUTDIRSNPS}/PreviewProjection_${2}.txt
# echo -e "[$(date "+%Y-%m-%d %T")] Parsing easySFS output ${INFILE}" >> ${PROGRESSLOG}
# for pop in ${POPS[@]}
# do
OUTFILE=${OUTDIRSNPS}/${pop}_${2}_PreviewProjection_Rformat.txt
echo "projection,snps" > ${OUTFILE}
grep -A1 "$pop$" ${INFILE} | tail -1 | \
sed 's/(//g' | sed 's/)//g' | sed 's/, /,/g' |  tr '\t' '\n' >> ${OUTFILE}
# done

echo -e "[$(date "+%Y-%m-%d %T")] Finish easySFS, perPop maxHet 0.75 for ${1}" >> ${PROGRESSLOG}

# perform plotting (redirect stderr to log as well)
# echo -e "[$(date "+%Y-%m-%d %T")] plotting ${INFILE}" >> ${PROGRESSLOG}
# cd ${OUTDIRSNPS}
# Rscript --vanilla ${PLOTSCRIPT} "${OUTDIRSNPS}" "${OUTDIRSNPS}/OptimalProjectionPlots/" ${2} >> ${LOG} 2>&1
}

###########################################################

## define variables 

HOMEDIR=/u/project/rwayne/jarobins/vaquita
WORKDIR=${HOMEDIR}/Phsi_neutral_sfs
VCFDIR=${WORKDIR}/Neutral_Regions/neutralVCFs
WORKSCRIPT=${WORKDIR}/get_neutral_SFS_pipeline/easySFS_a.py
PLOTSCRIPT=${WORKDIR}/get_neutral_SFS_pipeline/easySFS_function_determineOptimalProjection.R
IDX=$(printf %01d ${SGE_TASK_ID})
pop=Vaquita

# Vaquita reference genome
REFDIR=${HOMEDIR}/reference_genomes/vaquita_genome/mPhoSin1.pri.cur.20190723_rename
REFERENCE=${REFDIR}/mPhoSin1.pri.cur.20190723_rename.fasta
EXCLUDE_SAMPLE=("z0001663" "z0004380" "z0004393" "z0004394" "z0185383")
# POPS=("Vaquita")


OUTDIRSNPS=${WORKDIR}/SFS_neutral/Neutral/SNPs
OUTDIRINVS=${WORKDIR}/SFS_neutral/Neutral/Invariants
mkdir -p ${OUTDIRSNPS}
mkdir -p ${OUTDIRINVS}
mkdir -p ${WORKDIR}/SFS_neutral/Neutral/logs

POPFILE=${WORKDIR}/get_neutral_SFS_pipeline/pop_map.txt
maxHetFilter=0.75 # still need max het filter if using only PASS sites in filtered vcf, because of population specific het (popHet) could still exceed maxHet 

# additional arguments 
mydate=$(date "+%Y%m%d")
LOG=${WORKDIR}/SFS_neutral/Neutral/logs/easySFS_ProjectionPreview_onlyPASSpopHet75_${mydate}_${IDX}.log
PROGRESSLOG=${WORKDIR}/SFS_neutral/Neutral/logs/ProjectionPreview_${mydate}_${IDX}_progress.log
###########################################################
## main 
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; easySFS projection Final (onlyPASS, popHet=0.75)" 

cd ${OUTDIRSNPS}
gatk_select_variants_onlyPASS ${VCFDIR}/Neutral_sites_SFS_${IDX}.vcf.gz ${OUTDIRSNPS}/SNPs_neutral_for_SFS_${IDX}.vcf.gz 

easySFS_projection_popHet75 ${OUTDIRSNPS}/SNPs_neutral_for_SFS_${IDX}.vcf.gz "neutral_SNPS_${IDX}"

cd ${OUTDIRINVS}
gatk_select_invariants_onlyPASS ${VCFDIR}/Neutral_sites_SFS_${IDX}.vcf.gz ${OUTDIRINVS}/Invariants_neutral_for_SFS_${IDX}.vcf.gz

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} ${IDX} Done" 
echo -e "[$(date "+%Y-%m-%d %T")] job for ${JOB_ID} ${IDX} Done" >> ${PROGRESSLOG}

conda deactivate 
