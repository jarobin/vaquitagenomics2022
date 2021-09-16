#! /bin/bash
#$ -wd /u/project/rwayne/meixilin/fin_whale/analyses/SFS
#$ -l h_rt=20:00:00,h_data=4G,h_vmem=10G
#$ -o /u/project/rwayne/meixilin/fin_whale/analyses/reports/easySFS_ProjectionPreview_20200713.out.txt
#$ -e /u/project/rwayne/meixilin/fin_whale/analyses/reports/easySFS_ProjectionPreview_20200713.err.txt
#$ -m abe

# Usage: qsub easySFS_1_ProjectionPreview_onlyPASSpopHet75.sh
# This script will get and parse project preview for a given vcf file 
# @modification Mon Jul 13 21:44:56 2020
# @modification Final settings. Add per population maxHet filters at 0.75.
set -eo pipefail

###########################################################
## import packages 
source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

###########################################################
## def functions 
# select variants 
# $1 = input VCF name
# $2 = output VCF name
gatk_select_variants_onlyPASS(){
echo -e "[$(date "+%Y-%m-%d %T")] Start gatk3 onlyPASS, filter ${EXCLUDE_SAMPLE[@]} for ${1}" >> ${LOG}
gatk3 -Xmx4G -R $REFERENCE -T SelectVariants \
$(for ss in ${EXCLUDE_SAMPLE[@]}; do echo "--exclude_sample_name ${ss} "; done) \
--excludeFiltered \
--variant ${1} -o ${2} &>> ${LOG}
}
# perform the projection; parse the projection and plot the projection
# $1 = input VCF name
# $2 = output file prefix (OUTDIR previously defined)
easySFS_projection_popHet75() {
# generate easySFS preview 
# -a keep all snps 
# -v verbose 
echo -e "[$(date "+%Y-%m-%d %T")] Start easySFS, perPop maxHet 0.75 for ${1}" >> ${LOG}
python $WORKSCRIPT -i ${1} -maxHetFilter ${maxHetFilter} -p ${POPFILE} --preview -a -v > ${OUTDIR}/${2}.easySFS.projPreview.txt 2>> ${LOG}

# parse the preview (the population had to match presupplied variables)
INFILE=${OUTDIR}/${2}.easySFS.projPreview.txt
echo -e "[$(date "+%Y-%m-%d %T")] Parsing easySFS output ${INFILE}" >> ${LOG}
for pop in ${POPS[@]}
do
OUTFILE=${OUTDIR}/${pop}.${2}.easySFS.projPreview.R.format.txt
echo "projection,snps" > ${OUTFILE}
grep -A1 "$pop$" ${INFILE} | tail -1 | \
sed 's/(//g' | sed 's/)//g' | sed 's/, /,/g' |  tr '\t' '\n' >> ${OUTFILE}
done

# perform plotting (redirect stderr to log as well)
echo -e "[$(date "+%Y-%m-%d %T")] plotting ${INFILE}" >> ${LOG}
cd $OUTDIR
Rscript --vanilla ${PLOTSCRIPT} "${OUTDIR}" "${OUTDIR}/OptimalProjectionPlots/" ${2} >> ${LOG} 2>&1
}

###########################################################
## def variables 
WORKDIR=/u/project/rwayne/meixilin/fin_whale/analyses
VCFDIR=/u/project/rwayne/snigenda/finwhale/filteredvcf/all48/Minke_canonical_CDS
SCRATCHDIR=/u/scratch/m/meixilin/finwhale/filteredvcf/all48/Minke_canonical_CDS
WORKSCRIPT=${WORKDIR}/scripts/SFS/easySFS/easySFS_ab.py
PLOTSCRIPT=${WORKDIR}/scripts/SFS/easySFS/easySFS_function_determineOptimalProjection.R

# Minke whale reference genome 
REFERENCE=/u/project/rwayne/snigenda/finwhale/cetacean_genomes/minke_whale_genome/GCF_000493695.1_BalAcu1.0/GCF_000493695.1_BalAcu1.0_genomic.fasta
EXCLUDE_SAMPLE=("ENPCA01" "ENPCA09" "ENPOR12" "GOC010" "GOC080" "GOC111")
POPS=("ENP" "GOC")

OUTDIR=${WORKDIR}/SFS/all48/Minke_canonical_CDS/easySFS/onlyPASS_popHet0.75
mkdir -p ${OUTDIR}

POPFILE="${WORKDIR}/scripts/config/pop_map.txt"
maxHetFilter=0.75 # still need max het filter if using only PASS sites in filtered vcf, because of population specific het (popHet) could still exceed maxHet 

# additional arguments 
mydate=$(date "+%Y%m%d")
LOG=/u/project/rwayne/meixilin/fin_whale/analyses/reports/easySFS_ProjectionPreview_onlyPASSpopHet75_${mydate}.log

###########################################################
## main 
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; easySFS projection Final (onlyPASS, popHet=0.75)" 

cd ${OUTDIR}
gatk_select_variants_onlyPASS ${VCFDIR}/JointCalls_filterpassmiss_syn_all.vcf.gz ${VCFDIR}/JointCalls_filterpass_m6_syn_all.vcf.gz
gatk_select_variants_onlyPASS ${VCFDIR}/JointCalls_filterpassmiss_nonsyn_all.vcf.gz ${VCFDIR}/JointCalls_filterpass_m6_nonsyn_all.vcf.gz
gatk_select_variants_onlyPASS ${VCFDIR}/JointCalls_filterpassmiss_LOF_all.vcf.gz ${VCFDIR}/JointCalls_filterpass_m6_LOF_all.vcf.gz

easySFS_projection_popHet75 ${VCFDIR}/JointCalls_filterpass_m6_syn_all.vcf.gz "synonymous.allpop_onlyPASS.popHet.75"
easySFS_projection_popHet75 ${VCFDIR}/JointCalls_filterpass_m6_nonsyn_all.vcf.gz "nonsynonymous.allpop_onlyPASS.popHet.75"
easySFS_projection_popHet75 ${VCFDIR}/JointCalls_filterpass_m6_LOF_all.vcf.gz "lossoffunc.allpop_onlyPASS.popHet.75"

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" 
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${LOG}
conda deactivate 
