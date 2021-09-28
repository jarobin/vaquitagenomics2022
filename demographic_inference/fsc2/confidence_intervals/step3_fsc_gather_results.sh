#!/bin/bash
#$ -l h_data=5G,h_rt=10:00:00
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/fastsimcoal
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/step3_fsc_gather_results.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/step3_fsc_gather_results.err.txt
#$ -m abe

# @version      v0
# @script       bash step3_fsc_gather_results.sh <model> <submodel> <pop> <prefix>
# @description  Gather results from replicated runs
# Author: Meixi Lin (meixilin@ucla.edu), modified by Sergio Nigenda
# Date: Mon Feb 22 08:56:10 2021

###########################################################
## import packages
set -eo pipefail

###########################################################
## def functions

###########################################################
## def variables
MODEL=${1}
SUBMODEL=${2}
POP=${3}
PREFIX=${4}
# OGFILE=${5} # The original output bestlhoods file from real data
BOOTPREFIX=${SUBMODEL}.${POP}_${PREFIX}

TODAY=$(date "+%Y%m%d")


HOMEDIR=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/fastsimcoal
SUMDIR=${HOMEDIR}/fsc_bootstrap/resultsSummaries/${SUBMODEL}
WORKDIR=${HOMEDIR}/fsc_bootstrap/${MODEL}/${POP}/${SUBMODEL}/${BOOTPREFIX}_bootFSC_n1pv
INDIR=${HOMEDIR}/fsc_bootstrap/1DModels/boots_input_files


###########################################################
## main
mkdir -p ${SUMDIR}
cd ${SUMDIR}

rsync -a ${INDIR}/${BOOTPREFIX}.bestlhoods ./

OUTFILE=${SUMDIR}/${BOOTPREFIX}_bootFSC_n1pv_SummaryT_${TODAY}.csv
FINALFILE=${SUMDIR}/${BOOTPREFIX}_bootFSC_n1pv_Summary_${TODAY}.csv
OGFILE=${BOOTPREFIX}.bestlhoods

HEADER=$(head -n1 ${WORKDIR}/boot_1/${BOOTPREFIX}_boot/${BOOTPREFIX}_boot.bestlhoods)
echo -e "runNum\t${HEADER}" | tr "\\t" "," > ${OUTFILE}

# start with the original file
if [ -f "$OGFILE" ]; then
    results=$(grep -v [A-Z] ${OGFILE})
else
    echo 'ERROR: ${OGFILE} missing!'
    exit 1
fi
echo -e "0\t$results"| tr "\\t" "," >> ${OUTFILE}


for ii in {1..100}
do
    INFILE=${WORKDIR}/boot_${ii}/${BOOTPREFIX}_boot/${BOOTPREFIX}_boot.bestlhoods
    if [ -f "$INFILE" ]; then
        results=$(grep -v [A-Z] ${INFILE})
    else
        echo 'ERROR: ${INFILE} missing!'
        exit 1
    fi
    echo -e "$ii\t$results"| tr "\\t" "," >> ${OUTFILE}
done

############################################################
# add a column with DiffLhood = MaxObsLhood - MaxEstLhood
MAXOBSCOL=$(awk 'BEGIN {FS = ","}; {print NF}' ${OUTFILE} | head -1)
MAXESTCOL=$((${MAXOBSCOL}-1))


awk -v maxobs="${MAXOBSCOL}" -v maxest="${MAXESTCOL}" 'BEGIN {FS = ","}; {print $0 "," $maxobs-$maxest}' ${OUTFILE} > ${FINALFILE}
rm ${OUTFILE}


echo -e "[$(date "+%Y-%m-%d %T")] Done"

