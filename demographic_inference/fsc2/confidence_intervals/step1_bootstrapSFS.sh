#!/bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=10G
#$ -N SimulateSFS
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/step1_bootstrapSFS.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/step1_bootstrapSFS.err.txt
#$ -m abe
#$ -M snigenda

# @version 		v0
# @usage		qsub step1_param_bootstrap.sh <model> <submodel> <population> <prefix> # For 1D models this usage line will change depending on the model and the best run for such model. For example qsub step1_bootstrapSFS.sh 1DModels 1D.1Epoch Vaquita r61 
# @description	Generate bootstrapped SFS from given parameter estimations
# Author: Meixi Lin (meixilin@ucla.edu), modified by Sergio Nigenda
# Date: Sun Feb 21 23:59:39 2021

# Example of variable definitions and how they are taken from the usage or submitting line.
# MODEL='1DModels' # The general model, it could also be "2DModels" if you are using two population models
# SUBMODEL='1D.1Epoch' # The specific model from which the CI will be calculated (e.g. 1D.1Epoch, 1D.2Epoch, 1D.3Epoch)
# POP='Vaquita' # Population name
# PREFIX='r61' # This is the run that had the best log-likelihood during parameter estimation for which the CIs will be inferred. The run number with the best LL varies from model to model.

############################################################
## import packages
sleep $((RANDOM % 120))

set -eo pipefail

############################################################
## def variables
MODEL=${1}
SUBMODEL=${2}
POP=${3}
PREFIX=${4}
BOOTPREFIX=${SUBMODEL}.${POP}_${PREFIX}
NSIM=100 # Number of simulation that will be performed
TODAY=$(date "+%Y%m%d")

HOMEDIR=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/fastsimcoal
WORKDIR=${HOMEDIR}/fsc_bootstrap/${MODEL}/${POP}/${SUBMODEL}/${BOOTPREFIX}_bootSFS
mkdir -p ${WORKDIR}
LOG=${HOMEDIR}/fsc_bootstrap/step1_bootstrapSFS_${BOOTPREFIX}_${TODAY}.log

# the bootstrap file
BOOTFILE=${HOMEDIR}/fsc_bootstrap/1DModels/boots_input_files/${BOOTPREFIX}_boot.par
# COMMITID=$(git --git-dir="${HOMEDIR}/scripts/.git" --work-tree="${HOMEDIR}/scripts" rev-parse master)

# the software
FSC26=/u/project/rwayne/software/finwhale/fastsimcoal/fsc26_linux64/fsc26

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Using ${BOOTFILE}; Outputting to ${WORKDIR}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}; Using ${BOOTFILE}; Outputting to ${WORKDIR}" > ${LOG}

cd ${WORKDIR}

rsync -a ${BOOTFILE} ./

# -n  --numsims 1000      : number of simulations to perform. Also applies for parameter estimation
# -j  --jobs              : output one simulated or bootstrapped SFS per file
#                            in a separate directory for easier analysis
#                            (requires -d or -m and -s0 options)
# -m  --msfs              : computes minor site frequency spectrum
# -s  --dnatosnp 2000     : output DNA as SNP data, and specify maximum no. of SNPs
#                            to output (use 0 to output all SNPs).
# -x  --noarloutput       : does not generate Arlequin output
# -I  --inf               : generates DNA mutations according to an
#                            infinite site (IS) mutation model
# -q  --quiet             : minimal message output to console

${FSC26} -i ${BOOTPREFIX}_boot.par -n ${NSIM} -j -m -s0 -x -I -q &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" >> ${LOG}


