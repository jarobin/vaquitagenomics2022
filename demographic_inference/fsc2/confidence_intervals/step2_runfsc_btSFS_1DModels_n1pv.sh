#!/bin/bash
#$ -l h_data=5G,h_rt=20:00:00
#$ -pe shared 8
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/fastsimcoal
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/step2_runfsc_btSFS_1DModels_n1pv.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/step2_runfsc_btSFS_1DModels_n1pv.err.txt
#$ -m abe
#$ -t 1-1000

# @version 		v0
# @usage		qsub step2_runfsc_btSFS_1DModels_n1pv.sh <model> <submodel> <population> <prefix>
# @description	Estimate parameters again from bootstrapped SFS generated from the maximum likelihood parameter estimations (n1pv: Each bootstrap only estimate once, use best estimate from data as initial values)
# Author: Meixi Lin (meixilin@ucla.edu), modified by Sergio Nigenda
# Date: Sun Feb 21 23:59:39 2021

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

NSIM=1000000
NCORE=8
TODAY=$(date "+%Y%m%d")

HOMEDIR=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/fastsimcoal
INDIR=${HOMEDIR}/fsc_bootstrap/1DModels/boots_input_files
WORKDIR=${HOMEDIR}/fsc_bootstrap/${MODEL}/${POP}/${SUBMODEL}/${BOOTPREFIX}_bootFSC_n1pv/boot_${SGE_TASK_ID}
mkdir -p ${WORKDIR}
LOGDIR=${HOMEDIR}/fsc_bootstrap/logs/step2_runfsc_btSFS_1DModels_n1pv/${BOOTPREFIX}
mkdir -p ${LOGDIR}
LOG=${LOGDIR}/${BOOTPREFIX}.${SGE_TASK_ID}.n1pv_${TODAY}.log

# the bootstrap SFS file
BOOTSFS=${HOMEDIR}/fsc_bootstrap/${MODEL}/${POP}/${SUBMODEL}/${BOOTPREFIX}_bootSFS/${BOOTPREFIX}_boot/${BOOTPREFIX}_boot_${SGE_TASK_ID}/${BOOTPREFIX}_boot_MAFpop0.obs

# the software
FSC26=/u/project/rwayne/software/finwhale/fastsimcoal/fsc26_linux64/fsc26

############################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; Using ${BOOTSFS}; Outputting to ${WORKDIR}"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID}; Using ${BOOTSFS}; Outputting to ${WORKDIR}" > ${LOG}

cd ${WORKDIR}

rsync -a ${BOOTSFS} ./
rsync -a ${INDIR}/${BOOTPREFIX}_boot.est ./
rsync -a ${INDIR}/${BOOTPREFIX}_boot.tpl ./
rsync -a ${INDIR}/${BOOTPREFIX}_boot.pv ./


# -n  --numsims 1000      : number of simulations to perform. Also applies for parameter estimation
# -m  --msfs              : computes minor site frequency spectrum
# -M  --maxlhood          : perform parameter estimation by max lhood from SFS
#                           values between iterations
# -L  --numloops 20       : number of loops (ECM cycles) to perform during
#                           lhood maximization. Default is 20
# -c   --cores 1          : number of openMP threads for parameter estimation
#                           (default=1, max=numBatches, use 0 to let openMP choose optimal value)

${FSC26} -t ${BOOTPREFIX}_boot.tpl -e ${BOOTPREFIX}_boot.est --initValues ${BOOTPREFIX}_boot.pv -n ${NSIM} -m -M -L 60 -c ${NCORE} -q &>> ${LOG}
exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL"
    exit 1
fi

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done"
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID}.${SGE_TASK_ID} Done" >> ${LOG}

