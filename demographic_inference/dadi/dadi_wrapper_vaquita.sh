#! /bin/bash
#$ -N dadi_1D.Models
#$ -l h_rt=22:00:00,h_data=10G,h_vmem=15G
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/dadi
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/dadi_1D.Models.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/dadi_1D.Models.err.txt
#$ -m abe
#$ -M snigenda

# @version      v0
# @usage        qsub dadi_wrapper_vaquita.sh
# @description  wrapper for the demographic models to test 
# Author: Annabel Beichman and Meixi Lin

# @modification Wed Nov 11 2020
# @modification Author: Paulina Nunez (pnunez@lcg.unam.mx) and Sergio Nigenda 
# @modification convert for use in vaquita neutral regions framework
# NOTE: In this case the out.txt and err.txt are important 

###########################################################

# import packages --------------------------- 

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -eo pipefail

## Def variables -------------------------------

generaldir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs
mu=1.08E-08 # in units of mutations bp-1 generation-1 (Dornburg et al. 2012)
hetFilter=0.75
todaysdate=`date +%Y%m%d`
workdir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/dadi
sfsdir=${generaldir}/SFS_neutral/Neutral/SFS_projection_fsc
pop="Vaquita"

# Models to test --------------------------
scripts='1D.4Epochmod.dadi.py' # this will change depending on the model that is being run (i.e "1D.1Epoch.dadi.py", "1D.2Epoch.dadi.py", "1D.1Bottleneck.dadi.py", "1D.4Epochmod.dadi.py") Note: the 1D.1Bottleneck model is our 3-epoch model

## Main _______________________________________________________

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} dadi_wrapper_vaquita.sh" 

# get total sites from total sites file that was written out as part of the easySFS scripts
L=`grep $pop $sfsdir/$pop-[0-9][0-9].totalSiteCount.L.withMonomorphic.txt | awk '{print $2}'`
for script in $scripts
do
    model=${script%.dadi.py}
    echo "starting inference for $pop for model $model"
    outdir=$workdir/$pop/$model
        
    mkdir -p $outdir
    # carry out inference with 100 replicates that start with different p0 perturbed params:
    for i in {1..100}
    do
        echo "carrying out inference $i for model $model for pop $pop" 
        # [0-9] indicates that it's a number, but not specific about proj value
        python $workdir/$script --runNum $i --pop $pop --mu $mu --L $L --sfs ${sfsdir}/${pop}-[0-9][0-9].sfs --outdir $outdir
    done

    echo "concatenating results"
    # get the header 
    grep rundate -m1 $outdir/${pop}.dadi.inference.${model}.runNum.1.output > $outdir/${pop}.dadi.inference.${model}.all.output.concat.txt
    for i in {1..100}
    do
        grep rundate -A1 $outdir/${pop}.dadi.inference.${model}.runNum.${i}.output | tail -n1 >> $outdir/${pop}.dadi.inference.${model}.all.output.concat.txt
    done
done

############################ deactivate virtualenv ###############

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${JOB_ID} Done" 
conda deactivate # deactivate virtualenv
