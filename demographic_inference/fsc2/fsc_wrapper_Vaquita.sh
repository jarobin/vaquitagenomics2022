#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=5G
#$ -N fscWrapper
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/fscwrapper.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/fscwrapper.err.txt
#$ -M snigenda
#$ -t 1-100

# This is a wrapper that will run 100 fastsimcoal iterations for each population for any list of models
# Author: Annabel Beichman , modified by Paulina Nunez (pnunez@lcg.unam.mx) and Sergio Nigenda
# Usage: qsub fsc_wrapper_Vaquita.sh

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

sleep $((RANDOM % 120))

set -o pipefail


# Defined directories -------------

wd=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/Demography/fastsimcoal
infDir=$wd/fastsimcoal_inference # specify where inference is happening
genericDir=$wd/ModelsFiles # location of generic FSC models
sfsDir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/SFS_projection_fsc

mkdir -p ${infDir}

# Programs -------------------------

fsc=/u/project/rwayne/snigenda/programs/fsc26_linux64/fsc26

# Parameters -----------------------

models='1D.4Epoch.Vaquita'
pop="Vaquita"
cores=8 
ss="26"
version="test1"

########################################  MAIN  #############################################


for model in $models
do
		
	echo "starting $pop, $model"
	header=${model}

	# Copy generic files into directory and update 
	outdir=$infDir/$pop/$model/$version/run_${SGE_TASK_ID} 
	mkdir -p $outdir 

	cp $genericDir/$model.tpl $genericDir/$model.est $outdir # copy .est and .tpl files to outdir
	sed -i'' "s/SAMPLE_SIZE/$ss/g" $outdir/$model.tpl # sub in the sample size; note you need double quotes for variable to be expanded
	
	# Get sfs into inference directory and rename to match .est and .tpl files 
	cp $sfsDir/${pop}_MAFpop0.obs $outdir/${header}_MAFpop0.obs # copy your sfs into the directory where you'll be doing the fsc inference 
	cd $outdir
	$fsc -t ${header}.tpl -n 500000 -m -e ${header}.est -M -L 60 -c${cores} -q


done

