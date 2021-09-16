#! /bin/bash
#$ -cwd
#$ -l h_rt=15:00:00,h_data=15G
#$ -N easySFSProjetion_chr
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/SFS_projection.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/SFS_projection.err.txt
#$ -t 1-21
#$ -m abe
#$ -M snigenda


# This script runs SFS projection  per chromosomes
# Author: Meixi Lin, modified by Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub SFS_projection_chr.sh

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

#Define directories ---------------------------------

homedir=/u/project/rwayne/jarobins/vaquita
workdir=${homedir}/Phsi_neutral_sfs/get_neutral_SFS_pipeline
easySFS=${workdir}/easySFS_a.py
vcfdir=${homedir}/Phsi_neutral_sfs/SFS_neutral/Neutral/SNPs
IDX=$(printf %01d ${SGE_TASK_ID})


#Main ---------------

cd ${workdir}
popfile=${workdir}/pop_map_filtered.txt
# outdir=${homedir}/Phsi_neutral_sfs/SFS_neutral/Neutral/SFS_projection_Vaquita
outdir=${homedir}/Phsi_neutral_sfs/SFS_neutral/Neutral/SFS_projection_Vaquita/projection_24_indiv

mkdir -p $outdir
cd ${outdir}

# python $easySFS -i ${vcfdir}/SNPs_neutral_for_SFS_${IDX}.vcf.gz -p ${popfile} --proj 26 -a -f -v -o SNPs_easySFS_projection_${IDX} -maxHetFilter 0.75 

python $easySFS -i ${vcfdir}/SNPs_neutral_for_SFS_${IDX}.vcf.gz -p ${popfile} --proj 24 -a -f -v -o SNPs_easySFS_projection_${IDX} -maxHetFilter 0.75

conda deactivate

