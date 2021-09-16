#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=20G
#$ -N easySFSProjetion_monomorphic
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/SFS_monomorphic.err.txt
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/SFS_monomorphic.out.txt
#$ -t 17-18
#$ -m abe
#$ -M snigenda


# This script calculates monomorphic site.
# Author: Paulina Nunez (pnunez@lcg.unam.mx), modified by Sergio Nigenda
# Usage: qsub SFS_count_monomorphic_sites.sh

source /u/local/Modules/default/init/modules.sh
module load python/2.7

set -o pipefail

#Define directories ---------------------------------

homedir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs
workdir=${homedir}/get_neutral_SFS_pipeline
vcfdir=${homedir}/Neutral_Regions/neutralVCFs
#outdir=${homedir}/SFS_neutral/Neutral/SFS_projection_with_monomorphics
outdir=${homedir}/SFS_neutral/Neutral/SFS_projection_with_monomorphics/projection_24_indiv
monoscript=${workdir}/getMonomorphicProjectionCounts.1D.2DSFS.py
IDX=$(printf %01d ${SGE_TASK_ID})


#Main ---------------

popfile=${workdir}/pop_map_filtered.txt
# projections="26"
projections="24"

mkdir -p ${outdir}
cd ${outdir}

########## get counts of monomorphic sites to add to the SFSes ############

python $monoscript --vcf ${vcfdir}/Neutral_sites_SFS_${IDX}.vcf.gz --popMap ${popfile} --proj ${projections} --popIDs Vaquita --outdir ${outdir} --outPREFIX ${IDX}

