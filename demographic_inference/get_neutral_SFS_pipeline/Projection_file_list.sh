#! /bin/bash
#$ -cwd
#$ -l h_rt=00:20:00,h_data=15G
#$ -N easySFSProjection_join
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/SFS_projection_names.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/SFS_projection_names.err.txt
#$ -m abe
#$ -M snigenda


# This script concentrate files names of SFS projection per chromosomes
# Author: Paulina Nunez (pnunez@lcg.unam.mx)
# Usage: qsub SFS_projection_join.sh

source /u/project/rwayne/software/finwhale/miniconda2/etc/profile.d/conda.sh
conda activate gentools

set -o pipefail

#Define directories ---------------------------------

# workdir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/SFS_projection_fsc
# SFSdir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/SFS_projection_Vaquita

workdir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/SFS_projection_fsc/projection_24_indiv
SFSdir=/u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/SFS_neutral/Neutral/SFS_projection_Vaquita/projection_24_indiv

mkdir -p ${workdir}
mkdir -p ${SFSdir}
pop=Vaquita

#Main ---------------

 
cd ${workdir}

x=$(ls ${SFSdir}) ; printf "%s\n" "$x" > directories.txt
while read line; do
	echo -e ${SFSdir}"/"$line"/fastsimcoal2/"${pop}"_MAFpop0.obs" >> ${pop}"_SFS_projection_files.txt"
done < directories.txt

conda deactivate
