#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,mfree=32G,highp
#$ -N easySFSPreview
#$ -o /net/harris/vol1/home/beichman/vaquita/reports/SFS
#$ -e /net/harris/vol1/home/beichman/vaquita/reports/SFS
#$ -m abe
#$ -M ab08028

####### Easy SFS
# https://github.com/isaacovercast/easySFS
# install:
# git clone git@github.com:isaacovercast/easySFS.git
# cd easySFS
# chmod +x *.py
# easySFS.py
module load modules modules-init modules-gs # initialize modules 
module load python/2.7.13 # aha then have to load modules that only become avail once python loaded (weird)
module load numpy/1.15.1 matplotlib/2.2.3 pandas/0.24.2  scipy/1.2.3 # need to load in THIS ORDER (prereqs)

# pip install dadi # dadi-2.1.0.tar.gz 
# installs it to: /net/harris/vol1/home/beichman/.local/lib/python2.7/site-packages
# can see it with pip show dadi 

maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/VaquitaDFE/
scriptdir=$gitdir/generate_sfs
vcfdir=/net/harris/vol1/home/beichman/vaquita/vcfs/cds_vcfs

popFile=$scriptdir/popFile.Vaquita.noRelatives.txt

easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  # this is my modification
outdir=/net/harris/vol1/home/beichman/vaquita/sfs/easySFS/projection_preview
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# want it to just output that info and continue , not prompt yes/no.
# this vcf has all snps across all categories (cds, neutral, etc.) with 0.9 max no call frac (v. liberal)
# and has had all individuals removed that won't go into the SFS
# going to do the actual projection for each category of site


# must gunzip it first 
# gunzip gunzip missenseSNPsOnly_ExcludeRelatives_vaquita_20_simple_PASS_autos_variants.vcf.gz 
# gunzip synonymousSNPsOnly_ExcludeRelatives_vaquita_20_simple_PASS_autos_variants.vcf.gz
categories="synonymous missense"
for category in $categories
do
### before running the preview, gunzip the vcfs!!! 
vcf=${category}SNPsOnly_ExcludeRelatives_vaquita_20_simple_PASS_autos_variants.vcf # must be gunzipped

$easySFS -i $vcfdir/${vcf} -p $popFile --preview -a -v > $outdir/${category}.easySFS.projPreview.txt
done 

