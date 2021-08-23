#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,mfree=32G
#$ -N easySFSProjection
#$ -o /net/harris/vol1/home/beichman/vaquita/reports/SFS
#$ -e /net/harris/vol1/home/beichman/vaquita/reports/SFS
#$ -m abe
#$ -M ab08028
module load modules modules-init modules-gs # initialize modules 
module load python/2.7.13 # aha then have to load modules that only become avail once python loaded (weird)
module load numpy/1.15.1 numpy/1.16.6 matplotlib/2.2.3 pandas/0.24.2  scipy/1.2.3 # need to load in THIS ORDER (prereqs)

bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)



todaysdate=`date +%Y%m%d`

maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)

gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/VaquitaDFE/
scriptdir=$gitdir/generate_sfs
vcfdir=/net/harris/vol1/home/beichman/vaquita/vcfs/cds_vcfs

popFile=$scriptdir/popFile.Vaquita.noRelatives.txt

easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  

perPopHetFilter=0.75 # excess het filtering

populations="GoC"
projections="24" # haploids



################################################################################
##################################### coding sites #############################
################################################################################


categories="synonymous missense"
for category in $categories
do
outdir=/net/harris/vol1/home/beichman/vaquita/sfs/easySFS/$category/projection-${todaysdate}-hetFilter-${perPopHetFilter}
mkdir -p $outdir
### before running the preview, gunzip the vcfs!!! 
vcf=${category}SNPsOnly_ExcludeRelatives_vaquita_20_simple_PASS_autos_variants.vcf # must be gunzipped

$easySFS -i $vcfdir/${vcf} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $perPopHetFilter
done 

########## get counts of monomorphic cds sites -- still have to think about how to scale these for syn/mis ############
outdir=/net/harris/vol1/home/beichman/vaquita/sfs/easySFS/allSitesCDS/projection-${todaysdate}-hetFilter-${perPopHetFilter}
mkdir -p $outdir
cdsVCF=cdsAllSites_ExcludeRelatives_vaquita_20_simple_PASS_autos.vcf.gz
python $scriptdir/getMonomorphicProjectionCounts.1D.2DSFS.py --vcf $vcfdir/${cdsVCF} --popMap $popFile --proj $projections --popIDs GoC --outdir $outdir

