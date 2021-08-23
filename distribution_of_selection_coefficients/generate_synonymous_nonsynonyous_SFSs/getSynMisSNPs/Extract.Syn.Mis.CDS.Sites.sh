#! /bin/bash
#$ -l h_rt=50:00:00,mfree=15G
#$ -o /net/harris/vol1/home/beichman/vaquita/reports/GATK/
#$ -e /net/harris/vol1/home/beichman/vaquita/reports/GATK/
#$ -m bea
#$ -M annabel.beichman@gmail.com
#$ -N ExtractSynMis

########### extract syn mis and cds sites from JAR's bed files, and exclude relatives #####
module load modules modules-init modules-gs # initialize modules 
java18=/usr/lib/jvm/java-1.8.0/bin/java
 ## need to call java directly; is no longer a module 
GATK=/net/harris/vol1/home/beichman/bin/GATK_3.7/GenomeAnalysisTK.jar
module load htslib/1.9
module load samtools/1.9
module load picard/2.21.7
projectdir=/net/harris/vol1/home/beichman/vaquita
wd=$projectdir/vcfs/ 
mkdir -p $wd/cds_vcfs

refdir=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename
REFERENCE=$refdir/mPhoSin1.pri.cur.20190723_rename.fasta
# once: need to generate .fai and .dict files: 
#samtools faidx $refdir/$REFERENCE   # generates .fai
# and need GATK dict file:
#picard=/net/gs/vol3/software/modules-sw/picard/2.23.3/Linux/CentOS7/x86_64/picard.jar

##$java18 -jar -Xmx1g $picard CreateSequenceDictionary -R $REFERENCE # generates .dict
# autosomes only:
# note the "20" refers to the sample size (NOT the chromosome!). this is across all chrs.
snpVCF=vaquita_20_simple_PASS_autos_variants.vcf.gz # this has SNPs only 
allSitesVCF=vaquita_20_simple_PASS_autos.vcf.gz # this has all sites 

# bed files:
bedDir=/net/harris/vol1/home/beichman/vaquita/coordinate_bed_files/variant_effects
#cdsBED=vaquita_20_simple_PASS_autos_variants.vcf.gz_CDS.bed # DON'T USE THIS! THIS IS JUST SNPS. I WANT ALL CDS(var+invar) COORDS FROM GFF FILE:
cdsBED=AllCDS.Coords.FromGTF.0Based.bed
synBED=vaquita_20_simple_PASS_autos_variants.vcf.gz_SYN.bed
misBED=vaquita_20_simple_PASS_autos_variants.vcf.gz_MIS.bed

### 5 relatives to exclude:
# z0001663, z0004380, z0004393, z0004394, z0185383 # exclude these due to high relatedness
############## synonymous snps ############
echo "starting synonymous"
$java18 -jar -Xmx8g -Djava.io.tmpdir=$TMPDIR $GATK -T SelectVariants \
-R ${REFERENCE} \
-L $bedDir/$synBED \
-V $wd/$snpVCF \
-o $wd/cds_vcfs/synonymousSNPsOnly_ExcludeRelatives_${snpVCF} \
--exclude_sample_name z0001663 \
--exclude_sample_name z0004380 \
--exclude_sample_name z0004393 \
--exclude_sample_name z0004394 \
--exclude_sample_name z0185383


########### missense snps #####################
echo "starting missense"
$java18 -jar -Xmx8g -Djava.io.tmpdir=$TMPDIR $GATK -T SelectVariants \
-R ${REFERENCE} \
-L $bedDir/$misBED \
-V $wd/$snpVCF \
-o $wd/cds_vcfs/missenseSNPsOnly_ExcludeRelatives_${snpVCF} \
--exclude_sample_name z0001663 \
--exclude_sample_name z0004380 \
--exclude_sample_name z0004393 \
--exclude_sample_name z0004394 \
--exclude_sample_name z0185383

############ coding (CDS) all sites ###########
echo "starting CDS"
$java18 -jar -Xmx8g -Djava.io.tmpdir=$TMPDIR $GATK -T SelectVariants \
-R ${REFERENCE} \
-L $bedDir/$cdsBED \
-V $wd/$allSitesVCF \
-o $wd/cds_vcfs/cdsAllSites_ExcludeRelatives_${allSitesVCF} \
--exclude_sample_name z0001663 \
--exclude_sample_name z0004380 \
--exclude_sample_name z0004393 \
--exclude_sample_name z0004394 \
--exclude_sample_name z0185383
