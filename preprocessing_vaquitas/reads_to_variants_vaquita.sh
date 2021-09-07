### Vaquita reads-to-variants pipeline (example commands)

# Workflow
# - Per read set:
#   - Trim reads with BBDUK (as needed)
#   - FastqToSam
#   - MarkIlluminaAdapters
#   - Alignment
#   - Qualimap
# - Per sample:
#   - MarkDuplicates
#   - Indel realignment: RealignerTargetCreator
#   - Indel realignment: IndelRealigner
#   - Base quality score recalibration (BQSR):
#     - Generate raw variants for BQSR: bcftools mpileup/call/filter
#     - Concatenate VCF files
#     - Generate recalibration table: BaseRecalibrator
#     - Generate recalibrated bam file: PrintReads
#     - Generate post-recalibration table: BaseRecalibrator
#     - Generate recalibration plots: AnalyzeCovariates
#     - Qualimap and DepthOfCoverage
#   - HaplotypeCaller
# - All samples jointly:
#   - GenotypeGVCFs
#   - Trim unused alternate alleles: SelectVariants
#   - VariantAnnotator
#   - SnpEff
#   - SIFT
#   - Mask repeats (identified with RepeatMasker + Tandem Repeats Finder)
#   - Custom filtering
#   - Mask regions with excess SNP density
#   - Simplify, extract SNPs
#   - Concatenate VCF files


################################################################################
### Set variables

export REFERENCE=~/project/vaquita/reference/mPhoSin1.pri/mPhoSin1.pri.fasta
export REPEATMASK=~/project/vaquita/reference/mPhoSin1.pri/mPhoSin1.pri.fasta_repeats_TRF_RM.bed

export NAME=z0184984				# Sample name
export RGID=z0184984_A				# Read group ID, unique to each read set
export RGLB=MK2010					# Library type
export RGPU=HCVTCCCXY.7.z0184984	# Flowcell ID, lane, sample name
export RGCN=HudsonAlpha				# Sequencing center


################################################################################
### PER-READ SET PROCESSING

### Trim reads with BBDUK (as needed)

export TRIM=4
export FASTQ=z0184984_Psin_RunHSX4972.3_n0009047_R1.fastq.gz

bbduk.sh -Xmx8G threads=1 ftl=${TRIM} ordered in=${FASTQ} out=${FASTQ%.fastq.gz}_trim.fastq.gz


### FastqToSam

export FQ1=z0184984_Psin_RunHSX4972.3_n0009047_R1_trim.fastq.gz
export FQ2=z0184984_Psin_RunHSX4972.3_n0009067_R2_trim.fastq.gz

picard -Xmx40G FastqToSam \
FASTQ=${FQ1} \
FASTQ2=${FQ2} \
OUTPUT=${RGID}_FastqToSam.bam \
READ_GROUP_NAME=${RGID} \
SAMPLE_NAME=${NAME} \
LIBRARY_NAME=${RGLB} \
PLATFORM_UNIT=${RGPU} \
SEQUENCING_CENTER=${RGCN} \
PLATFORM=illumina \
TMP_DIR=./temp


### MarkIlluminaAdapters

picard -Xmx40G MarkIlluminaAdapters \
INPUT=${RGID}_FastqToSam.bam \
OUTPUT=${RGID}_MarkIlluminaAdapters.bam \
METRICS=${RGID}_MarkIlluminaAdapters_metrics.txt \
TMP_DIR=./temp


### Alignment (SamToFastq | bwa mem | MergeBamAlignment)

export NUMTHREADS=8

picard -Xmx40G SamToFastq \
INPUT=${RGID}_MarkIlluminaAdapters.bam \
FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=./temp \
| bwa mem -M -t ${NUMTHREADS} -p ${REFERENCE} /dev/stdin \
| picard -Xmx40G MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=${RGID}_FastqToSam.bam \
OUTPUT=${RGID}_Aligned.bam \
R=${REFERENCE} CREATE_INDEX=true \
ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=./temp 


### Qualimap

export NUMTHREADS=8

qualimap bamqc -bam ${RGID}_Aligned.bam -c -nt ${NUMTHREADS} --java-mem-size=40G


################################################################################
### PER-SAMPLE PROCESSING

### MarkDuplicates

picard -Xmx32G -Djava.io.tmpdir=./temp MarkDuplicates \
$(for i in *_Aligned.bam; do echo "INPUT=${i} "; done) \
OUTPUT=${NAME}_MarkDuplicates.bam \
METRICS_FILE=${NAME}_MarkDuplicates_metrics.txt \
MAX_RECORDS_IN_RAM=150000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
CREATE_INDEX=true \
TMP_DIR=./temp


### Indel realignment 1 (RealignerTargetCreator)

export NUMTHREADS=4

gatk3 -Xmx32G -Djava.io.tmpdir=./temp -T RealignerTargetCreator \
-nt ${NUMTHREADS} \
-R ${REFERENCE} \
-I ${NAME}_MarkDuplicates.bam \
-o ${NAME}_Realign.bam.intervals


### Indel realignment 2 (IndelRealigner)

gatk3 -Xmx32G -Djava.io.tmpdir=./temp -T IndelRealigner \
-R ${REFERENCE} \
-I ${NAME}_MarkDuplicates.bam \
-o ${NAME}_Realign.bam \
-targetIntervals ${NAME}_Realign.bam.intervals


### Generate raw variants for BQSR (per chromosome)

bcftools mpileup -f ${REFERENCE} -r ${CHR} --annotate INFO/ADF,INFO/ADR \
-Q 20 -q 30 --rf 2 --ff 256 --ff 1024 -Ou ${NAME}_Realign.bam \
| bcftools call -mv -Ou \
| bcftools filter -i 'N_ALT=1 && INFO/ADF[1]>0 && INFO/ADR[1]>0' \
-Oz -o ${NAME}_BQSR_${CHR}.vcf.gz

tabix -p vcf ${NAME}_BQSR_${CHR}.vcf.gz


### Concatenate VCF files

ls ${NAME}_BQSR_*.vcf.gz > ${NAME}_BQSR_vcf.list

bcftools concat -f ${NAME}_BQSR_vcf.list -Oz -o ${NAME}_BQSR.vcf.gz

tabix -p vcf ${NAME}_BQSR.vcf.gz


### BQSR: make recalibration table (BaseRecalibrator)

export NUMTHREADS=8

gatk3 -Xmx40g -Djava.io.tmpdir=./temp -T BaseRecalibrator \
-nct ${NUMTHREADS} \
-R ${REFERENCE} \
-I ${NAME}_Realign.bam \
-knownSites ${NAME}_BQSR.vcf.gz \
-o ${NAME}_BQSR_recal.table1


### BQSR: make recalibrated bam file (PrintReads)

export NUMTHREADS=8

gatk3 -Xmx40g -Djava.io.tmpdir=./temp -T PrintReads \
-nct ${NUMTHREADS} \
-R ${REFERENCE} \
-BQSR ${NAME}_BQSR_recal.table1 \
--disable_indel_quals \
-I ${NAME}_Realign.bam \
-o ${NAME}_Recal.bam


### BQSR: make post-recalibration table (BaseRecalibrator)

export NUMTHREADS=8

gatk3 -Xmx40g -Djava.io.tmpdir=./temp -T BaseRecalibrator \
-nct ${NUMTHREADS} \
-R ${REFERENCE} \
-I ${NAME}_Realign.bam \
-knownSites ${NAME}_BQSR.vcf.gz \
-BQSR ${NAME}_BQSR_recal.table1 \
-o ${NAME}_BQSR_recal.table2


### Make recalibration plots (AnalyzeCovariates)

gatk3 -Xmx16g -Djava.io.tmpdir=./temp -T AnalyzeCovariates \
-R ${REFERENCE} \
-before ${NAME}_BQSR_recal.table1 \
-after ${NAME}_BQSR_recal.table2 \
-plots ${NAME}_BQSR_recal.plots.pdf


### Qualimap and DepthOfCoverage

export NUMTHREADS=8

qualimap bamqc -bam ${NAME}_Recal.bam -c -nt ${NUMTHREADS} --java-mem-size=96G

gatk3 -Xmx96G -Djava.io.tmpdir=./temp -T DepthOfCoverage \
-R ${REFERENCE} \
-rf MappingQualityUnavailable -mmq 20 -mbq 20 -omitBaseOutput -omitIntervals \
-I ${NAME}_Recal.bam \
-o ${NAME}_Recal.bam


### Generate gVCF files (per chromosome)

export NUMTHREADS=4

gatk3 -Xmx40g -Djava.io.tmpdir=./temp -T HaplotypeCaller \
-nct ${NUMTHREADS} \
-R ${REFERENCE} \
-ERC BP_RESOLUTION \
-mbq 20 \
-out_mode EMIT_ALL_SITES \
-I ${NAME}_Recal.bam \
-L ${REFERENCE}.contiglist_${CHR}.list \
-o ${NAME}_${CHR}.g.vcf.gz


################################################################################
### JOINT VCF FILE PROCESSING

### Generate joint VCF files from per-sample gVCF files (per chromosome)

gatk3 -Xmx32g -Djava.io.tmpdir=./temp -T GenotypeGVCFs \
-R ${REFERENCE} \
-allSites \
-stand_call_conf 0 \
-L ${REFERENCE}.contiglist_${CHR}.list \
$(for j in *_${CHR}.g.vcf.gz; do echo "-V ${j} "; done) \
-o vaquita_20_${CHR}.vcf.gz


### Trim unused alternate alleles

gatk3 -Xmx32g -Djava.io.tmpdir=./temp -T SelectVariants \
-R ${REFERENCE} \
-trimAlternates \
-L ${REFERENCE}.contiglist_${CHR}.list \
-V vaquita_20_${CHR}.vcf.gz \
-o vaquita_20_${CHR}_TrimAlt.vcf.gz


### Add annotations

gatk3 -Xmx32g -Djava.io.tmpdir=./temp -T VariantAnnotator \
-R ${REFERENCE} \
-G StandardAnnotation \
-A VariantType \
-A AlleleBalance \
-L ${REFERENCE}.contiglist_${CHR}.list \
-V vaquita_20_${CHR}_TrimAlt.vcf.gz \
-o vaquita_20_${CHR}_TrimAlt_AddAnnot.vcf.gz


### SnpEff annotation
# Note: Following instructions in SnpEff manual, custom SnpEff annotation database created
# from GCF_008692025.1_mPhoSin1.pri_genomic.gtf.gz.

export SNPEFF=~/project/programs/snpEff/snpEff.jar
export DATABASE=mPhoSin1_RefSeq
export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot.vcf.gz
export OUT=${VCF%.vcf*}_snpEff.vcf.gz
export CSV=${OUT%.vcf.gz}_summary.csv

java -Xmx8g -jar ${SNPEFF} -v -nodownload -csvStats ${CSV} ${DATABASE} ${VCF} \
| bgzip > ${OUT}

tabix -p vcf ${OUT}


### SIFT annotation
# Note: Following instructions at https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB,
# custom SIFT database created from GCF_008692025.1_mPhoSin1.pri_genomic.gtf.gz with
# SIFT4G algorithm (https://github.com/rvaser/sift4g).

export SIFT=~/project/programs/sift/SIFT4G_Annotator_v2.4.jar
export DATABASE=~/project/programs/sift/databases/GCF_008692025.1_mPhoSin1.pri 
export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot_snpEff.vcf.gz
export LOG=${VCF%.vcf.gz}_SIFT.log

gzip -cd ${VCF} > ${VCF%.gz}

java -Xmx16g -jar ${SIFT} -c -t -d ${DATABASE} -i ${VCF%.gz} |& tee ${LOG}

awk '{if ($1~"^#"){ gsub("##SIFT_Threshold: 0.05", "##SIFT_Threshold=0.05"); print } \
else { gsub(" ", "_"); print }}' ${VCF%.vcf.gz}_SIFTpredictions.vcf \
| bgzip > ${VCF%.vcf.gz}_SIFT.vcf.gz

tabix -p vcf ${VCF%.vcf.gz}_SIFT.vcf.gz

rm ${VCF%.gz}
rm ${VCF%.vcf.gz}_SIFTpredictions.vcf


### Mask repeats (from RepeatMasker and TandemRepeatsFinder)

export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot_snpEff_SIFT.vcf.gz

gatk3 -Xmx16g -Djava.io.tmpdir=./temp \
-T VariantFiltration \
-R ${REFERENCE} \
-mask ${REPEATMASK} -maskName "FAIL_Repeat" \
-l ERROR \
-V ${VCF} \
-o ${VCF%.vcf.gz}_Mask.vcf.gz


### Custom filtering (see filterVCF_vaquita.py)

export SCRIPT=~/project/vaquita/scripts/filterVCF_vaquita.py
export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot_snpEff_SIFT_Mask.vcf.gz

python2.7 ${SCRIPT} ${VCF} | bgzip > ${VCF%.vcf.gz}_Filter.vcf.gz

tabix -p vcf ${VCF%.vcf.gz}_Filter.vcf.gz


### Apply SNP density mask (see generate_SNP_density_mask.sh)

export MASK=vaquita_20_highSNPdensity.bed
export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot_snpEff_SIFT_Mask_Filter.vcf.gz

gatk3 -Xmx16g -Djava.io.tmpdir=./temp \
-T VariantFiltration \
-R ${REFERENCE} \
-mask ${MASK} -maskName "FAIL_SNPdensity" \
-l ERROR \
-V ${VCF} \
-o ${VCF%.vcf.gz}_Mask2.vcf.gz


### Simplify & exclude sites failing filters, extract SNPs

export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot_snpEff_SIFT_Mask_Filter_Mask2.vcf.gz

bcftools annotate \
-x ^INFO/AC,INFO/AF,INFO/AN,INFO/ANN,INFO/LOF,INFO/NMD,INFO/SIFTINFO,INFO/VariantType,FORMAT \
-Ov ${VCF} \
| awk '{ if ( $0~/^#/ || $7=="PASS" ){print $0}}' \
| bgzip > vaquita_20_${CHR}_simple_PASS.vcf.gz

tabix -p vcf vaquita_20_${CHR}_simple_PASS.vcf.gz

zcat vaquita_20_${CHR}_simple_PASS.vcf.gz \
| awk '{ if ( $0~/^#/ || $8~"VariantType=SNP" ){print $0}}' \
| bgzip > vaquita_20_${CHR}_simple_PASS_variants.vcf.gz

tabix -p vcf vaquita_20_${CHR}_simple_PASS_variants.vcf.gz


### Concatenate

ls -v vaquita_20_chr{1..21}_simple_PASS.vcf.gz > vcf.list1
ls -v vaquita_20_chr{1..21}_simple_PASS_variants.vcf.gz > vcf.list2

export NUMTHREADS=8

bcftools concat -f vcf.list1 --threads ${NUMTHREADS} -Oz -o vaquita_20_simple_PASS_autos.vcf.gz
bcftools concat -f vcf.list2 --threads ${NUMTHREADS} -Oz -o vaquita_20_simple_PASS_autos_variants.vcf.gz

