# Cetacean reads to variants pipeline

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
# - Collect statistics on all samples
# - Per chromosome, per sample:
#   - HaplotypeCaller
#   - GenotypeGVCFs
#   - Trim unused alternate alleles: SelectVariants
#   - VariantAnnotator
#   - SnpEff
#   - SIFT
#   - Mask
#   - Custom filtering
#   - Simplify
# - Concatenate VCF files per sample


################################################################################
### Set variables

export REFERENCE=~/project/vaquita/reference/GCF_000331955.2_Oorc_1.1/GCF_000331955.2_Oorc_1.1_genomic.fa
export REPEATMASK=~/project/vaquita/reference/GCF_000331955.2_Oorc_1.1/GCF_000331955.2_Oorc_1.1_repeats_TRF_RM.bed

export NAME=SAMN01180276			# Sample name
export RGID=SRR574967				# Read set


################################################################################
### PER-READ SET PROCESSING

### Trim reads with BBDUK (as needed)

export TRIM=10
export FASTQ=SRR574967_1.fastq.gz

bbduk.sh -Xmx8G threads=1 ftl=${TRIM} ordered in=${FASTQ} out=${FASTQ%.fastq.gz}_trim.fastq.gz


### FastqToSam

export FQ1=SRR574967_1_trim.fastq.gz
export FQ2=SRR574967_2_trim.fastq.gz

picard -Xmx56G FastqToSam \
FASTQ=${FQ1} \
FASTQ2=${FQ2} \
OUTPUT=${RGID}_FastqToSam.bam \
READ_GROUP_NAME=${RGID} \
SAMPLE_NAME=${NAME} \
PLATFORM=illumina \
TMP_DIR=./temp 


### MarkIlluminaAdapters

picard -Xmx56G MarkIlluminaAdapters \
INPUT=${RGID}_FastqToSam.bam \
OUTPUT=${RGID}_MarkIlluminaAdapters.bam \
METRICS=${RGID}_MarkIlluminaAdapters_metrics.txt \
TMP_DIR=./temp 


### Alignment (SamToFastq | bwa mem | MergeBamAlignment)

export NUMTHREADS=8

picard -Xmx56G SamToFastq \
INPUT=${RGID}_MarkIlluminaAdapters.bam \
FASTQ=/dev/stdout CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR=./temp \
| bwa mem -M -t ${NUMTHREADS} -p ${REFERENCE} /dev/stdin \
| picard -Xmx56G MergeBamAlignment \
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

qualimap bamqc -bam ${RGID}_Aligned.bam -c -nt ${NUMTHREADS} --java-mem-size=56G


################################################################################
### PER-SAMPLE PROCESSING

### MarkDuplicates

picard -Xmx56G -Djava.io.tmpdir=./temp MarkDuplicates \
$(for i in *_Aligned.bam; do echo "INPUT=${i} "; done) \
OUTPUT=${NAME}_MarkDuplicates.bam \
METRICS_FILE=${NAME}_MarkDuplicates_metrics.txt \
MAX_RECORDS_IN_RAM=150000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
CREATE_INDEX=true \
TMP_DIR=./temp 


### Indel realignment 1 (RealignerTargetCreator)

export NUMTHREADS=4

gatk3 -Xmx56G -Djava.io.tmpdir=./temp -T RealignerTargetCreator \
-nt ${NUMTHREADS} \
-R ${REFERENCE} \
-I ${NAME}_MarkDuplicates.bam \
-o ${NAME}_Realign.bam.intervals 


### Indel realignment 2 (IndelRealigner)

gatk3 -Xmx56G -Djava.io.tmpdir=./temp -T IndelRealigner \
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


### Concatenate VCF files (and exclude sites with ambiguous IUPAC nucleotides present in some assemblies)

ls ${NAME}_BQSR_*.vcf.gz > ${NAME}_BQSR_vcf.list

bcftools concat -f ${NAME}_BQSR_vcf.list -Ov \
| awk '{if ($0~"^#") {print $0} else {if ($4!~"U|R|Y|S|W|K|M|B|D|H|V|N" && $5!~"U|R|Y|S|W|K|M|B|D|H|V|N") {print $0} } }' \
| bgzip > ${NAME}_BQSR.vcf.gz

tabix -p vcf ${NAME}_BQSR.vcf.gz


### BQSR: make recalibration table (BaseRecalibrator)

export NUMTHREADS=8

gatk3 -Xmx56g -Djava.io.tmpdir=./temp -T BaseRecalibrator \
-nct ${NUMTHREADS} \
-R ${REFERENCE} \
-I ${NAME}_Realign.bam \
-knownSites ${NAME}_BQSR.vcf.gz \
-o ${NAME}_BQSR_recal.table1 


### BQSR: make recalibrated bam file (PrintReads)

export NUMTHREADS=8

gatk3 -Xmx56g -Djava.io.tmpdir=./temp -T PrintReads \
-nct ${NUMTHREADS} \
-R ${REFERENCE} \
-BQSR ${NAME}_BQSR_recal.table1 \
--disable_indel_quals \
-I ${NAME}_Realign.bam \
-o ${NAME}_Recal.bam 


### BQSR: make post-recalibration table (BaseRecalibrator)

export NUMTHREADS=8

gatk3 -Xmx80g -Djava.io.tmpdir=./temp -T BaseRecalibrator \
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


################################################################################
### Get coverage and scaffold length stats from all samples, genome assemblies

### Get chr/scaff coverage stats and calculate coverage thresholds for later filtering

for i in */*Recal*stats ; do FILE=$(ls ${i}/*html) ; 
TOP=$(grep -n "Chromosome stats" ${FILE} | cut -d':' -f1)
BOTTOM=$(grep -n "summary section" ${FILE} | tail -n 1 | cut -d':' -f1)
head -n ${BOTTOM} ${FILE} | tail -n +${TOP} | grep "^<td>" \
| sed 's/<b>//g' | sed 's/<\/b>//g' | sed 's/<td>//g' | sed 's/<\/td>//g' \
| sed 's/\ /_/g' | sed 's/,//g' | xargs -n 5 > ${i}/chr_cov.txt
done

cat SAMN*/*summary | grep "^SAMN" \
| awk '{printf "%s\t%s\t%s\t%s\n", $1, $3, int($3/3), int(2*$3)}' \
> gatkDOC_SAMN_mean_min_max.txt


### Identify scaffolds with coverage below 0.75x or above 1.5x the mean

for NAME in SAMN* ; do 
MEAN=$(grep "mean coverageData" ${NAME}/*Recal*stats/genome_results.txt | sed 's/\ //g' | cut -d'=' -f2 | sed 's/X//g')
echo -e "${NAME}\t${MEAN}" | awk '{printf "%s\t%s\t%s\t%s\n", $1, $2, 0.75*$2, 1.5*$2}'
done > contig_qmapDOC_SAMN_mean_min_max.txt

while read -r NAME MEAN MIN MAX ; do
grep -v "^Name" ${NAME}/${NAME}_Recal_stats/chr_cov.txt \
| awk -v MIN=${MIN} -v MAX=${MAX} '{if($4<MIN || $4>MAX){print $1}}' > ${NAME}/${NAME}_contigs_exclude.list
done < contig_qmapDOC_SAMN_mean_min_max.txt


### Check percentage of reference contained in scaffolds above various length thresholds

for i in GCF* ; do 
TOTAL=$(awk '{sum+=$2}END{print sum}' ${i}/*sizes) ;
A=$(awk -v TOTAL=${TOTAL} -v MIN=1000000 '$2>=MIN{sum+=$2}END{print sum/TOTAL}' ${i}/*sizes)
B=$(awk -v TOTAL=${TOTAL} -v MIN=500000 '$2>=MIN{sum+=$2}END{print sum/TOTAL}' ${i}/*sizes)
C=$(awk -v TOTAL=${TOTAL} -v MIN=250000 '$2>=MIN{sum+=$2}END{print sum/TOTAL}' ${i}/*sizes)
D=$(awk -v TOTAL=${TOTAL} -v MIN=100000 '$2>=MIN{sum+=$2}END{print sum/TOTAL}' ${i}/*sizes)
echo -e "${i}\t${TOTAL}\t${A}\t${B}\t${C}\t${D}" ; done

# GCF_000331955.2_Oorc_1.1	2372919875	0.985199	0.992237	0.995273	0.996727
# GCF_000493695.1_BalAcu1.0	2431687698	0.955892	0.972092	0.981173	0.988513
# GCF_002288925.2_ASM228892v3	2362774659	0.988267	0.989761	0.990926	0.991665
# GCF_002837175.2_ASM283717v2	2512149402	0.943448	0.943448	0.944834	0.954461
# GCF_003031525.1_Neophocaena_asiaeorientalis_V1	2284628084	0.910193	0.957801	0.977691	0.988598
# GCF_003676395.1_ASM367639v1	2334467171	0.983916	0.986342	0.986767	0.988164
# GCF_005190385.1_NGI_Narwhal_1	2355574979	0.983677	0.985191	0.985579	0.985852
# GCF_006547405.1_ASM654740v1	2333870982	0.963146	0.974255	0.97756	0.981934
# GCF_009873245.2_mBalMus1.pri.v3	2374868943	0.998277	0.998277	0.998432	0.998626
# GCF_011762595.1_mTurTru1.mat.Y	2378522213	0.988323	0.989373	0.991375	0.995151


### Identify scaffolds <500 kb in length
# Note: >94% of all genome assemblies contained in scaffolds >=500 kb

for i in GCF* ; do 
OUT=${i}/${i}_genomic_contigs_shorter_than_500kb.list
awk '$2<500000{print $1}' ${i}/${i}_genomic.sizes > ${OUT}
done


################################################################################
### PER-CHROMOSOME PROCESSING

### Generate gVCF files

export NUMTHREADS=4

gatk3 -Xmx56g -Djava.io.tmpdir=./temp -T HaplotypeCaller \
-nct ${NUMTHREADS} \
-R ${REFERENCE} \
-ERC BP_RESOLUTION \
-mbq 20 \
-out_mode EMIT_ALL_SITES \
-I ${NAME}_Recal.bam \
-L ${REFERENCE}.contiglist_${CHR}.list \
-o ${NAME}_${CHR}.g.vcf.gz 


### Generate VCF file from gVCF files

gatk3 -Xmx32g -Djava.io.tmpdir=./temp -T GenotypeGVCFs \
-R ${REFERENCE} \
-allSites \
-stand_call_conf 0 \
-L ${REFERENCE}.contiglist_${CHR}.list \
$(for j in *_${CHR}.g.vcf.gz; do echo "-V ${j} "; done) \
-o ${NAME}_${CHR}.vcf.gz 


### Trim unused alternate alleles from VCF file

gatk3 -Xmx32g -Djava.io.tmpdir=./temp -T SelectVariants \
-R ${REFERENCE} \
-trimAlternates \
-L ${REFERENCE}.contiglist_${CHR}.list \
-V ${NAME}_${CHR}.vcf.gz \
-o ${NAME}_${CHR}_TrimAlt.vcf.gz 


### Add annotations to VCF file

gatk3 -Xmx32g -Djava.io.tmpdir=./temp -T VariantAnnotator \
-R ${REFERENCE} \
-G StandardAnnotation \
-A VariantType \
-A AlleleBalance \
-L ${REFERENCE}.contiglist_${CHR}.list \
-V ${NAME}_${CHR}_TrimAlt.vcf.gz \
-o ${NAME}_${CHR}_TrimAlt_AddAnnot.vcf.gz 


### SnpEff annotation

export SNPEFF=~/project/programs/snpEff/snpEff.jar
export DATABASE=GCF_000331955.2_Oorc_1.1
export VCF=${NAME}_${CHR}_TrimAlt_AddAnnot.vcf.gz
export CSV=${OUT%.vcf.gz}_summary.csv

java -Xmx16g -jar ${SNPEFF} -v -nodownload -csvStats ${CSV} ${DATABASE} ${VCF} \
| bgzip > ${VCF%.vcf.gz}_snpEff.vcf.gz

tabix -p vcf ${VCF%.vcf.gz}_snpEff.vcf.gz


### SIFT annotation

export SIFT=~/project/programs/sift/SIFT4G_Annotator_v2.4.jar
export DATABASE=~/project/programs/sift/databases/GCF_000331955.2_Oorc_1.1 
export VCF=${NAME}_${CHR}_TrimAlt_AddAnnot_snpEff.vcf.gz
export LOG=${VCF%.vcf.gz}_SIFT.log

gzip -cd ${VCF} > ${VCF%.gz}

java -Xmx16g -jar ${SIFT} -c -t -d ${DATABASE} -i ${VCF%.gz} |& tee ${LOG}

awk '{if ($1~"^#"){ gsub("##SIFT_Threshold: 0.05", "##SIFT_Threshold=0.05"); print } \
else { gsub(" ", "_"); print }}' ${VCF%.vcf.gz}_SIFTpredictions.vcf \
| bgzip > ${VCF%.vcf.gz}_SIFT.vcf.gz

tabix -p vcf ${VCF%.vcf.gz}_SIFT.vcf.gz

rm ${VCF%.gz}
rm ${VCF%.vcf.gz}_SIFTpredictions.vcf


### Mask
# - Exclude short chrs/scaffs and those with aberrant coverage (based on sample depth)
# - Mask repeats (from RepeatMasker and TandemRepeatsFinder)
# - .sizes is a two-column tab-delimited file containing chromosome name, length

export DATABASE=GCF_000331955.2_Oorc_1.1
export SIZES=${DATABASE}_genomic.sizes
export CONTIGMASK_A=${DATABASE}_genomic_contigs_shorter_than_500kb.list
export CONTIGMASK_B=${NAME}_contigs_exclude.list
export VCF=${NAME}_${CHR}_TrimAlt_AddAnnot_snpEff_SIFT.vcf.gz

cat ${CONTIGMASK_A} ${CONTIGMASK_B} | sort -V | uniq > excludeContig.list

while read -r CONTIG ; do
awk -v var=${CONTIG} '$1==var{printf "%s\t0\t%s\n", $1, $2}' ${SIZES}
done < excludeContig.list > excludeContig.bed

cat excludeContig.bed ${REPEATMASK} \
| sort -k1,1 -k2,2n \
| bedtools merge -i stdin \
| bedtools sort -g ${SIZES} -i stdin > ${NAME}_mask.bed

gatk3 -Xmx16g -Djava.io.tmpdir=./temp \
-T VariantFiltration \
-R ${REFERENCE} \
-mask ${NAME}_mask.bed -maskName "FAIL_Mask" \
-l ERROR \
-V ${VCF} \
-o ${VCF%.vcf.gz}_Mask.vcf.gz


### Custom filtering (see filterVCF_cetaceans.py)
# Min and max thresholds from gatkDOC_SAMN_mean_min_max.txt, calculated earlier

export MIN
export MAX

export SCRIPT=~/project/vaquita/scripts/filterVCF_cetaceans.py
export VCF=${NAME}_${CHR}_TrimAlt_AddAnnot_snpEff_SIFT_Mask.vcf.gz

python2.7 ${SCRIPT} ${VCF} ${MIN} ${MAX} | bgzip > ${VCF%.vcf.gz}_Filter.vcf.gz

tabix -p vcf ${VCF%.vcf.gz}_Filter.vcf.gz


### Simplify VCF files and exclude sites failing filters

export VCF=${NAME}_${CHR}_TrimAlt_AddAnnot_snpEff_SIFT_Mask_Filter.vcf.gz

bcftools annotate -i 'FILTER="PASS"' \
-x ^INFO/ANN,INFO/LOF,INFO/NMD,INFO/SIFTINFO,FORMAT \
-Oz -o ${VCF%.vcf.gz}_Simple_PASS.vcf.gz ${VCF}

tabix -p vcf ${VCF%.vcf.gz}_Simple_PASS.vcf.gz


################################################################################

### Concatenate

ls *_Simple_PASS.vcf.gz > ${NAME}_vcf.list_Simple

bcftools concat -f ${NAME}_vcf.list_Simple -Oz -o ${NAME}_Simple_PASS.vcf.gz

tabix -p vcf ${NAME}_Simple_PASS.vcf.gz

