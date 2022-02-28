### Identify regions with excess SNP density (>0.01) for masking

### Use bedtools to generate coordinates of 10 kb windows with 1 kb step size
# Note: .sizes is a two-column tab-delimited file containing chromosome name, length

export SIZES=~/project/vaquita/reference/mPhoSin1.pri.fasta.sizes
export WINS=${SIZES%.sizes}_10kbwins_1kbstep.bed

bedtools makewindows -g ${SIZES} -w 10000 -s 1000 > ${WINS}


### Count SNPs in windows (per chromosome)

export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot_Rename_snpEff_SIFT_Mask_Filter.vcf.gz

awk -v chr=${CHR} '$1==chr' ${WINS} > wins_${CHR}.temp

while read -r CHROM START END ; do
SNPCOUNT=$(tabix ${VCF} ${CHROM}:$((START+1))-${END} | awk '$7=="PASS" && $8~"VariantType=SNP"' | wc -l)
echo -e "${CHROM}\t$((START+1))\t${END}\t${SNPCOUNT}" >> ${VCF%.vcf.gz}_snpcount_${CHR}.bed
done < wins_${CHR}.temp

rm wins_${CHR}.temp


### Generate mask bed file with regions of excess SNP density (>0.01)

export VCF=vaquita_20_${CHR}_TrimAlt_AddAnnot_Rename_snpEff_SIFT_Mask_Filter.vcf.gz
export SIZES=~/project/vaquita/reference/mPhoSin1.pri.fasta.sizes

cat $(ls -v ${VCF%.vcf.gz}_snpcount_*.bed) \
| awk '$4/($3-$2+1)>0.01' \
| awk '{printf "%s\t%s\t%s\n", $1, $2-1, $3}' \
| bedtools merge -i stdin > temp.bed
bedtools sort -g ${SIZES} -i temp.bed > vaquita_20_highSNPdensity.bed


