### Runs of homozygosity (ROH)

### bcftools

# Generate list of samples, exclude related individuals for calculating allele frequencies
# in bcftools roh command later
VCF=vaquita_20_simple_PASS_autos_variants.vcf.gz
zcat ${VCF} | head -n 1000 | grep "^#" | tail -n 1 | cut -f10- | tr '\t' '\n' > samples.list
KEEP=samples_keep.list
grep -v "z0001663\|z0004380\|z0004393\|z0004394\|z0185383" samples.list > ${KEEP}

# Identify ROH
bcftools roh -e ${KEEP} -G 30 -Orz -o ${VCF}_bcftoolsROH.txt.gz ${VCF}

# Reformat output
zcat ${VCF}_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> ${VCF}_bcftoolsROH.txt

# Run with a pseudo-genome with entirely 0/0 genotypes
IN=vaquita_20_simple_PASS_autos_variants.vcf.gz
OUT=vaquita_20_simple_PASS_autos_variants_pseudohom.vcf.gz
KEEP=samples_keep.list

zcat ${IN} | head -n 1000 | grep "^#" > head.tmp
# Manually edit head.tmp to add a sample name to the last line ("pseudohom")
zcat ${IN} | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat head.tmp - | bgzip > ${OUT}
bcftools roh -e ${KEEP} -G 30 -Orz -o ${OUT}_bcftoolsROH.txt.gz ${OUT}
zcat ${OUT}_bcftoolsROH.txt.gz \
| awk -v s=pseudohom 'BEGIN{sum=0}{if ($2==s){sum+=$6}}END{printf "%s\t%s\n", s, sum}'
# pseudohom	ROH length: 2237397028

# Calculate Froh using max length calculated from pseudo-genome
DATA=vaquita_20_simple_PASS_autos_variants.vcf.gz_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/2237397028, num, sum/num}'
done < samples.list

# z0000703	0.0683867	64	2.39075e+06
# z0001649	0.0378107	40	2.11494e+06
# z0001654	0.0480418	51	2.10762e+06
# z0001660	0.0638652	58	2.46365e+06
# z0001663	0.0965794	68	3.17774e+06
# z0004379	0.0362249	44	1.84203e+06
# z0004380	0.0788742	76	2.32201e+06
# z0004381	0.0535957	60	1.99858e+06
# z0004382	0.0432992	44	2.20176e+06
# z0004390	0.0552752	62	1.99472e+06
# z0004393	0.0619033	76	1.8224e+06
# z0004394	0.0483124	53	2.03951e+06
# z0183496	0.0566541	59	2.14843e+06
# z0184983	0.0414674	44	2.10862e+06
# z0184984	0.0457263	56	1.82693e+06
# z0185383	0.0653769	62	2.35926e+06
# z0185384	0.0728995	66	2.47129e+06
# z0185385	0.0507922	55	2.06622e+06
# z0186934	0.0253311	34	1.66693e+06
# z0186935	0.0340173	48	1.58563e+06 


### vcftools

VCF=vaquita_20_simple_PASS_autos_variants.vcf.gz
for CHR in chr{1..21} ; do
vcftools --gzvcf ${VCF} --LROH --chr ${CHR} --out ${VCF}_vcftools_${CHR}
done

OUT=${VCF}_vcftools.LROH
cp ${VCF}_vcftools_chr1.LROH  ${OUT}
cat vaquita_20_simple_PASS_autos_variants.vcf.gz_vcftools_chr{2..21}.LROH | grep -v "^CHROM" >> ${OUT}


### plink

# Convert simplified, concatenated VCF file to plink format
export PLINK=~/project/programs/plink1.9/plink
export VCF=vaquita_20_simple_PASS_autos_variants.vcf.gz
${PLINK} --make-bed --keep-allele-order --double-id --chr-set 21 --vcf ${VCF} --out ${VCF%.vcf.gz}

# Identify ROH and reformat output
FILE=vaquita_20_simple_PASS_autos_variants
${PLINK} --bfile ${VCF%.vcf.gz} --out ${VCF%.vcf.gz}_ROH --homozyg
sed -i -E 's/^\ +//g' ${FILE}_ROH.hom
sed -i -E 's/\ +/\t/g' ${FILE}_ROH.hom 

