'''
Custom filtering of VCF files

Input: raw VCF
Output: filtered VCF prints to screen
- Sites failing filters are marked as FAIL_? or WARN_? in the 7th column
- Sites where REF is in [A,C,G,T] and ALT is in [A,C,G,T,.] go on to genotype filtering if AD and DP present in FORMAT
- Filtered out genotypes are changed to './.', all others reported
- Sites with non-reference alleles remaining after genotype filtering also filtered based on values in INFO column

Possible usage:

SCRIPT=filterVCF.py
python ${SCRIPT} myfile.vcf.gz [min] [max] | bgzip > myfile_filtered.vcf.gz
tabix -p vcf myfile_filtered.vcf.gz

'''

import sys
import gzip
import re

vcf_file = sys.argv[1]
VCF = gzip.open(vcf_file, 'rt')


# Min depth (1/3x mean)
minD=int(sys.argv[2])

# Max depth (2x mean)
maxD=int(sys.argv[3])


# Individual genotype filtering function
#  - Genotypes failing filters are set to missing (./.)
#  - Applies individual min and max depth filters
#  - Filters heterozygotes if the allele balance (REF/DP) is <20% or >80%
#  - Filters homozygotes if more than 10% of the alleles are different type
#  - 'sample' is the sample name
#  - 'GT_entry' is the entire genotype entry for that individual (typically GT:AD:DP:GQ)
#  - 'ADpos' is the position of the AD field in FORMAT (determined below)
#  - 'DPpos' is the position of the DP field in FORMAT (determined below)

def GTfilter(sample, GT_entry, ADpos, DPpos):
    if GT_entry[:1]=='.' : return GT_entry
    else:
        gt=GT_entry.split(':')
        nocall='./.:' + ':'.join(gt[1:])
        if gt[0] in ('0/0','0/1','1/1') and gt[DPpos]!='.':
            DP=int(gt[DPpos])
            if minD<=DP<=maxD:
                REF=float(gt[ADpos].split(',')[0])
                AB=float(REF/DP)
                if gt[0]=='0/0':
                    if AB>=0.9: return GT_entry
                    else: return nocall
                elif gt[0]=='0/1':
                    if 0.2<=AB<=0.8: return GT_entry
                    else: return nocall
                elif gt[0]=='1/1':
                    if AB<=0.1: return GT_entry
                    else: return nocall
                else: return nocall
            else: return nocall
        else: return nocall


# Get list of samples in VCF file
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break


# Go back to beginning of file
VCF.seek(0)


# Write pre-existing header lines & add new lines describing filters being applied
for line0 in VCF:
    if line0.startswith('#'):
        if line0.startswith('##FORMAT'):
            sys.stdout.write('##FILTER=<ID=FAIL_REF,Description="Reference allele not one of [A,C,G,T].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_ALT,Description="Alternate allele not one of [A,C,G,T,.].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noADi,Description="AD not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noDPi,Description="DP not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noGT,Description="No called genotypes remain after filtering.">\n')
            sys.stdout.write('##FILTER=<ID=WARN_missing,Description="Excess missingness (>25% of samples uncalled or set to missing).">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_QD,Description="QD < 4.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_FS,Description="FS > 60.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_MQ,Description="MQ < 40.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_MQRankSum,Description="MQRankSum < -12.5.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_ReadPosRankSum,Description="ReadPosRankSum < -8.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_SOR,Description="SOR > 3.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_excessHet,Description="Excess heterozygosity (>75% of genotypes are 0/1).">\n')
            sys.stdout.write(line0)
            break
        else: sys.stdout.write(line0)


# Go through VCF file line by line to apply filters
for line0 in VCF:
    if line0.startswith('#'):
        sys.stdout.write(line0)
        continue
    line=line0.strip().split('\t')

### Site filtering:
### Keep any filters that have already been applied
    filter=[]
    if line[6] not in ('.', 'PASS'):
        filter.append(line[6])

### Check REF allele
    if line[3] not in ['A','C','G','T']:       
        filter.append('FAIL_REF') 

### Check ALT allele
    if line[4] not in ['A','C','G','T','.']:
        filter.append('FAIL_ALT') 

### Access INFO field annotations
    if ';' in line[7]:
        INFO=line[7].split(';')
        d=dict(x.split('=') for x in INFO)
    else:
        INFO=line[7]
        if '=' in INFO:
            d={INFO.split('=')[0]:INFO.split('=')[1]}
        else: d={}

### Get the position of AD, DP in genotype fields
    if 'AD' in line[8]:
        ADpos=line[8].split(':').index('AD')
    else: filter.append('FAIL_noADi')

    if 'DP' in line[8]:
        DPpos=line[8].split(':').index('DP')
    else: filter.append('FAIL_noDPi')

### If any filters failed, write out line and continue
    if filter!=[]:
        sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:])) )
        continue

### Genotype filtering:
    GT_list=[]
    for i in range(0,len(samples)):
        GT=GTfilter(samples[i],line[i+9],ADpos,DPpos)
        GT_list.append(GT)

### Recalculate AC, AN, AF for INFO (after this step, modified INFO values will be output)
    REF=2*[x[:3] for x in GT_list].count('0/0') + [x[:3] for x in GT_list].count('0/1')
    ALT=2*[x[:3] for x in GT_list].count('1/1') + [x[:3] for x in GT_list].count('0/1')
    if REF+ALT==0:
        filter.append('FAIL_noGT')
        sys.stdout.write('%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:9]), '\t'.join(GT_list)) )
        continue    
    d['AC']=ALT
    d['AN']=REF+ALT
    d['AF']=round(float(ALT)/(float(REF)+float(ALT)), 4)

### Set VariantType, outputting sites with just hom. REF genotypes without further filtering
    if ALT==0:
        d['VariantType']='NO_VARIATION'   
        if filter==[]:
            filter.append('PASS')
        sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )
        continue
    elif REF==0:
        d['VariantType']='NO_VARIATION'
    else:
        d['VariantType']='SNP'

### Fail sites with poor variant metrics
    if 'QD' in d and float(d['QD']) < 4.0:
        filter.append('FAIL_QD')
    if 'FS' in d and float(d['FS']) > 60.0:
        filter.append('FAIL_FS')
    if 'MQ' in d and float(d['MQ']) < 40.0:
        filter.append('FAIL_MQ')
    if 'MQRankSum' in d and float(d['MQRankSum']) < -12.5:
        filter.append('FAIL_MQRankSum')
    if 'ReadPosRankSum' in d and float(d['ReadPosRankSum']) < -8.0:
        filter.append('FAIL_ReadPosRankSum')
    if 'SOR' in d and float(d['SOR']) > 3.0:
        filter.append('FAIL_SOR')

### Write out new line
    if filter==[]:
        filter.append('PASS')
    sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )


# Close files and exit
VCF.close()
exit()
