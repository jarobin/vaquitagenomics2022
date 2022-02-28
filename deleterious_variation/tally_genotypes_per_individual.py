# Script to tally genotype counts and number of genotypes called/not called 
# Usage example: python GTCountsPerInd.py chr1_filtered.vcf.gz

import sys
import gzip

# Open input file (gzipped VCF file)
filename = sys.argv[1]
VCF = gzip.open(filename, 'r')

# Get list of samples
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break

# Initialize counters
nocall=[0]*len(samples)
call=[0]*len(samples)
homR=[0]*len(samples)
het=[0]*len(samples)
homA=[0]*len(samples)
other=[0]*len(samples)

# Get counts
for line in VCF:
    if line.startswith('#'): continue
    line=line.strip().split('\t')
    if line[6] not in ('.','PASS'): continue
    for i in range(0,len(samples)):
        GT=line[i+9]
        if GT[:1]=='.': nocall[i]+=1
        else: 
            call[i]+=1
            if GT[:3]=='0/0': homR[i]+=1
            elif GT[:3]=='0/1': het[i]+=1
            elif GT[:3]=='1/1': homA[i]+=1
            else: other[i]+=1

# Write out final counts
output = open(filename + '_GTCountsPerInd.txt', 'w')

output.write('sample\tcall\tnocall\thomR\thet\thomA\tother\n')
for i in range(len(samples)):
    output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (samples[i], call[i], nocall[i], homR[i], het[i], homA[i], other[i]))

# Close files and exit  
output.close()
VCF.close()
exit()

