### make_CDS_bed_from_gtf.py
# Jacqueline Robinson 2021

# Generate a bed file from GTF input
# Output columns are: chr, start, end, transcript_id, strand
# Output is one line per CDS line in the GTF file
#
# Usage:
# SCRIPT=make_CDS_bed_from_gtf.py
# python ${SCRIPT} genes.gtf
# OR
# python ${SCRIPT} genes.gtf.gz

import sys
import gzip

filename=sys.argv[1]
if filename[-3:]==".gz":
    infile=gzip.open(filename, 'rt')
    outfile=open(filename[:-3]+'_CDS.bed', 'w')
else:
    infile=open(filename,'r')
    outfile=open(filename+'_CDS.bed', 'w')

feature_type="CDS"

for line in infile:
    if line.startswith('#'): continue
    if "unknown_transcript_1" in line: continue
    line=line.strip().split('\t')
    if line[2]==feature_type:
        chr=line[0]
        start=int(line[3])-1
        end=int(line[4])
        strand=line[6]
        attr=[x.strip() for x in line[8].split(';')]
        #gene_id=[x.split()[1].strip('\"') for x in attr if x.startswith('gene_id "')][0]
        transcript_id=[x.split()[1].strip('\"') for x in attr if x.startswith('transcript_id "')][0]
        outfile.write("%s\t%s\t%s\t%s\t%s\n" % (chr, start, end, transcript_id, strand))

outfile.close()
infile.close()

exit()

