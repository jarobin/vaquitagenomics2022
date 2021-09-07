### get_longest_transcript_per_gene.py
# Jacqueline Robinson 2021

# Extract longest transcript per gene from GTF file based on summed CDS length
# Also ensures CDS length is evenly divisible by 3
#
# Usage: 
# SCRIPT=get_longest_transcript_per_gene.py
# python ${SCRIPT} genes.gtf
# OR
# python ${SCRIPT} genes.gtf.gz

import sys
import gzip
import operator

filename=sys.argv[1]
if filename[-3:]==".gz":
    infile=gzip.open(filename, 'rt')
    outfile=open(filename[:-3]+'_longest_transcript_per_gene.txt', 'w')
else:
    infile=open(filename,'r')
    outfile=open(filename+'_longest_transcript_per_gene.txt', 'w')

feature_type="CDS"

# Generate a dictionary with one entry per gene
# Each gene entry is itself a dictionary with all transcripts for that gene and their summed CDS lengths
gene_dict={}

for line in infile:
    if line.startswith('#'): continue
    if "unknown_transcript_1" in line: continue
    line=line.strip().split('\t')
    if line[2]==feature_type:
        attr=[x.strip() for x in line[8].split(';')]
        gene_id=[x.split()[1].strip('\"') for x in attr if x.startswith('gene_id "')][0]
        transcript_id=[x.split()[1].strip('\"') for x in attr if x.startswith('transcript_id "')][0]
        feature_length=int(line[4])-int(line[3])+1
        if gene_id in gene_dict:
            if transcript_id in gene_dict[gene_id]:
                gene_dict[gene_id][transcript_id]+=feature_length
            else:
                gene_dict[gene_id][transcript_id]=feature_length
        else:
            new={transcript_id:feature_length}
            gene_dict[gene_id]=new

# Generate a new dictionary, only keeping genes with CDS length that is a multiple of 3
new_dict={}
for g in gene_dict:
    temp_dict={key:val for key,val in gene_dict[g].items() if val % 3 == 0}
    if len(temp_dict)>0: new_dict[g]=temp_dict

# Get the longest transcripts per gene
longest_transcripts=[max(new_dict[x].iteritems(), key=operator.itemgetter(1))[0] for x in new_dict] 

for transcript in longest_transcripts:
    outfile.write("%s\n" % transcript)

outfile.close()
infile.close()

exit()

