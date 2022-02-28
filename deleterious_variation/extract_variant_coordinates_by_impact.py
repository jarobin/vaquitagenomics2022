# From a VCF file, create bed files with coordinates of mutations in different classes:
# - CDS: all coding sites
# - LOF: loss of function (as determined by SnpEff or SIFT)
# - MIS: all missense
# - DEL: deleterious missense (as determined by SIFT)
# - TOL: tolerated missense (as determined by SIFT)
# - SYN: synonymous

# For sites with multiple annotations, the most deleterious annotation is prioritized:
# LOF > deleterious nonsynonymous > tolerated nonsynonymous > synonymous

import sys
import gzip

# Open input file
vcf_file = sys.argv[1]
VCF = gzip.open(vcf_file, 'rt')

# Variant types
# CDS: all coding sites
# LOF: loss of function
# MIS: all missense
# DEL: deleterious missense
# TOL: tolerated missense
# SYN: synonymous

# Open output files
my_cds=open(vcf_file+'_CDS.bed', 'w')
my_lof=open(vcf_file+'_LOF.bed', 'w')
my_mis=open(vcf_file+'_MIS.bed', 'w')
my_del=open(vcf_file+'_DEL.bed', 'w')
my_tol=open(vcf_file+'_TOL.bed', 'w')
my_syn=open(vcf_file+'_SYN.bed', 'w')

# Go line by line and check annotations, prioritizing by most deleterious impact
for line0 in VCF:
    if line0.startswith('#'): continue
    line=line0.strip().split('\t')
    if "LOF=" not in line[7] and "SIFTINFO=" not in line[7]: continue
    if "|protein_coding|" not in line[7] and "|CDS|" not in line[7]: continue
    my_cds.write("%s\t%s\t%s\n" % (line[0], int(line[1])-1, line[1]))
    if ';' in line[7]:
        INFO=line[7].split(';')
        d=dict(x.split('=') for x in INFO)
    else:
        INFO=line[7]
        d={INFO.split('=')[0]:INFO.split('=')[1]}
    if "LOF" in d.keys(): my_lof.write("%s\t%s\t%s\n" % (line[0], int(line[1])-1, line[1])) ; continue
    if "SIFTINFO" in d.keys():
        ss=d['SIFTINFO'].split(',')
        if any("|NONSYNONYMOUS|" in x for x in ss): 
            my_mis.write("%s\t%s\t%s\n" % (line[0], int(line[1])-1, line[1]))
            if any("|NONSYNONYMOUS|" in x and "|DELETERIOUS" in x for x in ss): my_del.write("%s\t%s\t%s\n" % (line[0], int(line[1])-1, line[1])) ; continue
            if any("|NONSYNONYMOUS|" in x and "|TOLERATED" in x for x in ss): my_tol.write("%s\t%s\t%s\n" % (line[0], int(line[1])-1, line[1])) ; continue
            continue
        if any("|SYNONYMOUS|" in x for x in ss): my_syn.write("%s\t%s\t%s\n" % (line[0], int(line[1])-1, line[1])) ; continue

# Close all files and exit
my_cds.close()
my_lof.close()
my_mis.close()
my_del.close()
my_tol.close()
my_syn.close()

VCF.close()
exit()

