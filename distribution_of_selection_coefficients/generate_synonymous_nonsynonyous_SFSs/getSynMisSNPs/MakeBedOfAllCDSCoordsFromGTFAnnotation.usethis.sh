########## make all sites cds coordinates:

wd=/net/harris/vol1/home/beichman/reference_genomes/vaquita/mPhoSin1.pri.cur.20190723_rename
annotation=GCF_008692025.1_mPhoSin1.pri_genomic_rename.gtf.gz
bedDir=/net/harris/vol1/home/beichman/vaquita/coordinate_bed_files/variant_effects
# want get exnoic sequence: 
# note that GTF is 1-based coordinates:
# have to subtract 1 from start b/c 0 based. don't ahve to subtract one from the end because while it is zero based (minus 1)
# it is also non-inclusive, so would want to add one. so that cancels out.
# count up exon regions:
# zcat $annotation | awk '{if($3=="CDS")print}' | wc -l # 630511 regions
# note no semicolon in the awk statement! 
zcat $wd/$annotation | awk 'BEGIN {OFS="\t"}; {if($3=="CDS") print $1,$4-1,$5}' > $bedDir/AllCDS.Coords.FromGTF.0Based.bed

wc -l $bedDir/AllCDS.Coords.FromGTF.0Based.bed #630511 -- matches good
