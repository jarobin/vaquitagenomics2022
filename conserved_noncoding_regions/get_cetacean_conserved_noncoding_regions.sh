# Get coordinates of conserved noncoding regions in cetacean genomes

# Required:
# - Mouse 60-way placental mammals conservation track from UCSC (see below)
# - Reciprocal best chain file for mouse genome to cetacean genome liftover
# - Genome annotation with CDS coordinates (GTF file)
# - Sizes file (2-column tab-delimited file with chromosomes and sizes for the cetacean genome)
# - Repeat coordinates (we use RepeatMasker + Tandem Repeats Finder coordinates, following UCSC)
# - Utilities: bedtools, liftOver (from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64)

# Mouse 60-way conservation track available from http://genome.ucsc.edu/cgi-bin/hgTables
# clade: mammal, genome: mouse, assembly: mm10
# group: all tracks, track: Conservation
# table: Placental El (phastConsElements60wayPlacental)
# region: genome
# gzipped BED output
# 4th column is lod=# (log odds), 5th column is transformed lod score (ranges from 0-1000)

# From UCSC: PhastCons is sensitive to "runs" of conserved sites, and is therefore effective for picking out conserved elements. 


# Keep columns 1-3 from bed file output by UCSC
cd ~/vaquita/liftover
zcat ../mm10_phastConsElements60wayPlacental.bed.gz \
| cut -f1-3 > mm10_phastConsElements60wayPlacental_coords.bed

# Select regions >=50 bp
IN=mm10_phastConsElements60wayPlacental_coords.bed
OUT=mm10_phastConsElements60wayPlacental_coords_min50bp.bed
awk '$3-$2>=50' ${IN} > ${OUT}

# Check summed lengths of regions
awk '{sum+=$3-$2}END{print sum}' mm10_phastConsElements60wayPlacental_coords.bed
# 127619787
awk '{sum+=$3-$2}END{print sum}' mm10_phastConsElements60wayPlacental_coords_min50bp.bed
# 59365871


# Lift coordinates from mouse genome to new genome

cd ~/vaquita/liftover
OLD=mm10_phastConsElements60wayPlacental_coords_min50bp.bed

# For loop: run on each cetacean genome
for i in GCF_* ; do
CHAIN=mm10.${i}.rbest.chain.gz
cp ${i}/${CHAIN} .
NEW=${OLD%.bed}_lifted_rbest_${i}.bed
UNMAPPED=${NEW%.bed}_unmapped.bed
liftOver ${OLD} ${CHAIN} ${NEW} ${UNMAPPED}
done


# Get CDS coordinates

cd ~/vaquita/cetaceans/reference
for i in GCF* ; do
cd ~/vaquita/cetaceans/reference/${i}
GTF=${i}_genomic.gtf.gz
SIZES=${i}_genomic.sizes
zcat ${GTF} \
| grep -v "^#" \
| awk '$3=="CDS"{printf "%s\t%s\t%s\n", $1, $4-1, $5}' \
| sort -k1,1 -k2,2n \
| bedtools merge -i stdin \
| bedtools sort -g ${SIZES} -i stdin > ${i}_genomic_CDS_all_flattened.bed
done


# Get coordinates of conserved regions that do not overlap with CDS or repeats

OLD=mm10_phastConsElements60wayPlacental_coords_min50bp.bed
cd ~/vaquita/liftover

for i in GCF_* ; do
CONSERVED=${OLD%.bed}_lifted_rbest_${i}.bed
CDS=~/vaquita/cetaceans/reference/${i}/${i}_genomic_CDS_all_flattened.bed
REPEATS=~/vaquita/cetaceans/reference/${i}/${i}_repeats_TRF_RM.bed
SIZES=~/vaquita/cetaceans/reference/${i}/${i}_genomic.sizes
# Exclude CDS
bedtools intersect -v -a ${CONSERVED} -b ${CDS} > ${i}_temp_a.bed
# Exclude repeats
bedtools intersect -v -a ${i}_temp_a.bed -b ${REPEATS} > ${i}_temp_b.bed
bedtools sort -g ${SIZES} -i ${i}_temp_b.bed > ${OLD%.bed}_lifted_rbest_${i}_noncoding_norepeats.bed
rm ${i}_temp_a.bed
rm ${i}_temp_b.bed
done

