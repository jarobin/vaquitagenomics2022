### Pipeline for identifying 0- and 4-fold degenerate sites
#
# Identifies and uses the longest transcript per gene that:
# - Has no errors from SnpEff (i.e. has a start and stop codon)
# - Has a total CDS length that is a multiple of 3
#   - Note: 'CDS' sequence begins with a start codon and ends with a stop codon. The 'exon' entries in the GTF file can contain UTR sequence.
#
# Input files: 
# - Reference genome (use the fasta files provided by RefSeq)
# - GTF file (use the GTF files provided by RefSeq)
# - SnpEff database generated with the same GTF file
#
# Output files:
# - Bed files with coordinates of zero-fold and four-fold degenerate sites
#
# Required software:
# - SnpEff
# - gffread
# - Python modules: HTSeq, Bio, itertools
#   - Note: Everything in Python 2, not Python 3
#
# Scripts used:
# - get_longest_transcript_per_gene.py
# - make_CDS_bed_from_gtf.py
# - Finding_SynNonSyn_Sites_nostopcodons_allfold.py (from Clare Marsden)
# - quick_remove_dups.py (from Clare Marsden)
#
# Usage: 
# export SNPEFF=~/project/programs/snpEff/snpEff.jar
# export GFFREAD=~/project/programs/gffread/gffread
# export NAME=GCF_011762595.1_mTurTru1.mat.Y
# export SNPEFF_DB=${NAME}
# export GTF_ORIG=/wynton/group/wall/jacqueline/vaquita/cetaceans/reference/${NAME}/${NAME}_genomic.gtf.gz
# export REFERENCE=/wynton/group/wall/jacqueline/vaquita/cetaceans/reference/${NAME}/${NAME}_genomic.fa
# ./get_zero_four_fold_coords.sh

set -eo pipefail

echo "Name: ${NAME}"
echo "GTF file: ${GTF_ORIG}"
echo "SnpEff database name: ${SNPEF_DB}"
echo "Reference fasta: ${REFERENCE}"

# Unzip original GTF file and create a temp GTF file (excludes lines with empty transcript_id)
echo "Unzipping GTF file and creating temp GTF file..."
TEMP_GTF=$(basename ${GTF_ORIG%.gz})
zcat ${GTF_ORIG} | grep -v 'transcript_id "";' > ${TEMP_GTF}

# Generate SnpEff log file to get ERROR/WARNING transcript IDs
echo "Getting list of failing transcripts from SnpEff dump..."
SNPEFF_LOG=${SNPEFF_DB}_snpEff_dump_debug_canonical.log
java -Xmx8g -jar ${SNPEFF} dump -v -txt -debug -nodownload -canon ${SNPEFF_DB} > /dev/null 2> ${SNPEFF_LOG}
FAIL_LIST=${SNPEFF_LOG%.log}_fail.list
grep "ERROR\|WARNING" ${SNPEFF_LOG} | tr ' ' '\n' | grep "XM_" | sed "s/'//g" | sort | uniq > ${FAIL_LIST}

# Use gffread to eliminate failing transcripts from unzipped GTF file
echo "Generating GTF file without failing transcripts..."
OUT=${TEMP_GTF%.gtf}_coding_noPseudo_noFail.gtf
${GFFREAD} -C --no-pseudo --nids ${FAIL_LIST} -T -o ${OUT} ${TEMP_GTF}

# Get list of longest transcripts per gene (which also have CDS length evenly divisible by 3)
echo "Getting list of longest transcript per gene..."
SCRIPT=get_longest_transcript_per_gene.py
IN=${TEMP_GTF%.gtf}_coding_noPseudo_noFail.gtf
python ${SCRIPT} ${IN}

# Use gffread to extract just the longest transcripts
echo "Generating GTF file with longest transcripts only..."
KEEP_LIST=${TEMP_GTF%.gtf}_coding_noPseudo_noFail.gtf_longest_transcript_per_gene.txt
IN=${TEMP_GTF%.gtf}_coding_noPseudo_noFail.gtf
OUT=${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf
${GFFREAD} --ids ${KEEP_LIST} -T -o ${OUT} ${IN}

# Convert to bed format
echo "Generating bed file of longest transcripts..."
SCRIPT=make_CDS_bed_from_gtf.py
IN=${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf
python ${SCRIPT} ${IN}

# Identify 0-,2-,3-,4-fold positions and remove duplicate sites
# Note: For duplicates, the first site is kept and all others encountered subsequently 
# with the same coordinate are discarded.
# Note: Codons containing any bases that are not in [A,C,G,T] are ignored.
echo "Running Finding_SynNonSyn_Sites_nostopcodons_allfold.py and quick_remove_dups.py..."
export SCRIPT1=Finding_SynNonSyn_Sites_nostopcodons_allfold.py
export SCRIPT2=quick_remove_dups.py
export IN=${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf_CDS.bed
export OUT=${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf_CDS.bed_fold.txt
python ${SCRIPT1} ${REFERENCE} ${IN} ${OUT}
python ${SCRIPT2} ${OUT}

# Generate bed files for 0- and 4-fold degenerate positions
echo "Generating bed files for zero- and four-fold degenerate positions..."
IN=nodups_${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf_CDS.bed_fold.txt
OUT_A=${TEMP_GTF%.gtf}_zerofold_degenerate.bed
OUT_B=${TEMP_GTF%.gtf}_fourfold_degenerate.bed
grep -v "^#" ${IN} | awk '$3=="Zero_Fold"{printf "%s\t%s\t%s\n", $1, $2-1, $2}' > ${OUT_A}
grep -v "^#" ${IN} | awk '$3=="Four_Fold"{printf "%s\t%s\t%s\n", $1, $2-1, $2}' > ${OUT_B}

# Remove intermediate files
rm ${TEMP_GTF}
rm ${SNPEFF_LOG}
rm ${FAIL_LIST}
rm ${TEMP_GTF%.gtf}_coding_noPseudo_noFail.gtf
rm ${TEMP_GTF%.gtf}_coding_noPseudo_noFail.gtf_longest_transcript_per_gene.txt
rm ${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf
rm ${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf_CDS.bed
rm ${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf_CDS.bed_fold.txt
rm nodups_${TEMP_GTF%.gtf}_coding_noPseudo_noFail_longest.gtf_CDS.bed_fold.txt

echo "SUCCESS"

