# SIFT database creation for cetacean genomes
#
# Run time is typically one to several days, depending on the genome
# Best to run on a single multi-core machine
#
# For instructions, see:
# - Database creation pipeline: https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB
# - SIFT4G: https://github.com/rvaser/sift4g
#   - NOTE: Used modified code for sift4g as detailed in https://github.com/rvaser/sift4g/issues/10
#
# Used faSplit for convenience (available from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
#
# Used uniref90.fasta protein database, downloaded May 1 2020 (109653977 sequences)

### Set variables, paths
# EXAMPLE: Beluga whale, Delphinapterus leucas
export PARENT_DIR=~/programs/sift/SIFT_databases/GCF_002288925.2_ASM228892v3_makeSiftDB
export ORG=delphinapterus_leucas
export ORG_VERSION=GCF_002288925.2_ASM228892v3
export FASTA=GCF_002288925.2_ASM228892v3_genomic.fna.gz
export GTF=GCF_002288925.2_ASM228892v3_genomic.gtf.gz
export URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/288/925/GCF_002288925.2_ASM228892v3

export SIFT4G_PATH=~/programs/sift/sift4g_mod/bin/sift4g
export PROTEIN_DB=~/programs/sift/SIFT_databases/uniref90.fasta

### Set up PARENT_DIR
mkdir ${PARENT_DIR}
cd ${PARENT_DIR}
mkdir chr-src
mkdir dbSNP
mkdir gene-annotation-src

### Get genome fasta and gene annotation
wget ${URL}/${FASTA}
wget ${URL}/${GTF}

mv ${GTF} gene-annotation-src/

### Split genome into individual fasta records, then gzip each
# faSplit available from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
gzip -cd ${FASTA} > ${FASTA%.fna.gz}.fa
faSplit byname ${FASTA%.fna.gz}.fa chr-src/
for i in chr-src/*.fa ; do gzip ${i} ; done
rm ${FASTA%.fna.gz}.fa

### Generate config file
cat > config <<EOF
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=2
MITO_GENETIC_CODE_TABLENAME=Vertebrate Mitochondrial

PARENT_DIR=${PARENT_DIR}
ORG=${ORG}
ORG_VERSION=${ORG_VERSION}

#Running SIFT 4G
SIFT4G_PATH=${SIFT4G_PATH}
PROTEIN_DB=${PROTEIN_DB}

# Sub-directories, don't need to change
GENE_DOWNLOAD_DEST=gene-annotation-src
CHR_DOWNLOAD_DEST=chr-src
LOGFILE=Log.txt
ZLOGFILE=Log2.txt
FASTA_DIR=fasta
SUBST_DIR=subst
ALIGN_DIR=SIFT_alignments
SIFT_SCORE_DIR=SIFT_predictions
SINGLE_REC_BY_CHR_DIR=singleRecords
SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores
DBSNP_DIR=dbSNP

# Doesn't need to change
FASTA_LOG=fasta.log
INVALID_LOG=invalid.log
PEPTIDE_LOG=peptide.log
ENS_PATTERN=ENS
SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord

EOF

### Run make-SIFT-db-all.pl
cd ~/programs/sift/scripts_to_build_SIFT_db/
CONFIG=${PARENT_DIR}/config
LOG=${PARENT_DIR}/makeSiftDB.log
perl ./make-SIFT-db-all.pl -config ${CONFIG} |& tee ${LOG}

### Archive
cd ~/programs/sift/SIFT_databases
cp -r ${ORG_VERSION}_makeSiftDB/${ORG_VERSION} .
tar -czf ${ORG_VERSION}_makeSiftDB.tar.gz ${ORG_VERSION}_makeSiftDB
