#! /bin/bash
#$ -wd /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs
#$ -l h_rt=23:00:00,h_data=10G,h_vmem=15G,highp
#$ -o /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/Blast_to_zebrafish.out.txt
#$ -e /u/project/rwayne/jarobins/vaquita/Phsi_neutral_sfs/reports/Blast_to_zebrafish.err.txt
#$ -pe shared 4
#$ -m abe

# Identifies regions that match to conserved regions (aligns to zebra fish genome)
# Adapted from Annabel Beichman's script (to analyze exomes) by Sergio Nigenda to analyze whole genome data.
# Usage: ./Identify_Neutral_Regions.sh [reference species]

source /u/local/Modules/default/init/modules.sh
module load bedtools

set -oe pipefail

########## Set variables, files and directories

REF=${1}

HOMEDIR=/u/project/rwayne/jarobins/vaquita
WORKDIR=${HOMEDIR}/Phsi_neutral_sfs/Neutral_Regions
OUTDIR=${WORKDIR}/neutralVCFs
REFDIR=${HOMEDIR}/reference_genomes/vaquita_genome/mPhoSin1.pri.cur.20190723_rename
REFERENCE=${REFDIR}/mPhoSin1.pri.cur.20190723_rename.fasta
ZFISHDIR=/u/project/rwayne/snigenda/finwhale/Other_genomes/zebraFish_genome
BLAST=/u/project/rwayne/software/CAPTURE/ncbi-blast-2.7.1+

# script to get gc content
getGC=/u/project/rwayne/snigenda/finwhale/scripts/Other/get_gc_content.pl


##### make directories were this information will be stored

mkdir -p ${WORKDIR}/get_fasta
mkdir -p ${WORKDIR}/GC_Content
mkdir -p ${WORKDIR}/zebra_fish


##### Logging

cd ${OUTDIR}

mkdir -p ./logs
mkdir -p ./temp

echo "[$(date "+%Y-%m-%d %T")] Start creating fasta file for ${REF} JOB_ID: ${JOB_ID}"
echo "The qsub input"
echo "${REF} ${JOB_ID}"


PROGRESSLOG=${OUTDIR}/logs/creates_fasta_and_blasting_progress.log
echo -e "[$(date "+%Y-%m-%d %T")] JOB_ID: ${JOB_ID}" > ${PROGRESSLOG}


########## Get fasta sequences

echo -e "[$(date "+%Y-%m-%d %T")] creating fasta files ... " >> ${PROGRESSLOG}
LOG=${OUTDIR}/logs/Step3_creating_fasta_and_blasting.log
date "+%Y-%m-%d %T" > ${LOG}

bedtools getfasta -fi ${REFERENCE} -bed ${WORKDIR}/get_fasta/min10kb_DistFromCDRs_noCpGIsland_noRepeat_mergedMaxDistance10_sorted.bed -fo ${WORKDIR}/get_fasta/min10kb_DistFromCDRs_noCpGIsland_noRepeat.fasta


########## Get GC content of each part of Fasta (exclude if >50%?) ##############

## for now not filtering on this; just generating it for interest. Not going to filter further since I already got rid of CpG islands.
perl ${getGC} ${WORKDIR}/get_fasta/min10kb_DistFromCDRs_noCpGIsland_noRepeat.fasta > ${WORKDIR}/GC_Content/min10kb_DistFromCDRs_noCpGIsland_noRepeat_GCcontent.txt

mv ${OUTDIR}/gc_out.txt ${WORKDIR}/GC_Content

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}


########## Blast against zebra fish genome to look for conservation #############

echo -e "[$(date "+%Y-%m-%d %T")]  Blasting ... " >> ${PROGRESSLOG}
date "+%Y-%m-%d %T" > ${LOG}

# do this once (I already downloaded the zebra fish genome and created a data base, so the lines are commented out):
cd ${ZFISHDIR}
#wget ftp://ftp.ensembl.org/pub/release-99/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz
#gunzip Danio_rerio.GRCz11.dna.toplevel.fa.gz
# ${BLAST}/bin/makeblastdb -in Danio_rerio.GRCz11.dna.toplevel.fa -out Dare_blastdb -dbtype nucl
${BLAST}/bin/blastn -query ${WORKDIR}/get_fasta/min10kb_DistFromCDRs_noCpGIsland_noRepeat.fasta -db ${ZFISHDIR}/Dare_blastdb -outfmt 7 -num_threads 4 > ${WORKDIR}/zebra_fish/neutralRegions_Blast_ZebraFish_blastn.out
# based on output, get regions with e-value < 1e-10 to exclude. You are getting their coordinates from their fasta name, not from the blast output
# so it is still 0-based even though blast output is 1-based.
cd ${WORKDIR}/zebra_fish
grep -v "#"  ${WORKDIR}/zebra_fish/neutralRegions_Blast_ZebraFish_blastn.out | awk '{if($11<1e-10)print $1}' | awk -F"[:-]" '{OFS="\t"; print $1,$2,$3}' | sort | uniq > ${WORKDIR}/zebra_fish/fish.matches.eval.1e-10.0based.bed

exitVal=${?}
if [ ${exitVal} -ne 0 ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL" >> ${PROGRESSLOG}
    exit 1
fi

date "+%Y-%m-%d %T" >> ${LOG}
echo -e "[$(date "+%Y-%m-%d %T")] Done" >> ${PROGRESSLOG}

