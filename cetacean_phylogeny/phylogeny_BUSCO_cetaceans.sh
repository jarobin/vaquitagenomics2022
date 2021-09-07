### Cetacean phylogeny workflow and notes
# - Run BUSCO on each genome assembly
# 	- v4.1.4
# 	- vertebrata_odb10.2020-09-10
# 		- number of species: 67, number of BUSCOs: 3354
# 	- cow + 11 cetaceans (not using N. phocaenoides)
# - Get BUSCOs found in all 12 genomes (n=2459)
# - Generate a fasta file for each BUSCO containing the 12 sequences (cat sequences from each species together into one file per busco)
# - For each busco fasta file: align (mafft) then trim (trimal)
# - Concatenate trimmed alignments into supermatrix with partition info (catsequences)
# - Generate trees with and without partitioning (iqtree v2.1.2 with ModelFinder, UFBoot, partitioning)
# 	- Input: 12 taxa with 167 partitions and 4421551 total sites (0% missing data)
# 	- Edge-linked-proportional partition model with separate substitution models and separate rates across sites
# 	- 1000 bootstrap reps
# 	- Results: 
#     - Partitioned model has better likelihood, AIC/BIC scores
# 	  - Partitioned tree is larger (branch lengths are longer), but looks proportionally the same as non-partitioned phylogeny)
# 	  - Consensus and maximum likelihood trees identical
# 	  - All branches have 100% bootstrap support and topology consistent with the phylogeny in McGowen et al., 2020
# 
### Requirements
# - BUSCO (v4)
# - mafft
# - trimal
# - catsequences (https://github.com/ChrisCreevey/catsequences)
# - iqtree (v2)


### Run BUSCO on each genome (here, busco_species_fasta.list has two columns: species_name fasta_filename)
export LINEAGE=vertebrata_odb10
export AUG_SPECIES=human
while read -r SPECIES FASTA ; do
export SPECIES
export FASTA
busco -f --offline --cpu 16 -m genome --augustus_species=${AUG_SPECIES} -i ${FASTA} -o ${SPECIES}.busco.${LINEAGE} -l ${LINEAGE}
done < busco_species_fasta.list


### When complete, move/copy all the run_{species_name} results folders to a single directory
### In directory with all vertebrata runs for cow+11 cetaceans (not using N. phocaenoides)...

### Get species count for each busco
cat run_*/full_table.tsv | grep -v "^#" | awk '$2=="Complete"{print $1}' | sort | uniq -c > busco_ids_complete_counts.txt
wc -l busco_ids_complete_counts.txt
# 3302 busco_ids_complete_counts.txt


### Get buscos in all 12 species
awk '$1==12' busco_ids_complete_counts.txt > busco_ids_complete_12.list
wc -l busco_ids_complete_12.list 
# 2459 busco_ids_complete_12.list


### Per busco, concatenate sequences from each species into one file
mkdir phylo
while read -r COUNT ID ; do
for DIR in run_* ; do
SPECIES=$(echo ${DIR} | sed 's/run_//g')
echo ">${SPECIES}" >> phylo/${ID}.fna
grep -v "^>" ${DIR}/busco_sequences/single_copy_busco_sequences/${ID}.fna >> phylo/${ID}.fna
done
done < busco_ids_complete_12.list 


### Align and trim each busco .fna
cd phylo
for FA in *.fna ; do
ALN=${FA/fna/aln}
mafft --thread 16 --auto ${FA} > ${ALN}
trimal -in ${ALN} -out ${ALN}.trimmed -automated1
done


### Concatenate trimmed alignments into supermatrix with partition info
CAT=~/project/busco_phylo/catsequences/catsequences
ls *trimmed > trimmed_alignments.list
${CAT} trimmed_alignments.list
# Reformat partition output file for iqtree
cat allseqs.partitions.txt | sed 's/^/DNA, /g' | sed 's/;//g' | sed -e 's/[[:space:]]\+/ /g' > allseqs.partitions.txt_1


### Non-partitioned tree inference with bootstrap
iqtree -s allseqs.fas -B 1000 -T AUTO -ntmax 32 --prefix T1


### Partitioned tree inference with model & partition finding and bootstrap
iqtree -s allseqs.fas -p allseqs.partitions.txt_1 -m MFP+MERGE -B 1000 -T AUTO -ntmax 32 --prefix T2


### IQ-TREE citation info
# To cite IQ-TREE please use:
# 
# B.Q. Minh, H.A. Schmidt, O. Chernomor, D. Schrempf, M.D. Woodhams, A. von Haeseler, R. Lanfear (2020) IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., 37:1530-1534. https://doi.org/10.1093/molbev/msaa015
# 
# To cite ModelFinder please use: 
# 
# Subha Kalyaanamoorthy, Bui Quang Minh, Thomas KF Wong, Arndt von Haeseler,
# and Lars S Jermiin (2017) ModelFinder: Fast model selection for
# accurate phylogenetic estimates. Nature Methods, 14:587–589.
# https://doi.org/10.1038/nmeth.4285
# 
# Since you used ultrafast bootstrap (UFBoot) please also cite: 
# 
# Diep Thi Hoang, Olga Chernomor, Arndt von Haeseler, Bui Quang Minh,
# and Le Sy Vinh (2018) UFBoot2: Improving the ultrafast bootstrap
# approximation. Mol. Biol. Evol., 35:518–522.
# https://doi.org/10.1093/molbev/msx281
# 
# Since you used partition models please also cite:
# 
# Olga Chernomor, Arndt von Haeseler, and Bui Quang Minh (2016)
# Terrace aware data structure for phylogenomic inference from
# supermatrices. Syst. Biol., 65:997-1008.
# https://doi.org/10.1093/sysbio/syw037
	
