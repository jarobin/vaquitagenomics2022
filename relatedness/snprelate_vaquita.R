### Vaquita SNPRelate

library(gdsfmt)
library(SNPRelate)
library(gplots)

# Table with sample names, sample year, sex, etc.
info=read.table("vaquita_sample_info.txt", header=T, sep="\t")

# Create GDS file
vcf="vaquita_20_simple_PASS_autos_variants.vcf.gz"
gf="vaquita_20_simple_PASS_autos_variants.gds"

snpgdsVCF2GDS(vcf, gf, method="biallelic.only")
# VCF Format ==> SNP GDS Format
# Method: exacting biallelic SNPs
# Number of samples: 20
# Parsing "vaquita_20_simple_PASS_autos_variants.vcf.gz" ...
# 	import 466872 variants.
# + genotype   { Bit2 20x466872, 2.2M } *
# Optimize the access efficiency ...
# Clean up the fragments of GDS file:
#     open the file 'vaquita_20_simple_PASS_autos_variants.gds' (4.7M)
#     # of fragments: 97
#     save to 'vaquita_20_simple_PASS_autos_variants.gds.tmp'
#     rename 'vaquita_20_simple_PASS_autos_variants.gds.tmp' (4.7M, reduced: 924B)
#     # of fragments: 20

# Rename samples
genofile <- snpgdsOpen(gf, readonly=FALSE)
new_samp_id=paste(info$SAMPLE, info$YEAR, sep="_")
add.gdsn(genofile, "sample.id", new_samp_id, replace=TRUE, compress="LZMA_RA", closezip=TRUE)
snpgdsClose(genofile)

# Prune by LD
gf="vaquita_20_simple_PASS_autos_variants.gds"
genofile <- snpgdsOpen(gf)

set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=.5)
snpset.id <- unlist(snpset)
# SNP pruning based on LD:
# Excluding 0 SNP on non-autosomes
# Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# Working space: 20 samples, 466,872 SNPs
#     using 1 (CPU) core
#     sliding window: 500,000 basepairs, Inf SNPs
#     |LD| threshold: 0.5
#     method: composite
# Chromosome 1: 13.96%, 4,852/34,755
# Chromosome 2: 13.80%, 4,505/32,653
# Chromosome 3: 8.73%, 5,218/59,743
# Chromosome 4: 13.58%, 3,576/26,340
# Chromosome 5: 13.43%, 3,490/25,994
# Chromosome 6: 14.53%, 2,920/20,095
# Chromosome 7: 14.01%, 2,993/21,366
# Chromosome 8: 13.83%, 2,885/20,862
# Chromosome 9: 9.02%, 3,258/36,116
# Chromosome 10: 14.37%, 2,675/18,618
# Chromosome 11: 13.43%, 2,803/20,866
# Chromosome 12: 13.99%, 2,391/17,092
# Chromosome 13: 13.67%, 2,430/17,770
# Chromosome 14: 13.46%, 2,437/18,107
# Chromosome 15: 13.49%, 2,398/17,776
# Chromosome 16: 13.81%, 2,286/16,554
# Chromosome 17: 14.40%, 2,071/14,377
# Chromosome 18: 12.05%, 2,007/16,661
# Chromosome 19: 13.82%, 1,606/11,625
# Chromosome 20: 13.57%, 1,623/11,959
# Chromosome 21: 13.20%, 996/7,543
# 59,420 markers are selected in total.

gf.pruned="vaquita_20_simple_PASS_autos_variants_pruned.gds"
snpgdsCreateGenoSet(gf, gf.pruned, snp.id=snpset.id)
# Create a GDS genotype file:
# The new dataset consists of 20 samples and 59420 SNPs
#     write sample.id
#     write snp.id
#     write snp.rs.id
#     write snp.position
#     write snp.chromosome
#     write snp.allele
# SNP genotypes are stored in SNP-major mode (Sample X SNP).

snpgdsClose(genofile)


### ANALYSIS

# Pruned dataset
gf="vaquita_20_simple_PASS_autos_variants_pruned.gds"
# Full dataset
#gf="vaquita_20_simple_PASS_autos_variants.gds"

genofile <- snpgdsOpen(gf)

# Heatmap with dendrogram and scale
ibs <- snpgdsIBS(genofile, num.thread=2)

mycols=colorRampPalette(c("firebrick1", "darkorange", "gold1", "white"))(20)

pdf(paste("plot_", gf, "_IBSheatmap2_20210611.pdf", sep=""), width=5, height=5, pointsize=10)
heatmap.2(ibs$ibs, scale = "none", col = mycols, trace = "none", density.info = "none", labRow=ibs$sample.id, labCol=ibs$sample.id, margins=c(8,8), key.xlab="", key.title="", keysize=1.1, breaks=seq(.8,1,by=.01), key.par=list(mar=c(4, 2.5, 2.5, 1.5)), extrafun=box())
dev.off()


# KING
ibd.homo <- snpgdsIBDKING(genofile, type="KING-homo")
dat1 <- snpgdsIBDSelection(ibd.homo)

#ibd.robust <- snpgdsIBDKING(genofile, type="KING-robust")
#dat2 <- snpgdsIBDSelection(ibd.robust)

pdf(paste("plot_", gf, "_KING-homo.pdf", sep=""), width=5, height=5, pointsize=10)
par(mar=c(4.5,4.5,1,1))
plot(dat1$k0, dat1$kinship, xlab="Pr(IBD=0)", ylab="Estimated Kinship Coefficient", main="", ylim=c(-.02,.25))
dev.off()


# Relatedness stats with pruned dataset
# 4 pairs with kinship coefficient around .25 (first-degree relatives, e.g. parent-offspring, full sib)
mean(dat1[which(dat1$kinship>.2),]$kinship)
# 0.2310487
# Next highest relatedness: 0.07919926 between z0004380 and z0004382
# Mean and standard error kinship coefficient excluding first-degree relationships
mean(dat1[which(dat1$kinship<.2),]$kinship)
# 0.007738109
se <- function(x) sqrt(var(x)/length(x))
se(dat1[which(dat1$kinship<.2),]$kinship)
# 0.00107724


# Relatedness stats with full dataset
# 4 pairs with kinship coefficient around .25 (first-degree relatives, e.g. parent-offspring, full sib)
mean(dat1[which(dat1$kinship>.2),]$kinship)
# 0.244342
# Next highest relatedness: 0.07817300 between z0004380 and z0004382
# Mean and standard error kinship coefficient excluding first-degree relationships
mean(dat1[which(dat1$kinship<.2),]$kinship)
# -0.009822925
se <- function(x) sqrt(var(x)/length(x))
se(dat1[which(dat1$kinship<.2),]$kinship)
# 0.001407837



