# Divergence and mutation rate estimation for vaquita

# Morin et al. (2020) used a mutation rate of 1.08e-8/site/gen, based on the estimated 
# nuclear mutation rate of 9.10e-10/site/year in odontocetes (Dornburg et al., 2012) and 
# a generation time in vaquitas of 11.9 years.

# Here we calculate a new estimate using genome-wide divergence

# Harbor porpoise (Phocoena phocoena) and finless porpoise (Neophocaena phocaenoides) 
# aligned to vaquita reference genome

# Minor filters applied:
# - Repeats masked
# - Autosomes only
# - Min depth 1/3x mean, max depth 2x mean
# - Min QD of 4
# - SNPs and monomorphic REF calls only


################################################################################
# Genotype stats

# Harbor porpoise
HOMREF1=1225065884
HOMALT1=4295333
HET1=1862503
TOTAL1=1231223720

# Heterozygosity
PI1=HET1/TOTAL1
# 0.001512725

# Divergence (number of alt alleles in hets and hom. alt divided by the total number of 
# alleles)
DIV1=(HET1+(2*HOMALT1))/(2*TOTAL1)
# 0.004245032


# Finless porpoise
HOMREF2=1159799038
HOMALT2=5641882
HET2=1102535
TOTAL2=1166543455

# Heterozygosity
PI2=HET2/TOTAL2
# 0.0009451298

# Divergence (number of alt alleles divided by the total number of alleles)
DIV2=(HET2+(2*HOMALT2))/(2*TOTAL2)
# 0.005308975


################################################################################
# Mutation rate estimation

# Generation time (in years)
G=11.9


# CHEHIDA ET AL. 2020 SPLIT TIME
Split=5.42e6

T=Split/G

# Assuming u=DIV/2T
# Here, we ignore the time to coalescence in the ancestral population, which is okay when 
# 2Tu >> pi

# Harbor porpoise
DIV1/(2*T)
# 4.660137e-09

# Finless porpoise
DIV2/(2*T)
# 5.828118e-09

# Assuming u=(DIV-PI)/2T
# Here, we incorporate time to coalescence in the ancestral population, and we assume 
# that current pi = ancestral pi

# Harbor porpoise
(DIV1-PI1)/(2*T)
# 2.999489e-09

# Finless porpoise
(DIV2-PI2)/(2*T)
# 4.790568e-09


# TIMETREE SPLIT TIMES
# Median times from TimeTree.org (two divergence times, one for each porpoise with vaquita)
Split1=7.33e6
Split2=8.15e6

T1=Split1/G
T2=Split1/G

# Assuming u=DIV/2T
# Here, we ignore the time to coalescence in the ancestral population, which is okay when 
# 2Tu >> pi

# Harbor porpoise
DIV1/(2*T1)
# 3.445831e-09

# Finless porpoise
DIV2/(2*T2)
# 4.309468e-09

# Assuming u=(DIV-PI)/2T
# Here, we incorporate time to coalescence in the ancestral population, and we assume 
# that current pi = ancestral pi

# Harbor porpoise
(DIV1-PI1)/(2*T1)
# 2.217903e-09

# Finless porpoise
(DIV2-PI2)/(2*T2)
# 3.542275e-09


################################################################################
# Plot mutation rates

# Rates
hp=c(4.66e-09, 3.00e-09, 3.45e-09, 2.22e-09)
fp=c(5.83e-09, 4.79e-09, 4.31e-09, 3.54e-09)

setwd("~/Vaquita/analysis/polarization_mutationrate")
pdf("vaquita_divergence_mutationrate.pdf", width=5, height=3, pointsize=10)

par(mfrow=c(1,2))
par(mar=c(5,4.5,1,1))

plot(c(rep(1, length(hp)), rep(2, length(fp))), c(hp, fp), xlim=c(0.5,2.5), 
ylim=c(1e-9, 7e-9), pch=c(1,16,2,17,1,16,2,17), xlab="", 
ylab=expression(paste(mu,"/site/generation", sep="")), axes=F)

axis(side=2)
axis(side=1, labels=F, at=c(1,2))
par(xpd=T)
text(c(1,2)+.25, rep(1e-10, 2), c("Harbor\nporpoise","Finless\nporpoise"), srt=45, adj=1)
box()
par(mar=c(5,1,1,1))
par(xpd=T)
plot(0,0, type="n", axes=F, xlab="", ylab="")
l=legend("left", legend=c("Chehida et al., 2020", "TimeTree (median)"), bty="n")

points(rep(l$text$x, 2)-c(.3,.1), rev(sort(rep(l$text$y, 2))), pch=c(1,16,2,17))

dev.off()

