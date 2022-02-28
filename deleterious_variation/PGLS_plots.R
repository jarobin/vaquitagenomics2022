# Deleterious variation in vaquita and cetacean genomes
# Phylogenetic generalized least squares (PGLS), plots

### Requirements
# - Species phylogeny with branch lengths
# - R packages: ape, geiger, nlme

# Load functions and packages for phylogenetic generalized least squares (PGLS)
`%nin%`=Negate(`%in%`)

library(ape)
library(geiger)
library(nlme)


# Get input data
setwd("~/Vaquita/analysis/GTcounts")
# Cetaceans, same-species ref
cs=read.table("cetaceans_samespecies_ref_GTcounts.txt", header=T, sep="\t")
# Cetaceans, blue whale ref
cd=read.table("cetaceans_bluewhale_ref_GTcounts.txt", header=T, sep="\t")
# Vaquitas, same-species ref
vs=read.table("vaquitas_samespecies_ref_joint_GTcounts.txt", header=T, sep="\t")
# Vaquitas, blue whale ref
vd=read.table("vaquitas_bluewhale_ref_GTcounts.txt", header=T, sep="\t")


# Assign common and scientific species names for same-species ref. dataset
sample=c(cs[which(cs$type=="ALL"),]$sample, "vaquita")

species=sample
species[grep("oorc", sample)]="orca"
species[grep("bacu", sample)]="minke whale"
species[grep("dleu", sample)]="beluga whale"
species[grep("pmac", sample)]="sperm whale"
species[grep("nasi", sample)]="Yangtze finless porpoise"
species[grep("npho", sample)]="Indo-Pacific finless porpoise"
species[grep("lobl", sample)]="Pacific white-sided dolphin"
species[grep("mmon", sample)]="narwhal"
species[grep("gmel", sample)]="long-finned pilot whale"
species[grep("bmus", sample)]="blue whale"
species[grep("ttru", sample)]="bottlenose dolphin"
species[grep("vaquita", sample)]="vaquita"

species_sci=sample
species_sci[grep("oorc", sample)]="orcinus_orca"
species_sci[grep("bacu", sample)]="balaenoptera_acutorostrata"
species_sci[grep("dleu", sample)]="delphinapterus_leucas"
species_sci[grep("pmac", sample)]="physeter_macrocephalus"
species_sci[grep("nasi", sample)]="neophocaena_asiaeorientalis"
species_sci[grep("npho", sample)]="neophocaena_phocaenoides"
species_sci[grep("lobl", sample)]="lagenorhynchus_obliquidens"
species_sci[grep("mmon", sample)]="monodon_monoceros"
species_sci[grep("gmel", sample)]="globicephala_melas"
species_sci[grep("bmus", sample)]="balaenoptera_musculus"
species_sci[grep("ttru", sample)]="tursiops_truncatus"
species_sci[grep("vaquita", sample)]="phocoena_sinus"


# Genome-wide het. (same-species reference)
gw_het=cs[which(cs$type=="ALL"),]$het/cs[which(cs$type=="ALL"),]$call
gw_het=c(gw_het, mean(vs[which(vs$type=="ALL"),]$het/vs[which(vs$type=="ALL"),]$call))


# Proportion of het. variants relative to synonymous het.
het_mis_syn=cs[which(cs$type=="MIS"),]$het/cs[which(cs$type=="SYN"),]$het
het_mis_syn=c(het_mis_syn, mean(vs[which(vs$type=="MIS"),]$het/vs[which(vs$type=="SYN"),]$het))

het_tol_syn=cs[which(cs$type=="TOL"),]$het/cs[which(cs$type=="SYN"),]$het
het_tol_syn=c(het_tol_syn, mean(vs[which(vs$type=="TOL"),]$het/vs[which(vs$type=="SYN"),]$het))

het_del_syn=cs[which(cs$type=="DEL"),]$het/cs[which(cs$type=="SYN"),]$het
het_del_syn=c(het_del_syn, mean(vs[which(vs$type=="DEL"),]$het/vs[which(vs$type=="SYN"),]$het))

het_lof_syn=cs[which(cs$type=="LOF"),]$het/cs[which(cs$type=="SYN"),]$het
het_lof_syn=c(het_lof_syn, mean(vs[which(vs$type=="LOF"),]$het/vs[which(vs$type=="SYN"),]$het))

df_het_p=data.frame(sample, species, species_sci, gw_het, het_mis_syn, het_tol_syn, het_del_syn, het_lof_syn)

rm(het_mis_syn)
rm(het_tol_syn)
rm(het_del_syn)
rm(het_lof_syn)


# Number of het. variants (normalize by mean number of calls across samples in exons/conserved noncoding regions)
calls_exon=c(cs[which(cs$type=="EXON"),]$call, mean(vs[which(vs$type=="EXON"),]$call))
mean_calls_exon=mean(calls_exon)

calls_cnc=c(cs[which(cs$type=="CNC"),]$call, mean(vs[which(vs$type=="CNC"),]$call))
mean_calls_cnc=mean(calls_cnc)

het_syn=mean_calls_exon*c(cs[which(cs$type=="SYN"),]$het/cs[which(cs$type=="EXON"),]$call, mean(vs[which(vs$type=="SYN"),]$het/vs[which(vs$type=="EXON"),]$call))

het_mis=mean_calls_exon*c(cs[which(cs$type=="MIS"),]$het/cs[which(cs$type=="EXON"),]$call, mean(vs[which(vs$type=="MIS"),]$het/vs[which(vs$type=="EXON"),]$call))

het_tol=mean_calls_exon*c(cs[which(cs$type=="TOL"),]$het/cs[which(cs$type=="EXON"),]$call, mean(vs[which(vs$type=="TOL"),]$het/vs[which(vs$type=="EXON"),]$call))

het_del=mean_calls_exon*c(cs[which(cs$type=="DEL"),]$het/cs[which(cs$type=="EXON"),]$call, mean(vs[which(vs$type=="DEL"),]$het/vs[which(vs$type=="EXON"),]$call))

het_lof=mean_calls_exon*c(cs[which(cs$type=="LOF"),]$het/cs[which(cs$type=="EXON"),]$call, mean(vs[which(vs$type=="LOF"),]$het/vs[which(vs$type=="EXON"),]$call))

het_cnc=mean_calls_cnc*c(cs[which(cs$type=="CNC"),]$het/cs[which(cs$type=="CNC"),]$call, mean(vs[which(vs$type=="CNC"),]$het/vs[which(vs$type=="CNC"),]$call))

df_het_n=data.frame(sample, species, species_sci, gw_het, het_syn, het_mis, het_tol, het_del, het_lof, het_cnc)

rm(het_syn)
rm(het_mis)
rm(het_tol)
rm(het_del)
rm(het_lof)
rm(het_cnc)


### Blue whale reference - for analysis of derived homozygotes

# Assign common and scientific species names for blue whale ref. dataset
sample=c(cd[which(cd$type=="ALL"),]$sample, "vaquita")

species=sample
species[grep("oorc", sample)]="orca"
#species[grep("bacu", sample)]="minke whale"
species[grep("dleu", sample)]="beluga whale"
#species[grep("pmac", sample)]="sperm whale"
species[grep("nasi", sample)]="Yangtze finless porpoise"
species[grep("npho", sample)]="Indo-Pacific finless porpoise"
species[grep("lobl", sample)]="Pacific white-sided dolphin"
species[grep("mmon", sample)]="narwhal"
species[grep("gmel", sample)]="long-finned pilot whale"
#species[grep("bmus", sample)]="blue whale"
species[grep("ttru", sample)]="bottlenose dolphin"
species[grep("vaquita", sample)]="vaquita"

species_sci=sample
species_sci[grep("oorc", sample)]="orcinus_orca"
#species_sci[grep("bacu", sample)]="balaenoptera_acutorostrata"
species_sci[grep("dleu", sample)]="delphinapterus_leucas"
#species_sci[grep("pmac", sample)]="physeter_macrocephalus"
species_sci[grep("nasi", sample)]="neophocaena_asiaeorientalis"
species_sci[grep("npho", sample)]="neophocaena_phocaenoides"
species_sci[grep("lobl", sample)]="lagenorhynchus_obliquidens"
species_sci[grep("mmon", sample)]="monodon_monoceros"
species_sci[grep("gmel", sample)]="globicephala_melas"
#species_sci[grep("bmus", sample)]="balaenoptera_musculus"
species_sci[grep("ttru", sample)]="tursiops_truncatus"
species_sci[grep("vaquita", sample)]="phocoena_sinus"


# Proportion of del. hom. variants relative to synonymous hom.
hom_mis_syn=cd[which(cd$type=="MIS"),]$homA/cd[which(cd$type=="SYN"),]$homA
hom_mis_syn=c(hom_mis_syn, mean(vd[which(vd$type=="MIS"),]$homA/vd[which(vd$type=="SYN"),]$homA))

hom_tol_syn=cd[which(cd$type=="TOL"),]$homA/cd[which(cd$type=="SYN"),]$homA
hom_tol_syn=c(hom_tol_syn, mean(vd[which(vd$type=="TOL"),]$homA/vd[which(vd$type=="SYN"),]$homA))

hom_del_syn=cd[which(cd$type=="DEL"),]$homA/cd[which(cd$type=="SYN"),]$homA
hom_del_syn=c(hom_del_syn, mean(vd[which(vd$type=="DEL"),]$homA/vd[which(vd$type=="SYN"),]$homA))

hom_lof_syn=cd[which(cd$type=="LOF"),]$homA/cd[which(cd$type=="SYN"),]$homA
hom_lof_syn=c(hom_lof_syn, mean(vd[which(vd$type=="LOF"),]$homA/vd[which(vd$type=="SYN"),]$homA))

gw_het=NULL
for(i in 1:length(sample)){gw_het[i]=df_het_n[which(df_het_n$sample==sample[i]),]$gw_het}

df_hom_p=data.frame(sample, species, species_sci, gw_het, hom_mis_syn, hom_tol_syn, hom_del_syn, hom_lof_syn)

rm(hom_mis_syn)
rm(hom_tol_syn)
rm(hom_del_syn)
rm(hom_lof_syn)


################################################################################
### PGLS

# Clear variables
rm(sample)
rm(species)
rm(species_sci)
rm(gw_het)

# Get tree
setwd("~/Vaquita/analysis/cetaceans/phylo")
mytree=read.tree("T2_ml.nwk")

# Drop cow from full tree
mytree2=drop.tip(mytree, "bos_taurus")

# Drop baleen whales, sperm whale for homozygote analyses (blue whale ref)
mytree3=drop.tip(mytree2, c("balaenoptera_musculus", "balaenoptera_acutorostrata", "physeter_macrocephalus") )

# Exclude I-P finless porpoise, which is not in the tree
df_het_p2=df_het_p[which(df_het_p$sample!="npho1"),]
df_het_n2=df_het_n[which(df_het_n$sample!="npho1"),]
df_hom_p2=df_hom_p[which(df_hom_p$sample!="npho1"),]

# Format data frame row names and check for consistency with the trees
row.names(df_het_p2)=df_het_p2$species_sci
row.names(df_het_n2)=df_het_n2$species_sci
row.names(df_hom_p2)=df_hom_p2$species_sci

name.check(mytree2, df_het_p2)
name.check(mytree2, df_het_n2)
name.check(mytree3, df_hom_p2)


### Regressions

### Het. proportion (normalized by syn) PGLS models (Fig. 2A,B)
pgls_het_del_syn=gls(het_del_syn ~ log10(gw_het), correlation=corBrownian(phy=mytree2, form=~species_sci), data=df_het_p2, method="ML")
pgls_het_lof_syn=gls(het_lof_syn ~ log10(gw_het), correlation=corBrownian(phy=mytree2, form=~species_sci), data=df_het_p2, method="ML")

# het. del/syn
as.data.frame(summary(pgls_het_del_syn)$tTable)
#                    Value  Std.Error    t-value    p-value
# (Intercept)   -0.1117674 0.12518016 -0.8928522 0.39518242
# log10(gw_het) -0.1186312 0.03853193 -3.0787769 0.01316595

# het. lof/syn
as.data.frame(summary(pgls_het_lof_syn)$tTable)
#                      Value   Std.Error    t-value     p-value
# (Intercept)   -0.005147690 0.006895077 -0.7465748 0.474370204
# log10(gw_het) -0.007215444 0.002122386 -3.3996847 0.007877262


### Het. number PGLS models (Fig. 2 C,D, Fig. S10)
pgls_het_del=gls(het_del ~ gw_het, correlation=corBrownian(phy=mytree2, form=~species_sci), data=df_het_n2, method="ML")
pgls_het_lof=gls(het_lof ~ gw_het, correlation=corBrownian(phy=mytree2, form=~species_sci), data=df_het_n2, method="ML")

pgls_het_syn=gls(het_syn ~ gw_het, correlation=corBrownian(phy=mytree2, form=~species_sci), data=df_het_n2, method="ML")
pgls_het_tol=gls(het_tol ~ gw_het, correlation=corBrownian(phy=mytree2, form=~species_sci), data=df_het_n2, method="ML")
pgls_het_cnc=gls(het_cnc ~ gw_het, correlation=corBrownian(phy=mytree2, form=~species_sci), data=df_het_n2, method="ML")

# het. del
as.data.frame(summary(pgls_het_del)$tTable)
#                    Value   Std.Error  t-value      p-value
# (Intercept)     459.0347    288.2721 1.592366 1.457664e-01
# gw_het      1450122.1058 152969.9115 9.479787 5.571798e-06

# het. lof
as.data.frame(summary(pgls_het_lof)$tTable)
#                   Value   Std.Error  t-value      p-value
# (Intercept)    37.43021    21.16528 1.768473 1.107730e-01
# gw_het      91523.32011 11231.23275 8.149000 1.909485e-05

# het. syn
as.data.frame(summary(pgls_het_syn)$tTable)
#                   Value   Std.Error   t-value      p-value
# (Intercept)    1127.309    702.9857  1.603601 1.432653e-01
# gw_het      6773829.411 373035.3085 18.158682 2.124516e-08

# het. tol
as.data.frame(summary(pgls_het_tol)$tTable)
#                   Value   Std.Error   t-value      p-value
# (Intercept)    1387.458    468.3604  2.962374 1.589810e-02
# gw_het      3819513.077 248532.7366 15.368249 9.136722e-08

# het. cnc
as.data.frame(summary(pgls_het_cnc)$tTable)
#                    Value   Std.Error   t-value      p-value
# (Intercept)     1487.149    672.7105  2.210682 5.438467e-02
# gw_het      15127844.743 356969.9574 42.378482 1.131568e-11


### Hom. proportion (normalized by syn) PGLS models (Fig. S8)
pgls_hom_del_syn=gls(hom_del_syn ~ log10(gw_het), correlation=corBrownian(phy=mytree3, form=~species_sci), data=df_hom_p2, method="ML")
pgls_hom_lof_syn=gls(hom_lof_syn ~ log10(gw_het), correlation=corBrownian(phy=mytree3, form=~species_sci), data=df_hom_p2, method="ML")
pgls_hom_tol_syn=gls(hom_tol_syn ~ log10(gw_het), correlation=corBrownian(phy=mytree3, form=~species_sci), data=df_hom_p2, method="ML")

# hom. del/syn
as.data.frame(summary(pgls_hom_del_syn)$tTable)
#                      Value   Std.Error   t-value      p-value
# (Intercept)    0.115688844 0.003821484 30.273276 8.620003e-08
# log10(gw_het) -0.005477762 0.001136965 -4.817881 2.946436e-03

# hom. lof/syn
as.data.frame(summary(pgls_hom_lof_syn)$tTable)
#                       Value    Std.Error   t-value      p-value
# (Intercept)    0.0116312434 0.0005065265 22.962753 4.469503e-07
# log10(gw_het) -0.0002690256 0.0001507014 -1.785157 1.244848e-01

# hom. tol/syn
as.data.frame(summary(pgls_hom_tol_syn)$tTable)
#                     Value   Std.Error   t-value      p-value
# (Intercept)   0.721860667 0.010632450 67.892220 6.869113e-10
# log10(gw_het) 0.003203002 0.003163358  1.012532 3.503598e-01


################################################################################
### Plots

# Function for axis ticks on log plots
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
 lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
 major.ticks <- pretty(lims,n=2)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
 major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
 labels <- sapply(major.ticks,function(i)
            as.expression(bquote(10^ .(i)))
          )
  axis(ax,at=major.ticks,labels=F)
  axis(ax,at=major.ticks,labels=labels,tick=F, line=-0.4)
 n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}


################################################################################
# Figure 2
pdf("vaquita_fig2_delvar.pdf", width=7.5, height=7, pointsize=12)

par(mfrow=c(2,2))
par(mar=c(2,4,2,1))

# A) Del het. proportions (normalized by syn)
plot(log10(df_het_p$gw_het), df_het_p$het_del_syn, type="n", xlim=c(-4.1,-2.4), ylim=c(0.145, 0.345), xlab="", ylab="", axes=F, xaxt="n")
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(ylab="Deleterious nonsynonymous/synonymous\nheterozygotes", line=2.15)
minor.ticks.axis(1,9,mn=-5,mx=-2)
abline(a=coef(pgls_het_del_syn)[1], b=coef(pgls_het_del_syn)[2], col="grey")
points(log10(df_het_p$gw_het), df_het_p$het_del_syn, pch=16)
text(log10(df_het_p$gw_het), df_het_p$het_del_syn, df_het_p$species)
box()


# B) LOF het. proportions (normalized by syn)
plot(log10(df_het_p$gw_het), df_het_p$het_lof_syn, type="n", xlim=c(-4.1,-2.4), ylim=c(0.01, 0.023), xlab="", ylab="", axes=F, xaxt="n")
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(ylab="LOF/synonymous\nheterozygotes", line=2.15)
minor.ticks.axis(1,9,mn=-5,mx=-2)
abline(a=coef(pgls_het_lof_syn)[1], b=coef(pgls_het_lof_syn)[2], col="grey")
points(log10(df_het_p$gw_het), df_het_p$het_lof_syn, pch=16)
text(log10(df_het_p$gw_het), df_het_p$het_lof_syn, df_het_p$species)
box()


# C) Del het. number
plot(df_het_n$gw_het, df_het_n$het_del, type="n", xlim=c(0,.0025), ylim=c(0, 4000), xlab="", ylab="", axes=F)
axis(side=1, line=0, labels=F, tick=T)
axis(side=1, line=-0.4, labels=T, tick=F, at=seq(0,0.002, by=0.001))
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="Number of deleterious\nnonsynonymous heterozygotes", line=2.15)
abline(a=coef(pgls_het_del)[1], b=coef(pgls_het_del)[2], col="grey")
points(df_het_n$gw_het, df_het_n$het_del, pch=16)
text(df_het_n$gw_het, df_het_n$het_del, df_het_n$species)
box()


# D) LOF het. number
plot(df_het_n$gw_het, df_het_n$het_lof, type="n", xlim=c(0,.0025), ylim=c(0, 250), xlab="", ylab="", axes=F)
axis(side=1, line=0, labels=F, tick=T)
axis(side=1, line=-0.4, labels=T, tick=F, at=seq(0,0.002, by=0.001))
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="Number of\nLOF heterozygotes", line=2.15)
abline(a=coef(pgls_het_lof)[1], b=coef(pgls_het_lof)[2], col="grey")
points(df_het_n$gw_het, df_het_n$het_lof, pch=16)
text(df_het_n$gw_het, df_het_n$het_lof, df_het_n$species)
box()

dev.off()


################################################################################
# Fig. S7: Vaquita individual variability

# Get sample year
info=read.table("~/Vaquita/analysis/sample_sex_year.txt", header=T, sep='\t')

# het proportion (vaquita reference)
het_tol_syn=vs[which(vs$type=="TOL"),]$het/vs[which(vs$type=="SYN"),]$het
het_del_syn=vs[which(vs$type=="DEL"),]$het/vs[which(vs$type=="SYN"),]$het
het_lof_syn=vs[which(vs$type=="LOF"),]$het/vs[which(vs$type=="SYN"),]$het


# hom proportion (blue whale reference)
hom_tol_syn=vd[which(vd$type=="TOL"),]$homA/vd[which(vd$type=="SYN"),]$homA
hom_del_syn=vd[which(vd$type=="DEL"),]$homA/vd[which(vd$type=="SYN"),]$homA
hom_lof_syn=vd[which(vd$type=="LOF"),]$homA/vd[which(vd$type=="SYN"),]$homA


# het number (vaquita reference)
het_syn=mean_calls_exon*(vs[which(vs$type=="SYN"),]$het/vs[which(vs$type=="EXON"),]$call)
het_tol=mean_calls_exon*(vs[which(vs$type=="TOL"),]$het/vs[which(vs$type=="EXON"),]$call)
het_del=mean_calls_exon*(vs[which(vs$type=="DEL"),]$het/vs[which(vs$type=="EXON"),]$call)
het_lof=mean_calls_exon*(vs[which(vs$type=="LOF"),]$het/vs[which(vs$type=="EXON"),]$call)
het_cnc=mean_calls_cnc*(vs[which(vs$type=="CNC"),]$het/vs[which(vs$type=="CNC"),]$call)

sample=vs[which(vs$type=="ALL"),]$sample
year=NULL
for (i in 1:length(sample)){year[i]=info[which(info$sample==sample[i]),]$year}

df_v=data.frame(sample, year, het_tol_syn, het_del_syn, het_lof_syn, hom_tol_syn, hom_del_syn, hom_lof_syn, het_syn, het_tol, het_del, het_lof, het_cnc)

pdf("vaquita_del_by_year.pdf", width=7.5, height=8, pointsize=12)

par(mfrow=c(4,3))
par(mar=c(3,6,1,1))

# Row 1: Heterozygote ratios
plot(df_v$year, df_v$het_tol_syn, type="n", xlab="", ylab="Tolerated\nnonsynonymous/synonymous\nheterozygotes")
abline(lm(df_v$het_tol_syn ~ df_v$year), col="grey")
points(df_v$year, df_v$het_tol_syn)
plot(df_v$year, df_v$het_del_syn, type="n", xlab="", ylab="Deleterious\nnonsynonymous/synonymous\nheterozygotes")
abline(lm(df_v$het_del_syn ~ df_v$year), col="grey")
points(df_v$year, df_v$het_del_syn)
plot(df_v$year, df_v$het_lof_syn, type="n", xlab="", ylab="LOF/synonymous\nheterozygotes")
abline(lm(df_v$het_lof_syn ~ df_v$year), col="grey")
points(df_v$year, df_v$het_lof_syn)

# Row 2: Homozygote ratios
plot(df_v$year, df_v$hom_tol_syn, type="n", xlab="", ylab="Tolerated\nnonsynonymous/synonymous\nhomozygotes")
abline(lm(df_v$hom_tol_syn ~ df_v$year), col="grey")
points(df_v$year, df_v$hom_tol_syn)
plot(df_v$year, df_v$hom_del_syn, type="n", xlab="", ylab="Deleterious\nnonsynonymous/synonymous\nhomozygotes")
abline(lm(df_v$hom_del_syn ~ df_v$year), col="grey")
points(df_v$year, df_v$hom_del_syn)
plot(df_v$year, df_v$hom_lof_syn, type="n", xlab="", ylab="LOF/synonymous\nhomozygotes")
abline(lm(df_v$hom_lof_syn ~ df_v$year), col="grey")
points(df_v$year, df_v$hom_lof_syn)

# Rows 3-4: Heterozygote numbers
plot(df_v$year, df_v$het_tol, type="n", xlab="", ylab="Number of tolerated\nnonsynonymous heterozygotes")
abline(lm(df_v$het_tol ~ df_v$year), col="grey")
points(df_v$year, df_v$het_tol)
plot(df_v$year, df_v$het_del, type="n", xlab="", ylab="Number of deleterious\nnonsynonymous heterozygotes")
abline(lm(df_v$het_del ~ df_v$year), col="grey")
points(df_v$year, df_v$het_del)
plot(df_v$year, df_v$het_lof, type="n", xlab="", ylab="Number of LOF\nheterozygotes")
abline(lm(df_v$het_lof ~ df_v$year), col="grey")
points(df_v$year, df_v$het_lof)
plot(df_v$year, df_v$het_syn, type="n", xlab="", ylab="Number of synonymous\nheterozygotes")
abline(lm(df_v$het_syn ~ df_v$year), col="grey")
points(df_v$year, df_v$het_syn)
plot(df_v$year, df_v$het_cnc, type="n", xlab="", ylab="Number of conserved\nnoncoding heterozygotes")
abline(lm(df_v$het_cnc ~ df_v$year), col="grey")
points(df_v$year, df_v$het_cnc)

dev.off()


################################################################################
# Fig. S8: Homozygote ratios

pdf("vaquita_cetacean_hom_del_syn.pdf", width=7.5, height=7, pointsize=12)

par(mfrow=c(2,2))
par(mar=c(2,4,2,1))

# Data: df_hom_p

# A) Deleterious
plot(log10(df_hom_p$gw_het), df_hom_p$het_del_syn, type="n", xlim=c(-4.1,-2.4), ylim=c(0.128, 0.14), xlab="", ylab="", axes=F, xaxt="n")
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="Deleterious nonsynonymous/synonymous\nhomozygotes", line=2.15)
minor.ticks.axis(1,9,mn=-5,mx=-2)
abline(a=coef(pgls_hom_del_syn)[1], b=coef(pgls_hom_del_syn)[2], col="grey")
points(log10(df_hom_p$gw_het), df_hom_p$hom_del_syn, pch=16)
text(log10(df_hom_p$gw_het), df_hom_p$hom_del_syn, df_hom_p$species)
box()

# B) LOF
plot(log10(df_hom_p$gw_het), df_hom_p$het_lof_syn, type="n", xlim=c(-4.1,-2.4), ylim=c(0.0122, 0.01275), xlab="", ylab="", axes=F, xaxt="n")
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="LOF/synonymous\nhomozygotes", line=2.15)
minor.ticks.axis(1,9,mn=-5,mx=-2)
abline(a=coef(pgls_hom_lof_syn)[1], b=coef(pgls_hom_lof_syn)[2], col="grey")
points(log10(df_hom_p$gw_het), df_hom_p$hom_lof_syn, pch=16)
text(log10(df_hom_p$gw_het), df_hom_p$hom_lof_syn, df_hom_p$species)
box()

# C) Tolerated
plot(log10(df_hom_p$gw_het), df_hom_p$het_tol_syn, type="n", xlim=c(-4.1,-2.4), ylim=c(0.695, 0.72), xlab="", ylab="", axes=F, xaxt="n")
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="Tolerated nonsynonymous/synonymous\nhomozygotes", line=2.15)
minor.ticks.axis(1,9,mn=-5,mx=-2)
abline(a=coef(pgls_hom_tol_syn)[1], b=coef(pgls_hom_tol_syn)[2], col="grey")
points(log10(df_hom_p$gw_het), df_hom_p$hom_tol_syn, pch=16)
text(log10(df_hom_p$gw_het), df_hom_p$hom_tol_syn, df_hom_p$species)
box()

dev.off()


################################################################################
# Fig. S10: Number of heterozygotes

pdf("vaquita_cetacean_het_del_supp.pdf", width=7.5, height=7, pointsize=12)

par(mfrow=c(2,2))
par(mar=c(2,4,2,1))

# Data: df_het_n

# A) Synonymous
plot(df_het_n$gw_het, df_het_n$het_syn, type="n", xlim=c(0,.0025), ylim=c(0, 18000), xlab="", ylab="", axes=F)
axis(side=1, line=0, labels=F, tick=T)
axis(side=1, line=-0.4, labels=T, tick=F, at=seq(0,0.002, by=0.001))
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="Number of synonymous\nheterozygotes", line=2.15)
abline(a=coef(pgls_het_syn)[1], b=coef(pgls_het_syn)[2], col="grey")
points(df_het_n$gw_het, df_het_n$het_syn, pch=16)
text(df_het_n$gw_het, df_het_n$het_syn, df_het_n$species)
box()

# B) Tolerated
plot(df_het_n$gw_het, df_het_n$het_tol, type="n", xlim=c(0,.0025), ylim=c(0, 11000), xlab="", ylab="", axes=F)
axis(side=1, line=0, labels=F, tick=T)
axis(side=1, line=-0.4, labels=T, tick=F, at=seq(0,0.002, by=0.001))
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="Number of tolerated\nnonsynonymous heterozygotes", line=2.15)
abline(a=coef(pgls_het_tol)[1], b=coef(pgls_het_tol)[2], col="grey")
points(df_het_n$gw_het, df_het_n$het_tol, pch=16)
text(df_het_n$gw_het, df_het_n$het_tol, df_het_n$species)
box()

# C) Conserved noncoding
plot(df_het_n$gw_het, df_het_n$het_cnc, type="n", xlim=c(0,.0025), ylim=c(0, 35000), xlab="", ylab="", axes=F)
axis(side=1, line=0, labels=F, tick=T)
axis(side=1, line=-0.4, labels=T, tick=F, at=seq(0,0.002, by=0.001))
axis(side=2, line=0, labels=F, tick=T)
axis(side=2, line=-0.4, labels=T, tick=F)
title(xlab="Genome-wide heterozygosity", line=-2)
title(ylab="Number of conserved\nnoncoding heterozygotes", line=2.15)
abline(a=coef(pgls_het_cnc)[1], b=coef(pgls_het_cnc)[2], col="grey")
points(df_het_n$gw_het, df_het_n$het_cnc, pch=16)
text(df_het_n$gw_het, df_het_n$het_cnc, df_het_n$species)
box()

dev.off()
