### Phylogenetic generalized least squares (PGLS)

### Requirements
# - Species phylogeny with branch lengths
# - R packages: ape, geiger, nlme

### PGLS in R

### Load packages
library(ape)
library(geiger)
library(nlme)

### Load data file (various loading and pre-formatting steps not shown)
# Objective: create a dataframe with the species names and data values 
# Include a column for common species name and a column for scientific species name.
# In this case, the data are genome-wide het. and numbers of het. deleterious mutations in
#  different categories (eg. LOF, deleterious nonsynonymous, etc.).
# Counts have been normalized to account for varying call rates by multiplying all 
$ proportions (i.e., heterozygosity values) by the mean call rate across all samples.

### Exclude Indo-Pacific finless porpoise because it's not in the tree
df3=df2[which(df2$sample!="npho1"),]
dd=df3[,3:9]

### Read in tree file
mytree=read.tree("T2_ml.nwk")

### Drop the outgroup from the tree and check that species names in dataframe match those in tree file
# Note: the check is between names in the tree file and the row names of the dataframe
row.names(dd)=df3$species_sci
mytree2=drop.tip(mytree, "bos_taurus")
name.check(mytree2, dd)
# OK

### Regression between genome-wide het. and LOF heterozygotes
lof_pglsModel=gls(lof_het ~ gw_het, correlation=corBrownian(phy=mytree2, form=~species_sci), data=dd, method="ML")

summary(lof_pglsModel)
# Generalized least squares fit by maximum likelihood
#   Model: lof_het ~ gw_het 
#   Data: dd 
#        AIC      BIC    logLik
#   93.92914 95.12283 -43.96457
# 
# Correlation Structure: corBrownian
#  Formula: ~species_sci 
#  Parameter estimate(s):
# numeric(0)
# 
# Coefficients:
#                Value Std.Error   t-value p-value
# (Intercept)    13.06    12.871  1.014718  0.3367
# gw_het      92952.75  6851.895 13.565991  0.0000
# 
#  Correlation: 
#        (Intr)
# gw_het -0.589
# 
# Standardized residuals:
#          Min           Q1          Med           Q3          Max 
# -0.979854318 -0.383300839  0.003588461  0.467718960  1.531499554 
# 
# Residual standard error: 18.28713 
# Degrees of freedom: 11 total; 9 residual

### Get full p-values
as.data.frame(summary(lof_pglsModel)$tTable)
#                   Value  Std.Error   t-value      p-value
# (Intercept)    13.06043   12.87099  1.014718 3.367411e-01
# gw_het      92952.75172 6851.89528 13.565991 2.691234e-07

### Plot and add regression line (including Indo-Pacific finless porpoise in the plot)
plot(df2$gw_het, df2$lof_het, type="n", xlim=c(0,.0025), ylim=c(0, 250), xlab="Genome-wide heterozygosity", ylab="Number of loss-of-function heterozygotes")
abline(a=coef(lof_pglsModel)[1], b=coef(lof_pglsModel)[2], col="grey")
points(df2$gw_het, df2$lof_het, pch=16)
text(df2$gw_het, df2$lof_het, df2$species)

### Generate predicted values from the model
predict(lof_pglsModel, data.frame(gw_het=1e-3))
# 106.0132

predict(lof_pglsModel, data.frame(gw_het=2e-3))
# 198.9659


### Repeat for other mutation classes (deleterious nonsynonymous, tolerated nonsynonymous, synonymous)
