`extract_variant_coordinates_by_impact.py` gets coordinates of mutations by class/impact (CDS, missense, synonymous, deleterious, tolerated, LOF) from annotated VCF files. The resulting bed files can then be used to extract the sites of interest from the VCF files to create subsets for further analysis.

`tally_genotypes_per_individual.py` gets counts of genotypes per sample from a VCF file.

`PGLS_plots.R` contains the steps for conducting phylogenetic generalized least squares (PGLS) regressions and generating plots of numbers/proportions of deleterious variants.
