See `vaquita_reads_to_variants.sh` and `cetaceans_reads_to_variants.sh` for the overall pipelines used to generate VCF files from fastq files. Note that an additional set of files using the blue whale as an outgroup reference genome was also generated, following the same pipelines.

`SIFT_database_creation.sh` contains the pipeline for building a custom SIFT database for a given reference genome. Once generated, the database can be used to annotate VCF files with the SIFT annotator, as done in the reads_to_variants scripts.

`SIFT_database_creation_submission_log.sh` contains run info for each SIFT database created.

`vaquita_filterVCF.py` and `vaquita_generate_SNP_density_mask.sh` are custom filtering scripts used for vaquita VCF files.

`cetaceans_filterVCF.py` is a custom filtering script used for non-vaquita cetacean VCF files.
