#!/usr/bin/env bash

# Download needed databases directly from ANNOVAR
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp35a humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_exome humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_ea humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
annotate_variation.pl -buildver hg19 -downdb -webfrom annovar intervar_20180118 humandb/



# Create a database for the latest COSMIC release
# Download http://www.openbioinformatics.org/annovar/download/prepare_annovar_user.pl and add to annovar folder
# Register at COSMIC
# Download:
# VCF/CosmicCodingMuts.vcf.gz
# VCF/CosmicNonCodingVariants.vcf.gz
# CosmicMutantExport.tsv.gz
# CosmicNCV.tsv.gz
# unzip all files and build database
# actual command to build the database for ANNOVAR
	# prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg19_cosmic87_coding.txt
	# prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.vcf > hg19_cosmic87_noncoding.txt
# Move both created files to the ~/annovar/humandb folder.
