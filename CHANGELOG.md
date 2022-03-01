# Changelog

## Release v4.0.0

* Integration of Mutect2 (GATK4) for tumorOnly and Panel variant calling
* Analysis of the TSO500 panel handels of the DNA and RNA part of the panel
* Calcualtion of various complex biomarkers like TMB, MSI, HRD,
* Calcualtion of more QC metrics, e.g. bioinformatic tumor cell content (purity)
* Various bugfixes
* Adjustments of the PDF report mainly structure and readability

## Release v3.1.0

* Updated ClinVar version to clinvar_20210123; re-run of ./setup.sh -t setup_tools from the MIRACUM-Pipe-docker repository recommended or manual installation of the updated annovar databases:

```bash
# annovar databases
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210123 humandb/
```

* Introduced duplicate removal for panel pipeline
* Implemented advanced qc and coverage metrics
* Included RNA Fusion detection for panels containing RNA data

## Minor Updates

* Update of several annovar databases; re-start of ./setup.sh -t setup_tools from the MIRACUM-Pipe-docker repository recommended or manual installation of the updated annovar databases:

```bash
# annovar databases
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp41a humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20200316 humandb/
```

## Release v3.0.0

* Interation of a "tumorOnly" Pipeline
* Implementation of a tNGS Pipeline

## Release v2.0.1

* minor bugfixes

## Release v2.0.0

* MIRACUM-Pipe extended to tNGS (tested with Illumina TruSight Tumor 170 panel)
* gnomAD database updated to gnomAD_genome v2.1.1
* UCSC SQL server is now used to annotate CNV regions per default
* Report extended and revised to include links to genome nexus and The Variant Interpretation for Cancer Consortium Meta-Knowledgebase
* stability improvments
* various bugfixes

## Release v1.2.0

* Added export of copy number alterations for later import into cBioPortal

## Release v1.1.0

* Included export into a maf file for later import into cBioPortal

## Minor Release v1.0.1

* Added support for more than 64GB of RAM
* Included txt2maf function as possible import format into cBioPortal

## Release v1.0.0

* Initial release of MIRACUM-Pipe v1.0.0
