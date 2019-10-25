# MIRACUM-Pipe
MIRACUM-Pipe incorporates tools for detecting single nucleotide variants (SNVs), insertions and deletions (InDels), loss of heterozygosity (LoH), copy number variations (CNVs) as well as quality statistics including the coverage on WES data. Various functional prediction and annotation databases are integrated to automatically annotate the identi-fied variants. The workflow is designed as a fully automated “one-click” solution from the raw sequencing files to the PDF report containing quality assessments, the identified and annotated variants (SNVs, InDels and LoH), copy number variations as well as a functional enrichment of the SNVs and CNVs respectively.

## Getting Started

This repo is intended to be run as docker (see [MIRACUM-Pipe-docker]()). Alternatively, you can pick a [Galaxy version]().

In some cases Docker or Galaxy might be inapropriate. Therefore, one can install this software on a [Debian 10](https://www.debian.org/) system.

### Prerequisites
We offer installation scripts. In order to setup up the main components execute the following scripts.

```bash
# install all os requirements
./debian/setup.sh

# install all R-packages
Rsctipt RScripts/install_packages.R

# install tools
./tools/install.sh
```

Further one needs to install some components manually. See [MIRACUM-Pipe-docker]() for more information.
You can also simply apply the [setup.sh]() file for this repo. Just download it into the root of your clone of this project.


Additional Files
	* CaptureRegions, e.g. Agilent SureSelect V5UTR according to the used CaptureKit
	* dbSNP vcf file corresbonding to the one used from ANNOVAR, e.g. snp150hg19.vcf.gz
	* FLAG genes (included, flag_genes.txt) (Shyr et al., 2014)
	* CancerGenes from OncoKB (needed)
	* Cancer hotspots from cancerhotspots.org (needed)
	* RVIS score (included, RVIS_score.txt) (Petrovski S. et al., 2013)
	* TARGET db (included, TARGET_db.txt) (https://software.broadinstitute.org/cancer/cga/target)
	* DGIdb from www.dgidb.org (needed)
	* Actionable alterations from OncoKB (needed)
	* Condel score from FannsDB (needed)
	* captured genes according to the used CaptureKit (Targets.txt supplied for Agilent SureSelect)
	* Cancer associated pathways (included, MTB_GeneSets.RData)
	* Consensus Pathways (included, Consensus.RData)
	* Reactome Pathways (included, Reactome.RData)
	* GO Pathways (included, GOGeneSets.RData)
	* HallmarksOfCancer (needed)

If available, the used versions are noted.

## Running MIRACUM-Pipe
Once you installed the system and setup the required tools and databases, you can simply run the pipeline by setting up patients and run `./miracum-pipe.sh`. For more documentation about this, see the documentation of [MIRACUM-Pipe-docker]().


### Adding new databases respectively update current databases
For the databases contained in ANNOVAR updating or installing new databases is very simple following the descriptions on the homepage of ANNOVAR. If the additional database is installed it has to be added in the yaml configuration file. Adding a database currently not supported by ANNOVAR is also possible but one has to add a few lines and the actual database to the corresponding part in the **filtering.R** script. See the documentation of the R scripts.

### Information about the included R-scripts
All neccessary informations about the included R-scripts can be found with the *?functionName* function in R.
The predefined report in Report.rnw contains in line 45 the author name which should be adjusted to the person generating the report.

## Example
For testing MIRACUM-Pipe the public available dataset from Texas (http://txcrb.org/data.html and https://www.nature.com/articles/sdata201610) can be used. We provide the PDF reports for the samples TCRBOA6 and TCRBOA7. Both samples were run with the default parameters.

## Authors

* Patrick Metzger
* Raphael Scheible
* Maria Hess
* Victor Jaravine
* Jochen Hochrein
* Martin Boeker
* Hauke Busch
* Geoffroy Andrieux
* Melanie Börries

## License
This work is licensed under [GNU Affero General Public License version 3](https://opensource.org/licenses/AGPL-3.0).

## Acknowledgments
We thank

* The Molecular Tumor Board Freburg Team
* The whole MIRACUM consortia
* The German Ministry of Education and Research (BMBF) for funding
* The developers from Control-FREEC for the code on CNV significance

