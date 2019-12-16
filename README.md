# MIRACUM-Pipe

MIRACUM-Pipe incorporates tools for detecting single nucleotide variants (SNVs), insertions and deletions (InDels), loss of heterozygosity (LoH), copy number variations (CNVs) as well as quality statistics including the coverage on WES data. Various functional prediction and annotation databases are integrated to automatically annotate the identi-fied variants. The workflow is designed as a fully automated “one-click” solution from the raw sequencing files to the PDF report containing quality assessments, the identified and annotated variants (SNVs, InDels and LoH), copy number variations as well as a functional enrichment of the SNVs and CNVs respectively.

## Getting Started

This repo is intended to be run as docker (see [MIRACUM-Pipe-docker](https://github.com/AG-Boerries/MIRACUM-Pipe-docker)). Alternatively, you can pick a [Galaxy version](https://github.com/AG-Boerries/MIRACUM-Pipe-galaxy).

In some cases Docker or Galaxy might be inapropriate. Therefore, one can install this software on a [Debian 10](https://www.debian.org/) system.

### Prerequisites

We offer installation scripts. In order to setup up the main components execute the following scripts.

```bash
# install all os requirements
./debian/setup.sh

# install all R-packages
Rscript RScripts/install_packages.R

# install tools
./tools/install.sh
```

Further one needs to install some components manually. See [MIRACUM-Pipe-docker](https://github.com/AG-Boerries/MIRACUM-Pipe-docker) for more information.
You can also simply apply the [setup.sh](https://github.com/AG-Boerries/MIRACUM-Pipe-docker/blob/master/setup.sh) file for this repo. Just download it into the root of your clone of this project.

## Running MIRACUM-Pipe

Once you installed the system and setup the required tools and databases, you can simply run the pipeline by setting up patients and run `./miracum_pipe.sh`. For more documentation about this, see the documentation of [MIRACUM-Pipe-docker](https://github.com/AG-Boerries/MIRACUM-Pipe-docker).

### Adding new databases respectively update current databases

For the databases contained in annovar updating or installing new databases is possible following the descriptions on the homepage of [annovar](http://annovar.openbioinformatics.org/en/latest/user-guide/). If you install additional databases they have to be added in the yaml configuration file, e.g. `conf/custom.yaml`. Be aware that the additionally added databases are not considered in the PDF Report.

### Information about the included R-scripts

All neccessary informations about the included R-scripts can be found with the *?functionName* function in R.

## Example

For testing MIRACUM-Pipe an example dataset is provided and downloaded via the **setup.sh** script. It contains the necessary capture region files (`V5UTR.bed`, `V5UTR_Targets.txt`), the sequence files (`germline_R{1/2}.fastq.gz`, `tumor_R{1/2}.fastq.gz`) and the `patient.yaml`. For the sake of runtime the fastq files were adjusted to contain only chromosome 12.

## Authors

* Patrick Metzger
* Raphael Scheible
* Maria Hess
* Victor Jaravine
* Martin Boeker
* Geoffroy Andrieux
* Melanie Börries

## License

This work is licensed under [GNU Affero General Public License version 3](https://opensource.org/licenses/AGPL-3.0).

## Acknowledgments

We thank

* The Molecular Tumor Board Freiburg Team
* The MIRACUM consortium
* The German Ministry of Education and Research (BMBF) for funding
* The developers from Control-FREEC for the code on CNV significance
