# Title     : install_packages
# Objective : installs all required R packages
# Created by: Raphael Scheible
# Created on: 10.07.2019
# Extended by: Patrick Metzger
# Extended on: 09.11.2021

options(repos=structure(c(CRAN="http://cloud.r-project.org")), timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase")
packages <- c(
    "foreach",
    "doMC",
    "openxlsx",
    "circlize",
    "knitr",
    "kableExtra",
    "OmicCircos",
    "rtracklayer",
    "org.Hs.eg.db",
    "Homo.sapiens",
    "RMySQL",
    "biomaRt",
    "Rsamtools",
    "gtrellis",
    "ComplexHeatmap",
    "VariantAnnotation",
    "BSgenome.Hsapiens.UCSC.hg19",
    "TxDb.Hsapiens.UCSC.hg19.knownGene",
    "gdata",
    "stringi",
    "tidyr",
    "ensembldb",
    "EnsDb.Hsapiens.v75",
    "GenomicRanges",
    "dplyr",
    "magrittr",
    "pracma",
    "getopt",
    "IRanges",
    "DNAcopy",
    "copynumber",
    "devtools",
    "optparse",
    "bedr",
    "stringr",
    "plyr"
)
BiocManager::install(pkgs = packages, update  = TRUE, ask = FALSE)
BiocManager::install(c("YAPSA", "SomaticSignatures"), update  = TRUE, ask = FALSE)
devtools::install_version("sequenza", version = "3.0.0", repos = "http://cran.us.r-project.org")
devtools::install_github('sztup/scarHRD', build_vignettes = FALSE)
