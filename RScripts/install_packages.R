# Title     : install_packages
# Objective : installs all required R packages
# Created by: Raphael Scheible
# Created on: 10.07.2019
# Extended by: Patrick Metzger
# Extended on: 14.05.2020

options(repos=structure(c(CRAN="http://cloud.r-project.org")))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.8")

# TODO package install in a more efficient way
packages <- c("foreach", "doMC", "openxlsx", "circlize", "knitr", "kableExtra", "OmicCircos", "rtracklayer", "org.Hs.eg.db",
              "Homo.sapiens", "RMySQL", "biomaRt", "Rsamtools", "gtrellis", "ComplexHeatmap",
              "YAPSA", "SomaticSignatures", "VariantAnnotation", "BSgenome.Hsapiens.UCSC.hg19",
              "TxDb.Hsapiens.UCSC.hg19.knownGene", "gdata", "stringi", "tidyr", "ensembldb",
              "EnsDb.Hsapiens.v75", "GenomicRanges", "dplyr", "magrittr", "pracma", "getopt", "IRanges", "DNAcopy")
BiocManager::install(pkgs = packages, update  = TRUE, ask = FALSE)