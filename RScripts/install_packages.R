# Title     : install_packages
# Objective : installs all required R packages
# Created by: Raphael Scheible
# Created on: 10.07.2019
# Extended by: Patrick Metzger
# Extended on: 14.05.2020

options(repos=structure(c(CRAN="http://cloud.r-project.org")), timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase")
packages <- c("foreach", "doMC", "openxlsx", "circlize", "knitr", "kableExtra", "OmicCircos", "rtracklayer", "org.Hs.eg.db",
              "Homo.sapiens", "RMySQL", "biomaRt", "Rsamtools", "gtrellis", "ComplexHeatmap",
              "VariantAnnotation", "BSgenome.Hsapiens.UCSC.hg19",
              "TxDb.Hsapiens.UCSC.hg19.knownGene", "gdata", "stringi", "tidyr", "ensembldb",
              "EnsDb.Hsapiens.v75", "GenomicRanges", "dplyr", "magrittr", "pracma", "getopt", "IRanges", "DNAcopy",
              "copynumber", "sequenza", "devtools")
BiocManager::install(pkgs = packages, update  = TRUE, ask = FALSE)

url <- "https://cran.r-project.org/src/contrib/Archive/lsei/lsei_1.2-0.1.tar.gz"
url2 <- "https://cran.r-project.org/src/contrib/Archive/MSIseq/MSIseq_1.0.0.tar.gz"
devtools::install_url(url)
devtools::install_url(url2)
BiocManager::install(c("YAPSA", "SomaticSignatures"), update  = FALSE, ask = FALSE)
devtools::install_github('sztup/scarHRD',build_vignettes = FALSE)