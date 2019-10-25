# Title     : install_packages
# Objective : installs all required R packages
# Created by: Raphael Scheible
# Created on: 10.07.19

options(repos=structure(c(CRAN="http://cran.us.r-project.org")))

install.packages(c( "foreach", "doMC"))
install.packages("openxlsx", dependencies=TRUE)
install.packages("RColorBrewer")
install.packages("RColorBrewer")
install.packages("circlize")
install.packages("kableExtra")
install.packages("knitr")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("OmicCircos")
BiocManager::install("rtracklayer")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("Homo.sapiens")

BiocManager::install("biomaRt")
BiocManager::install("Rsamtools")
BiocManager::install("gtrellis")
BiocManager::install("ComplexHeatmap")
BiocManager::install("YAPSA")
BiocManager::install("SomaticSignatures")
BiocManager::install("VariantAnnotation")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("gdata")
