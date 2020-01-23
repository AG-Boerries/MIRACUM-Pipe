# Title     : install_packages
# Objective : installs all required R packages
# Created by: Raphael Scheible
# Created on: 10.07.19
# Adjusted on: 10.12.19 by Patrick Metzger

options(repos=structure(c(CRAN="http://cloud.r-project.org")))

install.packages(c( "foreach", "doMC"))
install.packages("openxlsx")
install.packages("circlize")
install.packages("knitr")
install.packages("kableExtra")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.8")

BiocManager::install("OmicCircos")
BiocManager::install("rtracklayer")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("Homo.sapiens")
BiocManager::install("RMySQL")
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
BiocManager::install("data.table")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("RMySQL")
#BiocManager::install("GSA")
