############################
# Filter Variants Function #
############################

filtering <- function(snpfile, indelfile, snpefffile_snp, snpefffile_indel,
                      outfile, path_data, path_script, covered_region, mode = "T"){
  #' Filter Variants
  #'
  #' @description Filters the somatic SNPs and InDel for analysis
  #'
  #' @param snpfile dataframe. Table of SNVs
  #' @param indelfile dataframe. Table of InDels
  #' @param snpefffile_snp dataframe. Table with Results of SnpEff for SNVs
  #' @param snpefffile_indel dataframe. Table with Results of SnpEff for indels
  #' @param outfile dataframe. Table with results
  #' @param path_data string. Directory of the databases
  #' @param path_script string. Directory of required scripts
  #' @param sureselect string. Kind of sequencer
  #' @param mode string. Mode for filtering: "T", "N", "LOH"
  #'
  #' @return returns list of
  #' @return Table dataframe. Filtered table of mutations
  #' @return tmb numerical. Tumor mutational burden
  #'
  #' @note Please make sure that the columns of your input and your databases
  #' @note have the right names.
  #'
  #' @details The mutations are filtered to find the pathogenic mutations.
  #' @details First only the mutations passing all the quality filters are
  #' @details considered. Then the tumor mutational burden is calculated in
  #' @details Tumor mode ("T"). Afterwards the filter for the functionality of
  #' @details the genetic region is applied and the mutations with the desired
  #' @details exonic functions are chosen. Only rare mutations are kept.
  #' @details In the following some additional information is extraced of the
  #' @details table: in Normal and Tumor mode ("N", "T") Variant Allele
  #' @details Frequency (VAF), Readcounts and Zygosity, in LoH mode ("LOH") VAF
  #' @details and Readcounts both for Tumor and Normal.
  #' @details The full gene names are added to the dataframe.
  #' @details In Addition to that some database queries are done and finally
  #' @details the canonical aminoacid and base changes are added.
  #' @details As last step Condel prediction scores are written into another
  #' @details column.
  #' @details After ordering the columns the final dataframe is written into
  #' @details the xlsx file "outfile".
  require(org.Hs.eg.db)
  require(stringr)
  require(openxlsx)
  require(Rsamtools)
  require(gdata)

  source(paste(path_script, "filtering_tools.R", sep = "/"))

  # Read Data
  x.snp <- read.csv(snpfile, stringsAsFactors = FALSE)
  x.indel <- read.csv(indelfile, stringsAsFactors = FALSE)
  x.snp <- x.snp[- c(1 : 25),]
  x.indel <- x.indel[- c(1 : 25),]
  x <- rbind(x.snp, x.indel)

  # Quality Filter
  id.pass <- grep("PASS", x$Otherinfo)
  if (length(id.pass) > 0) {
    x <- x[id.pass,]
  } else {
    stop("No variant passed quality filter!")
  }

  if (mode == "T") {
    # TumorMutationBurden
    tmb <- tumbu(x, covered_region)
  } else {
    tmb <- 0
  }

  # Filter for function
  x <- filt(x, "intergenic")
  x <- filt(x, "intronic")
  x <- filt(x, "ncRNA_exonic")
  x <- filt(x, "UTR")
  # Filter for exonic function
  test <- as.character(x$ExonicFunc.refGene)
  syn.snv <- which(test == "synonymous SNV")
  if (length(syn.snv) > 0) {
    x <- x[- syn.snv,]
  }
  # Filter for rare mutations
  x <- rare(x)

  # Extract VAF, Readcounts (and Zygosity)
  x <- vrz(x, mode)

  if (dim(x)[1] != 0) {
    # Include GeneName
    x$Start <- as.numeric(as.character(x$Start))
    x$End <- as.numeric(as.character(x$End))
    x$Gene.refGene <- as.character(x$Gene.refGene)
    x <- gene_name(x)

    # Further Annotation
    # Database Queries
    x <- isflag(x, dbfile = paste(path_data, "flag_genes.txt", sep = "/"))
    x <- isogtsg(x, dbfile = paste(path_data, "CancerGenesList.txt",
    sep = "/"))
    # x <- ishs(x, paste(path_data, "hotspots_V2.txt", sep = "/"))
    # x <- isihs(x, paste(path_data, "hotspots_V2_indel.txt", sep = "/"))
    x <- ishs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- isihs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- rvis(x, paste(path_data, "RVIS_score.txt", sep = "/"))
    x <- trgt(x, paste(path_data, "TARGET_db.txt", sep = "/"))
    x <- dgidb(x, paste(path_data, "DGIdb_interactions.tsv", sep = "/"))
    x <- oncokb(x, paste(path_data, "allActionableVariants.txt", sep = "/"))
    if (dim(x)[1] != 0) {
      x <- snpeff(x, snpefffile_snp, snpefffile_indel)
      x.condel <- addCondel(x, paste(path_data, "fannsdb.tsv.gz", sep = "/"))

      if (mode == "N" | mode == "T") {
        ids <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
        "Gene.refGene", "GeneName", "ExonicFunc.refGene",
        "AAChange.refGene", "AAChange.SnpEff", "CChange.SnpEff",
        "gnomAD_exome_NFE", "ExAC_NFE",
        "esp6500siv2_ea", "EUR.sites.2015_08",
        "Variant_Allele_Frequency", "Variant_Reads",
        "Zygosity", "is_tumorsuppressor", "is_oncogene", "is_hotspot",
        "is_flag", "target", "DGIdb", "condel.label", "cosmic_coding",
        "CLNSIG", "CLINSIG.SnpEff", "InterVar_automated",
        "CADD_phred", "DANN_score", "SIFT_pred",
        "Polyphen2_HDIV_pred", "avsnp150", "rvis")
      } else if (mode == "LOH") {
        ids <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
        "Gene.refGene", "GeneName", "ExonicFunc.refGene",
        "AAChange.refGene", "AAChange.SnpEff", "CChange.SnpEff",
        "gnomAD_exome_NFE", "ExAC_NFE", "esp6500siv2_ea",
        "EUR.sites.2015_08", "VAF_Normal", "VAF_Tumor", "Count_Normal",
        "Count_Tumor", "is_tumorsuppressor", "is_oncogene", "is_hotspot",
        "is_flag", "target", "DGIdb", "condel.label", "cosmic_coding",
        "CLNSIG", "CLINSIG.SnpEff", "InterVar_automated",
        "CADD_phred", "DANN_score", "SIFT_pred", "Polyphen2_HDIV_pred",
        "avsnp150", "rvis")
      }
      idx <- match(ids, colnames(x.condel))
      tot <- seq(1, ncol(x.condel))
      idx2 <- setdiff(tot, idx)

      x <- x.condel[, c(idx, idx2)]
      write.xlsx(x, outfile, keepNA = FALSE, rowNames = FALSE, firstRow = TRUE)
      return(list(table = x, tmb = tmb))
    } else if (mode == "N" | mode == "T") {
      print("No SNVs passed filter!")
      write.xlsx(x, outfile, keepNA = FALSE, rowNames = FALSE, firstRow = TRUE)
      return(list(table = x, tmb = tmb))
    } else if (mode == "LOH") {
      print("No LOH passed filter!")
      write.xlsx(x, outfile, keepNA = FALSE, rowNames = FALSE, firstRow = TRUE)
      return(list(table = x, tmb = tmb))
    }
  }
}
