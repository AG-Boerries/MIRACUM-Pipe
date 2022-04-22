#### Filter Variants Function

filtering <- function(
  snpfile,
  indelfile,
  snpefffile_snp,
  snpefffile_indel,
  outfile,
  path_data,
  path_script,
  mode = "T",
  center = "Freiburg",
  id = id,
  protocol,
  sureselect,
  vaf = 5,
  min_var_count = 4,
  maf = 0.001,
  actionable_genes = NA,
  covered_exons = covered_exons,
  cov_t = 1,
  sureselect_type
  ) {
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
  #' @param sureselect string. Kind of sequencer; path to bed file
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
  require(GenomicRanges)

  source(paste(path_script, "filtering_tools.R", sep = "/"))

  # Read Data
  x.snp <- read.csv(
    file  = snpfile,
    stringsAsFactors = FALSE,
    comment.char = "#",
    header = TRUE
  )
  x.indel <- read.csv(
    file = indelfile,
    stringsAsFactors = FALSE,
    comment.char = "#",
    header = TRUE
  )
  x <- rbind(x.snp, x.indel)

  # ANNOVAR changed name of "Otherinfo" column in latest release to "Otherinfo1"
  colnames(x)[grep(pattern = "Otherinfo1", x = colnames(x))] <- "Otherinfo"

  # Filter for actionable genes in Germline
  if (protocol == "somaticGermline" & mode == "N") {
    x <- actionable(x, "Gene.refGene", actionable_genes)
  }

  # Quality Filter
  id.pass <- grep("PASS", x$Otherinfo)
  if (length(id.pass) > 0) {
    x <- x[id.pass, ]
  } else {
    stop("No variant passed quality filter!")
  }

  # Extract VAF, Readcounts (and Zygosity)
  x <- vrz(x, mode, protocol = protocol)

  # Filter for targets region in panels
  #if (protocol == "panelTumor") {
  x <- target_check(x, sureselect)
  #}

  # Calcualte covered region from bed file
  cov_region <- covered_region(sureselect = sureselect, mode = mode)

  # Additional Filter for VAF
  x <- exclude(x, vaf = vaf)

  # Remove Variants with Variant Read Count below 4
  x <- mrc(x = x, min_var_count = min_var_count)

  # snpEff
  x <- snpeff(x, snpefffile_snp, snpefffile_indel, protocol)

  # replace synonymous entries from refGene with snpEff consequence
  test <- as.character(x$ExonicFunc.refGene)
  syn.snv <- which(test == "synonymous SNV")
  if (length(syn.snv) > 0) {
    x$ExonicFunc.refGene[syn.snv] <- x[syn.snv, "Consequence_snpEff"]
  }

  # Filter for function
  x <- filt(x, "intergenic")
  x <- filt(x, "intronic")
  x <- filt(x, "ncRNA_exonic")
  x <- filt(x, "UTR")
  x <- filt(x, "upstream")
  x <- filt(x, "downstream")

  # TumorMutationBurden
  if (mode == "T") {
    if (protocol %in% c("somatic", "somaticGermline", "tumorOnly")) {
      id_ex <- which(x$Func.refGene == "exonic")
      x_coding <- x[id_ex, ]
      tmb <- tmb_ex(x_coding, covered_exons, mode = "T", cov_t)
    }
  } else {
    tmb <- list(tmb = NULL, exon_region = NULL, used_exon_region = NULL,
                number_used_mutations_tmb = NULL)
  }
  
  # Filter for exonic function
  test <- as.character(x$ExonicFunc.refGene)
  syn.snv <- which(test %in% c("synonymous SNV", "synonymous_variant", "unknown", "synonymous_variant;synonymous_variant", ""))
  if (length(syn.snv) > 0) {
    x <- x[- syn.snv,]
  }

  # Filter for synonymous_variants from snpEff
  test <- as.character(x$Consequence_snpEff)
  syn.snv <- which(test == "synonymous_variant")
  if (length(syn.snv) > 0) {
    x <- x[- syn.snv,]
  }

  # Filter for rare mutations
  x <- rare(x, maf)

  # TMB calculation for TSO500; only rare mutations; assumption: rare mutations are "somatic"
  if (sureselect_type == "TSO500") {
    if (mode == "T"){
      id_ex <- which(x$Func.refGene == "exonic")
      x_coding <- x[id_ex, ]
      tmb <- tmb_ex(x_coding, covered_exons, mode = "T", cov_t)
    } else {
      tmb <- list(tmb = NULL, exon_region = NULL, used_exon_region = NULL,
                  number_used_mutations_tmb = NULL)
    }
  }

  
  # TumorMutationBurden
  # if (mode == "T") {
  #   if (protocol %in% c("panelTumor")) {
  #    if (sureselect_type %in% c("TSO500")) {
  #      tmb <- tumbu(x, 1.27)
  #     } else {
  #       tmb <- tmb_ex(x, coveredExons, mode = "T", cov_t)
  #     }
  #   }
  # } else {
  #   tmb <- tmb <- list(tmb = NULL, exon_region = NULL, used_exon_region = NULL,
  #                      number_used_mutations_tmb = NULL)
  # }

  if (dim(x)[1] != 0) {
    # Include GeneName
    x$Start <- as.numeric(as.character(x$Start))
    x$End <- as.numeric(as.character(x$End))
    x$Gene.refGene <- as.character(x$Gene.refGene)
    x <- gene_name(x)

    # Further Annotation
    # Database Queries
    x <- isflag(x, dbfile = paste(path_data, "flag_genes.txt", sep = "/"))
    x <- isogtsg(x, dbfile = paste(path_data, "cancerGeneList.tsv",sep = "/"))
    x <- ishs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- isihs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- rvis(x, paste(path_data, "RVIS_score.txt", sep = "/"))
    x <- trgt(x, paste(path_data, "TARGET_db.txt", sep = "/"))
    x <- dgidb(x, paste(path_data, "DGIdb_interactions.tsv", sep = "/"))
    x <- oncokb(x, paste(path_data, "oncokb_biomarker_drug_associations.tsv", sep = "/"))
    x.condel <- addCondel(x, paste(path_data, "fannsdb.tsv.gz", sep = "/"))

    if (mode == "N" | mode == "T") {
      ids <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
      "Gene.refGene", "GeneName", "ExonicFunc.refGene",
      "AAChange.refGene", "AAChange", "CChange",
      "AF_popmax",
      "Variant_Allele_Frequency", "Variant_Reads",
      "Zygosity", "is_tumorsuppressor", "is_oncogene", "is_hotspot",
      "is_flag", "target", "DGIdb", "condel.label", "cosmic_coding",
      "CLNSIG", "CLINSIG", "InterVar_automated",
      "CADD_phred", "DANN_score", "SIFT_pred",
      "Polyphen2_HDIV_pred", "avsnp150", "rvis")

    } else if (mode == "LOH") {
      ids <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
      "Gene.refGene", "GeneName", "ExonicFunc.refGene",
      "AAChange.refGene", "AAChange", "CChange",
      "AF_popmax", "VAF_Normal", "VAF_Tumor", "Count_Normal",
      "Count_Tumor", "is_tumorsuppressor", "is_oncogene", "is_hotspot",
      "is_flag", "target", "DGIdb", "condel.label", "cosmic_coding",
      "CLNSIG", "CLINSIG", "InterVar_automated",
      "CADD_phred", "DANN_score", "SIFT_pred", "Polyphen2_HDIV_pred",
      "avsnp150", "rvis")
    }
    idx <- match(ids, colnames(x.condel))
    tot <- seq(1, ncol(x.condel))
    idx2 <- setdiff(tot, idx)

    x <- x.condel[, c(idx, idx2)]
    x$ExonicFunc.refGene[which(is.na(x$ExonicFunc.refGene) | x$ExonicFunc.refGene == "NA")] <- x$Func.refGene[which(is.na(x$ExonicFunc.refGene))]

    out.maf <- txt2maf(
      input = x,
      Center = center,
      refBuild = 'GRCh37',
      id = id,
      sep = "\t",
      idCol = NULL,
      Mutation_Status = mode,
      protocol = protocol,
      snv_vcf = snpefffile_snp,
      indel_vcf = snpefffile_indel
    )

    return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )

  } else if (mode == "N" | mode == "T") {
    print("No SNVs passed filter!")
    out.maf <- data.frame()
    msi <- NULL
    return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )

  } else if (mode == "LOH") {
    print("No LOH passed filter!")
    out.maf <- data.frame()
    return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )
  }
}

filtering_mutect2 <- function(
  snpfile,
  snpefffile,
  outfile,
  path_data,
  path_script,
  mode = "T",
  center = "Freiburg",
  id = id,
  protocol,
  sureselect,
  vaf = 5,
  min_var_count = 20,
  maf = 0.01,
  actionable_genes = NA,
  covered_exons = covered_exons,
  cov_t = 1,
  sureselect_type
  ) {
  #' Filter Variants
  #'
  #' @description Filters the somatic SNPs and InDel for analysis
  #'
  #' @param snpfile dataframe. Table of SNVs
  #' @param snpefffile_snp dataframe. Table with Results of SnpEff for SNVs
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
  require(stringr)
  require(openxlsx)
  require(Rsamtools)
  require(GenomicRanges)

  # Read Data
  x <- read.delim(file = snpfile, header = T, stringsAsFactors = F, comment.char = "#")
  
  # Filter for targeted region in tNGS
  x <- target_check(x, sureselect)

  # Calcualte covered region from bed file
  cov_region <- covered_region(sureselect = sureselect, mode = mode)
    
  id <- grep(pattern = "CHROM", x = x$Chr)
  if (length(id) > 0) {
    x <- x[-id, ]
  }
  # Quality Filter
  if (protocol ==  "tumorOnly" & sureselect_type == "V5UTR"){
    x <- x
  } else if (sureselect_type %in% c("V5UTR", "V6", "V6UTR")) {
    id.pass <- grep("PASS", x$Otherinfo10)
    if (length(id.pass) > 0) {
      x <- x[id.pass, ]
    } else {
      stop("No variant passed quality filter!")
    }
  } else {
    id.pass <- grep("PASS", x$Otherinfo10)
    if (length(id.pass) > 0) {
      x <- x[id.pass, ]
    } else {
      stop("No variant passed quality filter!")
    }
  }
  
  # Extract VAF, Readcounts (and Zygosity)
  x <- vrz_gatk(x = x, mode = mode, protocol = protocol)
  
  # VAF Filter
  x <- exclude_gatk(x, vaf = vaf/100)
  
  # Remove Variants with Variant Read Count below 4/20
  x <- mrc(x = x, min_var_count = min_var_count)
  # snpEFF Annotation
  x <- snpeff(
    x = x,
    sef_snp = snpefffile,
    protocol = protocol
  )
  
  # replace synonymous entries from refGene with snpEff consequence
  test <- as.character(x$ExonicFunc.refGene)
  syn.snv <- which(test == "synonymous SNV")
  if (length(syn.snv) > 0) {
    x$ExonicFunc.refGene[syn.snv] <- x[syn.snv, "Consequence_snpEff"]
  }
  
  # Filter for exonic function
  x <- filt(x, "intergenic")
  x <- filt(x, "intronic")
  x <- filt(x, "ncRNA_exonic")
  x <- filt(x, "UTR")
  x <- filt(x, "upstream")
  x <- filt(x, "downstream")
  
  # TMB calculation for WES
  if (protocol == "somaticGermline" | protocol == "somatic"){
    if (mode == "T") {
      id_ex <- which(x$Func.refGene == "exonic")
      x_coding <- x[id_ex, ]
      tmb <- tmb_ex(x_coding, covered_exons, mode = "T", cov_t)
    } else {
      tmb <- list(tmb = NULL, exon_region = NULL, used_exon_region = NULL,
                  number_used_mutations_tmb = NULL)
    }
  }
  
  # Filter for exonic function
  test <- as.character(x$ExonicFunc.refGene)
  syn.snv <- which(test %in% c("synonymous SNV", "synonymous_variant", "unknown", "synonymous_variant;synonymous_variant", ""))
  if (length(syn.snv) > 0) {
    x <- x[- syn.snv,]
  }
  
  # Filter for synonymous_variants from snpEff
  test <- as.character(x$Consequence_snpEff)
  syn.snv <- which(test == "synonymous_variant")
  if (length(syn.snv) > 0) {
    x <- x[- syn.snv,]
  }
  
  # Filter for rare mutations
  x <- rare(x, maf = maf)
  
  # TMB calculation for panel (TSO500) and tumorOnly WES; only rare mutations; assumption: rare mutations are "somatic"
  if (protocol == "panelTumor" | protocol == "tumorOnly") {
    if (mode == "T") {
      id_ex <- which(x$Func.refGene == "exonic")
      x_coding <- x[id_ex, ]
      tmb <- tmb_ex(x_coding, covered_exons, mode = "T", cov_t)
    } else {
      tmb <- list(tmb = NULL, exon_region = NULL, used_exon_region = NULL,
                  number_used_mutations_tmb = NULL)
    }
  }
  
  if (dim(x)[1] != 0) {
    # Include GeneName
    x$Start <- as.numeric(as.character(x$Start))
    x$End <- as.numeric(as.character(x$End))
    x$Gene.refGene <- as.character(x$Gene.refGene)
    x <- gene_name(x)
    
    # Further Annotation
    # Database Queries
    x <- isflag(x, dbfile = paste(path_data, "flag_genes.txt", sep = "/"))
    x <- isogtsg(x, dbfile = paste(path_data, "cancerGeneList.tsv", sep = "/"))
    x <- ishs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- isihs(x, paste(path_data, "hotspots_v2.xls", sep = "/"))
    x <- rvis(x, paste(path_data, "RVIS_score.txt", sep = "/"))
    x <- trgt(x, paste(path_data, "TARGET_db.txt", sep = "/"))
    x <- dgidb(x, paste(path_data, "DGIdb_interactions.tsv", sep = "/"))
    x <- oncokb(x, paste(path_data, "oncokb_biomarker_drug_associations.tsv", sep = "/"))
    
    x.condel <- addCondel(x, paste(path_data, "fannsdb.tsv.gz", sep = "/"))
    
    if (mode == "N" | mode == "T"){
      ids <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
               "Gene.refGene", "GeneName", "ExonicFunc.refGene",
               "AAChange.refGene", "AAChange", "CChange",
               "AF_popmax", 
               "Variant_Allele_Frequency", "Variant_Reads",
               "Zygosity", "is_tumorsuppressor", "is_oncogene", "is_hotspot",
               "is_flag", "target", "DGIdb", "condel.label", "cosmic_coding",
               "CLNSIG", "CLINSIG", "InterVar_automated",
               "CADD_phred", "DANN_score", "SIFT_pred",
               "Polyphen2_HDIV_pred", "avsnp150", "rvis", "REVEL_score", "Consequence_snpEff")
    } else if (mode == "LOH"){
      ids <- c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
               "Gene.refGene", "GeneName", "ExonicFunc.refGene",
               "AAChange.refGene", "AAChange", "CChange",
               "AF_popmax", "VAF_Normal", "VAF_Tumor", "Count_Normal",
               "Count_Tumor", "is_tumorsuppressor", "is_oncogene",
               "is_hotspot", "is_flag", "target", "DGIdb", "condel.label",
               "cosmic_coding", "CLNSIG", "CLINSIG",
               "InterVar_automated", "CADD_phred", "DANN_score",
               "SIFT_pred", "Polyphen2_HDIV_pred", "avsnp150", "rvis", "REVEL_score", "Consequence_snpEff")
    }
    idx <- match(ids, colnames(x.condel))
    tot <- seq(1, ncol(x.condel))
    idx2 <- setdiff(tot, idx)
    x <- x.condel[, c(idx, idx2)]    
    
    # MAF
    out.maf <- txt2maf_mutect2(
      input = x,
      protocol = protocol,
      snv_vcf = snpefffile,
      Center = center,
      refBuild = 'hg19',
      id = sample,
      sep = "\t",
      idCol = NULL,
      Mutation_Status = mode
    )
    
    return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )

  } else if (mode == "N" | mode == "T") {
    print("No SNVs passed filter!")
    out.maf <- data.frame()
   return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )

  } else if (mode == "LOH") {
    print("No LOH passed filter!")
    out.maf <- data.frame()
    return(
      list(
        table = x,
        tmb = tmb$tmb,
        exon_region = tmb$exon_region,
        maf = out.maf,
        covered_region = cov_region,
        used_exon_region = tmb$used_exon_region,
        number_used_mutations_tmb = tmb$number_used_mutations_tmb
      )
    )
  }
}
