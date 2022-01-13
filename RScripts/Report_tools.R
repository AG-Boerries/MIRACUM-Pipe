## key_Results
keys <- function(
  mut_sig,
  mutation_analysis_result,
  mutation_analysis_result_gd,
  cnv_analysis_results,
  filt_result_td,
  filt_result_gd,
  med_tmb,
  protocol,
  fusions
) {

  if (sureselect_type %in% c("V5UTR", "V6UTR", "V6")){
    sureselect_type <- paste("Agilent SureSelect", sureselect_type, sep = " ")
  } else if (sureselect_type %in% c("TSO500", "TruSight_Tumor")) {
    sureselect_type <- paste("Illumina", sureselect_type, sep = " ")
  }

  if (protocol == "somaticGermline" | protocol == "somatic") {
    if (mutation_analysis_result$msi < 10) {
      msi_helper <- "MSS"
    } else {
      msi_helper <- "MSI"
    }
    brca_helper <- which(mut_sig$output$Summary$Signature ==  "AC3")
    if (length(brca_helper) == 1 & mut_sig$output$Summary["AC3", 3] > 1.0) {
      brca_helper <- paste0(
        round(
          mut_sig$output$Summary["AC3", 3],
          digits = 1
        ),
        " (", round(mut_sig$output$Summary["AC3", 4], digits = 1),
        ";", round(mut_sig$output$Summary["AC3", 5], digits = 1) , ")"
      )
    } else {
      brca_helper <- "<1%"
    }
  }
  if (protocol == "panelTumor" | protocol == "tumorOnly") {
    if (mutation_analysis_result$msi < 20) {
      msi_helper <- "MSS"
    } else {
      msi_helper <- "MSI"
    }
    brca_helper <- which(mut_sig ==  "AC3")
    if (length(brca_helper) == 1 & mut_sig["AC3", 3]*100 > 1) {
      brca_helper <- paste0(round(mut_sig["AC3", 3]*100, digits = 1))
    } else {
      brca_helper <- "<1%"
    }
  }
  if (protocol == "panelTumor") {
    if (!is.null(fusions)) {
      fus_tmp <- fusions$Fusion_OV$Fusionen
      fus_tmp <- paste(fus_tmp, collapse = ", ")
    } else {
      fus_tmp <- "Keine"
    }
  }

  if (!is.null(med_tmb$med) & !is.null(med_tmb$sd)) {
    tmb_helper <- paste0(
      med_tmb$med, " (",
      as.numeric(med_tmb$sd[1], digits = 2),
      "-", as.numeric(med_tmb$sd[2], digits = 2),
      ")"
    )
  } else if (!is.null(med_tmb$med)) {
    tmb_helper <- med_tmb$med
  } else {
    tmb_helper <- "-"
  }
  if (protocol == "somaticGermline" | protocol == "somatic") {
    mut_tab1 <- data.frame(
      Eigenschaften = c(
        "Capture Kit",
        "Abgedeckte Region (total)",
        "Abgedeckte Region (exonisch)",
        paste0("Mutationslast (exonisch, VAF > ", vaf, "%)"),
        paste0("Mittlere TMB der Entität", " (", entity, ")"),
        paste0("Anzahl somatischer Mutationen inkl. LoH (VAF > ", vaf, "%)"),
        paste0("BRCAness (%) inkl. KI", " (VAF > ", vaf, "%)"),
        "Mikrosatelliten Status (Score)",
        "HRD-Score (HRD-LoH|TAI|LST)",
        "Purity (%)",
        "Ploidity",
        "Anzahl CN- Regionen",
        paste0("Anzahl seltener Keimbahnmutationen (VAF > ", vaf, "%)")
      ), Wert = c(
        as.character(sureselect_type),
        paste(round(x = as.numeric(filt_result_td$covered_region), digits = 2), "Mb", sep = ""),
        paste(round(x = as.numeric(filt_result_td$exon_region), digits = 2), "Mb", sep = ""),
        paste0(round(x = filt_result_td$tmb, digits = 2), "/Mb", " (", filt_result_td$number_used_mutations_tmb, "/", round(x = as.numeric(filt_result_td$used_exon_region), digits = 2), ")"),
        tmb_helper,
        as.character(round(x = sum(as.numeric(mutation_analysis_result$mut_tab[, 2])), digits = 0)),
        brca_helper,
        paste(msi_helper," (", mutation_analysis_result$msi, "%)", sep = ""),
        cnv_analysis_results$hrd$score,
        cnv_analysis_results$purity$purity,
        cnv_analysis_results$purity$ploidy,
        paste0(round(x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1], digits = 0), " Regionen"),
        as.character(round(x = sum(as.numeric(mutation_analysis_result_gd$mut_tab[, 3])), digits = 0))
      )
    )
  }
  if (protocol == "panelTumor") {
    mut_tab1 <- data.frame(
      Eigenschaften = c(
        "Capture Kit",
        "Abgedeckte Region (total)",
        "Abgedeckte Region (exonisch)",
        paste0("Mutationslast (exonisch, VAF > ", vaf, "%)"),
        paste0("Mittlere TMB der Entität", " (", entity, ")"),
        paste0("Anzahl der Mutationen (VAF > ", vaf, "%)"),
        paste0("BRCAness (%)", " (VAF > ", vaf, "%)"),
        "Mikrosatelliten Status (Score)",
        "HRD-Score (HRD-LoH|TAI|LST)",
        "Purity (%)",
        "Ploidity",
        "Anzahl CN- Regionen",
        "Fusionen"
      ), Wert = c(
        as.character(sureselect_type),
        paste(round(x = as.numeric(filt_result_td$covered_region), digits = 2), "Mb", sep = ""),
        paste(round(x = as.numeric(filt_result_td$exon_region), digits = 2), "Mb", sep = ""),
        paste0(round(x = filt_result_td$tmb, digits = 2), "/Mb", " (", filt_result_td$number_used_mutations_tmb, "/", round(x = as.numeric(filt_result_td$used_exon_region), digits = 2), ")"),
        tmb_helper,
        as.character(round(x = sum(as.numeric(mutation_analysis_result$mut_tab[, 3])), digits = 0)),
        brca_helper,
        paste(msi_helper," (", mutation_analysis_result$msi, "%)", sep = ""),
        cnv_analysis_results$hrd$score,
        cnv_analysis_results$purity$purity,
        cnv_analysis_results$purity$ploidy,
        paste0(round(x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1], digits = 0), " Regionen"),
        fus_tmp
      )
    )
  }
  return(mut_tab1)
}

keys2 <- function(mut_sig_ana, mutation_analysis_result, mutation_analysis_result_gd, mutation_analysis_result_10,
                 filt_result_td_10, cnv_analysis_results, filt_result_gd, med_tmb){
  # Check for BRCAness
  helper <- which(mut_sig_ana$output$Summary$Signature ==  "AC3")
  
  # Put Results of Studies together
  if (as.integer(dim(mutation_analysis_result$studies$Soratram)[1] + dim(mutation_analysis_result_gd$studies$Soratram)[1]) == 0) {
    studies3 <- "Soratram: negativ"
  } else {
    mut_st <- paste0(mutation_analysis_result$studies$Soratram$Symbol,
                     "(", mutation_analysis_result$studies$Soratram$AAChange, ")")
    mut_st <- paste0(mut_st, collapse = ", ")
    studies3 <- paste0("Soratram:", mut_st)
  }
  if(dim(mutation_analysis_result$studies$JDQ443A)[1] != 0){
    studies_jd <- "JDQ443A: positiv"
  } else {
    studies_jd <- "JDQ443A: negativ"
  }
  if (length(helper) == 1 & mut_sig_ana$output$Summary["AC3", 3] > 1.0) {
    studies2 <- "TopArt: BRCA positiv."
  } else {
    studies2 <- "TopArt: BRCA negativ."
  }
  
  if(cnv_analysis_results$hrd$sum >= 42) {
    help_hrd <- "++"
  } else if (cnv_analysis_results$hrd$sum >= 30) {
    help_hrd <- "+"
  } else {
    help_hrd <- "o"
  }
  
  if(cnv_analysis_results$purity$purity > 0.8) {
    help_tzg <- "++"
  } else if (cnv_analysis_results$purity$purity > 0.6) {
    help_tzg <- "+"
  } else if (cnv_analysis_results$purity$purity > 0.5) {
    help_tzg <- "o"
  } else {
    help_tzg <- "-"
  }
  
  if (!is.null(med_tmb$sd) & !is.null(med_tmb$med)) {
    if (filt_result_td$tmb < med_tmb$sd[1]) {
      help_tmb <- "--"
    } else if (filt_result_td$tmb > med_tmb$sd[2]) {
      help_tmb <- "++"
    } else if (filt_result_td$tmb < med_tmb$med) {
      help_tmb <- "-"
    } else if (filt_result_td$tmb > med_tmb$med) {
      help_tmb <- "+"
    } else {
      help_tmb <- "o"
    }
  } else if (!is.null(med_tmb$sd)) {
    if (filt_result_td$tmb < med_tmb$med) {
      help_tmb <- "-"
    } else if (filt_result_td$tmb > med_tmb$med) {
      help_tmb <- "+"
    }
  } else {
    help_tmb <- "o"
  }
  
  if (length(helper) == "0") {
    help_ac3 <- "--"
  } else if(mut_sig_ana$output$Summary["AC3", 4] <= 0) {
    help_ac3 <- "+"
  } else if(mut_sig_ana$output$Summary["AC3", 4] > 0) {
    help_ac3 <- "++"
  }
  if (dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1] < 100) {
    help_cnv <- "-"
  } else if (dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1] < 150) {
    help_cnv <- "o"
  } else if (dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1] < 200) {
    help_cnv <- "+"
  } else {
    help_cnv <- "++"
  }
  if (length(which(as.numeric(gsub(pattern = "%", replacement = "",
                                   x = filt_result_gd$table$Variant_Allele_Frequency)) > 10)) < 150) {
    help_germ <- "-"
  } else if (length(which(as.numeric(gsub(pattern = "%", replacement = "",
                                          x = filt_result_gd$table$Variant_Allele_Frequency)) > 10) < 200)) {
    help_germ <- "o"
  } else if (length(which(as.numeric(gsub(pattern = "%", replacement = "",
                                          x = filt_result_gd$table$Variant_Allele_Frequency)) > 10) < 300)) {
    help_germ <- "+"
  } else {
    help_germ <- "++"
  }
  
  if (length(grep(pattern = "positiv", x = studies3)) > 0 & length(grep(pattern = "positive", x = studies2)) > 0) {
    help_studies <- "+/+"
  } else if (length(grep(pattern = "positiv", x = studies3)) > 0) {
    help_studies <- "+/-"
  } else if (length(grep(pattern = "positiv", x = studies2)) > 0) {
    help_studies <- "-/+" 
  } else {
    help_studies <- "-/-"
  }
  if (length(grep(pattern = "positiv", x = studies_jd)) > 0) {
    help_studies <- paste0(help_studies, "/+")
  } else {
    help_studies <- paste0(help_studies, "/-")
  }
  
  if(is.null(mutation_analysis_result$msi)) {
    msi_help <- c("", "")
  } else {
    if (mutation_analysis_result$msi < 3.5) {
      msi_help <- c("MSS", "o")
    } else if (mutation_analysis_result$msi < 10) {
      msi_help <- c("MSI-L", "+")
    } else {
      msi_help <- c("MSI-H", "++")
    }
  }
  
  if (length(helper) == 1 & mut_sig_ana$output$Summary["AC3", 3] < 1.0) {helper <- c()}
  if (length(helper) == 1){
    mut_tab1 <- data.frame(Eigenschaften = c("Mutationslast", "Anzahl somatischer Mutationen", "BRCAness", "HRD-Score (LoH|LST|TAI)",
                                             "Anzahl CN- Regionen", "Anzahl seltener Keimbahnmutationen", "Studiencheck", "", "Tumorzellgehalt (Ploidität)", "MSI"),
                           Wert1 = c(paste0(round(x = filt_result_td_10$tmb, digits = 2),"/Mb"),
                                     as.character(round(x = sum(as.numeric(mutation_analysis_result_10$mut_tab[, 2])), digits = 0)),
                                     paste0(round(x = mut_sig_ana$output$Summary[helper, 3], digits = 2), " %"),
                                     cnv_analysis_results$hrd$score,
                                     paste0(round(x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1], digits = 0), " Regionen"),
                                     length(which(as.numeric(gsub(pattern = "%", replacement = "", x = filt_result_gd$table$Variant_Allele_Frequency)) > 10)),
                                     studies3, studies_jd, cnv_analysis_results$purity$purity, mutation_analysis_result$msi),
                           Wert2 = c(paste0(round(x = filt_result_td$tmb, digits = 2), "/Mb"),
                                     as.character(round(x = sum(as.numeric(mutation_analysis_result$mut_tab[, 2])), digits = 0)),
                                     "-",
                                     "", 
                                     "-",
                                     length(which(as.numeric(gsub(pattern = "%", replacement = "", x = filt_result_gd$table$Variant_Allele_Frequency)) > 5)),
                                     studies2, "", paste0("(", cnv_analysis_results$purity$ploidy, ")"), msi_help[1]),
                           Wert3 = c(help_tmb, "", help_ac3, help_hrd, help_cnv, help_germ, help_studies, "", help_tzg, msi_help[2])
    )
  } else {
    mut_tab1 <- data.frame(Eigenschaften = c("Mutationslast", "Anzahl somatischer Mutationen", "BRCAness", "HRD-Score (LoH|LST|TAI)",
                                             "Anzahl CNV- Regionen", "Anzahl seltener Keimbahnmutationen", "Studiencheck", "", "Tumorzellgehalt (Ploidität)", "MSI"),
                           Wert1 = c(paste0(round(x = filt_result_td_10$tmb, digits = 2), "/Mb"),
                                     round(x = sum(as.numeric(mutation_analysis_result_10$mut_tab[, 2])), digits = 0),
                                     "< 1.0%", cnv_analysis_results$hrd$score,
                                     paste0(round(x = dim(cnv_analysis_results$cnvs_annotated$CNVsAnnotated)[1], digits = 0), " Regionen"),
                                     length(which(as.numeric(gsub(pattern = "%", replacement = "", x = filt_result_gd$table$Variant_Allele_Frequency)) > 10)),
                                     studies3, studies_jd, cnv_analysis_results$purity$purity, mutation_analysis_result$msi),
                           Wert2 = c(paste0(round(x = filt_result_td$tmb, digits = 2), "/Mb"),
                                     as.character(round(x = sum(as.numeric(mutation_analysis_result$mut_tab[, 2])), digits = 0)),
                                     "-",
                                     "-",
                                     "-",
                                     length(which(as.numeric(gsub(pattern = "%", replacement = "", x = filt_result_gd$table$Variant_Allele_Frequency)) > 5)),
                                     studies2, "", paste0("(", cnv_analysis_results$purity$ploidy, ")"), msi_help[1]),
                           Wert3 = c(help_tmb, "", help_ac3, help_hrd, help_cnv, help_germ, help_studies, "", help_tzg, msi_help[2])
    )
  }
  colnames(mut_tab1) <- c("Eigenschaften", "Wert (VAF > 10%)", "Wert (VAF > 5%)", "Einordnung")
  
  if (!is.null(med_tmb$med) & !is.null(med_tmb$sd)) {
    tmb <- data.frame(Eigenschaft = "Mittlere TMB der Entität", Wert1 = entity,
                      Wert2 = paste0(med_tmb$med, " (", as.numeric(med_tmb$sd[1], digits = 2), "-", as.numeric(med_tmb$sd[2], digits = 2), ")"),
                      Wert3 = "")
    colnames(tmb) <- colnames(mut_tab1)
    mut_tab1 <- rbind(mut_tab1, tmb)
    mut_tab1 <- mut_tab1[c(1, 11, 2:10), ]
  } else if (!is.null(med_tmb$med)) {
    tmb <- data.frame(Eigenschaft = "Mittlere TMB der Entität", Wert1 = entity,
                      Wert2 = med_tmb$med,
                      Wert3 = "")
    colnames(tmb) <- colnames(mut_tab1)
    mut_tab1 <- rbind(mut_tab1, tmb)
    mut_tab1 <- mut_tab1[c(1, 11, 2:10), ]
  } else {
    tmb <- data.frame(Eigenschaft = "Mittlere TMB der Entität", Wert1 = "Nicht gelistet",
                      Wert2 = "-",
                      Wert3 = "")
    colnames(tmb) <- colnames(mut_tab1)
    mut_tab1 <- rbind(mut_tab1, tmb)
    mut_tab1 <- mut_tab1[c(1, 11, 2:10), ]
  }
  
  return(mut_tab1)
}

l_gen_nex <- function(df, type = "SNV") {
  if ("Gene.refGene" %in% colnames(df)) {
    df$Symbol <- df$Gene.refGene
  }
  if ( type == "SNV") {
    vec_gen_nex <- paste0(
      "{\\href{https://www.genomenexus.org/variant/",
      df$Chr, ":g.", df$Start, df$Ref, "\\%3E",
      df$Alt, "}{", df$Symbol, "}}"
    )
  } else if (type == "DEL") {
    vec_gen_nex <- paste0(
      "{\\href{https://www.genomenexus.org/variant/",
      df$Chr, ":g.", df$Start,
      df$Ref, "\\%3Edel}{",
      df$Symbol, "}}"
    )
  } else if (type == "INS") {
    vec_gen_nex <- paste0(
      "{\\href{https://www.genomenexus.org/variant/",
      df$Chr, ":g.", df$Start,
      "_", (df$Start + nchar(as.character(df$Alt))),
     "ins\\%3El", df$Alt, "}{", df$Symbol, "}}"
    )
  }
  return(vec_gen_nex)
}
meta <- function(df) {
  df$AAChange <- substr(x = df$AAChange, start = 3, stop = nchar(df$AAChange))
  hyper_refs <- paste0(
    "https://search.cancervariants.org/\\#", df$Symbol, "\\%20",
    df$AAChange
  )
  vec_meta <- paste0(
    "{\\href{", hyper_refs, "}{", df$AAChange, "}}"
  )
  vec_meta <- gsub(
    pattern = "_", replacement = "\\_",
    x = vec_meta, fixed = TRUE
  )
  return(vec_meta)
}
varsome <- function(df, mode) {
  baseURL <- "https://varsome.com/variant/hg19/"
  if (mode == "1") {
    completeVarsomeLink <- paste0(
      "{\\href{", baseURL, df$Chr, "\\%3A",df$Start,
      "\\%3A", df$Ref, "\\%3A", df$Alt, "}{Link}}"
    )
  } else if (mode == "2"){
    completeVarsomeLink <- paste0(
      "{\\href{", baseURL, df$Chr, "\\%3A",df$Start,
      "\\%3A", df$Ref, "\\%3A", df$Alt, "}{", df$ExonicFunc.refGene, "}}"
    )
  }
    completeVarsomeLink <- gsub(
      pattern = "-", replacement = "", fixed = TRUE,
      x = completeVarsomeLink
    )
  return(completeVarsomeLink)
}
acmg <- function(df) {
  df$InterVar <- gsub(
    pattern = "Pathogenic", replacement = 5,
    x = df$InterVar
  )
  df$InterVar <- gsub(
    pattern = "Likely pathogenic", replacement = 4,
    x = df$InterVar
  )
  df$InterVar <- gsub(
    pattern = "Uncertain significance", replacement = 3,
    x = df$InterVar
  )
  df$InterVar <- gsub(
    pattern = "Likely benign", replacement = 2,
    x = df$InterVar
  )
  df$InterVar <- gsub(
    pattern = "Benign", replacement = 1,
    x = df$InterVar
  )
  df$InterVar[is.na(df$InterVar)] <- "0"
  df$InterVar <- gsub(
    pattern = ".", replacement = 0,
    x = df$InterVar, fixed = TRUE
  )
  # ClinVar
  if ("CLNSIG" %in% colnames(df)) {
    df$CLNSIG_new <- df$CLNSIG
  } else {
    df$CLNSIG_new <- df$ClinVar
  }
  df$CLNSIG_new <- gsub(
    pattern = "\\x2c_other", fixed = TRUE,
    replacement = "", x = df$CLNSIG_new
  )
  df$CLNSIG_new <- gsub(
    pattern = "Pathogenic",
    replacement = 5, x = df$CLNSIG_new)
  df$CLNSIG_new <- gsub(
    pattern = "Likely_pathogenic",
    replacement = 4, x = df$CLNSIG_new)
  df$CLNSIG_new <- gsub(
    pattern = "Uncertain_significance",
    replacement = 3, x = df$CLNSIG_new)
  df$CLNSIG_new <- gsub(
    pattern = "Likely_benign",
    replacement = 2, x = df$CLNSIG_new)
  df$CLNSIG_new <- gsub(
    pattern = "Benign",
    replacement = 1, x = df$CLNSIG_new)
  df$CLNSIG_new <- gsub(
    pattern = "Conflicting_interpretations_of_pathogenicity",
    replacement = 4, x = df$CLNSIG_new)
  df$CLNSIG_new <- gsub(
    pattern = ".",
    replacement = 0, x = df$CLNSIG_new, fixed = TRUE)
  df$CLNSIG_new[is.na(df$CLNSIG_new)] <- 0
  df$CLNSIG <- gsub(
    pattern = "\\x2c_other", fixed = TRUE,
    replacement = "", x = df$CLNSIG)
  df$CLNSIG <- gsub(
    pattern = "Pathogenic",
    replacement = 5, x = df$CLNSIG)
  df$CLNSIG <- gsub(
    pattern = "Likely_pathogenic",
    replacement = 4, x = df$CLNSIG)
  df$CLNSIG <- gsub(
    pattern = "Uncertain_significance",
    replacement = 3, x = df$CLNSIG)
  df$CLNSIG <- gsub(
    pattern = "Likely_benign",
    replacement = 2, x = df$CLNSIG)
  df$CLNSIG <- gsub(
    pattern = "Benign",
    replacement = 1, x = df$CLNSIG)
  df$CLNSIG <- gsub(
    pattern = "Conflicting_interpretations_of_pathogenicity",
    replacement = "4*", x = df$CLNSIG)
  df$CLNSIG[is.na(df$CLNSIG)] <- "."
  # combined classification; use most severe
  df$Classification <- unlist(lapply(
    strsplit(x = pmax(df$InterVar,
    df$CLNSIG_new),
    split = "/"), function(x) max(x)
  ))
  df$InterVar <- gsub(
    pattern = 0, replacement = ".",
    x = df$InterVar, fixed = TRUE)
  vec_acmg <- paste0(df$InterVar, " | ", df$CLNSIG)
  return(cbind(vec_acmg, df$Classification))
}
revel <- function(df) {
  if("REVEL" %in% colnames(df)) {
    vec_rev <- round(as.numeric(df$REVEL), digits = 1)
  } else {
    vec_rev <- round(as.numeric(df$REVEL_score), digits = 1)
  }
  vec_rev[is.na(vec_rev)] <- "."
  vec_cat <- ifelse(test = vec_rev > 0.5, yes = "D", no = "N")
  vec_cat <- paste0(vec_cat, " (", vec_rev, ")")
  vec_cat[vec_rev == "."] <- "."
  
  return(list(revel = vec_rev, cat = vec_cat))
}

cosmic <- function(df) {
  if ("COSMIC" %in% colnames(df)) {
    vec_cos <- lapply(
      strsplit(x = as.character(df$COSMIC), split = "="),
      function(x) { return(unlist(x)[2]) }
    )
  } else {
    vec_cos <- lapply(
      strsplit(x = as.character(df$cosmic_coding), split = "="),
      function(x) { return(unlist(x)[2]) }
    )
  }
  vec_cos <- lapply(
    strsplit(x = as.character(vec_cos), split = ";"),
    function(x){return(unlist(x)[1])}
  )
  vec_cos <- gsub(pattern = ",", replacement = ", ", x = vec_cos)
  vec_cos <- gsub(pattern = "COSV", replacement = "", x = vec_cos)
  vec_cos[which(is.na(vec_cos))] <- "."
  return(vec_cos)
}
ex_func <- function(df) {
  if ("ExonicFunction" %in% colnames(df)) {
    vec_exf <- df$ExonicFunction
  } else {
    vec_exf <- df$ExonicFunc.refGene
  }
  vec_exf <- gsub(
    pattern = "nonsynonymous SNV",
    replacement = "nsSNV", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "stopgain",
    replacement = "SG", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "nonframeshift deletion",
    replacement = "nfsDel", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "frameshift deletion",
    replacement = "fsDel", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "nonframeshift insertion",
    replacement = "nfsIns", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "frameshift insertion",
    replacement = "fsIns", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "startloss",
    replacement = "SL", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "splice_region_variant&synonymous_variant",
    replacement = "splice", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "nonframeshift substitution",
    replacement = "nfsSub", x = vec_exf
  )
  vec_exf <- gsub(
    pattern = "missense_variant",
    replacement = "nsSNV", x = vec_exf
  )
  return(vec_exf)
}


highlight <- function(muts_tab, protocol) {
  # Select mutations in tumorsuppressors and oncogenes as well as potentially deleterious
  highlight <- muts_tab[which(
    muts_tab$is_oncogene == 1 | muts_tab$is_tumorsuppressor == 1 |
    muts_tab$InterVar_automated %in% c("Pathogenic", "Likely pathogenic") |
    muts_tab$CLNSIG %in% c(
      "Pathogenic",
      "Likely_pathogenic",
      "Conflicting_interpretations_of_pathogenicity"
    )
  ),
  c(
    "Gene.refGene",
    "AAChange",
    "Variant_Allele_Frequency",
    "InterVar_automated",
    "is_tumorsuppressor",
    "is_oncogene",
    "is_hotspot",
    "CLNSIG",
    "REVEL_score",
    "Chr",
    "Start",
    "Ref",
    "Alt"
  )]
  colnames(highlight) <- c(
    "Symbol",
    "AAChange",
    "VAF",
    "InterVar",
    "TSG",
    "OG",
    "HS",
    "CLNSIG",
    "REVEL",
    "Chr",
    "Start",
    "Ref",
    "Alt"
  )
   if (dim(highlight)[1] != 0) {
    ### Interactive links ###
    # Genome Nexus
    highlight$Gene.refGene_new <- highlight$Symbol
    if (length(which(highlight$Alt == "-" | highlight$Ref == "-")) > 0) {
      highlight$Gene.refGene_new[- which(
        highlight$Alt == "-" | highlight$Ref == "-"
      )] <- l_gen_nex(df = highlight[-which(
        highlight$Alt == "-" | highlight$Ref == "-"
      ), ], type = "SNV") 
      if (length(which(highlight$Alt == "-")) > 0){
        highlight$Gene.refGene_new[which(
          highlight$Alt == "-"
        )] <- l_gen_nex(df = highlight[which(
          highlight$Alt == "-"
        ), ], type = "DEL")
      }
      if (length(which(highlight$Ref == "-")) > 0) {
      highlight$Gene.refGene_new[which(
        highlight$Ref == "-"
      )] <- l_gen_nex(df = highlight[which(
        highlight$Ref == "-"
      )], type = "INS")
      }
    } else {
      highlight$Gene.refGene_new <- l_gen_nex(df = highlight, type = "SNV")

    }
    # Cancer Consortium Meta-Knowledgebase
    highlight$AAChange <- meta(df = highlight)

    # VarSome links
    highlight$Varsome <- varsome(df = highlight, mode = "1")

    # VAF
    highlight$VAF <- gsub(pattern = "%", replacement = "", x = highlight$VAF)
    if (protocol == "panelTumor") {
      highlight$VAF <- as.numeric(highlight$VAF)*100
    }

    # InterVar (ACMG)
    highlight$Classification <- acmg(df = highlight)[, 2]
    highlight$combineInterVarClinVar <- acmg(df = highlight)[, 1]

    # REVEL
    res_revel <- revel(df = highlight)
    highlight$REVEL <- res_revel$revel
    highlight$REVEL_cat <- res_revel$cat

    # Order table
    highlight <- highlight[order(
      as.numeric(highlight$VAF), decreasing = TRUE
    ), , drop = FALSE]
    highlight <- highlight[order(
      as.numeric(highlight$Classification), decreasing = TRUE
    ), , drop = FALSE]
    id_hs <- which(highlight$HS != 0)

    # combined ACMG classification
    highlight$Cat <- highlight$Classification
    highlight$Cat[duplicated(highlight$Classification)] <- "."
    highlight$Cat <- gsub(
      pattern = 5, replacement = "Pathogen", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 4, replacement = "Wahrsch. Pathogen", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 3, replacement = "VUS", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 2, replacement = "Wahrsch. Gutartig", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 1, replacement = "Gutartig", x = highlight$Cat
    )
    highlight$Cat <- gsub(
      pattern = 0, replacement = "Nicht Klassifiziert", x = highlight$Cat
    )

    # cancer genes
    highlight$Cancergene <- "."
    highlight$Cancergene[which(highlight$TSG == 1)] <- "TSG"
    highlight$Cancergene[which(highlight$OG == 1)] <- paste0(
      "OG",
      highlight$Cancergene[which(highlight$OG == 1)]
    )
    highlight$Cancergene <- gsub(
      pattern = "OGTSG",
      replacement = "both",
      x = highlight$Cancergene
    )

    # output
    highlight <- highlight[, c(
      "Cat",
      "combineInterVarClinVar",
      "REVEL_cat",
      "Gene.refGene_new",
      "AAChange",
      "VAF",
      "Cancergene",
      "Varsome"
    )]
    colnames(highlight) <- c(
      "Kategorie",
      "InterVar | ClinVar",
      "REVEL",
      "Gen",
      "AA-Austausch",
      "VAF [\\%]",
      "Cancergene",
      "VarSome"
    )

  } else {
    highlight <- data.frame()
    id_hs <- NULL
  }
    return(list(highlight = highlight, id_hs = id_hs))
}

highlight_detail <- function(muts_tab, Mode = "Tumor", protocol) {
 highlight <- muts_tab
  if (dim(muts_tab)[1] == 0) {
     highlight <- data.frame()
     id_hs <- c()
    } else {
      # Genome Nexus
      highlight$Gene.refGene <- unlist(
        lapply(
          strsplit(
            x = as.character(highlight$Gene.refGene), split = ";", fixed = TRUE
          ), function(x) { return(x[[1]]) }
        )
      )
      highlight$Gene.refGene_new <- highlight$Gene.refGene 
      if (length(which(highlight$Alt == "-" | highlight$Ref == "-")) > 0) {
        highlight$Gene.refGene_new[-which(
          highlight$Alt == "-" | highlight$Ref == "-"
        )] <- l_gen_nex(df = highlight[-which(
          highlight$Alt == "-" | highlight$Ref == "-"
        ), ], type = "SNV") 
        if (length(which(highlight$Alt == "-")) > 0) {
          highlight$Gene.refGene_new[which(
            highlight$Alt == "-"
          )] <- l_gen_nex(df = highlight[which(
            highlight$Alt == "-"
          ), ], type = "DEL") 
        }
        if (length(which(highlight$Ref == "-")) > 0) {
          highlight$Gene.refGene_new[which(
            highlight$Ref == "-"
          )] <- l_gen_nex(df = highlight[which(
            highlight$Ref == "-"
          ), ], type = "INS")
        }
      } else {
        highlight$Gene.refGene_new <- l_gen_nex(df = highlight, type = "SNV")
      }
      # Cancer Consortium Meta-Knowledgebase
      highlight$AAChange <- meta(df = highlight)

      # Function
      highlight$ExonicFunc.refGene <- ex_func(df = highlight)

      # VarSome links
      highlight$Varsome <- varsome(df = highlight, mode = "2")

      # VAF
      if (Mode %in% c("Tumor", "Germline")) {
        highlight$VAF <- gsub(
          pattern = "%", replacement = "",
          x = highlight$Variant_Allele_Frequency
        )
        if (protocol == "panelTumor") {
          highlight$VAF <- as.numeric(highlight$VAF)*100
        }
        highlight$VAF <- paste0(
          highlight$VAF, " (", highlight$Variant_Reads, ")"
        )
      } else if (Mode == "LoH") {
        highlight$VAF_Normal <- gsub(
          pattern = "%", replacement = "",
          x = highlight$VAF_Normal, fixed = TRUE
        )
        highlight$VAF_Tumor <- gsub(
          pattern = "%",replacement = "",
          x = highlight$VAF_Tumor, fixed = TRUE
        )
        highlight$VAF_normal <- paste0(
          highlight$VAF_Normal, " (", highlight$Count_Normal, ")"
        )
        highlight$VAF_tumor <- paste0(
          highlight$VAF_Tumor, " (", highlight$Count_Tumor, ")"
        )
      }
      # InterVar (ACMG)
      highlight$combineInterVarClinVar <- acmg(df = highlight)[, 1]
      highlight$Classification <- acmg(df = highlight)[, 2]

      # REVEL
      res_revel <- revel(df = highlight)
      highlight$REVEL <- res_revel$revel
      highlight$REVEL_cat <- res_revel$cat

      # COSMIC
      highlight$cosmic <- cosmic(df = highlight)

      # Order table
      if (Mode %in% c("Tumor", "Germline")) {
        highlight <- highlight[order(
          as.numeric(gsub(
            pattern = "%", replacement = "",
            x = highlight$Variant_Allele_Frequency
          )),
          decreasing = TRUE
        ), , drop = FALSE]
      } else if (Mode == "LoH") {
        highlight <- highlight[order(
          as.numeric(highlight$VAF_Tumor), decreasing = TRUE
        ), , drop = FALSE]
      }
      highlight <- highlight[order(
        as.numeric(highlight$Classification), decreasing = TRUE
      ), , drop = FALSE]
      id_hs <- which(highlight$is_hotspot != 0)

      # cancer genes
      highlight$Cancergene <- "."
      highlight$Cancergene[which(highlight$is_tumorsuppressor == 1)] <- "TSG"
      highlight$Cancergene[which(highlight$is_oncogene == 1)] <- "OG"
      highlight$Cancergene[which(
        highlight$is_oncogene == 1 & highlight$is_tumorsuppressor == 1
      )] <- "both"

      # Population frequency
      highlight$AF_nfe <- as.character(format(
        as.numeric(highlight$AF_nfe), scientific = TRUE, digits = 2
      ))
      highlight$AF_nfe[is.na(highlight$AF_nfe)] <- "."

      # output
      if (Mode == "Tumor") {
        muts_tab <- highlight[, c("Gene.refGene_new", "AAChange", "Varsome", "VAF", "AF_nfe", "combineInterVarClinVar", "REVEL_cat", "cosmic", "Cancergene")]
        colnames(muts_tab) <- c(
          "Gen",
          "AA-Austausch",
          "Funktion / VarSome",
          "VAF [\\%] (Coverage)",
          "MAF",
          "InterVar | ClinVar",
          "REVEL",
          "Cosmic",
          "Cancergene"
        )
      } else if (Mode == "LoH") {
        muts_tab <- highlight[, c("Gene.refGene_new", "AAChange", "Varsome", "VAF_tumor", "VAF_normal", "AF_nfe", "combineInterVarClinVar", "REVEL_cat", "cosmic", "Cancergene")]
        colnames(muts_tab) <- c(
          "Gen",
          "AA-Austausch",
          "Funktion / VarSome",
          "VAF [\\%] (Coverage) Tumor",
          "VAF [\\%] (Coverage) Keimbahn",
          "MAF",
          "InterVar | ClinVar",
          "REVEL",
         "Cosmic",
         "Cancergene"
        )
      } else if (Mode == "Germline") {
        muts_tab <- highlight[, c("Gene.refGene_new", "AAChange", "Varsome", "VAF", "AF_nfe", "combineInterVarClinVar", "REVEL_cat", "cosmic", "Cancergene")]
        colnames(muts_tab) <- c(
          "Gen",
          "AA-Austausch",
          "Funktion / VarSome",
          "VAF [\\%] (Coverage)",
          "MAF",
          "InterVar | ClinVar",
          "REVEL",
          "Cosmic",
          "Cancergene"
        )
      }
    }
    return(list(muts_tab = muts_tab, id_hs = id_hs))
 }
 
summary_quality <- function(stats, protocol) {
  if (protocol == "somaticGermline" | protocol == "somatic") {
    q_t <- c()
    if (round(x = sum(
      stats$cover$cov[[2]][, 2] * stats$cover$cov[[2]][, 5]
    ), digits = 2) < 80) {
      q_t[1] <- "Akzeptabel" } else {
        q_t[1] <- "Sehr gut"
      }
    if (round(stats$cover$perc[[2]][1], digits = 2) * 100 > 80) {
      q_t[2] <- "Sehr gut"
    } else { q_t[2] <- "Akzeptabel" }
    if (round(stats$cover$perc[[2]][2], digits = 2) * 100 > 80) {
      q_t[3] <- "Sehr gut"
    } else {q_t[3] <- "Akzeptabel" }
    if (as.numeric(stats$avreads$tin) > 100) {
      q_t[4] <- "Sehr gut"
    } else {q_t[4] <- "Akzeptabel" }
    if (as.numeric(
      stats$qc_check$gc_content[[2]]
    ) >= 40 && as.numeric(stats$qc_check$gc_content[[2]]) <= 60) {
      q_t[5] <- "Sehr gut"
    } else {q_t[5] <- "Akzeptabel"}
    if (round(stats$qc_check$mean_QC[[2]], digits = 2) > 30) {
      q_t[6] <- "Sehr gut"
    } else {q_t[6] <- "Akzeptabel" }
    qp_t <- c(
      round(
        x = sum(
          stats$cover$cov[[2]][, 2] * stats$cover$cov[[2]][, 5]
        ), digits = 2
      ), paste0(round(stats$cover$perc[[2]][1], digits = 2) * 100, "%"),
      paste0(round(stats$cover$perc[[2]][2], digits = 2) * 100, "%"),
      as.numeric(stats$avreads$tin),
      as.numeric(stats$qc_check$gc_content[[2]]),
      round(stats$qc_check$mean_QC[[2]], digits = 2)
    )
    names(qp_t) <- c(
      "Mittlere Coverage",
      "Coverage > 8",
      "Coverage > 40",
      "Insertlänge",
      "GC-Anteil",
      "Mittlere Qualität"
    )
    if (length(which(q_t != "Sehr gut") > 0)) {
      warn_t <- paste0(
        names(qp_t)[which(q_t != "Sehr gut")],
        ": ", qp_t[which(q_t != "Sehr gut")]
      )
      q_t1 <- paste0(warn_t, collapse = ", ")
    } else {
      q_t1 <- "Keine."
    }
    q_n <- c()
    if (round(
      x = sum(stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]),
      digits = 2
    ) < 80) {
      q_n[1] <- "Akzeptabel"} else {
        q_n[1] <- "Sehr gut" 
      }
    if (round(stats$cover$perc[[1]][1], digits = 2) * 100 > 80) {
      q_n[2] <- "Sehr gut"
    } else { q_n[2] <- "Akzeptabel"
    }
    if (round(stats$cover$perc[[1]][2], digits = 2) * 100 > 80) {
      q_n[3] <- "Sehr gut"
    } else { q_n[3] <- "Akzeptabel"
    }
    if (as.numeric(stats$avreads$gin) > 100) {
      q_n[4] <- "Sehr gut"
    } else { q_n[4] <- "Akzeptabel"
    }
    if (as.numeric(
      stats$qc_check$gc_content[[1]]
    ) >= 40 && as.numeric(
      stats$qc_check$gc_content[[1]]
    ) <= 60) {
      q_n[5] <- "Sehr gut"
    } else {
      q_n[5] <- "Akzeptabel"
    }
    if (round(stats$qc_check$mean_QC[[1]], digits = 2) > 30) {
      q_n[6] <- "Sehr gut"
    } else {q_n[6] <- "Akzeptabel"
    }
    qp_n <- c(
      round(
        x = sum(stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]),
        digits = 2
      ),
      paste0(round(stats$cover$perc[[1]][1], digits = 2) * 100 , "%"),
      paste0(round(stats$cover$perc[[1]][2], digits = 2) * 100, "%"),
      as.numeric(stats$avreads$gin),
      as.numeric(stats$qc_check$gc_content[[1]]),
      round(stats$qc_check$mean_QC[[1]], digits = 2)
    )
    names(qp_n) <- c(
      "Mittlere Coverage",
      "Coverage > 8",
      "Coverage > 40",
      "Insertlänge",
      "GC-Anteil",
      "Mittlere Qualität"
    )
    if (length(which(q_n != "Sehr gut") > 0)) {
      warn_n <- paste0(
        names(qp_n)[which(q_n != "Sehr gut")],
        ": ",
        qp_n[which(q_n != "Sehr gut")]
      )
      q_n1 <- paste0(warn_n, collapse = ", ")
    } else {
      q_n1 <- "Keine"
    }
    tab <- rbind(c("Tumor", q_t1), c("Keimbahn", q_n1))
    colnames(tab) <- c("Probe" , "Auff\"alligkeiten")
    return(tab)
  }
  if (protocol == "panelTumor" | protocol == "tumorOnly") {
    q_t <- c()
    if (round(x = sum(
      stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]
    ), digits = 2) < 150) {
      q_t[1] <- "Akzeptabel" } else {
        q_t[1] <- "Sehr gut"
      }
    if (round(stats$cover$perc[[1]][1], digits = 2) * 100 > 90) {
      q_t[2] <- "Sehr gut"
    } else { q_t[2] <- "Akzeptabel" }
    if (round(stats$cover$perc[[1]][2], digits = 2) * 100 > 90) {
      q_t[3] <- "Sehr gut"
    } else {q_t[3] <- "Akzeptabel" }
    if (as.numeric(stats$avreads$tin) > 100) {
      q_t[4] <- "Sehr gut"
    } else {q_t[4] <- "Akzeptabel" }
    if (as.numeric(
      stats$qc_check$gc_content[[1]]
    ) >= 40 && as.numeric(stats$qc_check$gc_content[[1]]) <= 60) {
      q_t[5] <- "Sehr gut"
    } else {q_t[5] <- "Akzeptabel"}
    if (round(stats$qc_check$mean_QC[[1]], digits = 2) > 30) {
      q_t[6] <- "Sehr gut"
    } else {q_t[6] <- "Akzeptabel" }
    qp_t <- c(
      round(
        x = sum(
          stats$cover$cov[[1]][, 2] * stats$cover$cov[[1]][, 5]
        ), digits = 2
      ), paste0(round(stats$cover$perc[[1]][1], digits = 2) * 100, "%"),
      paste0(round(stats$cover$perc[[1]][2], digits = 2) * 100, "%"),
      as.numeric(stats$avreads$tin),
      as.numeric(stats$qc_check$gc_content[[1]]),
      round(stats$qc_check$mean_QC[[1]], digits = 2)
    )
    names(qp_t) <- c(
      "Mittlere Coverage",
      "Coverage > 50",
      "Coverage > 150",
      "Insertlänge",
      "GC-Anteil",
      "Mittlere Qualität"
    )
    if (length(which(q_t != "Sehr gut") > 0)) {
      warn_t <- paste0(
        names(qp_t)[which(q_t != "Sehr gut")],
        ": ", qp_t[which(q_t != "Sehr gut")]
      )
      q_t1 <- paste0(warn_t, collapse = ", ")
    } else {
      q_t1 <- "Keine."
    }
    tab <- rbind(c("Tumor", q_t1))
    colnames(tab) <- c("Probe" , "Auff\"alligkeiten")
    return(tab)
  }
}

sum_muts <- function(tmp) {
  colnames(tmp) <- c(
    "Mutationstyp",
    "Anzahl",
    "Zygosit\"at",
    "Tumorsuppressoren",
    "Onkogene",
    "Hotspots"
  )
  tmp[c(1, 4), 3] <- "homozygot"
  tmp[c(2, 5), 3] <- "heterozygot"

  return(tmp)
}

cnv_cg <- function(gene_loci, type = "OG") {
  idl <- which(as.numeric(gene_loci$cn) < 1)
  idg <- which(as.numeric(gene_loci$cn) > 4)
  id_d <- which(is.na(as.numeric(gene_loci$cn)))
  if (type == "TSG") {
    colnames_df <- c("TSG", "CN", "CN-Typ", "Status")
  } else {
    colnames_df <- c("OG", "CN", "CN-Typ", "Status")
  }

  if(length(idl) > 0) {
    tmpl <- gene_loci[idl, c(4:6)]
    tmpl$Status <- "Loss"
    tmpl <- tmpl[order(as.numeric(tmpl$cn), decreasing = FALSE), ]
  } else {
    tmpl <- data.frame("", "", "", "")
  }
  colnames(tmpl) <- colnames_df

  if(length(idg) > 0) {
    tmpg <- gene_loci[idg, c(4:6)]
    tmpg$Status <- "Gain"
    tmpg <- tmpg[order(as.numeric(tmpg$cn), decreasing = TRUE), ]
  } else {
    tmpg <- data.frame("", "", "", "")
  }
  colnames(tmpg) <- colnames_df

  if (length(id_d) > 0){
    tmpd <- gene_loci[id_d, c(4:6)]
  } else {
    tmpd <- data.frame("", "", "")
  }
  colnames(tmpd) <- colnames_df[1:3]

  if((tmpg[1,1] != "") & (tmpl[1,1] != "")){
    id_l <- which(tmpl$`CN-Typ` < 1)
    id_g <- which(tmpg$`CN-Typ` > 7)
  }

  return(list(Gains = tmpg, Losses = tmpl, Different = tmpd))
}

cnv_panel <- function(cnv_results) {
  tmp <- cnv_analysis_results$out[, c("Gene", "CopyNumber", "Status", "Type", "Cancergene")]
  colnames(tmp) <- c("Gen", "Kopien", "Status", "CN-Typ", "Cancergene")
  return(tmp)
}

pathws_cnv <- function(df) {
  if (sum(unlist(lapply(df, function(x){return(dim(x)[1])}))) == 0) {
    return(NULL)
  } else {
    df <- lapply(
      1:length(df),
      function(id) {
        if(nrow(df[[id]])) {
          cbind(df[[id]], names(df)[id])
        }
      }
    )
    table_new <- lapply(df, function(x){
      if (!is.null(x)) {
        new_mat <- as.data.frame(
          matrix(
            nrow = length(unique(x[, 2])),
            ncol = 4,
            data = 0
          )
        )
        colnames(new_mat) <- c("Signalweg", "Status", "Kopien", "Gene")
        new_mat[, 1] <- x[1, 4]
        new_mat[, 3] <- unique(x[, 2])[
          order(unique(x[, 2]), decreasing = FALSE)
        ]
        for (i in 1:dim(new_mat)[1]) {
          new_mat[i, 4] <- paste(
            x[which(as.numeric(x[, 2]) == new_mat[i, 3]), 1],
            collapse = ", "
          )
          new_mat[i, 2] <- ifelse(
            test = as.numeric(new_mat[i, 3]) < 2,
            yes = "Loss",
            no = "Gain"
          )
        }
        return(new_mat)
      }
    })
    df <- data.frame(
      rbind(
        table_new[[1]],
        table_new[[2]],
        table_new[[3]],
        table_new[[4]],
        table_new[[5]]
      )
    )
    df[which(duplicated(df$Signalweg)), 1] <- "."
    df[which(df[, ] == "ddr"), 1] <- "DNA DamResp"
    df[which(df[, ] == "pam"), 1] <- "PI3K-AKT-mTOR"
    df[which(df[, ] == "rme"), 1] <- "RAF-MEK-ERK"
    df[which(df[, ] == "tyk"), 1] <- "Tyr Kin"
    df[which(df[, ] == "cec"), 1] <- "Cell Cycle"
    return(df)
  }
}

pthws_mut <- function(df, protocol) {
  id_to <- which(df$Pathway == "Topart")
  if (length(id_to) != 0) {
    df <- df[-c(id_to:dim(df)[1]), ]
  }
  if (class(df) == "list" | is.null(df)) {
    if (is.null(df)) {
      return(NULL)
    } else if (sum(lapply(df, function(x){return(dim(x)[1])})) == 0){
      return(NULL)
    }
  } else if (class(df) == "data.frame") {
    df <- df[, c(
      "Pathway",
      "Symbol",
      "AAChange",
      "ExonicFunc",
      "VAF",
      "Reads",
      "MAF",
      "InterVar_automated",
      "CLNSIG",
      "REVEL_score",
      "COSMIC"
    )]
    colnames(df) <- c(
      "Pathway",
      "Symbol",
      "AAChange",
      "ExonicFunction",
      "VAF",
      "Reads",
      "MAF",
      "InterVar",
      "ClinVar",
      "REVEL",
      "COSMIC"
    )
    # Pathway
    df$Pathway <- gsub(
      pattern = "DNA Damage Response",
      replacement = "DNA DamResp",
      x = df$Pathway
    )
    df$Pathway <- gsub(
      pattern = "Tyrosine Kinases",
      replacement = "Tyr Kin",
      x = df$Pathway
    )
    # AAChange
    df$AAChange <- substr(
      x = df$AAChange, start = 3, stop = nchar(df$AAChange)
    )
    # VAF
    df$VAF <- gsub(
      pattern = "%", replacement = "", x = df$VAF, fixed = TRUE
    )
    if (protocol == "panelTumor") {
      df$VAF <- as.numeric(df$VAF)*100
    }
    df$VAF <- paste0(df$VAF, " (", df$Reads, ")")
    # AAChange
    df$AAChange <- gsub(
      pattern = "_", replacement = "\\_", x = df$AAChange, fixed = TRUE
    )
    # COSMIC
    df$cosmic <- cosmic(df)
    # Function
    df$ExonicFunction <- ex_func(df)
    # REVEL
    res_revel <- revel(df)
    df$REVEL <- res_revel$revel
    df$REVEL_cat <- res_revel$cat
    # InterVar (ACMG)
    df$combineInterVarClinVar <- acmg(df = df)[, 1]
    # Population frequency
    df$MAF <- as.character(
      format(as.numeric(df$MAF), scientific = TRUE, digits = 2)
    )
    # output
    df <- df[, c(1:5, 7, 14, 13, 12)]
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    colnames(df) <- c(
      "Signalweg",
      "Gen",
      "AA-Austausch",
      "Funktion",
      "VAF [\\%] (Coverage)",
      "MAF",
      "InterVar | ClinVar",
      "REVEL",
      "Cosmic"
    )
    return(df)
  }
}

topart_mut <- function(df, protocol) {
  id_to <- which(df$Pathway == "Topart")
  if (length(id_to) == 0) {
    return(NULL)
  } else {
    df <- df[id_to:dim(df)[1], c(
      "Pathway",
      "Symbol",
      "AAChange",
      "ExonicFunc",
      "VAF",
      "Reads",
      "MAF",
      "InterVar_automated",
      "CLNSIG",
      "REVEL_score",
      "COSMIC"
    )]
    colnames(df) <- c(
      "Pathway",
      "Symbol",
      "AAChange",
      "ExonicFunction",
      "VAF",
      "Reads",
      "MAF",
      "InterVar",
      "ClinVar",
      "REVEL",
      "COSMIC"
    )
    # AAChange
    df$AAChange <- substr(
      x = df$AAChange, start = 3, stop = nchar(df$AAChange)
    )
    # VAF
    df$VAF <- gsub(pattern = "%", replacement = "", x = df$VAF, fixed = TRUE)
    if (protocol == "panelTumor") {
      df$VAF <- as.numeric(df$VAF)*100
    }
    df$VAF <- paste0(df$VAF, " (", df$Reads, ")")
    # AAChange
    df$AAChange <- gsub(
      pattern = "_", replacement = "\\_", x = df$AAChange, fixed = TRUE
    )
    # COSMIC
    df$cosmic <- cosmic(df)
    # Function
    df$ExonicFunction <- ex_func(df)
    # REVEL
    res_revel <- revel(df)
    df$REVEL <- res_revel$revel
    df$REVEL_cat <- res_revel$cat
    # InterVar (ACMG)
    df$combineInterVarClinVar <- acmg(df = df)[, 1]
    # Population frequency
    df$MAF <- as.character(format(
      as.numeric(df$MAF), scientific = TRUE, digits = 2
    ))

    # output
    df <- df[, c(1:5, 7, 14, 13, 12)]
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    dup <- which(duplicated(df))
    if (length(dup) > 0){
      df[dup, 4] <- paste0(df[dup, 4], ".")
    }
    colnames(df) <- c(
      "Signalweg",
      "Gen",
      "AA-Austausch",
      "Funktion",
      "VAF [\\%] (Coverage)",
      "MAF",
      "InterVar | ClinVar",
      "REVEL",
      "Cosmic"
    )
    return(df)
  }
}

med_tmb <- function(entity) {
  if (entity == "ACC") {
      med <- 0.82
      sd <- c(0.55, 1.86)
  } else if (entity == "BLCA") {
      med <- 5.26
      sd <- c(2.79, 9.36)
  } else if (entity == "BRCA") {
    med <- 2.63
    sd <- c(0.20, 290.80)
  } else if (entity == "CESC") {
      med <- 3.32
      sd <- c(2.00, 6.07)
  } else if (entity == "CHOL") {
      med <- 1.71
      sd <- c(1.25, 2.80)
  } else if (entity == "COAD") {
      med <- 3.68
      sd <- c(2.67, 5.71)
  } else if (entity == "ESCA") {
      med <- 4.18
      sd <- c(3.16, 6.05)
  } else if (entity == "HNSC") {
      med <- 3.34
      sd <- c(2.05, 5.51)
  } else if (entity == "KIRC") {
      med <- 1.63
      sd <- c(1.21, 2.13)
  } else if (entity == "KIRP") {
      med <- 2.13
      sd <- c(1.34, 2.89)
  } else if (entity == "LIHC") {
      med <- 2.89
      sd <- c(2.03, 3.95)
  } else if (entity == "LUAD") {
      med <- 6.17
      sd <- c(2.47, 12.66)
  } else if (entity == "LUSC") {
      med <- 7.13
      sd <- c(5.01, 10.47)
  } else if (entity == "MESO") {
      med <- 0.95
      sd <- c(0.66, 1.16)
  } else if (entity == "OV") {
      med <- 2.13
      sd <- c(1.50, 3.16)
  } else if (entity == "PAAD") {
      med <- 1.09
      sd <- c(0.76, 1.47)
  } else if (entity == "SKCM") {
      med <- 13.09
      sd <- c(5.91, 28.09)
  } else if (entity == "STAD") {
      med <- 3.5
      sd <- c(2.16, 7.45)
  } else if (entity == "UCEC") {
      med <- 2.63
      sd <- c(1.58, 19.47)
  } else if (entity == "UCS") {
      med <- 1.37
      sd <- c(1.21, 1.79)
  } else if (entity == "UVM") {
      med <- 0.34
      sd <- c(0.26, 0.45)
  } else {
      sd <- NULL
      if (entity == "STA") {
        med <- 3.3
      } else if (entity == "US") {
        med <- 2.6
      } else if (entity == "STR") {
        med <- 2.5
      } else if (entity %in% c(
        "STL", "STU", "BOS", "STS", "STRE",
        "STMPNST", "UPG", "ULS"
      )) {
        med <- 2.5
      } else if (entity == "STMS") {
        med <- 2.2
      } else if (entity == "SIG") {
        med <- 1.8
      } else if (entity == "SG") {
        med <- 1.8
      } else if (entity %in% c(
        "STC", "SSTSFT", "STFS", "UESS",
        "BCS", "STES", "STLS", "STDSRCT",
        "STSS", "STRSA"
      )) {
        med <- 1.7
      } else if (entity == "BC") {
        med <- 1.3
      } else if (entity == "STF") {
        med <- 0.9
      } else if (entity == "MRT") {
        med <- 0.2
      } else {
        med <- NULL
      }
  }
  return(list(med = med, sd = sd))
}
