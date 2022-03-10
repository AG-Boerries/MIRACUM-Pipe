mut_sig_wCI <- function(
  vcf_file = NULL,
  cutoff = 0.01,
  sample = NULL,
  sureselect_type,
  path_script,
  ref_genome,
  target_capture_cor_factors,
  path_output,
  sample_name,
  outfile_cbioportal,
  vaf
  ) {
  library(getopt)
  library(dplyr)
  library(magrittr)
  library(rtracklayer)
  library(Rsamtools)
  library(YAPSA)
  library(IRanges)
  library(data.table)

  load(target_capture_cor_factors)

  if (sureselect_type == "V6"){
    targetCapture <- "Agilent6withoutUTRs"
  } else if (sureselect_type == "V6UTR"){
    targetCapture <- "Agilent6withUTRs"
  } else if (sureselect_type == "V5UTR"){
    targetCapture <- "Agilent5withUTRs"
  } else if (sureselect_type %in% names(targetCapture_cor_factors)){
    targetCapture <- sureselect_type
  } else {
      error("Unsupported Capture Region Kit.")
  }
  
  print(ref_genome)
  refGenome_path <- file.path(ref_genome)# hg19
  refGenome <- FaFile(refGenome_path)
  print(refGenome)

  commonCutoff <- cutoff

  mostInformativeCombo <- "Valid_norm"
  sequencingType <- "WES"
  wordLength <- 3

  # read in the snv data
  vcf_like_df <- read.delim(
    file = vcf_file,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")
  )
  id <- which(vcf_like_df$FILTER == "PASS")
  vcf_like_df <- vcf_like_df[id, ]
  vcf_like_df$FREQ <- as.numeric(
    gsub(
      pattern = "%", replacement = "", fixed = TRUE,
      unlist(
        lapply(
          strsplit(vcf_like_df$TUMOR, split = ":", fixed = TRUE),
          function(i) i[6])
      )
    )
  )
  id <- which(vcf_like_df$FREQ >= vaf)
  if(length(id) > 0) {
    vcf_like_df <- vcf_like_df[id, ]
  }
  number_of_SNVs <- nrow(vcf_like_df)

  # load data for signatures and cutoffs from the package
  data(sigs)
  data(cutoffs)
  CosmicValid_absCutoffVector <- cutoffCosmicValid_abs_df[6, ]
  CosmicValid_normCutoffVector <- cutoffCosmicValid_rel_df[6, ]
  CosmicValid_genCutoffVector <- rep(commonCutoff, ncol(AlexCosmicValid_sig_df))
  CosmicArtif_absCutoffVector <- cutoffCosmicArtif_abs_df[6, ]
  CosmicArtif_normCutoffVector <- cutoffCosmicArtif_rel_df[6, ]
  CosmicArtif_genCutoffVector <- rep(commonCutoff, ncol(AlexCosmicArtif_sig_df))
  
  # filter for onTarget if necessary
  targetCapture_gr <- NULL
  if(sequencingType == "WES" & !is.null(targetCapture_gr)){
    vcf_like_gr <- makeGRangesFromDataFrame(
      vcf_like_df, start.field = "POS", end.field = "POS", 
      keep.extra.columns = T)
    filtered_vcf_like_gr <- subsetByOverlaps(vcf_like_gr, targetCapture_gr)
    filtered_vcf_like_df <- as.data.frame(filtered_vcf_like_gr)
    filtered_vcf_like_df <- filtered_vcf_like_df[, -c(3,4,5)]
    names(filtered_vcf_like_df) = names(vcf_like_df)
  } else if(sequencingType == "WES") {
    filtered_vcf_like_df <- vcf_like_df
  }
  filtered_number_of_SNVs <- nrow(filtered_vcf_like_df)
  
  # create the mutational catalog
  mutation_catalogue_list <- create_mutation_catalogue_from_df(
    this_df = filtered_vcf_like_df,
    this_seqnames.field = "CHROM",
    this_refGenome = refGenome,
    this_wordLength = wordLength,
    this_rownames = rownames(AlexCosmicValid_sig_df),
    this_verbose = 0)
  mutation_catalogue_df <- as.data.frame(mutation_catalogue_list$matrix)
  
  # correct for triplet content
    cor_list <- targetCapture_cor_factors[[targetCapture]]
    corrected_catalogue_df <- normalizeMotifs_otherRownames(mutation_catalogue_df, cor_list$rel_cor)
  
  # run analysis of mutational signatures
  CosmicValid_abs_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df = corrected_catalogue_df,
    in_signatures_df = AlexCosmicValid_sig_df,
    in_cutoff_vector = CosmicValid_absCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind = AlexCosmicValid_sigInd_df)
  CosmicValid_norm_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df = corrected_catalogue_df,
    in_signatures_df = AlexCosmicValid_sig_df,
    in_cutoff_vector = CosmicValid_normCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind = AlexCosmicValid_sigInd_df)
  CosmicValid_gen_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df = corrected_catalogue_df,
    in_signatures_df = AlexCosmicValid_sig_df,
    in_cutoff_vector = CosmicValid_genCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind = AlexCosmicValid_sigInd_df)
  CosmicArtif_abs_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df = corrected_catalogue_df,
    in_signatures_df = AlexCosmicArtif_sig_df,
    in_cutoff_vector = CosmicArtif_absCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind = AlexCosmicArtif_sigInd_df)
  CosmicArtif_norm_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df = corrected_catalogue_df,
    in_signatures_df = AlexCosmicArtif_sig_df,
    in_cutoff_vector = CosmicArtif_normCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind = AlexCosmicArtif_sigInd_df)
  CosmicArtif_gen_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df = corrected_catalogue_df,
    in_signatures_df = AlexCosmicArtif_sig_df,
    in_cutoff_vector = CosmicArtif_genCutoffVector,
    in_filename = NULL,
    in_method = "abs",
    in_sig_ind = AlexCosmicArtif_sigInd_df)
  
  # Merge results of different decompositions
  exposures_list <- list(Valid_gen = CosmicValid_gen_LCDlist$exposures,
                         Valid_abs = CosmicValid_abs_LCDlist$exposures,
                         Valid_norm = CosmicValid_norm_LCDlist$exposures,
                         Artif_gen = CosmicArtif_gen_LCDlist$exposures,
                         Artif_abs = CosmicArtif_abs_LCDlist$exposures,
                         Artif_norm = CosmicArtif_norm_LCDlist$exposures)
  name_list <- names(exposures_list)
  exposures_list <- lapply(seq_along(exposures_list), FUN = function(current_ind){
    current_df <- exposures_list[[current_ind]]
    names(current_df) <- name_list[current_ind]
    return(current_df)
  })
  names(exposures_list) <- name_list
  my_exposures_df <- merge_exposures(exposures_list, AlexCosmicValid_sig_df)
  my_normExposures_df <- normalize_df_per_dim(my_exposures_df, in_dimension = 2)
  my_sigInd_df <- do.call("rbind",list(
    CosmicValid_gen_LCDlist$out_sig_ind[,c(1:4)],
    CosmicValid_abs_LCDlist$out_sig_ind[,c(1:4)],
    CosmicValid_norm_LCDlist$out_sig_ind[,c(1:4)],
    CosmicArtif_gen_LCDlist$out_sig_ind[,c(1:4)],
    CosmicArtif_abs_LCDlist$out_sig_ind[,c(1:4)],
    CosmicArtif_norm_LCDlist$out_sig_ind[,c(1:4)]))
  my_sigInd_df <- my_sigInd_df[!duplicated(my_sigInd_df$sig),]
  if(length(grep("unknown|unassigned",my_sigInd_df$sig))>0){
    choice_ind <- which(my_sigInd_df$sig %in% c("unknown","unassigned"))
    my_sigInd_df$index[choice_ind] <- 
      seq(from=max(my_sigInd_df$index)+1,
          to=max(my_sigInd_df$index)+length(choice_ind), by=1)
  }
  my_sigInd_df <- my_sigInd_df[order(my_sigInd_df$index),]
  
  mostInformativeVector <- rep(NA, ncol(my_exposures_df))
  names(mostInformativeVector) <- names(exposures_list)
  mostInformativeVector[mostInformativeCombo] <- "most informative"
  annotation_df <- data.frame(
    signatureSet = gsub("_.*$", "", names(exposures_list)),
    algorithm = rep(c("nnls", "YAPSA", "YAPSA"), 2),
    trainingSet = rep(c(NA, "abs. exposures", 
                        "norm. exposures"), 2),
    mostInformative = mostInformativeVector
  )
  annotation_col <- list(signatureSet = c("Valid" = "blue", "Artif" = "red"),
                         algorithm = c("nnls" = "yellow",
                                       # "deconstructSigs" = "purple",
                                       "YAPSA" = "green"),
                         trainingSet = c("abs. exposures" = "black",
                                         "norm. exposures" = "grey40"),
                         mostInformative = c("most informative" = "black"))
  my_labels <- paste(my_sigInd_df$sig, my_sigInd_df$process, sep = " ")
  temp_exposures_df <- my_exposures_df
  rownames(temp_exposures_df) <- my_labels
  temp_sigInd_df <- my_sigInd_df
  temp_sigInd_df$sig <- my_labels
  subgroups_df <- data.frame(
    PID = names(exposures_list), subgroup = "test", sum = 1, compl_sum = 1, 
    index = seq_along(exposures_list), col = "#FF0000FF")
  # Compute confidence intervals
  reduced_names <- names(exposures_list)
  complete_df_list <- 
    lapply(reduced_names, function(current_condition){
      current_exposures_df <- exposures_list[[current_condition]]
      current_sigs <- rownames(current_exposures_df)
      current_signatures_df <- AlexCosmicArtif_sig_df[, current_sigs, drop = FALSE]
      suppressWarnings(variateExp(
        in_catalogue_df = corrected_catalogue_df,
        in_sig_df = current_signatures_df,
        in_exposures_df = current_exposures_df,
        in_sigLevel = 0.025, in_delta = 0.4))
    })
  complete_df_list <- lapply(complete_df_list, function(current_complete_df){
    current_complete_df$norm_exposure <- current_complete_df$exposure / 
      current_complete_df[which(current_complete_df$sig == "total"), "exposure"]
    current_complete_df$norm_lower <- current_complete_df$norm_exposure * 
      current_complete_df$relLower
    current_complete_df$norm_upper <- current_complete_df$norm_exposure * 
      current_complete_df$relUpper
    return(current_complete_df)
  })
  complete_df <- do.call(rbind, complete_df_list)

  sign_list <- list(CosmicValid_gen_LCDlist, CosmicValid_abs_LCDlist, CosmicValid_norm_LCDlist)

  cutoffPerc <- cutoff * 100

  out <- data.frame(
    Signature = sign_list[[which(!is.na(annotation_df$mostInformative))]]$out_sig_ind_df[, 1],
    Process = sign_list[[which(!is.na(annotation_df$mostInformative))]]$out_sig_ind_df[, 4],
    Percentage = sign_list[[which(!is.na(annotation_df$mostInformative))]]$norm_exposures*100,
    CI_lb = complete_df_list[[which(!is.na(annotation_df$mostInformative))]]$norm_lower[1:length(sign_list[[which(!is.na(annotation_df$mostInformative))]]$out_sig_ind_df[, 1])]*100,
    CI_ub = complete_df_list[[which(!is.na(annotation_df$mostInformative))]]$norm_upper[1:length(sign_list[[which(!is.na(annotation_df$mostInformative))]]$out_sig_ind_df[, 1])]*100)

  output <- list("Mutation_Signature_Catalog" = mutation_catalogue_df,
                 "Signatures_Identified" = sign_list[[which(!is.na(annotation_df$mostInformative))]]$out_sig_ind_df[, c(1:2, 4)],
                 "Normalized_Exposures" = my_normExposures_df[, rownames(annotation_df)[which(!is.na(annotation_df$mostInformative))], drop = FALSE],
                 "Confidence_Intervals" = complete_df,
                 Summary = out)
  colnames(output$Mutation_Signature_Catalog) <- c(sample)
  colnames(output$Normalized_Exposures) <- c(sample)
  colnames(output$Summary)[3] <- c(sample)

  # Is 0 in the CI in any?
  CI <- length(which(output$Confidence_Intervals$lower[which(output$Confidence_Intervals$sig == "AC3")] < 0))
  mutsig2cbioportal_wes(signatures = out, all_signatures = AlexCosmicValid_sigInd_df, id = sample_name, outfile_cbioportal = outfile_cbioportal)

  return(list(output = output, CI = CI))
}
