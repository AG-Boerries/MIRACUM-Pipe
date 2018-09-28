mutation_signature_analysis <- function(vcf_file = NULL, cutoff = 0.01,
                                      sample_name = NULL, only_coding = FALSE,
                                      path_data){
  #' Mutation Signature Analysis
  #'
  #' @description Mutation Signature Analysis adopted from YAPSA package
  #' @description by Daniel Huebschmann et al.
  #' 
  #' @param vcf_file string. Filename of Input
  #' @param cutoff numerical. Signature is kept when over cutoff (default: 0.01)
  #' @param sample_name string. Name of the sample (default: NULL)
  #' @param only_coding logical. Use only mutations in coding regions
  #' @param (default: F)
  #' @param path_data string. Path to data
  #'
  #' @return returns list of
  #' @return CosmicValid_cutoffGen_LCDlist dataframe. List of cutoff genes
  #' @return mutationCataloge dataframe. Cataloge of mutations
  #' 
  #' @details For details see the Documentation of the YAPSA package.
  require("YAPSA")
  require("SomaticSignatures")
  require("VariantAnnotation")
  require("BSgenome.Hsapiens.UCSC.hg19")
  require("TxDb.Hsapiens.UCSC.hg19.knownGene")
  require("openxlsx")
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  ## Loading Siganture Info
  Alex_signatures_path <- paste(path_data, "signatures.txt", sep = "/")
  AlexInitialArtif_sig_df <- read.csv(Alex_signatures_path, header=TRUE,
                                      sep="\t")
  Alex_rownames <- paste(AlexInitialArtif_sig_df[, 1],
                         AlexInitialArtif_sig_df[, 2], sep=" ")
  select_ind <- grep("Signature", names(AlexInitialArtif_sig_df))
  AlexInitialArtif_sig_df <- AlexInitialArtif_sig_df[, select_ind]
  number_of_Alex_sigs <- dim(AlexInitialArtif_sig_df)[2]
  names(AlexInitialArtif_sig_df) <- gsub("Signature\\.","A",
                                         names(AlexInitialArtif_sig_df))
  rownames(AlexInitialArtif_sig_df) <- Alex_rownames
  AlexInitialValid_sig_df <- AlexInitialArtif_sig_df[,grep("^A[0-9]+",
                             names(AlexInitialArtif_sig_df))]
  number_of_Alex_validated_sigs <- dim(AlexInitialValid_sig_df)[2]
  Alex_COSMIC_signatures_path <- paste(path_data,
                                       "signatures_probabilities.txt",
                                       sep = "/")
  AlexCosmicValid_sig_df <- read.csv(Alex_COSMIC_signatures_path, header=TRUE,
                                     sep="\t")
  Alex_COSMIC_rownames <- paste(AlexCosmicValid_sig_df[, 1],
                                AlexCosmicValid_sig_df[, 2], sep = " ")
  COSMIC_select_ind <- grep("Signature", names(AlexCosmicValid_sig_df))
  AlexCosmicValid_sig_df <- AlexCosmicValid_sig_df[, COSMIC_select_ind]
  number_of_Alex_COSMIC_sigs <- dim(AlexCosmicValid_sig_df)[2]
  names(AlexCosmicValid_sig_df) <- gsub("Signature\\.","AC",
                                        names(AlexCosmicValid_sig_df))
  rownames(AlexCosmicValid_sig_df) <- Alex_COSMIC_rownames
  COSMIC_order_ind <- match(Alex_rownames,Alex_COSMIC_rownames)
  AlexCosmicValid_sig_df <- AlexCosmicValid_sig_df[COSMIC_order_ind, ]
  
  signature_colour_vector <- c("darkgreen", "green", "pink", "goldenrod",
                               "lightblue", "blue", "orangered", "yellow",
                               "orange", "brown", "purple", "red",
                               "darkblue", "magenta", "maroon", "yellowgreen",
                               "violet", "lightgreen", "sienna4", "deeppink",
                               "darkorchid", "seagreen", "grey10", "grey30",
                               "grey50", "grey70", "grey90")
  bio_process_vector <- c("spontaneous deamination", "spontaneous deamination",
                          "APOBEC", "BRCA1_2", "Smoking", "unknown",
                          "defect DNA MMR", "UV light exposure", "unknown",
                          "IG hypermutation", "POL E mutations",
                          "temozolomide", "unknown", "APOBEC", "unknown",
                          "unknown", "unknown", "unknown", "unknown",
                          "unknown", "unknown", "unknown",
                          "nonvalidated", "nonvalidated", "nonvalidated",
                          "nonvalidated", "nonvalidated")
  AlexInitialArtif_sigInd_df <- data.frame(
                                sig = colnames(AlexInitialArtif_sig_df))
  AlexInitialArtif_sigInd_df$index <- seq_len(dim(AlexInitialArtif_sigInd_df)[1])
  AlexInitialArtif_sigInd_df$colour <- signature_colour_vector
  AlexInitialArtif_sigInd_df$process <- bio_process_vector
  
  COSMIC_signature_colour_vector <- c("green", "pink", "goldenrod",
                                      "lightblue", "blue", "orangered",
                                      "yellow", "orange", "brown", "purple",
                                      "red", "darkblue", "magenta", "maroon",
                                      "yellowgreen", "violet", "lightgreen",
                                      "sienna4", "deeppink", "darkorchid",
                                      "seagreen", "grey", "darkgrey",
                                      "black", "yellow4", "coral2", "chocolate2",
                                      "navyblue", "plum", "springgreen")
  COSMIC_bio_process_vector <- c("spontaneous deamination", "APOBEC",
                                 "defect DNA DSB repair hom. recomb.",
                                 "tobacco mutatgens, benzo(a)pyrene",
                                 "unknown",
                                 "defect DNA MMR, found in MSI tumors",
                                 "UV light exposure", "unknown",
                                 "POL eta and SHM", "altered POL E",
                                 "alkylating agents, temozolomide",
                                 "unknown", "APOBEC", "unknown",
                                 "defect DNA MMR", "unknown", "unknown",
                                 "unknown", "unknown",
                                 "associated w. small indels at repeats",
                                 "unknown", "aristocholic acid","unknown",
                                 "aflatoxin", "unknown", "defect DNA MMR",
                                 "unknown", "unknown", "tobacco chewing",
                                 "unknown")
  AlexCosmicValid_sigInd_df <- data.frame(
    sig = colnames(AlexCosmicValid_sig_df))
  AlexCosmicValid_sigInd_df$index <- seq_len(dim(AlexCosmicValid_sigInd_df)[1])
  AlexCosmicValid_sigInd_df$colour <- COSMIC_signature_colour_vector
  AlexCosmicValid_sigInd_df$process <- COSMIC_bio_process_vector
  
  current_sig_df <- AlexCosmicValid_sig_df
  current_sigInd_df <- AlexCosmicValid_sigInd_df
  
  chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22", "chrX", "chrY")
  ########
  # MAIN #
  ########
  
  if(is.null(sample_name)){
    sam_na <- unlist(strsplit(vcf_file, "/", fixed = TRUE))
    sample_name <- unlist(strsplit(sam_na[length(sam_na)], ".", fixed = T))[1]
  }
  
  ## load vcf Files
  sample.vcf <- readVcf(vcf_file,"hg19")
  
  ## create data.frame from mutations
  if(only_coding == TRUE){
    mutations.coding <- predictCoding(sample.vcf, txdb, seqSource = Hsapiens)
    mutations <- data.frame(CHROM = as.character(seqnames(mutations.coding)),
                            POS = start(mutations.coding),
                            REF = as.character(mutations.coding$REF),
                            ALT = as.character(unlist(mutations.coding$ALT)),
                            Type = as.character(mutations.coding$CONSEQUENCE),
                            FILT = as.character(mutations.coding$FILTER))
    mutations <- mutations[!mutations$Type == "synonymous", ]
    mutations$PID <- sample_name
    mutations$SUBGROUP <- sample_name
  } else {  
    mutations <- data.frame(CHROM = as.character(seqnames(sample.vcf)),
                            POS = start(sample.vcf),
                            REF = as.character(ref(sample.vcf)),
                            ALT = as.character(unlist(alt(sample.vcf))),
                            FILT = as.character(filt(sample.vcf)))
    mutations$PID <- sample_name
    mutations$SUBGROUP <- sample_name
  }
  
  idx <- which (mutations$CHROM %in% chromosomes)
  mutations <- mutations[idx, ]
  idx <- which (mutations$FILT == "PASS")
  mutations <- mutations[idx, ]
  
  
  df <- mutations[(mutations$REF %in% DNA_BASES
                   & mutations$ALT %in% DNA_BASES), ]

  MutCat_List <- create_mutation_catalogue_from_df(this_df = df,
                 this_refGenome = BSgenome.Hsapiens.UCSC.hg19,
                 this_start.field = "POS", this_seqnames.field = "CHROM",
                 this_end.field = "POS", this_wordLength = 3,
                 this_PID.field = "PID", this_subgroup.field = "subgroup")
  mutCat_df <- as.data.frame(MutCat_List$matrix)
  
  ## LCD with general cutoff
  my_cutoff <- cutoff #signature is kept if it represents at least xx% of all snv's
  general_cutoff_vector <- rep(my_cutoff, dim(current_sig_df)[2])
  CosmicValid_cutoffGen_LCDlist <- LCD_complex_cutoff(
    in_mutation_catalogue_df = mutCat_df, in_signatures_df = current_sig_df,
    in_cutoff_vector = general_cutoff_vector,
    in_sig_ind_df = current_sigInd_df)

  cutoffPerc <- cutoff * 100
  
  out <- data.frame(
    Signature = CosmicValid_cutoffGen_LCDlist$out_sig_ind_df$sig,
    Process = CosmicValid_cutoffGen_LCDlist$out_sig_ind_df$process,
    Percentage = CosmicValid_cutoffGen_LCDlist$norm_exposures[,1])
  
  output <- list("Mutation_Signature_Catalog" = mutCat_df,
        "Signatures_Identified" = CosmicValid_cutoffGen_LCDlist$out_sig_ind_df,
        "Normalized_Exposures" = CosmicValid_cutoffGen_LCDlist$norm_exposures,
        Summary = out)
  write.xlsx(output, paste0(sample_name, "_Mutation_Signature_cutoff_",
             cutoffPerc, "Percent.xlsx"), rowNames = T, firstRow = T,
             headerStyle = createStyle(textDecoration = 'bold'))

  return(list(CosmicValid_cutoffGen_LCDlist = CosmicValid_cutoffGen_LCDlist,
              mutationCataloge = mutCat_df))
}