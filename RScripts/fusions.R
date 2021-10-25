#### Fusions
fusions_ana <- function(fus_file, path_data, fusion_genes = NULL) {
  fusions <- read.delim(fus_file, header = TRUE)
  if(dim(fusions)[1] == 0) {
    fusion_table <- NULL
    flt_fus_ov <- NULL
    fus_38 <- NULL
  } else {
    # Prepare data for the fusion plots
    fus_file_38 <- paste0(
      substr(
        fus_file,
        start = 1,
        stop = nchar(fus_file) - 8),
        "txt"
      )
    path_agfus <- substr(fus_file, start = 1, stop = nchar(fus_file) - 43)
    ag_fus <- paste0(
      "agfusion batch -f ",
      fus_file_38,
      " -a fusioncatcher -o ",
      path_agfus,
      "/plots -db ",
      path_data,
      "/agfusion.homo_sapiens.95.db --recolor \"Pkinase_Tyr;red\" --recolor \"Pkinase;red\" "
    )
    #system(ag_fus)

    fusions_38 <- read.delim(fus_file_38, header = TRUE)
    fus_38 <- fusions_38[, c(1:2, 9:10)]
    colnames(fus_38) <- c("Gen1", "Gen2", "Bruch1", "Bruch2")

    sep <- strsplit(x = as.character(fus_38$Bruch1), split = ":")
    fus_38$Chr1 <- unlist(lapply(sep, function(x){return( x[[1]])}))
    fus_38$Pos1 <- unlist(lapply(sep, function(x){return( x[[2]])}))
    sep <- strsplit(x = as.character(fus_38$Bruch2), split = ":")
    fus_38$Chr2 <- unlist(lapply(sep, function(x){return( x[[1]])}))
    fus_38$Pos2 <- unlist(lapply(sep, function(x){return( x[[2]])}))

    fus_38$paths <- paste0(path_agfus, "/plots/", fus_38$Gen1,
                            "-", fus_38$Pos1, "_", fus_38$Gen2, "-",
                            fus_38$Pos2)
    fus_38$file <- NA
    for (i in 1:dim(fus_38)[1]){
      files <- list.files(path = fus_38$paths[i], pattern = ".png$")
      if( length(files) > 0){
        fus_38$file[i] <- files[-grep(pattern = "exon", x = files)]
      }
    }

    fus_38 <- fus_38[, c(1:2, 9:10)]

    ftd_fus <- fusions[, c(1:3, 9:10, 15:16, 5:6)]
    colnames(ftd_fus) <- c("Gen1", "Gen2", "Beschreibung", "Bruch1", "Bruch2",
                         "Sequenz", "Effekt", "Supporting", "Split Reads")
    # Check for detected Genes only
    if (!is.null(fusion_genes)) {
      fus_genes <- read.delim(
        file = fusion_genes,
        header = F
      )$V1
      id1 <- which(ftd_fus$Gen1 %in% fus_genes)
      id2 <- which(ftd_fus$Gen2 %in% fus_genes)
      id <- union(id1, id2)
      ftd_fus <- ftd_fus[id, ]
    }

    # Check for False Positives
    if (dim(ftd_fus)[1] == 0) {
      cat("No Fusions detected.")
      fusion_table <- data.frame()
      flt_fus_ov <- NULL
      fus_38 <- NULL
    } else {
      beschr <- strsplit(x = as.character(ftd_fus$Beschreibung), split = ",",
                       fixed = TRUE)
      v_high <- c(
        "1000genome", "banned", "bodymap2", "cacg", "conjoing", "cortex",
        "cta", "ctb", "ctc", "ctd", "distance1000bp", "distance100kbp",
        "distance10kbp", "duplicates", "ensemble_fully_overlapping",
        "ensemble_partially_overlapping", "ensemble_same_strand_overlapping",
        "fragments", "gtex", "hpa", "mt", "no_protein", "pair_pseudo_genes",
        "paralogs", "readthrough", "refseq_fully_overlapping",
        "refseq_same_strand_overlapping", "rp11", "rp", "rrna", "similar_reads",
        "similar_symbols", "ucsc_fully_overlapping",
        "ucsc_same_strand_overlapping", "ucsc_partially_overlapping",
        "short_repeats", "long_repeats", "short_distance"
      )
      id_fp <- which(
        unlist(
          lapply(beschr, function(x){any(unlist(x) %in% v_high)})
        )
      )
      ftd_fus$FalsePositive <- "Nein"
      ftd_fus$FalsePositive[id_fp] <- "Ja"
      # Check for Known Fusions
      idx <- grep(pattern = "known", x = ftd_fus$Beschreibung, value = FALSE)
      ftd_fus$Bekannt <- "Nein"
      ftd_fus$Bekannt[idx] <- "Ja"
      gen12 <- paste0(ftd_fus$Gen1, "--", ftd_fus$Gen2)
      gen21 <- paste0(ftd_fus$Gen2, "--", ftd_fus$Gen1)
      idr <- which(gen12 %in% gen21)
      ftd_fus$Reziprok <- "Nein"
      ftd_fus$Reziprok[idr] <- "Ja"
      ftd_fus$Readthrough <- "Nein"
      ftd_fus$Readthrough[grep(
        pattern = "readthrough", x = as.character(ftd_fus$Beschreibung))]
      ftd_fus$Cosmic <- "Nein"
      ftd_fus$Cosmic[grep(pattern = "cosmic", ftd_fus$Beschreibung)] <- "Ja"

      flt_fus_ov <- data.frame(Fusionen = gen12[which(!duplicated(gen12))])
      flt_fus_ov$Bekannt <- ftd_fus$Bekannt[which(!duplicated(gen12))]
      flt_fus_ov$Reziprok <- ftd_fus$Reziprok[which(!duplicated(gen12))]
      flt_fus_ov$Readthrough <- ftd_fus$Readthrough[which(!duplicated(gen12))]
      flt_fus_ov$FalsePositive <- ftd_fus$FalsePositive[which(!duplicated(gen12))]
      flt_fus_ov$Cosmic <- ftd_fus$Cosmic[which(!duplicated(gen12))]
      if(length(which(flt_fus_ov$FalsePositive == "Ja")) > 0){
        fps <- flt_fus_ov$Fusionen[which(flt_fus_ov$FalsePositive == "Ja")]
        flt_fus_ov <- flt_fus_ov[-which(flt_fus_ov$FalsePositive == "Ja"), -c(5)]
      } else {
        flt_fus_ov <- flt_fus_ov[, -c(5)]
      }

      # keep only in-frame and out-of-frame fusions
      keep <- which(ftd_fus$Effekt %in% c("in-frame", "out-of-frame"))
      ftd_fus <- ftd_fus[keep,]
      fusions_ids <- paste(ftd_fus$Gen1, ftd_fus$Gen2, sep =  "--")
      flt_fus_ov <- flt_fus_ov[which(flt_fus_ov$Fusionen %in% fusions_ids),]

      # Additionally set a cutoff on split reads
      keep <- ftd_fus[,9] > 3
      ftd_fus <- ftd_fus[keep,]
      fusions_ids <- paste(ftd_fus$Gen1, ftd_fus$Gen2, sep =  "--")
      flt_fus_ov <- flt_fus_ov[which(flt_fus_ov$Fusionen %in% fusions_ids),]

      if (dim(ftd_fus)[1] == 0) {
        cat("No in-frame/out-of-frame Fusions detected.")
        fusion_table <- data.frame()
        flt_fus_ov <- NULL
        fus_38 <- NULL
      } else {
        fusion_table <- ftd_fus[
          which(
            paste(
              ftd_fus$Gen1, ftd_fus$Gen2 ,sep =  "--"
            )
            %in% flt_fus_ov$Fusionen
          ),
          c(1, 4, 2, 5:7, 8:9)
        ]
        fus_38 <- fus_38[
          which(
            paste(
              fus_38$Gen1, fus_38$Gen2, sep = "--"
            )
            %in% flt_fus_ov$Fusionen),
        ]
      }
    }
  }

  return(list(Fusion_OV = flt_fus_ov, Table = fusion_table, Plots = fus_38))
}

fusions2cbioportal <- function(fusions, sample, outfile_fusions_cbioportal){
  #' export fusions to cbioportal
  #' @param fusions. fusions_ana output (list)

  # cbioportal help
  # Hugo_Symbol: A HUGO gene symbol.
  # Entrez_Gene_Id: A Entrez Gene identifier.
  # Center: The sequencing center.
  # Tumor_Sample_Barcode: This is the sample ID.
  # Fusion: A description of the fusion, e.g., "TMPRSS2-ERG fusion".
  # DNA_support: Fusion detected from DNA sequence data, "yes" or "no".
  # RNA_support: Fusion detected from RNA sequence data, "yes" or "no".
  # Method: Fusion detected algorithm/tool.
  # Frame: "in-frame" or "frameshift".
  # Fusion_Status (OPTIONAL): An assessment of the mutation type (i.e., "SOMATIC", "GERMLINE", "UNKNOWN", or empty)
  ## example
  #Hugo_Symbol<TAB>Entrez_Gene_Id<TAB>Center<TAB>Tumor_Sample_Barcode<TAB>Fusion<TAB>DNA_support<TAB>RNA_support<TAB>Method<TAB>Frame>
  #ALK<TAB>238<TAB>center.edu<TAB>SAMPLE_ID_1<TAB>Fusion<TAB>unknown<TAB>yes<TAB>unknown<TAB>in-frame
  #ALK<TAB>238<TAB>center.edu<TAB>SAMPLE_ID_2<TAB>Fusion<TAB>unknown<TAB>yes<TAB>unknown<TAB>in-frame
  #RET<TAB>5979<TAB>center.edu<TAB>SAMPLE_ID_3<TAB>Fusion<TAB>unknown<TAB>yes<TAB>unknown<TAB>in-frame

  fusion_table <- fusions$Table[, c(1:4, 6)]
  fusion_table$Effekt <- as.character(fusion_table$Effekt)
  fusion_table$Effekt[which(fusion_table$Effekt == "out-of-frame")] <- "frameshift"

  fusions.cbioportal <- data.frame(
    "Hugo_Symbol" = c(
      as.character(fusion_table$Gen1),
      as.character(fusion_table$Gen2)
    ),
    "Entrez_Gene_Id" = c(
      unlist(
        lapply(
          mget(as.character(fusion_table$Gen1), org.Hs.egSYMBOL2EG),
          function(x) x[1]
        )
      ),
      unlist(
        lapply(
          mget(as.character(fusion_table$Gen2), org.Hs.egSYMBOL2EG),
          function(x) x[1]
        )
      )
    ),
    "Center" = "Freiburg",
    "Tumor_Sample_Barcode" = paste(as.character(sample),"TD",sep = "_"),
    "Fusion" = c(
      paste0(fusion_table$Gen1, "-", fusion_table$Gen2),
      paste0(fusion_table$Gen1, "-", fusion_table$Gen2)
    ),
    "DNA_support" = "unkown",
    "RNA_support" = "yes",
    "Method" = "fusioncatcher",
    "Frame" = c(fusion_table$Effekt, fusion_table$Effekt),
    "Fusion_Status" = "UNKNOWN"
  )

  write.table(
    x = fusions.cbioportal,
    file = outfile_fusions_cbioportal,
    append = F,
    quote = F,
    sep = "\t",
    col.names = T,
    row.names = F)
}
