###################
# Helper Function for Mutation Analysis #
find_indel <- function(list){
  #' Indel Finder
  #'
  #' @description Find Indels in list
  #'
  #' @param list dataframe. List of mutations
  #'
  #' @return id indexvector. List of indices
  #'
  #' @details Given a list of mutations find the indels by checking the
  #' @details reference and alternative bases. 
  id.ref <- grep("-", list$Ref)
  id.alt <- grep("-", list$Alt)
  id <- unique( c(id.ref, id.alt) )
  return(id)
}

div <- function(x_s, x_l, no_loh){
  #' Mutation separation
  #'
  #' @description Separate mutions
  #'
  #' @param x_s dataframe. List of somatic mutations
  #' @param x_l dataframe. List of LoH mutations
  #' @param no_loh logical. Logical describing existence of LoH mutations
  #'
  #' @return list of
  #' @return x_s_snp dataframe. List of somatic SNVs
  #' @return x_s_indel dataframe. List of somatic InDels
  #' @return x_l_snp dataframe. List of LoH SNVs
  #' @return x_l_indel dataframe. List of LoH InDels
  #' @return no_loh logical. Describing existence of LoH mutations
  #' @return no_indel_somatic logical. Describing existence of somatic Indels
  #' @return no_snp logical. Describing existence of SNVs
  #' @return no_indel_loh logical. Describing existence of LoH Indels
  #'
  #' @details Build separate dataframes for different type of mutations.
  #' @details Split somatic and LoH mutations in SNVs and InDels. Return also
  #' @details some logicals for existence of mutations at all.
  no_indel_somatic <- FALSE
  indel_s <- find_indel(x_s)
  if (length(indel_s) > 0){
    x_s_snp <- x_s[-indel_s, ]
    x_s_indel <- x_s[indel_s, ]
  } else {
    x_s_snp <- x_s
    x_s_indel <- data.frame()
    cat("No Indels in Somatic!\n")
    no_indel_somatic <- TRUE
  }
  if (dim(x_s_snp)[1] > 0){
    no_snp <- FALSE
  } else {
    no_snp <- TRUE
    }

  no_indel_loh <- FALSE
  if (!no_loh){
    indel_l <- find_indel(x_l)
    if (length(indel_l) == 0){
      no_indel_loh <- TRUE
      x_l_snp <- x_l
      x_l_indel <- data.frame()
      cat("No Indels in LOH!\n")
      } else {
      x_l_snp <- x_l[-indel_l, ]
      x_l_indel <- x_l[indel_l, ]
      }
  } else {
    x_l_snp <- data.frame()
    x_l_indel <- data.frame()
  }
  return(list(x_s_snp = x_s_snp, x_s_indel = x_s_indel, x_l_snp = x_l_snp,
              x_l_indel = x_l_indel, no_loh = no_loh,
              no_indel_somatic = no_indel_somatic, no_snp = no_snp,
              no_indel_loh = no_indel_loh))
}

mut_tab <- function(x_s_snp, x_s_indel, x_l_snp, x_l_indel){
  #' Mutation Table
  #'
  #' @description Build Mutation Table
  #'
  #' @param x_s_snp dataframe. List of somatic SNVs
  #' @param x_s_indel dataframe. List of somatic InDels
  #' @param x_l_snp dataframe. List of LoH SNVs
  #' @param x_l_indel dataframe. List of LoH InDels
  #'
  #' @return muta_tab dataframe. Summary table of mutations
  #'
  #' @details A summary table is build that shows the number of mutations
  #' @details separating SNVs/InDels, somatic/LoH and Zygosity.
  #' @details The table is also written in MutationTable.txt.
  muta_tab <- matrix(data = "-", nrow = 6, ncol = 6)
  colnames(muta_tab) <- c("Mutationtype", "Number of exonic", "Zygosity",
                               "Tumorsuppressor", "Oncogene", "Hotspot")
  muta_tab[, 1] <- c("somatic SNV", "somatic SNV", "LoH SNV",
                         "somatic InDel", "somatic InDel", "LoH InDel")

  muta_tab[1, 2] <- sum(x_s_snp$Zygosity == "hom")
  muta_tab[2, 2] <- sum(x_s_snp$Zygosity == "het")
  muta_tab[3, 2] <- dim(x_l_snp)[1]
  muta_tab[4, 2] <- sum(x_s_indel$Zygosity == "hom")
  muta_tab[5, 2] <- sum(x_s_indel$Zygosity == "het")
  muta_tab[6, 2] <- dim(x_l_indel)[1]

  muta_tab[c(1:2, 4:5), 3] <- rep(c("homozygous", "heterozygous"), times = 2)

  muta_tab[1, 4] <- sum(x_s_snp$is_tumorsuppressor == 1
                        & x_s_snp$Zygosity == "hom")
  muta_tab[2, 4] <- sum(x_s_snp$is_tumorsuppressor == 1
                        & x_s_snp$Zygosity == "het")
  muta_tab[3, 4] <- sum(x_l_snp$is_tumorsuppressor == 1)
  muta_tab[4, 4] <- sum(x_s_indel$is_tumorsuppressor == 1
                        & x_s_indel$Zygosity == "hom")
  muta_tab[5, 4] <- sum(x_s_indel$is_tumorsuppressor == 1
                        & x_s_indel$Zygosity == "het")
  muta_tab[6, 4] <- sum(x_l_indel$is_tumorsuppressor == 1)

  muta_tab[1, 5] <- sum(x_s_snp$is_oncogene == 1 & x_s_snp$Zygosity == "hom")
  muta_tab[2, 5] <- sum(x_s_snp$is_oncogene == 1 & x_s_snp$Zygosity == "het")
  muta_tab[3, 5] <- sum(x_l_snp$is_oncogene == 1)
  muta_tab[4, 5] <- sum(x_s_indel$is_oncogene == 1
                        & x_s_indel$Zygosity == "hom")
  muta_tab[5, 5] <- sum(x_s_indel$is_oncogene == 1
                        & x_s_indel$Zygosity == "het")
  muta_tab[6, 5] <- sum(x_l_indel$is_oncogene == 1)

  muta_tab[1, 6] <- sum(x_s_snp$is_hotspot != 0 & x_s_snp$Zygosity == "hom")
  muta_tab[2, 6] <- sum(x_s_snp$is_hotspot != 0 & x_s_snp$Zygosity == "het")
  muta_tab[3, 6] <- sum(x_l_snp$is_hotspot != 0 )
  muta_tab[4, 6] <- sum(x_s_indel$is_hotspot != 0
                        & x_s_indel$Zygosity == "hom")
  muta_tab[5, 6] <- sum(x_s_indel$is_hotspot != 0
                        & x_s_indel$Zygosity == "het")
  muta_tab[6, 6] <- sum(x_l_indel$is_hotspot != 0 )

  write.table(x = muta_tab, file = "MutationTable.txt", quote = F,
              sep = "\t", row.names = F, col.names = T)
  print(muta_tab)
 return(muta_tab = muta_tab)
}

mut_stats <- function(x_s, x_l = NULL, tumbu) {
  #' Mutation Statistics
  #'
  #' @description Print Number of Somatic Mutations
  #'
  #' @param x_s dataframe. List of somatic mutations
  #' @param x_l dataframe. List of LoH mutations
  #' @param tumbu numerical. Tumor mutational burden
  #'
  #' @return list of
  #' @return tot_mut numerical. Total number of mutations
  #' @return som_mut numerical. Number of all somatic mutations
  #' @return loh_mut numerical. Number of all LoH mutations
  #'
  #' @details Statistical number are calculated and printed. Furthermore
  #' @details the results are printed in mutationsStats.txt.
  sink(file = "mutationsStats.txt", append = T, split = T)
  print(paste(dim(x_s)[1], "somatic mutations", sep = " ")) 
  if (is.null(x_l)){
    print("0 LoH")
    dim_x_l <- 0
  } else {
    print(paste(dim(x_l)[1], "LoH", sep = " "))
    dim_x_l <- dim(x_l)[1]
    }
  print(paste("total number of mutations:", dim(x_s)[1] + dim_x_l,
              sep = " "))
  print(paste("mutational burden:", tumbu, sep = " "))
  sink()
  totalmutationnumber <- dim(x_s)[1] + dim_x_l
  somaticmutations <- dim(x_s)[1]
  lohmutations <- dim_x_l
  return(list(tot_mut <- totalmutationnumber, som_mut <- somaticmutations,
              loh_mut <- lohmutations))
}

tables <- function(x_s, x_l = NULL){
  #' Create Tables
  #'
  #' @description Write Tables for Tumorsuppressors/Oncogenes, all somatic
  #' @description mutations, LoH mutations
  #'
  #' @param x_s dataframe. List of somatic mutations
  #' @param x_l dataframe. List of LoH mutations
  #'
  #' @return list of
  #' @return ts_og_table dataframe. List of mutations in tumorsuppressors and
  #' @return oncogenes
  #' @return sm_table dataframe. List if somatic mutations
  #' @return lm_table dataframe. List of LoH mutations
  #'
  #' @details The Tables TumorSuppressor-OncogeneTable.xlsx,
  #' @details somaticMutations.xlsx and lohMutations.xlsx are generated and
  #' @details stored as well as returned.
  ts_og_table <- x_s[x_s$is_tumorsuppressor == 1 |
                       x_s$is_oncogene == 1,
                     c("Gene.refGene", "GeneName",
                       "ExonicFunc.refGene", "AAChange.SnpEff",
                       "Variant_Allele_Frequency",
                       "Zygosity", "Variant_Reads", "is_tumorsuppressor",
                       "is_oncogene", "is_hotspot", "target",
                       "gnomAD_exome_NFE", "CADD13_PHRED", "condel.label",
                       "CLINSIG.SnpEff", "cosmic84_coding"), drop = FALSE]
  write.xlsx(ts_og_table, file = "TumorSuppressor-OncogeneTable.xlsx",
             quote = F, sep = "\t", row.names = F, col.names = T)

  sm_table <- x_s[, c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                     "AAChange.SnpEff", "Variant_Allele_Frequency", "Zygosity",
                     "Variant_Reads", "is_tumorsuppressor", "is_oncogene",
                     "is_hotspot", "target", "gnomAD_exome_NFE",
                     "CADD13_PHRED", "condel.label", "CLINSIG.SnpEff",
                     "cosmic84_coding"), drop = FALSE]
  write.xlsx(sm_table, file = "somaticMutations.xlsx",
             quote = F, sep = "\t", row.names = F, col.names = T)
  if (!is.null(x_l)){
    lm_table <- x_l[, c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                       "AAChange.SnpEff", "VAF_Normal", "VAF_Tumor",
                       "Count_Normal", "Count_Tumor", "is_tumorsuppressor",
                       "is_oncogene", "is_hotspot", "target",
                       "gnomAD_exome_NFE", "CADD13_PHRED", "condel.label",
                       "CLINSIG.SnpEff", "cosmic84_coding"), drop = FALSE]
    write.xlsx(lm_table, file = "lohMutations.xlsx",
               quote = F, sep = "\t", row.names = F, col.names = T)
  }else{
    lm_table <- data.frame()
    
  }
  return(list(ts_og_table = ts_og_table, sm_table = sm_table,
              lm_table = lm_table))
}

get_mapping_matrix <- function(annovar_table, row_index)
{
  if (length(row_index)==0) return(NULL)
  col_index <- match(c("Chr", "Start", "Gene.refGene"), colnames(annovar_table))
  new_table <- annovar_table[row_index, col_index]
  new_table[] <- lapply(new_table, as.character)
  new_table$Start <- as.numeric(new_table$Start)
  return(new_table)
}

add_default_value <- function(mapping_matrix)
{
  if (is.null(mapping_matrix)) return(mapping_matrix)
  new_matrix <- cbind(mapping_matrix, rep(1, nrow(mapping_matrix)))
  colnames(new_matrix)[ncol(new_matrix)] <- "Value"
  return(new_matrix)	
}

rename_chr <- function(mapping_matrix)
{
  new_matrix <- mapping_matrix
  new_matrix[,1] <- gsub("chr", "", new_matrix[,1])
  return(new_matrix)	
}

duplicate_first_raw <- function(mapping_matrix)
{
  if (is.null(mapping_matrix)) return(NULL)
  new_matrix <- mapping_matrix
  new_matrix <- rbind(new_matrix[1,], new_matrix)
  new_matrix[1,c(4:ncol(new_matrix))] <- 0
  return(new_matrix)
}

circos_colors <- function(x_s_snp = NULL, x_s_indel = NULL, x_l_snp = NULL,
                          x_l_indel = NULL, no_loh,
                          no_indel_somatic, no_snp, no_indel_loh){
  #' Circos Colors
  #'
  #' @description Prepare List and colors for Circosplot
  #'
  #' @param x_s_snp dataframe. List of somatic SNVs
  #' @param x_s_indel dataframe. List of somatic InDels
  #' @param x_l_snp dataframe. List of LoH SNVs
  #' @param x_l_indel dataframe. List of LoH InDels
  #' @param no_loh logical. Describing existence of LoH mutations
  #' @param no_indel_somatic logical. Describing existence of somatic Indels
  #' @param no_snp logical. Describing existence of SNVs
  #' @param no_indel_loh logical. Describing existence of LoH Indels  #'
  #' 
  #' @return list of
  #' @return map_mat matrix. Mutationmatrix to be plotted
  #' @return circoscolors vector of strings. Colors for Circosplot
  #'
  #' @details For the Circosplot it is important to know the number of
  #' @details different types of mutation (Somatic SNV, Somatic InDel, LoH Snp,
  #' @details LoH InDel). So the corresponding colors can be chosen.
  #' @details Additionaly there is a matrix needed that includes all the
  #' @details mutations.
  if (no_indel_somatic == FALSE & no_snp == FALSE & no_loh == FALSE
     & no_indel_loh == FALSE){
    list1 <- list(x_s_snp, x_s_indel, x_l_snp, x_l_indel)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_s_indel), 1:nrow(x_l_snp),
                 1:nrow(x_l_indel))
    circoscolors <- c("#FF0000CC", "#008000CC", "#00FFFFCC", "#8000FFCC")
  }
  if (no_indel_somatic == FALSE & no_snp == FALSE & no_loh == FALSE
     & no_indel_loh == TRUE){
    list1 <- list(x_s_snp, x_s_indel, x_l_snp)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_s_indel), 1:nrow(x_l_snp))
    circoscolors <- c("#FF0000CC", "#008000CC", "#00FFFFCC")
  }
  if (no_indel_somatic == FALSE & no_snp == FALSE & no_loh == TRUE){
    list1 <- list(x_s_snp, x_s_indel)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_s_indel))
    circoscolors <- c("#FF0000CC", "#008000CC")
  }
  if (no_indel_somatic == FALSE & no_snp == TRUE & no_loh == FALSE
     & no_indel_loh == TRUE){
    list1 <- list(x_s_indel, x_l_snp)
    idxs <- list(1:nrow(x_s_indel), 1:nrow(x_l_snp))
    circoscolors <- c("#008000CC", "#00FFFFCC")
  }
  if (no_indel_somatic == FALSE & no_snp == TRUE & no_loh == TRUE){
    list1 <- list(x_s_indel)
    idxs <- list(1:nrow(x_s_indel))
    circoscolors <- c("#008000CC")
  }
  if (no_indel_somatic == TRUE & no_snp == FALSE & no_loh == FALSE
     & no_indel_loh == FALSE){
    list1 <- list(x_s_snp, x_l_snp, x_l_indel)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_l_snp), 1:nrow(x_l_indel))
    circoscolors <- c("#FF0000CC", "#00FFFFCC", "#8000FFCC")
  }
  if (no_indel_somatic == TRUE & no_snp == FALSE & no_loh == FALSE
     & no_indel_loh == TRUE){
    list1 <- list(x_s_snp, x_l_snp)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_l_snp))
    circoscolors <- c("#FF0000CC", "#00FFFFCC")
  }
  if (no_indel_somatic == TRUE & no_snp == FALSE & no_loh == TRUE){
    list1 <- list(x_s_snp)
    idxs <- list(1:nrow(x_s_snp))
    circoscolors <- c("#FF0000CC")
  }

  mapping_matrices <- lapply(seq(1:length(list1)),
         function(i) return(get_mapping_matrix(list1[[i]], idxs[[i]])))
  oc_matrices <- lapply(mapping_matrices, add_default_value)
  oc_matrices <- lapply(oc_matrices, rename_chr)
  # fix bug omnic circos
  oc_matrices <- lapply(oc_matrices, duplicate_first_raw)

  return(list(map_mat = oc_matrices, circoscolors = circoscolors))
}

omicCircosUni <- function(listOfMap, label = NULL, minR, outfile,
                          circosColors = NULL) {
  #' omic Circos Uni
  #'
  #' @description Create the Circosplot
  #'
  #' @param listOfMap matrix. Mutationmatrix to be plotted
  #' @param minR numerical. Minimum radius
  #' @param outfile string. Name of output file
  #' @param circosColors vector of strings. Colors for Circosplot
  #' 
  #' @details This function plots the human genome on a circle.
  #' @details The mutations are then arranged by location on smaller
  #' @details concentric circle. Each mutation type gets an extra circle.
  #' @details If there is no mutation of a mutation type, the circle is
  #' @details excluded and there are less circles. The plot is stored in
  #' @details the given output file.
 
   # Human chromosomes
  data(UCSC.hg19.chr)
  ref <- UCSC.hg19.chr
  ref[,1] <- gsub("chr", "", ref[,1])
  db <- segAnglePo(ref, seg = as.character(unique(ref[,1])))
  colors <- rainbow(24, alpha = 0.8)
  
  # Parameters
  labelR <- 350
  chrR <- labelR - 25
  if(is.null(circosColors)){
    circosColors <- rainbow(length(listOfMap), alpha = 0.8)
  }
  circosW <- floor((chrR - minR) / length(listOfMap))
  circosR <- chrR - 1.5 * circosW
  
  # Plot
  pdf(outfile)
  par(mar=c(2, 2, 2, 2))
  plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, xlab = "", ylab = "",
       main = "")
  circos(R = chrR, cir = db, type = "chr", col = colors, print.chr.lab = TRUE,
         W = 2, scale = TRUE, lwd=1.5)
  if (!is.null(label)){
    circos(R = labelR, cir = db, W = 20, mapping = label, type = "label",
           col = "black", side = "out", cex = 0.4, lwd=1.5)
  }
  for (i in 1:length(listOfMap)) {
    if (!is.null(listOfMap[[i]])) {
      circos(R = circosR, cir = db, W = circosW, mapping = listOfMap[[i]],
             type = "b", col = circosColors[i], col.v = 4, lwd = 1.5)
      circosR <- circosR - 1.5 *circosW
    }
    # Label
    if (length(circosColors) == 4){
      text(0,75, "Somatic SNV", adj = 0, col = "#FF0000CC")
      text(0,50, "Somatic InDel", adj = 0, col = "#008000CC")
      text(0,25, "LoH SNV", adj = 0, col = "#00FFFFCC")
      text(0,0, "LoH InDel", adj = 0, col = "#8000FFCC")
    } else{
      l <- length(circosColors)
      l <- l*25
      if ("#FF0000CC" %in% circosColors){
        text(0,l, "Somatic SNV", adj = 0, col = "#FF0000CC")
        l <- l-25
      }
      if ("#008000CC" %in% circosColors){
        text(0,l, "Somatic InDel", adj = 0, col = "#008000CC")
        l <- l-25
      }
      if ("#00FFFFCC" %in% circosColors){
        text(0,l, "LoH SNV", adj = 0, col = "#00FFFFCC")
        l <- l-25
      }
      if ("#8000FFCC" %in% circosColors){
        text(0,l, "LoH InDel", adj = 0, col = "#8000FFCC")
      }
    }
  }
  dev.off()	
}

write_all_mut <- function(x_s, x_l = NULL){
  #' Write all Mutations
  #'
  #' @description Write all Mutations in a xlsx File
  #'
  #' @param x_s dataframe. List of somatic mutations
  #' @param x_l dataframe. List of LoH mutations
  #' 
  #' @return list of
  #' @return all_muts dataframe. Table of all mutations
  #' @return mut vector of strings. List of mutated genes
  #' 
  #' @details A table with all mutations (somatic and LoH) is saved in
  #' @detials All_Mutations_Somatic.xlsx.
  mutations_somatic <- as.character(x_s$Gene.refGene)
  mutations_somatic <- unique(mutations_somatic)
  if (!is.null(x_l)){
    mutations_loh <- as.character(x_l$Gene.refGene)
    mutations_loh <- unique(mutations_loh)
  } else {
    mutations_loh <- c()
  }

  mut <- unique(c(mutations_somatic, mutations_loh))

  tmp1 <- x_s[, c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                  "Variant_Allele_Frequency", "Variant_Reads",
                  "AAChange.SnpEff", "is_tumorsuppressor",
                  "is_oncogene", "is_hotspot", "target", "gnomAD_exome_NFE",
                  "CADD13_PHRED", "condel.label", "CLINSIG.SnpEff",
                  "cosmic84_coding"),
              drop = FALSE]
  if (!is.null(x_l)){
    tmp2 <- x_l[, c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                    "VAF_Tumor", "Count_Tumor", "AAChange.SnpEff",
                    "is_tumorsuppressor", "is_oncogene",
                    "is_hotspot", "target", "gnomAD_exome_NFE", "CADD13_PHRED",
                    "condel.label", "CLINSIG.SnpEff", "cosmic84_coding"),
                drop = FALSE]
    }
  col_names <- c("Symbol", "GeneName", "ExonicFunc", "VAF", "Reads",
                 "AAChange", "TSG", "OG", "HS", "target", "MAF", "CADD13",
                 "Condel", "CLINSIG", "COSMIC84")
  colnames(tmp1) <- col_names
  if (!is.null(x_l)) {
    colnames(tmp2) <- col_names
  }
  if (!is.null(x_l)) {
    all_mutations <- rbind(tmp1, tmp2)
    } else {
      all_mutations <- tmp1
    }
  write.xlsx(all_mutations, file = "All_Mutations_Somatic_LoH.xlsx",
             keepNA = FALSE, rowNames = FALSE, firstRow = TRUE)
  return(list(all_muts = all_mutations, mut = mut))
}

prep_pwa <- function(targets, mut){
  #' Preparation for Pathway Analysis
  #'
  #' @description Preparation for Pathway Analysis
  #'
  #' @param targets dataframe. Table of sequenced targets
  #' @param mut vector of strings. List of mutated genes
  #' 
  #' @return list of
  #' @return de_genes dataframe. List of mutated genes
  #' @return universe dataframe. List of targeted genes
  #' 
  #' @details For a pathway analysis the mutated genes and the universe of all
  #' @details eventually mutated genes must be known. Therefor a list of
  #' @details targeted genes is extracted and then the mutated genes are found.
  xx <- as.list(org.Hs.egSYMBOL)
  y <- unlist(xx)
  y2 <- unique(y)

  t <- as.character(targets$V1)

  id <- t %in% y2

  t2 <- t[id]

  mut.entrez <- as.character(mget(mut, revmap(org.Hs.egSYMBOL),
                                  ifnotfound = NA))
  t2.entrez <- as.character(mget(t2, revmap(org.Hs.egSYMBOL), ifnotfound = NA))
  return(list(de_genes = mut.entrez, universe = t2.entrez))
}

hyperG <- function(geneSets, DEgenes, universe, org.library, cutoff = 0.1,
                   mincount = 2, parallel = T, adj.P.Val = F, set.size = NULL){
  #' Functional Analysis
  #'
  #' @description Functinal pathway Analysis with hyper geometric testing
  #'
  #' @param geneSets dataframe. Table of pathways with corresponding genes
  #' @param DEgenes dataframe. List of mutated genes
  #' @param universe dataframe. List of all possibly mutated genes.
  #' @param org.library dataframe. Genome wide annotation
  #' @param cutoff numerical. Cutoff for p-value (default: 0.1)
  #' @param mincount numerical. Minimal number of counts
  #' @param parallel logical. Calculate parallel or not (default: T)
  #' @param adj.P.Val logical. Use adjusted p-value instead (default: F)
  #' @param set.size vector. Define minimal and maximal pathway size (default: NULL)
  #' 
  #' @return results dataframe. Ordered list of pathways with p-value.
  #' 
  #' @details Given a list of pathways the goal of this function is to find by
  #' @details a hypergeomtrical test those who are mostly influenced by a
  #' @details number of mutated genes.
  library(foreach)
  library(doMC)
  if (parallel){
    registerDoMC(cores = detectCores())
    cores = detectCores()
  } else {
    cores = 1
  }
  if (!is.null(set.size)){
    print("Set Size Limits")
    idx <- lapply(geneSets, function(x) {
      length(x) <= set.size[2] & length(x) >= set.size[1]
    }
    )
    geneSets <- geneSets[unlist(idx)]
  }
  results <- mclapply(1:length(geneSets), function(i){
    results <- matrix(data = NA, ncol = 8, nrow = 1)
    colnames(results) <- c("Term", "Count", "Size", "p-value", "adj.P.Val",
                           "odds ratio", "Entrez", "Symbol")
    geneSet <- intersect(universe, geneSets[[i]])
    a <- length(intersect(DEgenes, geneSet))
    b <- length(setdiff(DEgenes, intersect(DEgenes, geneSet)))
    c <- length(setdiff(geneSet, intersect(DEgenes, geneSet)))
    d <- length(setdiff(universe, DEgenes)) - c
    contigency.matrix <- cbind(c(a, b), c(c, d))
    res <- fisher.test(contigency.matrix, alternative = "greater")
    results[1, "Term"] <- names(geneSets)[i]
    results[1, "Count"] <- a
    results[1, "Size"] <- length(geneSets[[i]])
    results[1, "p-value"] <- res$p.value
    results[1, "odds ratio"] <- res$estimate[[1]]
    # find genes annotated in the consensus term
    if(a > 0){
      genes <- intersect(DEgenes, geneSet)
      eid <- genes
      eid <- eid[order(eid)]
      results[1, "Entrez"] <- paste(eid, collapse = "|")
    }
    return(results)
  }
  , mc.cores = cores)
  
  results <- as.data.frame(do.call(rbind, results))
  for(i in c(2, 3, 4, 5)){
    results[, i] <- as.numeric(as.character(results[, i]))
  }
  
  if(nrow(results) != 1){
    results <- results[order(results[, "p-value"],decreasing = FALSE), ]
    results[, "adj.P.Val"] <- p.adjust(results[, "p-value"], "BH")
    if(adj.P.Val){
      results <- as.data.frame(subset(results, results[, "adj.P.Val"] <= cutoff))
    }else{
      results <- as.data.frame(subset(results, results[, "p-value"] <= cutoff))
    }
    results <- as.data.frame(subset(results, results[, "Count"] >= mincount))
  }else results <- as.data.frame(results)
  
  org.symb <- gsub(".db", "SYMBOL", org.library)
  # find genes 
  results$Symbol <- sapply(results$Entrez, function(x){
    y <- unlist(strsplit(as.character(x), "|", fixed = T))
    syms <- paste(unlist(lapply(mget(y, eval(parse(text = org.symb)),
                                     ifnotfound = NA), function(x) x[1])),
                  collapse = "|")
  })
  return(results)
}


get_terms <- function(dataset, outfile, mut.entrez, t2.entrez){
  #' Pathway Analysis
  #'
  #' @description Pathway Analysis with hyper geometric testing
  #'
  #' @param dataset dataframe. Table of pathways with corresponding genes
  #' @param outfile string. Name of output file
  #' @param mut.entrez dataframe. List of mutated genes
  #' @param t2.entrez dataframe. List of all possibly mutated genes.
  #' 
  #' @return results dataframe. Ordered list of pathway terms with p-value.
  #' 
  #' @details In this function the functional pathway analysis of hyperG is
  #' @details used and a corresponding xlsx file is written.
  ds_test <- hyperG(geneSets = dataset, DEgenes = mut.entrez,
                    universe = t2.entrez, org.library = "org.Hs.eg.db",
                    cutoff = 0.05, mincount = 3, parallel = T, adj.P.Val = F)
  if (nrow(ds_test) > 0){
    write.xlsx(ds_test, outfile, keepNA = FALSE, rowNames = FALSE,
               firstRow = TRUE)
  }
  return(ds_res = ds_test)
}

write_mtb_genesets <- function(mulist, mtb.genesets, outfile_mtb_geneset){
  #' Important Pathways
  #' 
  #' @description Find mutations in important pathways
  #'
  #' @param mulist dataframe. Table of mutated genes
  #' @param mtb.genesets dataframe. List of important pathways
  #' @param outfile_mtb_geneset string. Name of output file
  #' 
  #' @return ch_mat matrix. Result matrix
  #' 
  #' @details For each mutated gene is checked, whether they are in one of five
  #' @details important pathways: "PI3K-AKT-mTOR", "RAF-MEK-ERK",
  #' @details "DNA Damage Response", "Cell Cycle", "Tyrosine Kinases".
  check_matrix <- matrix(0, nrow = length(mulist), ncol = length(mtb.genesets))
  rownames(check_matrix) <- mulist
  colnames(check_matrix) <- names(mtb.genesets)

  for (i in 1:length(mulist)){
    for (j in 1:length(mtb.genesets)){
      if (mulist[i] %in% mtb.genesets[[j]]){
        check_matrix[i, j] <- 1
      }
    }
  }
  write.xlsx(check_matrix, outfile_mtb_genesets, rowNames = TRUE,
             firstRow = TRUE)
  return(ch_mat = check_matrix)
}

imp_pws <- function(ch_mat, all_muts){
  #' Important Pathways
  #' 
  #' @description Sort mutations into important pathways
  #'
  #' @param ch_mat matrix. Result matrix
  #' @param all_muts dataframe. Table of mutations
  #' 
  #' @return important_pws dataframe. Table of mutations in important pathways
  #' 
  #' @details A Table for the important Pathways is built that contains all the
  #' @details mutations found belonging to them. The results are writte in file
  #' @details "MTB_Pathways-Gene.txt". 
  sink(file = "MTB_Pathways-Gene.txt", append = T, split = T)
  print("PI3K-AKT-mTOR:")
  tmp <- rownames(ch_mat)[ch_mat[, 3] == 1]
  print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  print("RAF-MEK-ERK:")
  tmp <- rownames(ch_mat)[ch_mat[, 4] == 1]
  print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  print("DNA Damage Response:")
  tmp <- rownames(ch_mat)[ch_mat[, 2] == 1]
  print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  print("Cell Cycle:")
  tmp <- rownames(ch_mat)[ch_mat[, 1] == 1]
  print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  print("Tyrosine Kinases:")
  tmp <- rownames(ch_mat)[ch_mat[, 5] == 1]
  print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  sink()
  #print("PI3K-AKT-mTOR:")
  pi3k <- rownames(ch_mat)[ch_mat[, 3] == 1]
  if (length(pi3k) != 0){
    pi3k_genes <-  all_muts[match(pi3k, all_muts$Symbol), ]
    pi3k_genes$Pathway <- "."
    pi3k_genes[1, "Pathway"] <- "PI3K-AKT-mTOR"
  } else {
    pi3k_genes <- c()
    }
  #print("RAF-MEK-ERK:")
  raf <- rownames(ch_mat)[ch_mat[, 4] == 1]
  if (length(raf) != 0){
    raf_genes <- all_muts[match(raf, all_muts$Symbol), ]
    raf_genes$Pathway <- "."
    raf_genes[1, "Pathway"] <- "RAF-MEK-ERK"
  } else {
    raf_genes <- c()
    }
  #print("DNA Damage Response:")
  dna_damage <- rownames(ch_mat)[ch_mat[, 2] == 1]
  if (length(dna_damage) != 0){
    dna_damage_genes <- all_muts[match(dna_damage, all_muts$Symbol), ]
    dna_damage_genes$Pathway <- "."
    dna_damage_genes[1, "Pathway"] <- "DNA Damage Response"
  } else {
    dna_damage_genes <- c()
    }
  #print("Cell Cycle:")
  cell_cycle <- rownames(ch_mat)[ch_mat[, 1] == 1]
  if (length(cell_cycle) != 0){
    cell_cycle_genes <- all_muts[match(cell_cycle, all_muts$Symbol), ]
    cell_cycle_genes$Pathway <- "."
    cell_cycle_genes[1, "Pathway"] <- "Cell Cycle"
  } else {
    cell_cycle_genes <- c()
    }
  #print("Tyrosine Kinases:")
  tyrosine <- rownames(ch_mat)[ch_mat[, 5] == 1]
  if (length(tyrosine) != 0){
    tyrosine_genes <- all_muts[match(tyrosine, all_muts$Symbol), ]
    tyrosine_genes$Pathway <- "."
    tyrosine_genes[1, "Pathway"] <- "Tyrosine Kinases"
  } else {
    tyrosine_genes <- c()
    }

  important_pathways <- rbind(pi3k_genes, raf_genes, dna_damage_genes,
                              cell_cycle_genes, tyrosine_genes)
  return(important_pws = important_pathways)
}
