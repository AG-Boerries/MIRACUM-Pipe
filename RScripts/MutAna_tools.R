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

find_indel_2 <- function(list){
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
  id.ref <- which(nchar(as.character(list$Ref)) > 1)
  id.alt <- which(nchar(as.character(list$Alt)) > 1)
  id <- unique( c(id.ref, id.alt) )
  return(id)
}

# div <- function(x_s, x_l, no_loh) {
#   #' Mutation separation
#   #'
#   #' @description Separate mutions
#   #'
#   #' @param x_s dataframe. List of somatic mutations
#   #' @param x_l dataframe. List of LoH mutations
#   #' @param no_loh logical. Logical describing existence of LoH mutations
#   #'
#   #' @return list of
#   #' @return x_s_snp dataframe. List of somatic SNVs
#   #' @return x_s_indel dataframe. List of somatic InDels
#   #' @return x_l_snp dataframe. List of LoH SNVs
#   #' @return x_l_indel dataframe. List of LoH InDels
#   #' @return no_loh logical. Describing existence of LoH mutations
#   #' @return no_indel_somatic logical. Describing existence of somatic Indels
#   #' @return no_snp logical. Describing existence of SNVs
#   #' @return no_indel_loh logical. Describing existence of LoH Indels
#   #'
#   #' @details Build separate dataframes for different type of mutations.
#   #' @details Split somatic and LoH mutations in SNVs and InDels. Return also
#   #' @details some logicals for existence of mutations at all.

#   no_indel_somatic <- FALSE
#   # if (protocol == "Tumor_Only" & manifest == "V5UTR") {
#   #  indel_s <- find_indel_2(x_s)
#   # } else {
#   #  indel_s <- find_indel(x_s)
#   # }
#   indel_s <- find_indel(x_s)
#   if (length(indel_s) > 0) {
#     x_s_snp <- x_s[-indel_s, ]
#     x_s_indel <- x_s[indel_s, ]
#   } else {
#     x_s_snp <- x_s
#     x_s_indel <- data.frame()
#     cat("No Indels in Somatic!\n")
#     no_indel_somatic <- TRUE
#   }

#   if (dim(x_s_snp)[1] > 0) {
#     no_snp <- FALSE
#   } else {
#     no_snp <- TRUE
#   }

#   no_indel_loh <- FALSE
#   no_snp_loh <- FALSE
#   if (!no_loh) {
#     indel_l <- find_indel(x_l)
#     if (length(indel_l) == 0) {
#       no_indel_loh <- TRUE
#       x_l_snp <- x_l
#       x_l_indel <- data.frame()
#       cat("No Indels in LOH!\n")
#     } else if (length(indel_l) == dim(x_l)[1]) {
#       no_snp_loh <- TRUE
#       x_l_snp <- data.frame()
#       x_l_indel <- x_l
#     } else {
#       x_l_snp <- x_l[-indel_l, ]
#       x_l_indel <- x_l[indel_l, ]
#     }
#   } else {
#     x_l_snp <- data.frame()
#     x_l_indel <- data.frame()
#   }
#   return(list(
#     x_s_snp = x_s_snp, x_s_indel = x_s_indel, x_l_snp = x_l_snp,
#     x_l_indel = x_l_indel, no_loh = no_loh,
#     no_indel_somatic = no_indel_somatic, no_snp = no_snp,
#     no_indel_loh = no_indel_loh, no_snp_loh = no_snp_loh
#   ))
# }

div <- function(x_s, x_l, no_loh, protocol, sureselect_type) {
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
  if (protocol == "Tumor_Only" & sureselect_type == "V5UTR") {
    indel_s <- find_indel_2(x_s)
  } else {
    indel_s <- find_indel(x_s)
  }
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
  no_snp_loh <- FALSE
  if (!no_loh){
    indel_l <- find_indel(x_l)
    if (length(indel_l) == 0){
      no_indel_loh <- TRUE
      x_l_snp <- x_l
      x_l_indel <- data.frame()
      cat("No Indels in LOH!\n")
      } else if (length(indel_l) == dim(x_l)[1]) {
      no_snp_loh <- TRUE
      x_l_snp <- data.frame()
      x_l_indel <- x_l
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
              no_indel_loh = no_indel_loh, no_snp_loh = no_snp_loh))
}

mut_tab <- function(x_s_snp, x_s_indel, x_l_snp, x_l_indel, protocol) {
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
  if (protocol == "somaticGermline" | protocol == "somatic") {
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
  
    return(muta_tab = muta_tab)
  } else {
    muta_tab <- matrix(data = '-', nrow = 4, ncol = 6)
    colnames(muta_tab) <- c("Mutationtype", "Zygosity", "Number", "Tumorsuppressor", "Oncogene", "Hotspot")
    muta_tab[,1] <- c("SNV", "SNV", "InDel", "InDel")
    
    muta_tab[c(1, 3), 2] <- "homozygous"
    muta_tab[c(2, 4), 2] <- "heterozygous"
    
    muta_tab[1, 3] <- sum(x_s_snp$Zygosity == "hom")
    muta_tab[2, 3] <- sum(x_s_snp$Zygosity == "het")
    muta_tab[3, 3] <- sum(x_s_indel$Zygosity == "hom")
    muta_tab[4, 3] <- sum(x_s_indel$Zygosity == "het")
    
    muta_tab[1, 4] <- sum(x_s_snp$is_tumorsuppressor == 1
                              & x_s_snp$Zygosity == "hom")
    muta_tab[2, 4] <- sum(x_s_snp$is_tumorsuppressor == 1
                          & x_s_snp$Zygosity == "het")
    muta_tab[3, 4] <- sum(x_s_indel$is_tumorsuppressor == 1
                          & x_s_indel$Zygosity == "hom")
    muta_tab[4, 4] <- sum(x_s_indel$is_tumorsuppressor == 1
                          & x_s_indel$Zygosity == "het")
    
    muta_tab[1, 5] <- sum(x_s_snp$is_oncogene == 1
                          & x_s_snp$Zygosity == "hom")
    muta_tab[2, 5] <- sum(x_s_snp$is_oncogene == 1
                          & x_s_snp$Zygosity == "het")
    muta_tab[3, 5] <- sum(x_s_indel$is_oncogene == 1
                          & x_s_indel$Zygosity == "hom")
    muta_tab[4, 5] <- sum(x_s_indel$is_oncogene == 1
                          & x_s_indel$Zygosity == "het")
    
    muta_tab[1, 6] <- sum(x_s_snp$is_hotspot != 0
                          & x_s_snp$Zygosity == "hom")
    muta_tab[2, 6] <- sum(x_s_snp$is_hotspot != 0
                          & x_s_snp$Zygosity == "het")
    muta_tab[3, 6] <- sum(x_s_indel$is_hotspot != 0
                          & x_s_indel$Zygosity == "hom")
    muta_tab[4, 6] <- sum(x_s_indel$is_hotspot != 0
                          & x_s_indel$Zygosity == "het")

    return(muta_tab = muta_tab)
  }
  print(muta_tab)
}

mut_stats <- function(x_s, x_l = NULL, tumbu, protocol) {
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
  if (protocol == "somaticGermline" | protocol == "somatic"){
    print(paste(dim(x_s)[1], "somatic mutations", sep = " "))
  } else {
    print(paste(dim(x_s)[1], "mutations", sep = " "))
  }
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
  totalmutationnumber <- dim(x_s)[1] + dim_x_l
  somaticmutations <- dim(x_s)[1]
  lohmutations <- dim_x_l
  return(list(tot_mut <- totalmutationnumber, som_mut <- somaticmutations,
              loh_mut <- lohmutations))
}

tables <- function(x_s, x_l = NULL, protocol) {
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
  #' @return sm_table dataframe. List of somatic mutations
  #' @return lm_table dataframe. List of LoH mutations
  #'
  #' @details The Tables TumorSuppressor-OncogeneTable,
  #' @details somaticMutations and lohMutations are generated and
  #' @details stored.

  col_names <- c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                 "AAChange", "Variant_Allele_Frequency", "Zygosity",
                 "Variant_Reads", "is_tumorsuppressor", "is_oncogene",
                 "is_hotspot", "target", "AF_popmax",
                 "CADD_phred", "condel.label", "REVEL_score", "CLNSIG", "InterVar_automated",
                 "cosmic_coding", "Chr", "Start", "Ref",
                 "Alt")

  ts_og_table <- data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(ts_og_table) <- col_names
  sm_table <- data.frame(matrix(ncol = length(col_names), nrow = 0))
  colnames(sm_table) <- col_names

  if (!is.null(x_s) && dim(x_s)[1]) {
    ts_og_table <- x_s[
      x_s$is_tumorsuppressor == 1 |
      x_s$is_oncogene == 1,
      col_names,
      drop = FALSE
    ]

    sm_table <- x_s[, col_names, drop = FALSE]
  }

  if (!is.null(x_l)){
    lm_table <- x_l[, c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                       "AAChange", "VAF_Normal", "VAF_Tumor",
                       "Count_Normal", "Count_Tumor", "is_tumorsuppressor",
                       "is_oncogene", "is_hotspot", "target",
                       "AF_popmax", "CADD_phred", "condel.label", "REVEL_score",
                       "CLNSIG", "InterVar_automated", "cosmic_coding", "Chr", "Start", "Ref",
                       "Alt"), drop = FALSE]

  } else {
    lm_table <- data.frame()
  }

  return(list(ts_og_table = ts_og_table, sm_table = sm_table,
              lm_table = lm_table))
}

get_mapping_matrix <- function(annovar_table, row_index) {
  if (length(row_index) == 0) return(NULL)
  col_index <- match(c("Chr", "Start", "Gene.refGene"), colnames(annovar_table))
  new_table <- annovar_table[row_index, col_index]
  new_table[] <- lapply(new_table, as.character)
  new_table$Start <- as.numeric(new_table$Start)
  return(new_table)
}

add_default_value <- function(mapping_matrix) {
  if (is.null(mapping_matrix)) return(mapping_matrix)
  new_matrix <- cbind(mapping_matrix, rep(1, nrow(mapping_matrix)))
  colnames(new_matrix)[ncol(new_matrix)] <- "Value"
  return(new_matrix)	
}

rename_chr <- function(mapping_matrix) {
  new_matrix <- mapping_matrix
  new_matrix[,1] <- gsub("chr", "", new_matrix[,1])
  return(new_matrix)	
}

duplicate_first_raw <- function(mapping_matrix) {
  if (is.null(mapping_matrix)) return(NULL)
  new_matrix <- mapping_matrix
  new_matrix <- rbind(new_matrix[1,], new_matrix)
  new_matrix[1,c(4:ncol(new_matrix))] <- 0
  return(new_matrix)
}

circos_colors <- function(
  x_s_snp = NULL,
  x_s_indel = NULL,
  x_l_snp = NULL,
  x_l_indel = NULL,
  no_loh,
  no_indel_somatic,
  no_snp,
  no_indel_loh,
  no_snp_loh
) {
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
  # print(no_loh)
  # print(no_indel_somatic)
  # print(no_snp)
  # print(no_indel_loh)
  # print(no_snp_loh)
  if (no_indel_somatic == FALSE & no_snp == FALSE & no_loh == FALSE & no_indel_loh == FALSE & no_snp_loh == FALSE){
    list1 <- list(x_s_snp, x_s_indel, x_l_snp, x_l_indel)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_s_indel), 1:nrow(x_l_snp), 1:nrow(x_l_indel))
    circoscolors <- c("#FF0000CC", "#008000CC", "#00FFFFCC", "#8000FFCC")
  } else if (no_indel_somatic == FALSE & no_snp == FALSE & no_loh == FALSE & no_indel_loh == FALSE & no_snp_loh == TRUE){
    list1 <- list(x_s_snp, x_s_indel, x_l_indel)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_s_indel), 1:nrow(x_l_indel))
    circoscolors <- c("#FF0000CC", "#008000CC", "#8000FFCC")
  } else if (no_indel_somatic == FALSE & no_snp == FALSE & no_loh == FALSE & no_indel_loh == TRUE){
    list1 <- list(x_s_snp, x_s_indel, x_l_snp)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_s_indel), 1:nrow(x_l_snp))
    circoscolors <- c("#FF0000CC", "#008000CC", "#00FFFFCC")
  } else if (no_indel_somatic == FALSE & no_snp == FALSE & no_loh == TRUE){
    list1 <- list(x_s_snp, x_s_indel)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_s_indel))
    circoscolors <- c("#FF0000CC", "#008000CC")
  } else if (no_indel_somatic == FALSE & no_snp == TRUE & no_loh == FALSE & no_indel_loh == TRUE){
    list1 <- list(x_s_indel, x_l_snp)
    idxs <- list(1:nrow(x_s_indel), 1:nrow(x_l_snp))
    circoscolors <- c("#008000CC", "#00FFFFCC")
  } else if (no_indel_somatic == FALSE & no_snp == TRUE & no_loh == TRUE){
    list1 <- list(x_s_indel)
    idxs <- list(1:nrow(x_s_indel))
    circoscolors <- c("#008000CC")
  } else if (no_indel_somatic == TRUE & no_snp == FALSE & no_loh == FALSE & no_indel_loh == FALSE){
    list1 <- list(x_s_snp, x_l_snp, x_l_indel)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_l_snp), 1:nrow(x_l_indel))
    circoscolors <- c("#FF0000CC", "#00FFFFCC", "#8000FFCC")
  } else if (no_indel_somatic == TRUE & no_snp == FALSE & no_loh == FALSE & no_indel_loh == TRUE){
    list1 <- list(x_s_snp, x_l_snp)
    idxs <- list(1:nrow(x_s_snp), 1:nrow(x_l_snp))
    circoscolors <- c("#FF0000CC", "#00FFFFCC")
  } else if (no_indel_somatic == TRUE & no_snp == FALSE & no_loh == TRUE){
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

omicCircosUni <- function(
  listOfMap,
  label = NULL,
  minR,
  outfile,
  circosColors = NULL,
  protocol,
  sureselect,
  sureselect_type
) {
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
  
  if (protocol == "panelTumor" & !(sureselect_type %in% c("V6","V5UTR","V6UTR"))) {
    tg <- read.delim(file = sureselect, header = FALSE)

    hili <- as.data.frame(matrix(NA, nrow = nrow(tg), ncol = 7))
    hili$V7 <- hili$V8 <- "#fff68f"
    hili$V1 <- 50
    hili$V2 <- 250
    hili$V3 <- hili$V5 <- tg$V1
    hili$V4 <- tg$V2
    hili$V6 <- tg$V3
    hili$V5 <- gsub("chr", "", hili$V5)
    hili$V3 <- gsub("chr", "", hili$V3)
    hili <- as.matrix(hili)
  }

  # Plot
  pdf(outfile)
  par(mar=c(2, 2, 2, 2))
  plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, xlab = "", ylab = "",
       main = "")
  circos(R = chrR, cir = db, type = "chr", col = colors, print.chr.lab = TRUE,
         W = 2, scale = TRUE, lwd=1.5)
  if (protocol == "panelTumor" & !(sureselect_type %in% c("V6","V5UTR","V6UTR"))){
    for (i in 1:dim(hili)[1]){
          circos(R=chrR, cir=db, W=40, mapping=hili[i, ], type = "hl", lwd=1.5)
        }
  }
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
  }
    # Label
  if (protocol == "somaticGermline" | protocol == "somatic") {
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
  } else {
    if (length(circosColors) == 4){
      text(0,75, "SNV", adj = 0, col = "#FF0000CC")
      text(0,50, "InDel", adj = 0, col = "#008000CC")
      text(0,25, "LoH SNV", adj = 0, col = "#00FFFFCC")
      text(0,0, "LoH InDel", adj = 0, col = "#8000FFCC")
    } else{
      l <- length(circosColors)
      l <- l*25
      if ("#FF0000CC" %in% circosColors){
        text(0,l, "SNV", adj = 0, col = "#FF0000CC")
        l <- l-25
      }
      if ("#008000CC" %in% circosColors){
        text(0,l, "InDel", adj = 0, col = "#008000CC")
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

# omicCircosFus2 <- function(
#   listOfMap,
#   fusions,
#   label = NULL,
#   minR,
#   outfile,
#   circosColors = NULL,
#   protocol,
#   sureselect
# ) {
#   #' omic Circos with Fusions
#   #'
#   #' @description Create the Circosplot
#   #'
#   #' @param listOfMap matrix. Mutationmatrix to be plotted
#   #' @param minR numerical. Minimum radius
#   #' @param outfile string. Name of output file
#   #' @param circosColors vector of strings. Colors for Circosplot
#   #' @param protocol string. Name of the analysis protocol
#   #' @param mode string. Name of the capture kit.
#   #'
#   #' @details This function plots the human genome on a circle.
#   #' @details The mutations are then arranged by location on smaller
#   #' @details concentric circle. Each mutation type gets an extra circle.
#   #' @details If there is no mutation of a mutation type, the circle is
#   #' @details excluded and there are less circles. The plot is stored in
#   #' @details the given output file.
#   #' @details For smaller panels the captured regions are highlighted in
#   #' @details the circosplot. That is controlled by protocol and mode. 


#   # Human chromosomes
#   data(UCSC.hg19.chr)
#   ref <- UCSC.hg19.chr
#   ref[,1] <- gsub("chr", "", ref[,1])
#   db <- segAnglePo(ref, seg = as.character(unique(ref[,1])))
#   colors <- rainbow(24, alpha = 0.8)
  
#   # Parameters
#   labelR <- 350
#   chrR <- labelR - 25
#   if(is.null(circosColors)){
#     circosColors <- rainbow(length(listOfMap), alpha = 0.8)
#   }
#   circosW <- floor((chrR - minR) / length(listOfMap))
#   circosR <- chrR - 1.5 * circosW + 50
  
#   # Add highlighted area for targeted regions
#   if (protocol == "panelTumor"){
#     tg <- read.delim(file = sureselect, header = FALSE)
    
#     hili <- as.data.frame(matrix(NA, nrow = nrow(tg), ncol = 7))
#     hili$V7 <- hili$V8 <- "#fff68f"
#     hili$V1 <- 50
#     hili$V2 <- 250
#     hili$V3 <- hili$V5 <- tg$V1
#     hili$V4 <- tg$V2
#     hili$V6 <- tg$V3
#     hili$V5 <- gsub("chr", "", hili$V5)
#     hili$V3 <- gsub("chr", "", hili$V3)
#     hili <- as.matrix(hili)
#   }
  
#   # Plot
#   pdf(outfile)
#   par(mar=c(2, 2, 2, 2))
#   plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, xlab = "", ylab = "",
#        main = "")
#   circos(R = chrR, cir = db, type = "chr", col = colors, print.chr.lab = TRUE,
#          W = 2, scale = TRUE, lwd=1.5)
#   if(protocol == "panelTumor"){
#       for (i in 1:dim(hili)[1]){
#         circos(R=chrR, cir=db, W=40, mapping=hili[i, ], type = "hl", lwd=1.5)
#     }
#   }
#   if (!is.null(label)){
#     circos(R = labelR, cir = db, W = 20, mapping = label, type = "label",
#            col = "black", side = "out", cex = 0.4, lwd=1.5)
#   }
#   for (i in 1:length(listOfMap)) {
#     if (!is.null(listOfMap[[i]])) {
#       circos(R = circosR, cir = db, W = circosW, mapping = listOfMap[[i]],
#              type = "b", col = circosColors[i], col.v = 4, lwd = 1.5)
#       circosR <- circosR - 1.5 *circosW
#     }
#   }
#   colors_fus <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
#                   "#D55E00", "#CC79A7", "#999999")
#   fusions2 <- fusions[, c(2:3, 1, 5:6, 4)]
#   fusions2 <- as.data.frame(fusions2)
#   circos(R = 75, cir = db, W = circosW, mapping = fusions2,
#          type = "link", lwd = 2, col = colors_fus)
#   # Label
#   if (length(circosColors) == 4){
#     text(0,75, "SNV", adj = 0, col = "#FF0000CC")
#     text(0,50, "InDel", adj = 0, col = "#008000CC")
#     text(0,25, "LoH SNV", adj = 0, col = "#00FFFFCC")
#     text(0,0, "LoH InDel", adj = 0, col = "#8000FFCC")
#   } else{
#     l <- length(circosColors)
#     l <- l*25
#     if ("#FF0000CC" %in% circosColors){
#       text(0,l, "SNV", adj = 0, col = "#FF0000CC")
#       l <- l-25
#     }
#     if ("#008000CC" %in% circosColors){
#       text(0,l, "InDel", adj = 0, col = "#008000CC")
#       l <- l-25
#     }
#     if ("#00FFFFCC" %in% circosColors){
#       text(0,l, "LoH SNV", adj = 0, col = "#00FFFFCC")
#       l <- l-25
#     }
#     if ("#8000FFCC" %in% circosColors){
#       text(0,l, "LoH InDel", adj = 0, col = "#8000FFCC")
#     }
#   }
#   dev.off()	
# }

omicCircosFus2 <- function(
  listOfMap,
  fusions,
  label = NULL,
  minR,
  outfile,
  circosColors = NULL,
  mode = "V6",
  protocol = "Tumor_Normal",
  path_data,
  trgt
) {
  #' omic Circos with Fusions
  #'
  #' @description Create the Circosplot
  #'
  #' @param listOfMap matrix. Mutationmatrix to be plotted
  #' @param minR numerical. Minimum radius
  #' @param outfile string. Name of output file
  #' @param circosColors vector of strings. Colors for Circosplot
  #' @param protocol string. Name of the analysis protocol
  #' @param mode string. Name of the capture kit.
  #' 
  #' @details This function plots the human genome on a circle.
  #' @details The mutations are then arranged by location on smaller
  #' @details concentric circle. Each mutation type gets an extra circle.
  #' @details If there is no mutation of a mutation type, the circle is
  #' @details excluded and there are less circles. The plot is stored in
  #' @details the given output file.
  #' @details For smaller panels the captured regions are highlighted in
  #' @details the circosplot. That is controlled by protocol and mode. 
  require(OmicCircos)

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
  circosR <- chrR - 1.5 * circosW + 50
  
  # Add highlighted area for targeted regions
  if (protocol != "Tumor_Normal" & mode %in% c("TruSight_Tumor", "TSO500", "TruSight_Amplicon", "Patho")){
    # Targeted Area Highlighting in CircosPlot 
    tg <- read.delim(file = trgt, header = FALSE)
    if (mode == "Patho"){
      tg <- tg[-c(1:2), ]
    }
    hili <- as.data.frame(matrix(NA, nrow = nrow(tg), ncol = 7))
    hili$V7 <- hili$V8 <- "#fff68f"
    hili$V1 <- 75
    hili$V2 <- 300
    hili$V3 <- hili$V5 <- tg$V1
    hili$V4 <- tg$V2
    hili$V6 <- tg$V3
    hili$V5 <- gsub("chr", "", hili$V5)
    hili$V3 <- gsub("chr", "", hili$V3)
    hili <- as.matrix(hili)
    if (mode == "TruSight_Amplicon"){
      hili <- hili[-c(213),]
    }
  }
  
  
  # Plot
  pdf(outfile)
  par(mar=c(2, 2, 2, 2))
  plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, xlab = "", ylab = "",
       main = "")
  circos(R = chrR, cir = db, type = "chr", col = colors, print.chr.lab = TRUE,
         W = 2, scale = TRUE, lwd=1.5)
  if(protocol != "Tumor_Normal" & mode %in% c("TruSight_Tumor", "TSO500",
                                              "TruSight_Amplicon", "Patho")){
      for (i in 1:dim(hili)[1]){
        circos(R=chrR, cir=db, W=40, mapping=hili[i, ], type = "hl", lwd=1.5)
    }
  }
  if (!is.null(label)){
    circos(R = labelR, cir = db, W = 20, mapping = label, type = "label",
           col = "black", side = "out", cex = 0.4, lwd=1.5)
  }
  for (i in 1:length(listOfMap)) {
    if (!is.null(listOfMap[[i]]) & !is.na(listOfMap[[i]])) {
      circos(R = circosR, cir = db, W = circosW, mapping = listOfMap[[i]],
             type = "b", col = circosColors[i], col.v = 4, lwd = 1.5)
      circosR <- circosR - 1.5 *circosW
    }
  }
  colors_fus <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                  "#D55E00", "#CC79A7", "#999999")
  fusions2 <- fusions[, c(2:3, 1, 5:6, 4)]
  fusions2 <- as.data.frame(fusions2)
  circos(R = 75, cir = db, W = circosW, mapping = fusions2,
         type = "link", lwd = 2, col = colors_fus)
  # Label
  if (length(circosColors) == 4){
    text(0,75, "SNV", adj = 0, col = "#FF0000CC")
    text(0,50, "InDel", adj = 0, col = "#008000CC")
    text(0,25, "LoH SNV", adj = 0, col = "#00FFFFCC")
    text(0,0, "LoH InDel", adj = 0, col = "#8000FFCC")
  } else{
    l <- length(circosColors)
    l <- l*25
    if ("#FF0000CC" %in% circosColors){
      text(0,l, "SNV", adj = 0, col = "#FF0000CC")
      l <- l-25
    }
    if ("#008000CC" %in% circosColors){
      text(0,l, "InDel", adj = 0, col = "#008000CC")
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
  #' @details A table with all mutations (somatic and LoH) is saved.
  col_names <- c("Symbol", "GeneName", "ExonicFunc", "VAF", "Reads",
                 "AAChange", "TSG", "OG", "HS", "target", "MAF", "CADD",
                 "Condel", "REVEL_score", "CLNSIG", "InterVar_automated", "COSMIC")
  all_mutations <- data.frame(matrix(ncol = length(col_names), nrow = 0))

  if (!is.null(x_s) && dim(x_s)[1]) {
    mutations_somatic <- as.character(x_s$Gene.refGene)
    mutations_somatic <- unique(mutations_somatic)

    somatic <- x_s[, c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                       "Variant_Allele_Frequency", "Variant_Reads",
                       "AAChange", "is_tumorsuppressor",
                       "is_oncogene", "is_hotspot", "target", "AF_popmax",
                       "CADD_phred", "condel.label", "REVEL_score", "CLNSIG", "InterVar_automated",
                       "cosmic_coding"),
                  drop = FALSE]

      colnames(somatic) <- col_names
      all_mutations <- rbind(all_mutations, somatic)
  } else {
    mutations_somatic <- c()
  }

  if (!is.null(x_l)){
    mutations_loh <- as.character(x_l$Gene.refGene)
    mutations_loh <- unique(mutations_loh)

    loh <- x_l[, c("Gene.refGene", "GeneName", "ExonicFunc.refGene",
                    "VAF_Tumor", "Count_Tumor", "AAChange",
                    "is_tumorsuppressor", "is_oncogene",
                    "is_hotspot", "target", "AF_popmax", "CADD_phred",
                    "condel.label", "REVEL_score", "CLNSIG", "InterVar_automated", "cosmic_coding"),
                drop = FALSE]

    colnames(loh) <- col_names
    all_mutations <- rbind(all_mutations, loh)
  } else {
    mutations_loh <- c()
  }

  mut <- unique(c(mutations_somatic, mutations_loh))

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
    if (!is.null(outfile)){
      write.xlsx(ds_test, outfile, keepNA = FALSE, rowNames = FALSE,
                firstRow = TRUE)
    }
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
  #' @details mutations found belonging to them.
  #print("PI3K-AKT-mTOR:")
  #tmp <- rownames(ch_mat)[ch_mat[, 3] == 1]
  #print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  #print("RAF-MEK-ERK:")
  #tmp <- rownames(ch_mat)[ch_mat[, 4] == 1]
  #print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  #print("DNA Damage Response:")
  #tmp <- rownames(ch_mat)[ch_mat[, 2] == 1]
  #print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  #print("Cell Cycle:")
  #tmp <- rownames(ch_mat)[ch_mat[, 1] == 1]
  #print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  #print("Tyrosine Kinases:")
  #tmp <- rownames(ch_mat)[ch_mat[, 5] == 1]
  #print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
  #print("TopArt:")
  #tmp <- rownames(ch_mat)[ch_mat[, 6] == 1]
  #print(all_muts[match(tmp, all_muts$Symbol), c("Symbol", "ExonicFunc")])
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
  #print("TopArt:")
  topart <- rownames(ch_mat)[ch_mat[, 6] == 1]
  if (length(topart) != 0){
    topart_genes <- all_muts[match(topart, all_muts$Symbol), ]
    topart_genes$Pathway <- "."
    topart_genes[1, "Pathway"] <- "Topart"
  } else {
    topart_genes <- c()
    }

  important_pathways <- rbind(pi3k_genes, raf_genes, dna_damage_genes,
                              cell_cycle_genes, tyrosine_genes, topart_genes)
  return(important_pws = important_pathways)
}

get_status <- function(table, inf_tab_snv, inf_tab_indel) {
  require(limma)
  
  if(dim(table)[1] > 0) {
    aa <- table$AAChange
    aa <- strsplit(x = aa, split = ";", fixed = TRUE)
    aa <- unlist(lapply(aa, function(x){return(x[1])}))
    id <- grep(pattern = "delins", x = aa)
    
    aa_short = c("H", "Q", "P", "R", "L", "D", "E", "A", "G", "V", "Y", "S", "C", "W", "F", "N", "K", "T", "I", "M", "fs", "X")
    aa_long = c("His", "Gln", "Pro", "Arg", "Leu", "Asp", "Glu", "Ala", "Gly", "Val", "Tyr", "Ser", "Cys", "Trp", "Phe", "Asn", "Lys", "Thr", "Ile", "Met", "fs", "X")
    names(aa_short) <- aa_long
    aa = gsub(aa, pattern = '*',replacement = 'X',fixed = T)
    aa.num = as.numeric(gsub("[^\\d]+", "", aa, perl=TRUE))
    aa = unlist(lapply(strsplit(aa , split = '.', fixed = T), function(s) s[2]))
    aa.split = strsplit(aa, split = "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=T)
    if (length(id) > 0) {
      for (i in 1:length(id)) {
        coord <- paste0(substr(x = aa.split[id][[i]][2], start = 1,
                               stop = nchar(aa.split[id][[i]][2]) - 3),
                        aa.split[id][[i]][3])
        delins <- substr(x = aa.split[id][[i]][4], start = 1, stop = 3)
        aa.split[id[i]] <- paste0(coord, delins)
      }
      ref <- unlist(lapply(aa.split[-id], function(x) {return(x[1])}))
      pos <- unlist(lapply(aa.split[-id], function(x) {return(x[2])}))
      alt <- unlist(lapply(aa.split[-id], function(x) {return(x[3])}))
      ref <- aa_short[ref]
      alt <- aa_short[alt]
      aa.short <- rep(NA, times = length(aa.split))
      aa.short[-id] = paste0(ref, pos, alt)
      aa.short[id] <- aa.split[id]
      aa.short <- unlist(aa.short)
    } else {
      ref <- aa_short[unlist(lapply(aa.split, function(x) {return(x[1])}))]
      pos <- unlist(lapply(aa.split, function(x) {return(x[2])}))
      alt <- aa_short[unlist(lapply(aa.split, function(x) {return(x[3])}))]
      aa.short = unlist(paste0(ref, pos, alt))
    }
    proteinChange = paste0("p.", aa.short)
    
    identifier <- table$Gene.refGene
    inf_tab_indel$ANNOTATION_control <- inf_tab_indel$CLASSIFICATION
    reference <- rbind(inf_tab_snv[, c("GENE", "AA_CHANGE", "ANNOTATION_control")],
                       inf_tab_indel[, c("GENE", "AA_CHANGE", "ANNOTATION_control")])
    id <- grep(pattern = ",", x = reference$GENE, ignore.case = FALSE)
    # Get up to date Symbols instead of Aliases
    reference$GENE_new <- alias2SymbolTable(alias = reference$GENE, species = "Hs")
    if (length(id) > 1) {
      newgenes <- strsplit(x = reference$GENE[id], split = ",", fixed = TRUE)
      newgenes <- lapply(newgenes, function(x){
        return(alias2SymbolTable(alias = x,species = "Hs"))
        })
      reference$GENE_new[id] <- lapply(newgenes, function(x){return(paste0(x, collapse = ";"))})
    }
    identifier_ref <- reference$GENE_new
    
    ids <- list()
    for (i in 1:length(identifier)) {
      ids[i] <- grep(pattern = identifier[i], x = identifier_ref)[1]
    }
    ids <- unlist(ids)
    table$Classification <- NA
    table$Classification <- reference$ANNOTATION_control[ids]
  }
  return(table)
}

msi <- function(msi_file) {
  if (file.exists(msi_file)) {
    msi <- read.table(
      file = msi_file,
      header = TRUE,
      stringsAsFactors = FALSE
    )
    msi_score <- msi[1, 3]
  } else {
    msi_score <- NULL
  }
  return(msi = msi_score)
}