#####################
# Mutation Analysis #


mutation_analysis <- function(
  loh,
  somatic,
  tumbu,
  outfile_circos,
  path_data,
  path_script,
  path_output,
  targets_txt,
  protocol,
  sureselect,
  sureselect_type
) {
  #' Mutation Analysis
  #'
  #' @description Analysis and Annotation of Mutations
  #'
  #' @param loh dataframe. Table of LoH mutations
  #' @param somatic dataframe. Table of somatic mutations
  #' @param tumbu numerical. Tumor mutational burden
  #' @param outfile_circos string. Name of output file
  #' @param outfile_go string. Name of output file
  #' @param outfile_reactome string. Name of output file
  #' @param outfile_consensus string. Name of output file
  #' @param outfile_hallmarks string. Name of output file
  #' @param outfile_mtb_genesets string. Name of output file
  #' @param path_data string. Directory of the databases
  #' @param path_script string. Directory of required scripts
  #' @param path_output string. Output directory
  #' @param targets_txt string. Path of the targets.txt file provided by the Capture Kit manufracturer
  #' @param protocol string. Type of protocol used for analyses
  #' @param sureselect string. Path to capture region bed file
  #'
  #' @return returns list of
  #' @return ts_og dataframe. Table of mutations in tumorsuppressor and
  #' @return oncogenes
  #' @return go dataframe. Table of GO pathways
  #' @return reactome dataframe. Table of Reactome pathways
  #' @return consensus dataframe. Table of Consensus pathways
  #' @return hallmarks dataframe. Table of Hallmarks of Cancer pathways
  #' @return mtb_genesets dataframe. Table of genes in important pathways
  #' @return tot_mut numerical. Total number of mutations
  #' @return som_mut numerical. Number of all somatic mutations
  #' @return loh_mut numerical. Number of all LoH mutations
  #' @return important_pathways dataframe. Table of mutated genes in important
  #' @return pathways
  #' @return som_mut_tab dataframe. Table of somatic mutations
  #' @return table_loh_mutations dataframe. Table of LoH mutations
  #' @return all_mutations dataframe. Table of all mutations
  #'
  #' @note The following files are produced by mutation_analysis:
  #' @note - "MutationTable.txt"
  #' @note - "mutationStats.txt"
  #' @note - "Tumorsuppressor-OncogeneTable.xlsx"
  #' @note - "somaticMutations.xlsx"
  #' @note - "lohMutations.xlsx"
  #' @note - Circosplot (pdf file)
  #' @note - "All_Mutations_Somatic.xlsx"
  #' @note - xlsx file for each of the four pathway datasets
  #' @note - xlsx file with the predefined important pathways
  #'
  #' @details The given mutations are analysed and an overview is generated.
  #' @details The mutations are classified and written in several tables
  #' @details including a table for cancer genes. A functional pathway analysis
  #' @details performed to get an overview of the mutation's impact.
  require(org.Hs.eg.db)
  require(openxlsx)
  require(Homo.sapiens)
  require(stringr)

  load(paste(path_data, "hallmarksOfCancer_GeneSets.RData", sep = "/"))
  load(paste(path_data, "MTB_Genesets.rda", sep = "/"))
  targets <- read.table(targets_txt)
  
  source(paste(path_script, "MutAna_tools.R", sep = "/"))

  if (protocol == "somaticGermline" | protocol == "somatic"){
    x_loh <- loh
  } else {
    x_loh <- data.frame()
  }
  x_somatic <- somatic
# Take only exonic mutations
  x_somatic$Func.refGene <- as.character(x_somatic$Func.refGene)
  id_ex <- which(x_somatic$Func.refGene == "exonic")
  x_s <- x_somatic[id_ex, ]
  no_loh <- FALSE
  if (is.null(x_loh) | protocol == "panelTumor" | protocol == "tumorOnly"){
    no_loh <- TRUE
  }
  if (!no_loh){
    x_loh$Func.refGene <- as.character(x_loh$Func.refGene)
    id_ex <- which(x_loh$Func.refGene == "exonic")
    if (length(id_ex) == 0) {
      x_l <- NULL
      no_loh <- TRUE
    } else {
      x_l <- x_loh[id_ex, ]
    }
  }
# split lists into sublists
  sub_lst <- div(x_s, x_l, no_loh)
# Write/print summary of mutationsnumber, diverse tables
  mutation_table <- mut_tab(sub_lst$x_s_snp, sub_lst$x_s_indel, sub_lst$x_l_snp,
                            sub_lst$x_l_indel, protocol)
  if (!no_loh) {
    mut_stats_res <- mut_stats(x_s = x_s, x_l = x_l, tumbu = tumbu, protocol = protocol)
    tbl <- tables(x_s = x_s, x_l = x_l, protocol = protocol)
    all_mut <- write_all_mut(x_s = x_s, x_l = x_l)
  } else {
    mut_stats_res <- mut_stats(x_s = x_s, tumbu = tumbu, protocol = protocol)
    tbl <- tables(x_s = x_s, protocol = protocol)
    all_mut <- write_all_mut(x_s = x_s)
  }

# produce circos plot
  if(sub_lst$no_snp & sub_lst$no_loh & sub_lst$no_indel_somatic) {
    cat("No mutations to draw in CircosPlot!")
    result_go <- NULL
    result_re <- NULL
    result_co <- NULL
    result_hm <- NULL
    
    # Check for important pathways
    check_mat <- NULL
    importantpws <- NULL
  
  } else {
    cc <- circos_colors(x_s_snp = sub_lst$x_s_snp, x_s_indel = sub_lst$x_s_indel,
                        x_l_snp = sub_lst$x_l_snp,
                        x_l_indel = sub_lst$x_l_indel, no_loh = sub_lst$no_loh,
                        no_indel_somatic = sub_lst$no_indel_somatic,
                        no_snp = sub_lst$no_snp,
                        no_indel_loh = sub_lst$no_indel_loh,
                        no_snp_loh = sub_lst$no_snp_loh)
    omicCircosUni(listOfMap = as.matrix(cc$map_mat), label = NULL, 125,
                  outfile_circos, circosColors = as.vector(cc$circoscolors),
                  protocol = protocol, sureselect = sureselect, sureselect_type)

    if (protocol == "somaticGermline" | protocol == "somatic"){
      # Pathway-Analysis
      prep <- prep_pwa(targets, all_mut$mut)
      result_hm <- get_terms(dataset = hallmarksOfCancer, mut.entrez = prep$de_genes, t2.entrez = prep$universe, outfile = NULL)
    } else {
      result_go <- NULL
      result_re <- NULL
      result_co <- NULL
      result_hm <- NULL
    }

  # Check for important pathways
    check_mat <- write_mtb_genesets(all_mut$mut, mtb.genesets)
    importantpws <- imp_pws(check_mat, all_mut$all_muts)
  }
  msi <- msi(msi_file)
  return(list(mut_tab = mutation_table,
              ts_og = tbl$ts_og_table,
              go = NULL, reactome = NULL,
              consensus = NULL, hallmarks = result_hm,
              mtb_genesets = check_mat,
              tot_mut = mut_stats_res$tot_mut,
              som_mut = mut_stats_res$som_mut,
              loh_mut = mut_stats_res$loh_mut,
              important_pathways = importantpws,
              som_mut_tab = tbl$sm_table,
              table_loh_mutations = tbl$lm_table,
              all_mutations = all_mut$all_muts,
              studies = studies
              msi = msi))
}
