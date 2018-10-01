#' CNV Analysis
#'
#' @description Analyse the
#'
#' @param ratio_file string. Filename for ratio
#' @param cnvs_file string. Filename for CNVs
#' @param cnv_pvalue_txt string. Filename of CNV regions file
#' @param outfile_plot string. Filename of output file
#' @param outfile_ideogram string. Filename of ideogram 
#' @param outfile string. Filename for table of CNV regions
#' @param outfile_onco string. Filename of Oncogene list
#' @param outfile_tumorsuppressors string. Filename of Tumorsuppressor list
#' @param outfile_loss string. Filename of genes with loss
#' @param outfile_gain string. Filename of genes with gain
#' @param outfile_dna_damage string. Filename of Output file
#' @param path_data string. Directory of the databases
#' @param path_script string. Directory of required scripts
#' 
#' @return returns list of
#' @return cnvs_annotated dataframe. Table of annotated CNVs
#' @return cnv_analysis_results dataframe. Results of CNV analysis
#' @return out dataframe. Table of genes with CNVs
#' @return ddr dataframe. Genes with CNVs in the DNA Damage Response pathway
#' 
#' @details Analysis of the genomic regions with CNVs.
cnv_analysis <- function(ratio_file, cnvs_file, cnv_pvalue_txt, outfile_plot,
                         outfile_ideogram, outfile, outfile_onco,
                         outfile_tumorsuppressors, outfile_loss, outfile_gain,
                         outfile_dna_damage, path_data, path_script){

  load(paste(path_data, "GOGeneSets.RData", sep = "/"))
  targets <- read.table(paste(path_data, "Targets.txt", sep = "/"))
  
  assess_significance(ratio_file = ratio_file, cnvs_file = cnvs_file,
                      outfile = cnv_pvalue_txt)
  
  make_cnv_graph(ratio_file = ratio_file_plot, ploidity = "2", outfile_plot = cnv_plot,
                outfile_ideogram = cnv_ideogram_plot)
  
  cnvs_annotated <- cnv_annotation(cnv_pvalue_txt = cnv_pvalue_txt,
                                   outfile = outfile,
                                   outfile_onco = outfile_onco,
                                   outfile_tumorsuppressors =
                                     outfile_tumorsuppressors,
                                   dbfile =  paste(path_data,
                                                   "CancerGenesList.txt",
                                                   sep = "/"),
                                   path_data = path_data, path_script
                                   = path_script)
  
  cnv_analysis_results <- cnv_processing(cnv_file = cnv_annot_out,
                                         targets = targets,
                                         outfile_loss = outfile_loss,
                                         outfile_gain = outfile_gain,
                                         go.bp = go.bp,
                                         path_data = path_data,
                                         path_script = path_script)
  
  out <- cnv_table(cnvs = cnvs_annotated$CNVsAnnotated)
  ddr <- cnv_dna_damage(input = out, outfile_dna_damage = outfile_dna_damage,
                        db = paste(path_data, "DNA_Damage_Response.txt",
                                   sep = "/"))
return(list(cnvs_annotated = cnvs_annotated, cnv_analysis_results
       = cnv_analysis_results, out = out, ddr = ddr))
  }