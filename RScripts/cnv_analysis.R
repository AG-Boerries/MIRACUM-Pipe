cnv_analysis <- function(
  ratio_file,
  cnvs_file,
  cpn_file,
  cnv_pvalue_txt,
  outfile_ideogram,
  path_data,
  path_script,
  targets_txt,
  ampl_genes_txt,
  id,
  outfile_cbioportal,
  outfile_seg,
  protocol,
  sureselect_type,
  purity_file,
  hrd_file,
  gender,
  ucsc_server,
  cnv_region_annotation
) {
#' CNV Analysis
#'
#' @description Analyse the
#'
#' @param ratio_file string. Filename for ratio
#' @param cnvs_file string. Filename for CNVs
#' @param cpn_file string. Filename for CPN
#' @param cnv_pvalue_txt string. Filename of CNV regions file
#' @param outfile_ideogram string. Filename of ideogram
#' @param path_data string. Directory of the databases
#' @param path_script string. Directory of required scripts
#' @param targets_txt string. Path of the targets.txt file provided by the Capture Kit manufracturer
#' @param id string. Sample ID
#' @param outfile_cbioportal string. Name of the outputfile
#' @param outfile_seg string. Name of the outputfile
#' @param protocol string. Type of protocol used for analyses
#'
#' @return returns list of
#' @return cnvs_annotated dataframe. Table of annotated CNVs
#' @return cnv_analysis_results dataframe. Results of CNV analysis
#' @return out dataframe. Table of genes with CNVs
#' @return ddr dataframe. Genes with CNVs in the DNA Damage Response pathway
#'
#' @details Analysis of the genomic regions with CNVs.
  load(paste(path_data, "hallmarksOfCancer_GeneSets.RData", sep = "/"))
  targets <- read.table(targets_txt)

  pvalue_txt <- assess_significance(
    ratio_file = ratio_file,
    cnvs_file = cnvs_file
  )

  make_cnv_ideo_sig(
    ratio_file = pvalue_txt,
    outfile_ideogram = cnv_ideogram_plot
  )

  cnvs_annotated <- cnv_annotation(
    cnv_pvalue_txt = pvalue_txt,
    dbfile =  paste(path_data, "cancerGeneList.tsv", sep = "/"),
    path_data = path_data,
    path_script = path_script,
    ucsc_server = ucsc_server,
    cnv_region_annotation = cnv_region_annotation
  )

  cnv_analysis_results <- cnv_processing(
    cnv_file = cnvs_annotated$CNVsAnnotated,
    targets = targets,
    db = hallmarksOfCancer,
    path_data = path_data,
    path_script = path_script
  )
  if (sureselect_type == "TSO500") {
    ampl_genes <- read.delim(
      file = ampl_genes_txt,
      header = F
    )$V1
  } else {
    ampl_genes <- NULL
  }
  out <- cnv_table(
    cnvs = cnvs_annotated$CNVsAnnotated,
    ampl_genes = ampl_genes
  )

  ddr <- cnv_pathways(
    input = out,
    db = paste(path_data,"DNA_Damage_Response.txt",sep = "/")
  )
  pam <- cnv_pathways(
    input = out,
    db = paste(path_data, "PI3K_AKT_mTOR.txt", sep = "/")
  )
  rme <- cnv_pathways(
    input = out,
    db = paste(path_data, "RAF_MEK_ERK.txt", sep = "/")
  )
  tyk <- cnv_pathways(
    input = out,
    db = paste(path_data, "Tyrosine_Kinases.txt", sep = "/")
  )
  cec <- cnv_pathways(
    input = out,
    db = paste(path_data, "Cell_Cycle.txt", sep = "/")
  )
  impa <- list(
    ddr = ddr,
    pam = pam,
    rme = rme,
    tyk = tyk,
    cec = cec
  )

  if (dim(out)[1] != 0) {
    cnvs2cbioportal(out, id, outfile_cbioportal, gender = gender, ampl_genes = ampl_genes)
    freec2seg(cnvs_file, cpn_file, id, outfile_seg)
  }
  
  print("HRD")
  hrd <- hrd_extr(hrd_file)
  print("Purity")
  pur <- purity_extr(purity_file)
  
  if (sureselect_type == "TSO500") {
    if(dim(out)[1]!=0){
      db <- read.delim(
        paste(
          path_data, "cancerGeneList.tsv", sep = "/"
        ),
        header = T, sep = "\t", colClasses = "character"
      )
      out$is_tumorsuppressor <- 0
      ts <- which(db$Is.Tumor.Suppressor.Gene == "Yes")
      idx <- which(as.character(out$Gene) %in% db$Hugo.Symbol[ts])
      out$is_tumorsuppressor[idx] <- 1
  
      out$is_oncogene <- 0
      og <- which(db$Is.Oncogene == "Yes")
      idx <- which (as.character(out$Gene) %in% db$Hugo.Symbol[og])
      out$is_oncogene[idx] <- 1
      out$Cancergene <- ""
      out$Cancergene[which(out$is_tumorsuppressor == 1)] <- "TSG"
      out$Cancergene[which(out$is_oncogene == 1)] <- paste0(
        "OG", out$Cancergene[which(out$is_oncogene == 1)])
      
      type <- get_type(Oncogenes = out, sureselect_type = sureselect_type)
      
      # merge outputs
      out <- merge(x = out, y = type$gene_loci, by.x = c("Gene", "CopyNumber"),
                   by.y = c("gene_name", "cn"))
      
    } else {
      type <- NULL
    }
    return(list(cnvs_annotated = cnvs_annotated, cnv_analysis_results
                = cnv_analysis_results, out = out, impa = impa,
                gene_loci = type$gene_loci, hrd = hrd, purity = pur))
  }
  
  type <- get_type(
    Oncogenes = cnvs_annotated$CNVOncogenes,
    Tumorsuppressor = cnvs_annotated$CNVTumorSuppressors,
    CNVsAnnotated = cnvs_annotated$CNVsAnnotated,
    sureselect_type = sureselect_type
  )
	
  
  return(
    list(
      cnvs_annotated = cnvs_annotated,
      cnv_analysis_results = cnv_analysis_results,
      out = out,
      impa = impa,
      gene_loci_onc = type$gene_loci_onc,
      gene_loci_tsg = type$gene_loci_tsg,
      hrd = hrd, purity = pur
    )
  )
}
