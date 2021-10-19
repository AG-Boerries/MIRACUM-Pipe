####################
# Variant Analysis #
####################

# Packages
library(OmicCircos)
library(foreach)
library(doMC)
library(rtracklayer)
library(openxlsx)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(RMySQL)
library(biomaRt)
library(RColorBrewer)
library(stringr)
library(Rsamtools)
library(gtrellis)
library(circlize)
library(ComplexHeatmap)
library(YAPSA)
library(SomaticSignatures)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(gdata)
library(gtools)
library(stringi)
library(tidyr)
library(dplyr)
library(magrittr)
library(pracma)
library(MSIseq)

args <- commandArgs()

###########
# General #

sample <- paste(args[6:7], collapse = "_")
protocol <- args[6]
id <- args[7]
germline <- args[8]
tumor <- args[9]
path_output <- paste(args[10], "Analyses/", sep = "/")
path_input <- paste(args[10], "WES/", sep = "/")
path_script <- args[11]
path_data <- args[12]
targets_txt <- args[13]
covered_region <- args[14]
author <- args[15]
center <- args[16]
bed_file <- args[17]
sureselect_type <- args[18]
ref_genome <- args[19]
targetCapture_cor_factors <- args[20]
vaf <- as.numeric(args[21])*100
min_var_count <- as.numeric(args[22])
maf_cutoff <- as.numeric(args[23])
actionable_genes <- args[24]
coveredExons <- args[25]

#############
# Functions #
# Mutation Analyses #
print("Load functions.")
source(paste(path_script, "filtering.R", sep = "/"))
source(paste(path_script, "filtering_tools.R", sep = "/"))
source(paste(path_script, "mutationAnalysis.R", sep = "/"))
source(paste(path_script, "MutAna_tools.R", sep = "/"))
source(paste(path_script, "stats.R", sep = "/"))
source(paste(path_script, "stats_tools.R", sep = "/"))
# Copy Number Analyses #
source(paste(path_script, "cnv_analysis.R", sep = "/"))
source(paste(path_script, "cnv_ana_tools.R", sep = "/"))
# Mutation Signature Analysis #
source(paste(path_script, "mutationSignatureAnalysisFunction.R", sep = "/"))
source(paste(path_script, "Mut_sig_tools.R", sep = "/"))
source(paste(path_script, "mut_sig_ana.R", sep = "/"))

##################
# FASTQC Reports #

if( protocol == "somaticGermline" | protocol == "somatic"){
  tumor <- paste0(path_input, strsplit(x = tumor, split = ".", fixed=T)[[1]][1], "_fastqc/Images/per_base_quality.png")
  germline <- paste0(path_input, strsplit(x = germline, split = ".", fixed=T)[[1]][1], "_fastqc/Images/per_base_quality.png")

  tumor_bsqr <- paste0(path_input, sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png")
  germline_bsqr <- paste0(path_input, sample, "_gd_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png")
}
if (protocol == "panelTumor"){
  tumor <- paste0(path_input, strsplit(x = tumor, split = ".", fixed=T)[[1]][1], "_fastqc/Images/per_base_quality.png")
  tumor_bsqr <- paste0(path_input, sample, "_td_output.sort.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png")
}
if (protocol == "tumorOnly"){
  tumor <- paste0(path_input, strsplit(x = tumor, split = ".", fixed=T)[[1]][1], "_fastqc/Images/per_base_quality.png")
  tumor_bsqr <- paste0(path_input, sample, "_td_output.sort.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png")
}

#####################
# Mutation Analyses
#####################
# DEFINE FILES
print("Preparations.")
if (protocol == "somatic" | protocol == "somaticGermline"){
  if (protocol == "somaticGermline"){
    # GERMLINE NORMAL
    snp_file_gd <- paste0(path_input, sample, "_vc.output.snp.Germline.hc.NORMAL.avinput.hg19_multianno.csv")
    indel_file_gd <- paste0(path_input, sample, "_vc.output.indel.Germline.hc.NORMAL.avinput.hg19_multianno.csv")
    filter_out_gd <- paste0(path_output, sample, "_vc_Germline_NORMAL.xlsx")
    snpefffile_snp_gd <- paste0(path_input, sample, "_vc.output.snp.NORMAL.SnpEff.vcf")
    snpefffile_indel_gd <- paste0(path_input, sample, "_vc.output.indel.NORMAL.SnpEff.vcf")
    outfile_circos_gd <- paste0(path_output, sample, "_GD_circos.pdf")
  }
  # SOMATIC TUMOR
  snp_file_td <- paste0(path_input, sample, "_vc.output.snp.Somatic.hc.TUMOR.avinput.hg19_multianno.csv")
  indel_file_td <- paste0(path_input, sample, "_vc.output.indel.Somatic.hc.TUMOR.avinput.hg19_multianno.csv")
  snpefffile_snp <- paste0(path_input, sample, "_vc.output.snp.Somatic.SnpEff.vcf")
  snpefffile_indel <- paste0(path_input, sample, "_vc.output.indel.Somatic.SnpEff.vcf")
  filter_out_td <- paste0(path_output, sample, "_VC_Somatic_TUMOR.xlsx")
  # LOH
  snp_file_loh <- paste0(path_input, sample, "_vc.output.snp.LOH.hc.avinput.hg19_multianno.csv")
  indel_file_loh <- paste0(path_input, sample, "_vc.output.indel.LOH.hc.avinput.hg19_multianno.csv")
  loh_out <- paste0(path_output, sample, "_vc_LOH.xlsx")
  snpefffile_snp_loh <- paste0(path_input, sample, "_vc.output.snp.LOH.SnpEff.vcf")
  snpefffile_indel_loh <- paste0(path_input, sample, "_vc.output.indel.LOH.SnpEff.vcf")
  # Results
  outfile_circos <- paste0(path_output, sample, "_TD_circos.pdf")
  coverage_out <- paste0(path_output, sample, "_coverage.pdf")
  # MAFs
  maf_gd <- paste0(path_output, sample, "_Germline.maf")
  maf_td <- paste0(path_output, sample, "_Somatic.maf")
  maf_loh <- paste0(path_output, sample, "_LoH.maf")
  maf_complete <- paste0(path_output, sample, ".maf")
}
if (protocol == "panelTumor" | protocol == "tumorOnly"){
  # TUMOR
  snp_file_td <- paste0(path_input, sample, "_vc.output.snp.Sample1.avinput.hg19_multianno.csv")
  indel_file_td <- paste0(path_input, sample, "_vc.output.indel.Sample1.avinput.hg19_multianno.csv")
  snpefffile_snp_td <- paste0(path_input, sample, "_vc.output.snp.SnpEff.vcf")
  snpefffile_indel_td <- paste0(path_input, sample, "_vc.output.indel.SnpEff.vcf")
  filter_out_td <- paste0(path_output, sample, "_VC_TUMOR.xlsx")
  # Results
  outfile_circos <- paste0(path_output, sample, "_TD_circos.pdf")
  coverage_out <- paste0(path_output, sample, "_coverage.pdf")
  # MAFs
  maf_td <- paste0(path_output, sample, "_TUMOR.maf")
  maf_complete <- paste0(path_output, sample, ".maf")
}

if (protocol == "somatic" | protocol == "somaticGermline"){
  # SOMATIC TUMOR
  print("Filtering for Tumor.")
  filt_result_td <- filtering(snpfile = snp_file_td, indelfile = indel_file_td,
                              snpefffile_snp = snpefffile_snp,
                              snpefffile_indel = snpefffile_indel,
                              outfile = filter_out_td, outfile_maf = maf_td,
                              path_data = path_data,
                              path_script = path_script, covered_region = covered_region,
                              mode ="T", center = center, id = id,
                              protocol = protocol, sureselect = bed_file,
                              vaf = vaf, min_var_count = min_var_count,
                              maf = maf_cutoff,
                              coveredExons = coveredExons, cov_t = 1,
                              sureselect_type = sureselect_type)

  if (protocol == "somaticGermline"){
    # GERMLINE NORMAL
    print("Filtering for Germline.")
    filt_result_gd <- filtering(snpfile = snp_file_gd, indelfile = indel_file_gd,
                                snpefffile_snp = snpefffile_snp_gd,
                                snpefffile_indel = snpefffile_indel_gd,
                                outfile =  filter_out_gd, outfile_maf = maf_gd,
                                path_data = path_data,
                                path_script = path_script, covered_region = NULL,
                                mode = "N", center = center, id = id,
                                protocol = protocol, sureselect = bed_file,
                                vaf = vaf, min_var_count = min_var_count,
                                maf = maf_cutoff, actionable_genes = actionable_genes,
                                coveredExons = coveredExons, cov_t = 1,
                                sureselect_type = sureselect_type)
    
    mutation_analysis_result_gd <- mutation_analysis(loh = NULL,
                                                     somatic = filt_result_gd$table,
                                                     tumbu = NULL,
                                                     outfile_circos = outfile_circos_gd,
                                                     path_data = path_data,
                                                     path_script = path_script,
                                                     targets_txt = targets_txt,
                                                     protocol = "panelTumor",
                                                     sureselect = bed_file,
                                                     sureselect_type = sureselect_type)
  }
  
  # LOH
  print("Filtering for LoH.")
  filt_result_loh <- filtering(snpfile = snp_file_loh, indelfile = indel_file_loh,
                               snpefffile_snp = snpefffile_snp_loh,
                               snpefffile_indel = snpefffile_indel_loh,
                               outfile = loh_out, outfile_maf = maf_loh,
                               path_data = path_data, path_script = path_script,
                               covered_region = NULL, mode = "LOH", center = center, id = id,
                               protocol = protocol, sureselect = bed_file,
                               vaf = vaf, min_var_count = min_var_count,
                               maf = maf_cutoff,
                               coveredExons = coveredExons, cov_t = 1,
                               sureselect_type = sureselect_type)
  # Analyses
  print("Variant Analyses.")
  mutation_analysis_result <- mutation_analysis(loh = filt_result_loh$table,
                                                somatic = filt_result_td$table,
                                                tumbu = filt_result_td$tmb,
                                                outfile_circos = outfile_circos,
                                                path_data = path_data,
                                                path_script = path_script,
                                                targets_txt = targets_txt,
                                                protocol = protocol, sureselect = bed_file,
                                                sureselect_type = sureselect_type)
}
if (protocol == "panelTumor" | protocol == "tumorOnly"){
  print("Filtering for Tumor.")
  filt_result_td <- filtering(snpfile = snp_file_td, indelfile = indel_file_td,
                              snpefffile_snp = snpefffile_snp_td,
                              snpefffile_indel = snpefffile_indel_td,
                              outfile = filter_out_td, outfile_maf = maf_td,
                              path_data = path_data,
                              path_script = path_script, covered_region = covered_region,
                              mode ="T", center = center, id = id,
                              protocol = protocol, sureselect = bed_file,
                              vaf = vaf, min_var_count = min_var_count,
                              maf = maf_cutoff,
                              coveredExons = coveredExons, cov_t = 1,
                              sureselect_type = sureselect_type)

  filt_result_loh <- list(table = NULL, tmb = NULL)

  # Analyses
  print("Variant Analyses.")
  mutation_analysis_result <- mutation_analysis(loh = filt_result_loh$table,
                                                somatic = filt_result_td$table,
                                                tumbu = filt_result_td$tmb,
                                                outfile_circos = outfile_circos,
                                                path_data = path_data,
                                                path_script = path_script,
                                                targets_txt = targets_txt,
                                                protocol = protocol, sureselect = bed_file,
                                                sureselect_type = sureselect_type)
}
  

# Combine MAF files to obtain one complete maf per patient
if (protocol == "somaticGermline" | protocol == "somatic"){
  if (protocol == "somaticGermline"){
      maf_comb <- smartbind(filt_result_td$maf,filt_result_gd$maf, filt_result_loh$maf)
      write.table(x = maf_comb, file = maf_complete , append = F, quote = F, sep = '\t', col.names = T, row.names = F)
  } else {
    maf_comb <- smartbind(filt_result_td$maf, filt_result_loh$maf)
    write.table(x = maf_comb, file = maf_complete , append = F, quote = F, sep = '\t', col.names = T, row.names = F)
  }
}
#} else {
#  maf_comb <- maf_td
#  write.table(x = maf_comb, file = maf_complete , append = F, quote = F, sep = '\t', col.names = T, row.names = F)
#}

###
# STATISTICS
## Input Files
print("Statistics.")
if (protocol == "somaticGermline" | protocol == "somatic"){
  stats_td <- paste0(path_input, sample, "_td_stats.txt")
  stats_gd <- paste0(path_input, sample, "_gd_stats.txt")
  ## Analysis
  stats <- stats(path = path_input, outfile_pdf = coverage_out,
                 stats_td = stats_td, stats_gd = stats_gd, protocol = protocol)
}
if (protocol == "panelTumor" | protocol == "tumorOnly"){
  stats_td <- paste0(path_input, sample, "_td_stats.txt")
  stats <- stats(path = path_input, outfile_pdf = coverage_out,
                 stats_td = stats_td, protocol = protocol)
}

########################
# Copy Number Analysis #
print("CNV Analyses.")
if (protocol == "somaticGermline" | protocol == "somatic"){
  ## Input/Output Files
  ratio_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_ratio.txt")
  cnvs_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_CNVs")
  cpn_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_sample.cpn")

  # Results
  cnv_pvalue_txt <- paste0(cnvs_file, ".p.value.txt")
  # cnv_plot <- paste0(path_output, sample, "_CNV_Plot.pdf")
  cnv_ideogram_plot <- paste0(path_output, sample,"_CNV_Plot_Ideogram.pdf")
  # cnv_annot_out <- paste0(path_output, sample, "_CNV.xlsx")
  # outfile_onco <- paste0(path_output, sample, "_CNV_Oncogene.xlsx")
  # outfile_ts <- paste0(path_output, sample, "_CNV_TumorSuppressors.xlsx")
  # outfile_loss <- paste0(path_output, sample, "_CNV_loss_GO.xlsx")
  # outfile_gain <- paste0(path_output, sample, "_CNV_gain3_GO.xlsx")
  # outfile_dna_damage <- paste0(path_output, sample, "_CNV_dna_damage.xlsx")
  outfile_cnvs_cbioportal <- paste0(path_output, sample, "_CNV_cbioportal.txt")
  outfile_cnvs_seg <- paste0(path_output, sample, "_CNV.seg")

  #   cnv_analysis_results <- cnv_analysis(ratio_file = ratio_file, cnvs_file = cnvs_file, cnv_pvalue_txt = cnv_pvalue_txt,
  #                                      outfile_plot = cnv_plot,
  #                                      outfile_ideogram = cnv_ideogram_plot,
  #                                      outfile = cnv_annot_out,
  #                                      outfile_onco = outfile_onco,
  #                                      outfile_tumorsuppressors = outfile_ts,
  #                                      outfile_loss = outfile_loss,
  #                                      outfile_gain = outfile_gain,
  #                                      outfile_dna_damage = outfile_dna_damage,
  #                                      path_data = path_data,
  #                                      path_script = path_script,
  #                                      targets_txt = targets_txt,
  #                                      outfile_cbioportal = outfile_cnvs_cbioportal,
  #                                      id = id,
  #                                      protocol = protocol)

  cnv_analysis_results <- cnv_analysis(ratio_file = ratio_file,
                                       cnvs_file = cnvs_file,
                                       cpn_file = cpn_file,
                                       cnv_pvalue_txt = cnv_pvalue_txt,
                                       outfile_ideogram = cnv_ideogram_plot,
                                       path_data = path_data,
                                       path_script = path_script,
                                       targets_txt = targets_txt,
                                       outfile_cbioportal = outfile_cnvs_cbioportal,
                                       outfile_seg = outfile_cnvs_seg,
                                       id = id,
                                       protocol = protocol)
}

if (protocol == "tumorOnly"){
  ## Input/Output Files
  ratio_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.rmdup.realigned.fixed.recal.bam_ratio.txt")
  cnvs_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.rmdup.realigned.fixed.recal.bam_CNVs")
  cpn_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.rmdup.realigned.fixed.recal.bam_sample.cpn")

  # Results
  cnv_pvalue_txt <- paste0(cnvs_file, ".p.value.txt")
  cnv_ideogram_plot <- paste0(path_output, sample,"_CNV_Plot_Ideogram.pdf")
  outfile_cnvs_cbioportal <- paste0(path_output, sample, "_CNV_cbioportal.txt")
  outfile_cnvs_seg <- paste0(path_output, sample, "_CNV.seg")

  cnv_analysis_results <- cnv_analysis(ratio_file = ratio_file,
                                       cnvs_file = cnvs_file,
                                       cpn_file = cpn_file,
                                       cnv_pvalue_txt = cnv_pvalue_txt,
                                       outfile_ideogram = cnv_ideogram_plot,
                                       path_data = path_data,
                                       path_script = path_script,
                                       targets_txt = targets_txt,
                                       outfile_cbioportal = outfile_cnvs_cbioportal,
                                       outfile_seg = outfile_cnvs_seg,
                                       id = id,
                                       protocol = protocol)
}
# # TODO CNV Panel
# if (protocol == "panelTumor"){
#   ## Input Files
#   cnvsFile <- paste0(path_input,"CNV/", sample, "_td_output.sort.realigned.fixed.cnr")
#   ## Output Files
#   outfile <- paste0(path_output, sample, "_CNV.xlsx")
#   outfile_ts_og <- paste0(path_output,sample, "_CNV_TSG_OG.xlsx")
#   outfile_ideogram <- paste0(path_output, sample, "_ideogram.pdf")
#   cnvs <- cnv_panel(input_file = cnvsFile, outfile = outfile,
#                     outfile_ts_og = outfile_ts_og, outfile_ideogram = outfile_ideogram,
#                     path_data = path_data, sureselect = sureselect, targets_txt = targets_txt, protocol = protocol)
# }


###############################
# Mutation Signature Analysis #
print("Mutation Signature Analysis.")
if( protocol == "somaticGermline" | protocol == "somatic"){
  somaticVCF <- paste0(path_input, sample,"_vc.output.snp.Somatic.hc.fpfilter.vcf")
  outfile_mutsig_cbioportal <- paste0(path_output, sample, "_mutsig_cbioportal")
  mut_sig_ana <- mut_sig_wCI(vcf_file = somaticVCF, cutoff = 0.01, sample = sample, sureselect_type = sureselect_type, path_script = path_script, ref_genome = ref_genome, targetCapture_cor_factors, path_output = path_output, sample_name = id, outfile_cbioportal = outfile_mutsig_cbioportal)
} else {
  vcf <- paste0(path_input, sample,"_vc.output.snp.fpfilter.vcf")
  outfile_mutsig_cbioportal <- paste0(path_output, sample, "_mutsig_cbioportal")
  mut_sig_analysis <- mutation_signature_analysis(vcf_file = vcf,
                                                  cutoff = 0.01,
                                                  sample_name = NULL,
                                                  only_coding = FALSE,
                                                  path_data = path_data,
                                                  path_output = path_output,
                                                  outfile_cbioportal = outfile_mutsig_cbioportal)
}

# Write Excel File
if (protocol == "somaticGermline"){
  output <- list(Somatic_Mutations = filt_result_td$table,
                 LoH_Mutations = filt_result_loh$table,
                 Germline_Mutations = filt_result_gd$table,
                 Mutation_Signatures = mut_sig_ana$output$Summary,
                 CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated)
  write.xlsx(x = output, file = paste0(path_output, "/", sample, "_results.xlsx"))
} else if (protocol == "somatic"){
  output <- list(Somatic_Mutations = filt_result_td$table,
                 LoH_Mutations = filt_result_loh$table,
                 Mutation_Signatures = mut_sig_ana$output$Summary,
                 CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated)
  write.xlsx(x = output, file = paste0(path_output, "/", sample, "_results.xlsx"))
  } else if (protocol == "panelTumor"){
    output <- list(Somatic_Mutations = filt_result_td$table)
    write.xlsx(x = output, file = paste0(path_output, "/", sample, "_results.xlsx"))
} else {
  output <- list(Somatic_Mutations = filt_result_td$table,
                 CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated)
  write.xlsx(x = output, file = paste0(path_output, "/", sample, "_results.xlsx"))
}

save.image(file = "MTB.RData")
