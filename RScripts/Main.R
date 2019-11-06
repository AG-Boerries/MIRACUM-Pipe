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

args <- commandArgs()

###########
# General #

sample <- paste(args[6:7], collapse = "_")
tumor <- args[8]
germline <- args[9]
path_output <- paste(args[10], "Analysis/", sep = "/")
path_input <- paste(args[10], "WES/", sep = "/")
path_script <- args[11]
path_data <- args[12]
targets_txt <- args[13]
covered_region <- args[14]
author <- args[15]

#############
# Functions #
# Mutation Analysis #
source(paste(path_script, "filtering.R", sep = "/"))
source(paste(path_script, "filtering_tools.R", sep = "/"))
source(paste(path_script, "mutationAnalysis.R", sep = "/"))
source(paste(path_script, "MutAna_tools.R", sep = "/"))
source(paste(path_script, "stats.R", sep = "/"))
source(paste(path_script, "stats_tools.R", sep = "/"))
# Copy Number Analysis #
source(paste(path_script, "cnv_analysis.R", sep = "/"))
source(paste(path_script, "cnv_ana_tools.R", sep = "/"))
# Mutation Signature Analysis #
source(paste(path_script, "mutationSignatureAnalysisFunction.R", sep = "/"))

##################
# FASTQC Reports #

tumor <- paste0(path_input, strsplit(x = tumor, split = ".", fixed=T)[[1]][1], "_fastqc/Images/per_base_quality.png")
germline <- paste0(path_input, strsplit(x = germline, split = ".", fixed=T)[[1]][1], "_fastqc/Images/per_base_quality.png")

tumor_bsqr <- paste0(path_input, sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png")
germline_bsqr <- paste0(path_input, sample, "_gd_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png")

#####################
# Mutation Analysis
#####################
# FILES
if (args[6] == "somaticGermline"){
  # GERMLINE NORMAL
  snp_file_gd <- paste0(path_input, sample, "_vc.output.snp.Germline.hc.NORMAL.avinput.hg19_multianno.csv")
  indel_file_gd <- paste0(path_input, sample, "_vc.output.indel.Germline.hc.NORMAL.avinput.hg19_multianno.csv")
  filter_out_gd <- paste0(sample, "_VC_Germline_NORMAL.xlsx")
  snpefffile_snp_gd <- paste0(path_input, sample, "_vc.output.snp.NORMAL.SnpEff.vcf")
  snpefffile_indel_gd <- paste0(path_input, sample, "_vc.output.indel.NORMAL.SnpEff.vcf")
}
# SOMATIC TUMOR
snp_file_td <- paste0(path_input, sample, "_vc.output.snp.Somatic.hc.TUMOR.avinput.hg19_multianno.csv")
indel_file_td <- paste0(path_input, sample, "_vc.output.indel.Somatic.hc.TUMOR.avinput.hg19_multianno.csv")
snpefffile_snp <- paste0(path_input, sample, "_vc.output.snp.Somatic.SnpEff.vcf")
snpefffile_indel <- paste0(path_input, sample,
                          "_vc.output.indel.Somatic.SnpEff.vcf")
filter_out_td <- paste0(path_output, sample, "_VC_Somatic_TUMOR.xlsx")
# LOH
snp_file_loh <- paste0(path_input, sample, "_vc.output.snp.LOH.hc.avinput.hg19_multianno.csv")
indel_file_loh <- paste0(path_input, sample, "_vc.output.indel.LOH.hc.avinput.hg19_multianno.csv")
loh_out <- paste0(path_output, sample, "_VC_LOH.xlsx")
snpefffile_snp_loh <- paste0(path_input, sample, "_vc.output.snp.LOH.SnpEff.vcf")
snpefffile_indel_loh <- paste0(path_input, sample, "_vc.output.indel.LOH.SnpEff.vcf")
# Results
outfile_circos <- paste0(path_output, sample, "_TD_circos.pdf")
outfile_go <- paste0(path_output, sample, "_TD_hyperGTest_GO.xlsx")
outfile_reactome <- paste0(path_output, sample, "_TD_hyperGTest_Reactome.xlsx")
outfile_consensus <- paste0(path_output, sample, "_TD_hyperGTest_Consensus.xlsx")
outfile_hallmarks <- paste0(path_output, sample, "_TD_hyperGTest_Hallmarks.xlsx")
outfile_mtb_genesets <- paste0(path_output, sample, "_TD_Genesets.xlsx")
coverage_out <- paste0(path_output, sample, "_coverage.pdf")

if (args[6] == "somaticGermline"){
  # GERMLINE NORMAL
  filt_result_gd <- filtering(snpfile = snp_file_gd, indelfile =  indel_file_gd,
                            snpefffile_snp = snpefffile_snp_gd,
                            snpefffile_indel = snpefffile_indel_gd,
                            outfile =  filter_out_gd, path_data = path_data,
                            path_script = path_script, covered_region = NULL,
                            mode = "N")
}
# SOMATIC TUMOR
filt_result_td <- filtering(snpfile = snp_file_td, indelfile = indel_file_td,
                          snpefffile_snp = snpefffile_snp,
                          snpefffile_indel = snpefffile_indel,
                          outfile = filter_out_td, path_data = path_data,
                          path_script = path_script, covered_region = covered_region,
                          mode ="T")
# LOH
filt_result_loh <- filtering(snpfile = snp_file_loh, indelfile = indel_file_loh,
                           snpefffile_snp = snpefffile_snp_loh,
                           snpefffile_indel = snpefffile_indel_loh,
                           outfile = loh_out,
                           path_data = path_data, path_script = path_script,
                           covered_region = NULL, mode = "LOH")
# Analysis
mutation_analysis_result <- mutation_analysis(loh = filt_result_loh$table,
                                              somatic = filt_result_td$table,
                                              tumbu = filt_result_td$tmb,
                                              outfile_circos = outfile_circos,
                                              outfile_go = outfile_go,
                                              outfile_reactome = outfile_reactome,
                                              outfile_consensus = outfile_consensus,
                                              outfile_hallmarks = outfile_hallmarks,
                                              outfile_mtb_genesets =
                                              outfile_mtb_genesets,
                                              path_data = path_data,
                                              path_script = path_script,
                                              targets_txt = targets_txt)

###
# STATISTICS
## Input Files
stats_td <- paste0(path_input, sample, "_gd_stats.txt")
stats_gd <- paste0(path_input, sample, "_td_stats.txt")
## Analysis
stats <- stats(path = path_input, outfile_pdf = coverage_out,
               stats_td = stats_td, stats_gd = stats_gd)

########################
# Copy Number Analysis #

## Input/Output Files
ratio_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_ratio.txt")
cnvs_file <- paste0(path_input, "CNV/", sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_CNVs")

# Results
cnv_pvalue_txt <- paste(path_output, sample, "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_CNVs",
                      ".p.value.txt", sep = "")
cnv_plot <- paste0(path_output, sample, "_CNV_Plot.pdf")
cnv_ideogram_plot <- paste0(path_output, sample,"_CNV_Plot_Ideogram.pdf")
cnv_annot_out <- paste0(path_output, sample, "_CNV.xlsx")
outfile_onco <- paste0(path_output, sample, "_CNV_Oncogene.xlsx")
outfile_ts <- paste0(path_output, sample, "_CNV_TumorSuppressors.xlsx")
outfile_loss <- paste0(path_output, sample, "_CNV_loss_GO.xlsx")
outfile_gain <- paste0(path_output, sample, "_CNV_gain3_GO.xlsx")
outfile_dna_damage <- paste0(path_output, sample, "_CNV_dna_damage.xlsx")

cnv_analysis_results <- cnv_analysis(ratio_file = ratio_file, cnvs_file
                                     = cnvs_file, cnv_pvalue_txt = cnv_pvalue_txt,
                                     outfile_plot = cnv_plot,
                                     outfile_ideogram = cnv_ideogram_plot,
                                     outfile = cnv_annot_out,
                                     outfile_onco = outfile_onco,
                                     outfile_tumorsuppressors = outfile_ts,
                                     outfile_loss = outfile_loss,
                                     outfile_gain = outfile_gain,
                                     outfile_dna_damage = outfile_dna_damage,
                                     path_data = path_data,
                                     path_script = path_script,
                                     targets_txt = targets_txt)
###############################
# Mutation Signature Analysis #

somaticVCF <- paste0(path_input, sample,
                     "_vc.output.snp.Somatic.hc.fpfilter.vcf")

mut_sig_analysis <- mutation_signature_analysis(vcf_file = somaticVCF,
                                            cutoff = 0.01,
                                            sample_name = NULL,
                                            only_coding = FALSE,
                                            path_data = path_data,
                                            path_output = path_output)

save.image(file = "WES.RData")
