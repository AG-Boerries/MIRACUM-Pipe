#####################
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
library(scarHRD)

args <- commandArgs()

#############
# Parameter #

sample <- paste(args[6:7], collapse = "_")
protocol <- args[6]
id <- args[7]
germline <- args[8]
tumor <- args[9]
base_path <- args[10]
path_output <- paste(base_path, "Analyses/", sep = "/")
path_input <- paste(base_path, "WES/", sep = "/")
path_script <- args[11]
path_data <- args[12]
targets_txt <- args[13]
#covered_region <- args[14]
author <- args[15]
center <- args[16]
bed_file <- args[17]
sureselect_type <- args[18]
ref_genome <- args[19]
target_capture_cor_factors <- args[20]
vaf <- as.numeric(args[21])*100
min_var_count <- as.numeric(args[22])
maf_cutoff <- as.numeric(args[23])
actionable_genes <- args[24]
covered_exons <- args[25]
entity <- as.character(args[26])
gender <- as.character(args[27])
fusion_genes <- args[28]
ampl_genes_txt <- args[29]

print(ref_genome)

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

if (protocol == "somaticGermline" | protocol == "somatic") {
  tumor <- paste0(
    path_input, strsplit(x = tumor, split = ".", fixed = T)[[1]][1],
    "_fastqc/Images/per_base_quality.png"
  )
  germline <- paste0(
    path_input, strsplit(x = germline, split = ".", fixed = T)[[1]][1],
    "_fastqc/Images/per_base_quality.png"
  )
  tumor_bsqr <- paste0(
    path_input, sample,
    "_td_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png"
  )
  germline_bsqr <- paste0(
    path_input, sample,
    "_gd_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png"
  )
}
if (protocol == "panelTumor") {
  tumor <- paste0(
    path_input, strsplit(x = tumor, split = ".", fixed = T)[[1]][1],
    "_fastqc/Images/per_base_quality.png"
  )
  tumor_bsqr <- paste0(
    path_input, sample,
    "_td_output.sort.filtered.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png"
  )
}
if (protocol == "tumorOnly") {
  tumor <- paste0(
    path_input, strsplit(x = tumor, split = ".", fixed = T)[[1]][1],
    "_fastqc/Images/per_base_quality.png"
  )
  tumor_bsqr <- paste0(
    path_input, sample,
    "_td_output.sort.rmdup.realigned.fixed.recal_fastqc/Images/per_base_quality.png"
  )
}

# DEFINE FILES
print("Preparations.")
if (protocol == "somatic" | protocol == "somaticGermline") {
  if (protocol == "somaticGermline") {
    # GERMLINE NORMAL
    snp_file_gd <- paste0(
      path_input, sample,
      "_vc.output.snp.Germline.hc.NORMAL.avinput.hg19_multianno.csv"
    )
    indel_file_gd <- paste0(
      path_input, sample,
      "_vc.output.indel.Germline.hc.NORMAL.avinput.hg19_multianno.csv"
    )
    filter_out_gd <- paste0(
      path_output, sample,
      "_vc_Germline_NORMAL.xlsx"
    )
    snpefffile_snp_gd <- paste0(
      path_input, sample,
      "_vc.output.snp.NORMAL.SnpEff.vcf"
    )
    snpefffile_indel_gd <- paste0(
      path_input, sample,
      "_vc.output.indel.NORMAL.SnpEff.vcf"
    )
    outfile_circos_gd <- paste0(
      path_output, sample,
     "_GD_circos.pdf"
    )
  }
  # SOMATIC TUMOR
  snp_file_td <- paste0(
    path_input, sample,
    "_vc.output.snp.Somatic.hc.TUMOR.avinput.hg19_multianno.csv"
  )
  indel_file_td <- paste0(
    path_input, sample,
    "_vc.output.indel.Somatic.hc.TUMOR.avinput.hg19_multianno.csv"
  )
  snpefffile_snp <- paste0(
    path_input, sample,
    "_vc.output.snp.Somatic.SnpEff.vcf"
  )
  snpefffile_indel <- paste0(
    path_input, sample,
    "_vc.output.indel.Somatic.SnpEff.vcf"
  )
  filter_out_td <- paste0(
    path_output, sample,
    "_VC_Somatic_TUMOR.xlsx"
  )
  msi_file <- paste0(
     path_input, "/", sample,
     "_vc_MSI"
  )
  # LOH
  snp_file_loh <- paste0(
    path_input, sample,
    "_vc.output.snp.LOH.hc.avinput.hg19_multianno.csv"
  )
  indel_file_loh <- paste0(path_input, sample,
  "_vc.output.indel.LOH.hc.avinput.hg19_multianno.csv")
  loh_out <- paste0(
    path_output, sample,
    "_vc_LOH.xlsx"
  )
  snpefffile_snp_loh <- paste0(
    path_input, sample,
    "_vc.output.snp.LOH.SnpEff.vcf"
  )
  snpefffile_indel_loh <- paste0(
    path_input, sample,
   "_vc.output.indel.LOH.SnpEff.vcf"
  )
  # Results
  outfile_circos <- paste0(path_output, sample, "_TD_circos.pdf")
  coverage_out <- paste0(path_output, sample, "_coverage.pdf")
  # MAFs
  maf_gd <- paste0(path_output, sample, "_Germline.maf")
  maf_td <- paste0(path_output, sample, "_Somatic.maf")
  maf_loh <- paste0(path_output, sample, "_LoH.maf")
  maf_complete <- paste0(path_output, sample, ".maf")
  # Biomarker
  biomarker_out <- paste0(path_output, sample, "_Biomarker.txt")
}
if (protocol == "panelTumor" | protocol == "tumorOnly") {
  # TUMOR
  snp_file_td <- paste0(
    path_input, sample,
    "_vc.output.snp.Sample1.avinput.hg19_multianno.csv"
  )
  indel_file_td <- paste0(
    path_input, sample,
    "_vc.output.indel.Sample1.avinput.hg19_multianno.csv"
  )
  snpefffile_snp_td <- paste0(
    path_input, sample,
    "_vc.output.snp.SnpEff.vcf"
  )
  snpefffile_indel_td <- paste0(
    path_input, sample,
    "_vc.output.indel.SnpEff.vcf"
  )
  filter_out_td <- paste0(
    path_output, sample,
    "_VC_TUMOR.xlsx"
  )
  mutect2_vcf <- paste0(
    path_input, "/", sample,
    "_td_gatk4_mutect2_filtered.hg19_multianno.txt"
  )
  mutect2_snpEff_vcf <- paste0(
    path_input, "/", sample,
    "_td_gatk4_mutect2_filtered_SnpEff.vcf"
  )
  msi_file <- paste0(
    path_input, "/", sample,
    "_td_MSI"
  )

  # Results
  outfile_circos <- paste0(path_output, sample, "_TD_circos.pdf")
  outfile_circos_mutect2 <- paste0(
    path_output, sample, "_TD_circos_mutect2.pdf"
  )
  coverage_out <- paste0(path_output, sample, "_coverage.pdf")
  # MAFs
  maf_td <- paste0(path_output, sample, "_TUMOR.maf")
  maf_complete <- paste0(path_output, sample, ".maf")
  maf_td_mutect2 <- paste0(
    path_output, "/", sample, "_Mutect2_TUMOR.maf"
  )
  # Biomarker
  biomarker_out <- paste0(path_output, sample, "_Biomarker.txt")
}


##############
# STATISTICS #
## Input Files
print("Statistics.")
if (protocol == "somaticGermline" | protocol == "somatic") {
  stats_td <- paste0(path_input, sample, "_td_stats.txt")
  stats_gd <- paste0(path_input, sample, "_gd_stats.txt")
  ## Analysis
  stats <- stats_func(
    path = path_input,
    outfile_pdf = coverage_out,
    stats_td = stats_td,
    stats_gd = stats_gd,
    protocol = protocol,
    sureselect_type = sureselect_type
  )
}
if (protocol == "panelTumor" | protocol == "tumorOnly") {
  stats_td <- paste0(
    path_input,
    sample,
    "_td_stats.txt"
  )
  stats <- stats_func(
    path = path_input,
    outfile_pdf = coverage_out,
    stats_td = stats_td,
    protocol = protocol,
    sureselect_type = sureselect_type
  )
}

#####################
# Mutation Analyses
#####################

if (protocol == "somatic" | protocol == "somaticGermline") {
  # SOMATIC TUMOR
  print("Filtering for Tumor.")
  filt_result_td <- filtering(
    snpfile = snp_file_td,
    indelfile = indel_file_td,
    snpefffile_snp = snpefffile_snp,
    snpefffile_indel = snpefffile_indel,
    outfile = filter_out_td,
    outfile_maf = maf_td,
    path_data = path_data,
    path_script = path_script,
    mode = "T",
    center = center,
    id = id,
    protocol = protocol,
    sureselect = bed_file,
    vaf = vaf,
    min_var_count = min_var_count,
    maf = maf_cutoff,
    covered_exons = covered_exons,
    cov_t = stats$cover_exons$perc[[2]][1],
    sureselect_type = sureselect_type
  )

  # LOH
  print("Filtering for LoH.")
  filt_result_loh <- filtering(
    snpfile = snp_file_loh,
    indelfile = indel_file_loh,
    snpefffile_snp = snpefffile_snp_loh,
    snpefffile_indel = snpefffile_indel_loh,
    outfile = loh_out,
    outfile_maf = maf_loh,
    path_data = path_data,
    path_script = path_script,
    mode = "LOH",
    center = center,
    id = id,
    protocol = protocol,
    sureselect = bed_file,
    vaf = vaf,
    min_var_count = min_var_count,
    maf = maf_cutoff,
    covered_exons = covered_exons,
    cov_t = 1,
    sureselect_type = sureselect_type
  )

  if (protocol == "somaticGermline") {
    # GERMLINE NORMAL
    print("Filtering for Germline.")
    filt_result_gd <- filtering(
      snpfile = snp_file_gd,
      indelfile = indel_file_gd,
      snpefffile_snp = snpefffile_snp_gd,
      snpefffile_indel = snpefffile_indel_gd,
      outfile =  filter_out_gd,
      outfile_maf = maf_gd,
      path_data = path_data,
      path_script = path_script,
      mode = "N",
      center = center,
      id = id,
      protocol = protocol,
      sureselect = bed_file,
      vaf = vaf,
      min_var_count = min_var_count,
      maf = maf_cutoff,
      actionable_genes = NA,
      covered_exons = covered_exons,
      cov_t = stats$cover_exons$perc[[1]][1],
      sureselect_type = sureselect_type
    )
    loh_correction <- loh_correction(
      filt_loh = filt_result_loh,
      filt_gd = filt_result_gd,
      protocol = protocol,
      vaf = vaf
    )
    filt_result_gd <- loh_correction$gd
    filt_result_loh <- loh_correction$loh

    mutation_analysis_result_gd <- mutation_analysis(
      loh = NULL,
      somatic = filt_result_gd$table,
      tumbu = NULL,
      outfile_circos = outfile_circos_gd,
      path_data = path_data,
      path_script = path_script,
      targets_txt = targets_txt,
      protocol = "panelTumor",
      sureselect = bed_file,
      sureselect_type = sureselect_type,
      msi_file = msi_file
    )
  }

  # Analyses
  print("Variant Analyses.")
  mutation_analysis_result <- mutation_analysis(
    loh = filt_result_loh$table,
    somatic = filt_result_td$table,
    tumbu = filt_result_td$tmb,
    outfile_circos = outfile_circos,
    path_data = path_data,
    path_script = path_script,
    targets_txt = targets_txt,
    protocol = protocol,
    sureselect = bed_file,
    sureselect_type = sureselect_type,
    msi_file = msi_file
  )
}
if (protocol == "panelTumor" | protocol == "tumorOnly") {
  print("Filtering for Tumor.")
  filt_result_td <- filtering(
    snpfile = snp_file_td,
    indelfile = indel_file_td,
    snpefffile_snp = snpefffile_snp_td,
    snpefffile_indel = snpefffile_indel_td,
    outfile = filter_out_td,
    outfile_maf = maf_td,
    path_data = path_data,
    path_script = path_script,
    mode = "T",
    center = center,
    id = id,
    protocol = protocol,
    sureselect = bed_file,
    vaf = vaf,
    min_var_count = min_var_count,
    maf = maf_cutoff,
    covered_exons = covered_exons,
    cov_t = stats$cover_exons$perc[[1]][1],
    sureselect_type = sureselect_type
  )
  # GATK4 Mutect2
  filt_result_td_mutect2 <- filtering_mutect2(
    snpfile = mutect2_vcf,
    snpefffile = mutect2_snpEff_vcf,
    outfile_maf = maf_td_mutect2,
    id = id,
    path_data = path_data,
    path_script = path_script,
    mode = "T",
    center = center,
    protocol = protocol,
    sureselect = bed_file,
    vaf = vaf,
    min_var_count = min_var_count,
    maf = maf_cutoff,
    covered_exons = covered_exons,
    cov_t = stats$cover_exons$perc[[1]][1],
    sureselect_type = sureselect_type
  )

  filt_result_loh <- list(table = NULL, tmb = NULL)

  # Analyses
  print("Variant Analyses.")
  mutation_analysis_result <- mutation_analysis(
    loh = filt_result_loh$table,
    somatic = filt_result_td$table,
    tumbu = filt_result_td$tmb,
    outfile_circos = outfile_circos,
    path_data = path_data,
    path_script = path_script,
    targets_txt = targets_txt,
    protocol = protocol,
    sureselect = bed_file,
    sureselect_type = sureselect_type,
    msi_file = msi_file
  )
  # GATK4 Mutect2
  mutation_analysis_result_mutect2 <- mutation_analysis(
    loh = filt_result_loh$table,
    somatic = filt_result_td_mutect2$table,
    tumbu = filt_result_td_mutect2$tmb,
    outfile_circos = outfile_circos_mutect2,
    path_data = path_data,
    path_script = path_script,
    targets_txt = targets_txt,
    protocol = protocol,
    sureselect = bed_file,
    sureselect_type = sureselect_type,
    msi_file = msi_file
  )
}


# Combine MAF files to obtain one complete maf per patient
if (protocol == "somaticGermline" | protocol == "somatic") {
  if (protocol == "somaticGermline") {
      maf_comb <- smartbind(
        filt_result_td$maf,
        filt_result_gd$maf,
        filt_result_loh$maf
      )
      write.table(
        x = maf_comb,
        file = maf_complete,
        append = F,
        quote = F,
        sep = "\t",
        col.names = T,
        row.names = F
      )
  } else {
    maf_comb <- smartbind(
      filt_result_td$maf,
      filt_result_loh$maf
    )
    write.table(
      x = maf_comb,
      file = maf_complete,
      append = F,
      quote = F,
      sep = "\t",
      col.names = T,
      row.names = F
    )
  }
} else {
  maf_comb <- filt_result_td_mutect2$maf
  write.table(
    x = maf_comb,
    file = maf_complete,
    append = F,
    quote = F,
    sep = "\t",
    col.names = T,
    row.names = F
  )
}

########################
# Copy Number Analysis #
print("CNV Analyses.")
if (protocol == "somaticGermline" | protocol == "somatic") {
  ## Input/Output Files
  ratio_file <- paste0(
    path_input, "CNV/", sample,
    "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_ratio.txt"
  )
  cnvs_file <- paste0(
    path_input, "CNV/", sample,
    "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_CNVs"
  )
  cpn_file <- paste0(
    path_input, "CNV/", sample,
    "_td_output.sort.filtered.rmdup.realigned.fixed.recal.bam_sample.cpn"
  )

  # HRD/Purity
  purity_file <- paste0(
    path_input, "/",
    sample,
    "_sequenza/",
    sample,
    "_alternative_solutions.txt"
  )
  hrd_file <- paste0(
    path_input,
    "/",
    sample,
    "_HRD.txt"
  )

  # Results
  cnv_pvalue_txt <- paste0(cnvs_file, ".p.value.txt")
  cnv_ideogram_plot <- paste0(path_output, sample,"_CNV_Plot_Ideogram.pdf")
  outfile_cnvs_cbioportal <- paste0(path_output, sample, "_CNV_cbioportal.txt")
  outfile_cnvs_seg <- paste0(path_output, sample, "_CNV.seg")

  cnv_analysis_results <- cnv_analysis(
    ratio_file = ratio_file,
    cnvs_file = cnvs_file,
    cpn_file = cpn_file,
    cnv_pvalue_txt = cnv_pvalue_txt,
    outfile_ideogram = cnv_ideogram_plot,
    path_data = path_data,
    path_script = path_script,
    targets_txt = targets_txt,
    ampl_genes_txt = ampl_genes_txt,
    outfile_cbioportal = outfile_cnvs_cbioportal,
    outfile_seg = outfile_cnvs_seg,
    id = id,
    protocol = protocol,
    sureselect_type = sureselect_type,
    gender = gender,
    purity_file = purity_file,
    hrd_file = hrd_file
  )
print("End CNVs")
}

if (protocol == "tumorOnly" | protocol == "panelTumor") {
  ## Input/Output Files
  ratio_file <- paste0(
    path_input, "CNV/", sample,
    "_td_output.sort.rmdup.realigned.fixed.recal.bam_ratio.txt"
  )
  cnvs_file <- paste0(
    path_input, "CNV/", sample,
    "_td_output.sort.rmdup.realigned.fixed.recal.bam_CNVs"
  )
  cpn_file <- paste0(
    path_input, "CNV/", sample,
    "_td_output.sort.rmdup.realigned.fixed.recal.bam_sample.cpn"
  )
  # HRD/Purity
  purity_file <- paste0(
    path_input, "/",
    sample,
    "_sequenza/",
    sample,
    "_alternative_solutions.txt"
  )
  hrd_file <- paste0(
    path_input,
    "/",
    sample,
    "_HRD.txt"
  )

  # Results
  cnv_pvalue_txt <- paste0(cnvs_file, ".p.value.txt")
  cnv_ideogram_plot <- paste0(path_output, sample, "_CNV_Plot_Ideogram.pdf")
  outfile_cnvs_cbioportal <- paste0(path_output, sample, "_CNV_cbioportal.txt")
  outfile_cnvs_seg <- paste0(path_output, sample, "_CNV.seg")

  cnv_analysis_results <- cnv_analysis(
    ratio_file = ratio_file,
    cnvs_file = cnvs_file,
    cpn_file = cpn_file,
    cnv_pvalue_txt = cnv_pvalue_txt,
    outfile_ideogram = cnv_ideogram_plot,
    path_data = path_data,
    path_script = path_script,
    targets_txt = targets_txt,
    ampl_genes_txt = ampl_genes_txt,
    outfile_cbioportal = outfile_cnvs_cbioportal,
    outfile_seg = outfile_cnvs_seg,
    id = id,
    protocol = protocol,
    sureselect_type = sureselect_type,
    gender = gender,
    purity_file = purity_file,
    hrd_file = hrd_file
  )
}

###############################
# Mutation Signature Analysis #
print("Mutation Signature Analysis.")
if (protocol == "somaticGermline" | protocol == "somatic") {
  somatic_vcf <- paste0(
    path_input, sample,
    "_vc.output.snp.Somatic.hc.fpfilter.vcf"
  )
  outfile_mutsig_cbioportal <- paste0(path_output, sample, "_mutsig_cbioportal")
  mut_sig_ana <- mut_sig_wCI(
    vcf_file = somatic_vcf,
    cutoff = 0.01,
    sample = sample,
    sureselect_type = sureselect_type,
    path_script = path_script,
    ref_genome = ref_genome,
    target_capture_cor_factors = target_capture_cor_factors,
    path_output = path_output,
    sample_name = id,
    outfile_cbioportal = outfile_mutsig_cbioportal,
    vaf = vaf
  )
} else {
  if (protocol == "panelTumor") {
    vcf <- paste0(path_input, sample, "_td_gatk4_mutect2_filtered.vcf")
  } else {
    #vcf <- paste0(path_input, sample, "_vc.output.snp.fpfilter.vcf")
    vcf <- paste0(path_input, sample, "_td_gatk4_mutect2_filtered.vcf")
  }
  outfile_mutsig_cbioportal <- paste0(path_output, sample, "_mutsig_cbioportal")
  mut_sig_analysis <- mutation_signature_analysis(
    vcf_file = vcf,
    cutoff = 0.01,
    sample_name = NULL,
    only_coding = FALSE,
    path_data = path_data,
    path_output = path_output,
    outfile_cbioportal = outfile_mutsig_cbioportal
  )
  mut_sig <- data.frame(Signature = mut_sig_analysis$CosmicValid_cutoffGen_LCDlist$out_sig_ind_df$sig,
                        Process = mut_sig_analysis$CosmicValid_cutoffGen_LCDlist$out_sig_ind_df$process,
                        Frequency = mut_sig_analysis$CosmicValid_cutoffGen_LCDlist$norm_exposures[,1])
  rownames(mut_sig) <- mut_sig$Signature
}

# Write Excel File
print("Write Excel Table.")
if (protocol == "somaticGermline") {
  print("1 - somaticGermline")
  output <- list(
    Somatic_Mutations = filt_result_td$table,
    LoH_Mutations = filt_result_loh$table,
    Germline_Mutations = filt_result_gd$table,
    Mutation_Signatures = mut_sig_ana$output$Summary,
    CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated
  )
  #save(output, file = "excel_table_output.RData")
  write.xlsx(
    x = output,
    file = paste0(path_output, "/", sample, "_results.xlsx"),
    overwrite = TRUE,
    na.string = "."
    )
} else if (protocol == "somatic") {
  print("2 - somatic")
  output <- list(
    Somatic_Mutations = filt_result_td$table,
    LoH_Mutations = filt_result_loh$table,
    Mutation_Signatures = mut_sig_ana$output$Summary,
    CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated
  )
  write.xlsx(
    x = output,
    file = paste0(path_output, "/", sample, "_results.xlsx"),
    overwrite = TRUE,
    na.string = "."
  )
  } else if (protocol == "panelTumor") {
    print("3 - panelTumor")
    output <- list(
      Mutations = filt_result_td_mutect2$table,
      CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated,
      Mutation_Signatures = mut_sig$Summary
    )
    write.xlsx(
      x = output,
      file = paste0(path_output, "/", sample, "_results.xlsx"),
      overwrite = TRUE,
      na.string = "."
    )
} else {
  print("4 - tumorOnly")
  output <- list(
    Mutations = filt_result_td_mutect2$table,
    CopyNumberVariations = cnv_analysis_results$cnvs_annotated$CNVsAnnotated,
    Mutation_Signatures = mut_sig$Summary
  )
  write.xlsx(
    x = output,
    file = paste0(path_output, "/", sample, "_results.xlsx"),
    overwrite = TRUE,
    na.string = "."
  )
}

# Export TMB, MSI, HRD and BRCAness as text file for cBioPortal
print("cBioPortal Export.")
if(protocol == "panelTumor" & sureselect_type == "TSO500") {
  if (mutation_analysis_result_mutect2$msi < 20) {
    msi_helper <- "Non-MSI-H"
  } else {
    msi_helper <- "Instable"
  }
  brca_helper <- which(mut_sig ==  "AC3")
  if (length(brca_helper) == 1 & mut_sig["AC3", 3] * 100 > 1) {
    brca_helper <- paste0(round(mut_sig["AC3", 3] * 100, digits = 1))
  } else {
    brca_helper <- "<1%"
  }
  biomarker <- data.frame(
    Tumor_Sample_Barcode = paste(as.character(id),"TD",sep = "_"),
    MSI_SCORE = mutation_analysis_result_mutect2$msi,
    MSI_TYPE = msi_helper,
    CVR_TMB_SCORE = filt_result_td_mutect2$tmb,
    BRCAness = brca_helper,
    HRD = cnv_analysis_results$hrd$sum,
    Purity = cnv_analysis_results$purity$purity,
    Ploidity = cnv_analysis_results$purity$ploidy
  )
  write.table(x = biomarker, file = biomarker_out , append = F, quote = F,
              sep = '\t', col.names = T, row.names = F)
}
if (protocol == "somaticGermline" | protocol == "somatic") {
  if (mutation_analysis_result$msi < 10) {
    msi_helper <- "Non-MSI-H"
  } else {
    msi_helper <- "Instable"
  }
  brca_helper <- which(mut_sig_ana$output$Summary$Signature ==  "AC3")
  if (length(brca_helper) == 1 & mut_sig_ana$output$Summary["AC3", 3] > 1.0) {
    brca_helper <- paste0(
      round(
        mut_sig_ana$output$Summary["AC3", 3], digits = 1
      ), " (", round(
        mut_sig_ana$output$Summary["AC3", 4], digits = 1
      ), ";", round(
        mut_sig_ana$output$Summary["AC3", 5], digits = 1
      ), ")"
    )
  } else {
    brca_helper <- "<1%"
  }

  biomarker <- data.frame(
    Tumor_Sample_Barcode = paste(as.character(id), "TD", sep = "_"),
    MSI_SCORE = mutation_analysis_result$msi,
    MSI_TYPE = msi_helper,
    CVR_TMB_SCORE = filt_result_td$tmb,
    BRCAness = brca_helper,
    HRD = cnv_analysis_results$hrd$sum,
    Purity = cnv_analysis_results$purity$purity,
    Ploidity = cnv_analysis_results$purity$ploidy
  )
  write.table(x = biomarker, file = biomarker_out, append = F, quote = F,
              sep = "\t", col.names = TRUE, row.names = F)
}
if (protocol == "tumorOnly") {
  if (mutation_analysis_result_mutect2$msi < 10) {
    msi_helper <- "Non-MSI-H"
  } else {
    msi_helper <- "Instable"
  }
  brca_helper <- which(mut_sig ==  "AC3")
  if (length(brca_helper) == 1 & mut_sig["AC3", 3]*100 > 1) {
    brca_helper <- paste0(round(mut_sig["AC3", 3]*100, digits = 1))
  } else {
    brca_helper <- "<1%"
  }

  biomarker <- data.frame(
    Tumor_Sample_Barcode = paste(as.character(id), "TD", sep = "_"),
    MSI_SCORE = mutation_analysis_result_mutect2$msi,
    MSI_TYPE = msi_helper,
    CVR_TMB_SCORE = filt_result_td$tmb,
    BRCAness = brca_helper,
    HRD = cnv_analysis_results$hrd$sum,
    Purity = cnv_analysis_results$purity$purity,
    Ploidity = cnv_analysis_results$purity$ploidy
  )
  write.table(
    x = biomarker,
    file = biomarker_out,
    append = F,
    quote = F,
    sep = "\t",
    col.names = T,
    row.names = F
  )
}

#######################
####### FUSIONS #######
#######################
print("Fusion Analysis.")
fus_file <- paste(
  path_input,
  "../RNA/fusioncatcher/final-list_candidate-fusion-genes.hg19.txt",
  sep = "/"
)
outfile_fusions_cbioportal <- paste0(
  path_output, "/", sample, "_fusions_cbioportal.txt"
)
if (file.exists(fus_file)) {
  source(paste(path_script, "fusions.R", sep = "/"))
  fusions <- fusions_ana(
    fus_file = fus_file,
    path_data = path_data,
    fusion_genes = fusion_genes
  )
  cat(dim(fusions$Table))
  # print Circosplot with Fusions
  if (!is.null(fusions$Table)) {
    if (dim(fusions$Table)[1] != 0) {
      if (length(which(duplicated(fusions$Table[, c(1,3)]))) != 0) {
        fus_tab <- fusions$Table[-which(
          duplicated(fusions$Table[, c(1,3)])
        ), c(1:4)]
      } else {
        fus_tab <- fusions$Table[, c(1:5) ]
      }
      id_1 <- grep(pattern = "not-converted", x = fus_tab[, 2])
      id_2 <- grep(pattern = "not-converted", x = fus_tab[, 4])
      if(length(union(id_1, id_2)) > 0){
        fus_tab <- fus_tab[-union(id_1, id_2), ]
      }
      # Prepare Fusionstable for CircosPlot
      sep <- strsplit(x = as.character(fus_tab$Bruch1), split = ":")
      fus_tab$Chr1 <- unlist(lapply(sep, function(x){return( x[[1]])}))
      fus_tab$Pos1 <- unlist(lapply(sep, function(x){return( x[[2]])}))
      sep <- strsplit(x = as.character(fus_tab$Bruch2), split = ":")
      fus_tab$Chr2 <- unlist(lapply(sep, function(x){return( x[[1]])}))
      fus_tab$Pos2 <- unlist(lapply(sep, function(x){return( x[[2]])}))
      id_genes <- union(
        which(
          nchar(fus_tab$Chr1) > 2
        ),
        which(
          nchar(fus_tab$Chr2) > 2
        ))
      if (length(id_genes) > 0){
        fus_tab <- fus_tab[-id_genes, c(
          "Gen1",
          "Chr1",
          "Pos1",
          "Gen2",
          "Chr2",
          "Pos2"
        )]
      } else {
        fus_tab <- fus_tab[, c("Gen1", "Chr1", "Pos1", "Gen2", "Chr2", "Pos2")]
      }
      if (sureselect_type == "TSO500") {
        sub <- div(
          filt_result_td_mutect2$table,
          NULL,
          TRUE,
          protocol = protocol,
          sureselect_type = "V5UTR"
        )
      } else {
        sub <- div(
          filt_result_td$table,
          NULL,
          TRUE,
          protocol = protocol,
          sureselect_type = "V5UTR"
        )
      }

      cc <- circos_colors(
        x_s_snp = sub$x_s_snp,
        x_s_indel = sub$x_s_indel,
        x_l_snp = sub$x_l_snp,
        x_l_indel = sub$x_l_indel, no_loh = sub$no_loh,
        no_indel_somatic = sub$no_indel_somatic,
        no_snp = sub$no_snp,
        no_indel_loh = sub$no_indel_loh,
        no_snp_loh = sub$no_snp_loh
      )
      omicCircosFus2(
        listOfMap = as.matrix(cc$map_mat),
        fusions = fus_tab,
        label = NULL, minR = 125,
        outfile = outfile_circos_mutect2,
        circosColors = as.vector(cc$circoscolors),
        mode = sureselect_type,
        protocol = "tumorOnly",
        path_data = path_data,
        trgt = bed_file
      )
      if (any(is.na(fusions$Plots$file))) {
        fus_38_files <- paste(
          fusions$Plots$paths[!is.na(fusions$Plots$file)],
          fusions$Plots$file[!is.na(fusions$Plots$file)], sep = "/"
        )
        for (i in 1:length(fus_38_files)) {
          system(paste0("cp ", fus_38_files[i], " ", path_output))
        }
      } else {
        fus_38_files <- paste(
          fusions$Plots$paths,
          fusions$Plots$file,
          sep = "/"
        )
        for (i in 1:length(fus_38_files)) {
          system(paste0("cp ", fus_38_files[i], " ", path_output))
        }
      }
      fusions2cbioportal(fusions, sample, outfile_fusions_cbioportal)
    } else {
      fusions <- NULL
    }
  } else {
    fusions <- NULL
  } 
} else {
  fusions <- NULL
}
save.image(file = "MTB.RData")

# Prepare Report
print("Report Preparation.")
source(paste(path_script, "Report_tools.R", sep = "/"))

if (protocol == "somaticGermline"){

  tmb_med <- med_tmb(as.character(entity))

  key_results <- keys(
    mut_sig = mut_sig_ana,
    mutation_analysis_result = mutation_analysis_result,
    mutation_analysis_result_gd = mutation_analysis_result_gd,
    filt_result_td = filt_result_td,
    cnv_analysis_results = cnv_analysis_results,
    filt_result_gd = filt_result_gd,
    med_tmb = tmb_med,
    protocol = protocol,
    fusions = fusions
  )

  highlight_table <- highlight(
    muts_tab = mutation_analysis_result$som_mut_tab,
    protocol = protocol
  )

  sq_tab <- summary_quality(stats = stats, protocol = protocol)

  som_muts <- sum_muts(
    tmp = mutation_analysis_result$mut_tab
  )
  sum_mut_cg <- highlight_detail(
    muts_tab = mutation_analysis_result$ts_og,
    Mode = "Tumor",
    protocol = protocol
  )
  if (length(which(
      mutation_analysis_result$table_loh_mutations$is_oncogene == 1 |
      mutation_analysis_result$table_loh_mutations$is_tumorsuppressor == 1
    )) > 0) {
    sum_loh_cg <- highlight_detail(
      muts_tab = mutation_analysis_result$table_loh_mutations[
        which(
          mutation_analysis_result$table_loh_mutations$is_oncogene == 1 |
          mutation_analysis_result$table_loh_mutations$is_tumorsuppressor == 1
        ),
      ],
      Mode = "LoH",
      protocol = protocol
    )
  } else {
    sum_loh_cg <- NULL
  }
  cnvs_og <- cnv_cg(
    gene_loci = cnv_analysis_results$gene_loci_onc,
    type = "OG"
  )
  cnvs_tsg <- cnv_cg(
    gene_loci = cnv_analysis_results$gene_loci_tsg,
    type = "TSG"
  )
  som_mut_pthw <- pthws_mut(
    df = mutation_analysis_result$important_pathways,
    protocol = protocol
  )
  som_mut_topart <- topart_mut(
    df = mutation_analysis_result$important_pathways,
    protocol = protocol
  )
  cnvs_pthws <- pathws_cnv(df = cnv_analysis_results$impa)

  germ_mut_cg <- highlight_detail(
    muts_tab = filt_result_gd$table,
    Mode = "Germline",
    protocol = protocol
  )
  germ_mut_pthw <- pthws_mut(
    df = mutation_analysis_result_gd$important_pathways,
    protocol = protocol
  )

  som_all <- highlight_detail(
    muts_tab = mutation_analysis_result$som_mut_tab,
    Mode = "Tumor",
    protocol = protocol
  )
  germ_all <- highlight_detail(
    muts_tab = mutation_analysis_result_gd$som_mut_tab,
    Mode = "Tumor",
    protocol = protocol
  )
  loh_all <- highlight_detail(
    muts_tab = mutation_analysis_result$table_loh_mutations,
    Mode = "LoH",
    protocol = protocol
  )

  save.image("Report.RData")
}
if (protocol == "panelTumor" | protocol == "tumorOnly") {

  tmb_med <- med_tmb(as.character(entity))
  key_results <- keys(
    mut_sig = mut_sig,
    mutation_analysis_result = mutation_analysis_result_mutect2,
    mutation_analysis_result_gd = NULL,
    filt_result_td = filt_result_td_mutect2,
    cnv_analysis_results = cnv_analysis_results,
    filt_result_gd = NULL,
    med_tmb = tmb_med,
    protocol = protocol,
    fusions = fusions
  )

  highlight_table <- highlight(
    muts_tab = mutation_analysis_result_mutect2$som_mut_tab, protocol = protocol
  )

  sq_tab <- summary_quality(stats = stats, protocol = protocol)
  #
  som_muts <- mutation_analysis_result_mutect2$mut_tab
  sum_mut_cg <- highlight_detail(
    muts_tab = mutation_analysis_result_mutect2$ts_og,
    Mode = "Tumor",
    protocol = protocol
  )
  sum_loh_cg <- NULL

  if (protocol == "panelTumor") {
    cnvs <- cnv_panel(cnv_results = cnv_analysis_results$out)
  }
  if (protocol == "tumorOnly") {
    cnvs_og <- cnv_cg(
      gene_loci = cnv_analysis_results$gene_loci_onc,
      type = "OG"
    )
    cnvs_tsg <- cnv_cg(
      gene_loci = cnv_analysis_results$gene_loci_tsg,
      type = "TSG"
    )
  }
  som_mut_pthw <- pthws_mut(
    df = mutation_analysis_result_mutect2$important_pathways,
    protocol = protocol
  )
  som_mut_topart <- topart_mut(
    df = mutation_analysis_result_mutect2$important_pathways,
    protocol = protocol
  )
  cnvs_pthws <- pathws_cnv(df = cnv_analysis_results$impa)

  germ_mut_cg <- NULL
  germ_mut_pthw <- NULL

  som_all <- highlight_detail(
    muts_tab = mutation_analysis_result_mutect2$som_mut_tab,
    Mode = "Tumor",
    protocol = protocol
  )
  germ_all <- NULL
  loh_all <- NULL

  save.image("Report.RData")
}
