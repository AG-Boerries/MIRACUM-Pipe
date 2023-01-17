#!/usr/bin/env bash

# common settings
readonly CFG_AUTHOR=$(get_config_value common.author "${PARAM_DIR_PATIENT}")
readonly CFG_CENTER=$(get_config_value common.center "${PARAM_DIR_PATIENT}")
readonly CFG_PROTOCOL=$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")
#readonly CFG_SEX=$(get_config_value sex "${PARAM_DIR_PATIENT}")
readonly CFG_ENTITY=$(get_config_value common.entity "${PARAM_DIR_PATIENT}")

# temporary folder
readonly DIR_TMP="$(get_config_value common.dirTmp "${PARAM_DIR_PATIENT}")/${PARAM_DIR_PATIENT}"

readonly CFG_FILE_TUMOR_R1=$(get_config_value common.files.tumor_R1 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_TUMOR_R2=$(get_config_value common.files.tumor_R2 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE_R1=$(get_config_value common.files.germline_R1 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE_R2=$(get_config_value common.files.germline_R2 "${PARAM_DIR_PATIENT}")

readonly CFG_PANEL_FILE_TUMOR=$(get_config_value common.files.panel.tumor "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_FILE_NUMBER=$(get_config_value common.files.panel.numberOfFiles "${PARAM_DIR_PATIENT}")

readonly CFG_FOLDER_RNA=$(get_config_value common.RNA.folder "${PARAM_DIR_PATIENT}")

# folder containing patient output
readonly DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
readonly DIR_WES="${DIR_TARGET}/WES"
readonly DIR_ANALYSES="${DIR_TARGET}/Analyses"
readonly DIR_RNA="${DIR_TARGET}/RNA"
readonly DIR_FUSIONS="${DIR_RNA}/fusioncatcher"

# end paths

# ucsc mysql server
readonly CFG_UCSC_SERVER=$(get_config_value common.ucscServer "${PARAM_DIR_PATIENT}")

# CNV annotation
readonly CFG_CNV_ANNOTATION=$(get_config_value common.cnvAnnotation "${PARAM_DIR_PATIENT}")

## Genome
readonly FILE_GENOME="${DIR_REF}/genome/$(get_config_value reference.genome "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_LENGTH="${DIR_CHROMOSOMES}/$(get_config_value reference.length "${PARAM_DIR_PATIENT}")"
readonly HRD_REF_WIG="${DIR_DATABASE}/$(get_config_value reference.hrdRef "${PARAM_DIR_PATIENT}")"

# depending on measurement machine
## SureSelect (Capture Kit)
readonly CFG_REFERENCE_CAPTUREREGIONS="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.captureRegions "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTUREGENES="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.captureGenes "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_COVEREDREGION="$(get_config_value reference.sequencing.coveredRegion "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTUREREGIONNAME="$(get_config_value reference.sequencing.captureRegionName "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTURECORFACTORS="${DIR_DATABASE}/$(get_config_value reference.sequencing.captureCorFactors "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_COVERED_EXONS="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.coveredExons "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_ACTIONABLEGENES="${DIR_DATABASE}/$(get_config_value reference.sequencing.actionableGenes "${PARAM_DIR_PATIENT}")"

# database for known variants
## dbSNP vcf File
readonly CFG_REFERENCE_DBSNP="${DIR_DBSNP}/$(get_config_value reference.dbSNP "${PARAM_DIR_PATIENT}")"
# END variables

# if no parallel computation, set processes to 1
if [[ -z ${PARALLEL_PROCESSES} ]]; then
  readonly TMP_PROCESSES=1
else
  readonly TMP_PROCESSES="${PARALLEL_PROCESSES}"
fi

# take cpucores/2
readonly CFG_COMMON_CPUCORES=$(($(get_config_value common.cpucores "${PARAM_DIR_PATIENT}")/${TMP_PROCESSES}))

# take memory/2
readonly tmp_memory=$(get_config_value common.memory "${PARAM_DIR_PATIENT}")
readonly CFG_COMMON_MEMORY="$(("${tmp_memory//[^0-9.]/}"/${TMP_PROCESSES}))${tmp_memory//[^a-zA-Z]/}"

# general parameters
readonly CFG_GENERAL_MINBASEQUAL=$(get_config_value general.minBaseQual "${PARAM_DIR_PATIENT}")
readonly CFG_GENERAL_MAFCUTOFF=$(get_config_value general.maf_cutoff "${PARAM_DIR_PATIENT}")
readonly CFG_GENERAL_MINVAF=$(get_config_value general.minVAF "${PARAM_DIR_PATIENT}")

# sametools mpileup
readonly CFG_SAMTOOLS_MPILEUP_MINMQ=$(get_config_value wes.samtools.mpileup.minMQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_ADJUSTMQ=$(get_config_value wes.samtools.mpileup.adjustMQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_MAXDEPTH=$(get_config_value wes.samtools.mpileup.maxDepth "${PARAM_DIR_PATIENT}")

# mpileup flags
readonly CFG_SAMTOOLS_MPILEUP_ILLUMINA=$(get_config_value wes.samtools.mpileup.illumina "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_COUNTORPHANS=$(get_config_value wes.samtools.mpileup.countOrphans "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_NOBAQ=$(get_config_value wes.samtools.mpileup.noBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_REDOBAQ=$(get_config_value wes.samtools.mpileup.redoBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_EXCLUDERG=$(get_config_value wes.samtools.mpileup.excludeRG "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_POSITION=$(get_config_value wes.samtools.mpileup.position "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_REGION=$(get_config_value wes.samtools.mpileup.region "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_IGNORERG=$(get_config_value wes.samtools.mpileup.ignoreRG "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_IGNOREOVERLAPS=$(get_config_value wes.samtools.mpileup.ignoreOverlaps "${PARAM_DIR_PATIENT}")

CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_ILLUMINA,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -6"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_COUNTORPHANS,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -A"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_NOBAQ,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -B"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_REDOBAQ,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -E"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_EXCLUDERG,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -G"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_POSITION,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -l"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_REGION,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -r"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_IGNORERG,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -R"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_IGNOREOVERLAPS,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -x"

# VarScan somatic
readonly CFG_VARSCAN_SOMATIC_MINCOVERAGE=$(get_config_value wes.varscan.somatic.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_TUMORPURITY=$(get_config_value wes.varscan.somatic.tumorPurity "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_MINFREQFORHOM=$(get_config_value wes.varscan.somatic.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_MINCOVERAGENORMAL=$(get_config_value wes.varscan.somatic.minCoverageNormal "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_MINCOVERAGETUMOR=$(get_config_value wes.varscan.somatic.minCoverageTumor "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_NORMALPURITY=$(get_config_value wes.varscan.somatic.normalPurity "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_PVALUE=$(get_config_value wes.varscan.somatic.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_SOMATICPVALUE=$(get_config_value wes.varscan.somatic.somaticPValue "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_STRANDFILTER=$(get_config_value wes.varscan.somatic.strandFilter "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_VALIDATION=$(get_config_value wes.varscan.somatic.validation "${PARAM_DIR_PATIENT}")

# Varscan processSomatic
readonly CFG_VARSCAN_PROCESSSOMATIC_MAXNORMALFREQ=$(get_config_value wes.varscan.processSomatic.maxNormalFreq "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PROCESSSOMATIC_PVALUE=$(get_config_value wes.varscan.processSomatic.pValue "${PARAM_DIR_PATIENT}")

# VarScan fpfilter
readonly CFG_VARSCAN_FPFILTER_MINVARCOUNT=$(get_config_value wes.varscan.fpfilter.minVarCount "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARCOUNTLC=$(get_config_value wes.varscan.fpfilter.minVarCountLC "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXSOMATICP=$(get_config_value wes.varscan.fpfilter.maxSomaticP "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXSOMATICPDEPTH=$(get_config_value wes.varscan.fpfilter.maxSomaticPDepth "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFREADPOS=$(get_config_value wes.varscan.fpfilter.minRefReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARREADPOS=$(get_config_value wes.varscan.fpfilter.minVarReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFDIST3=$(get_config_value wes.varscan.fpfilter.minRefDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARDIST3=$(get_config_value wes.varscan.fpfilter.minVarDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINSTRANDEDNESS=$(get_config_value wes.varscan.fpfilter.minStrandedness "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINSTRANDREADS=$(get_config_value wes.varscan.fpfilter.minStrandReads "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXBASEQUALDIFF=$(get_config_value wes.varscan.fpfilter.maxBasequalDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFAVGRL=$(get_config_value wes.varscan.fpfilter.minRefAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARAVGRL=$(get_config_value wes.varscan.fpfilter.minVarAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXRLDIFF=$(get_config_value wes.varscan.fpfilter.maxRlDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXREFMMQS=$(get_config_value wes.varscan.fpfilter.maxRefMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXVARMMQS=$(get_config_value wes.varscan.fpfilter.maxVarMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINMMQSDIFF=$(get_config_value wes.varscan.fpfilter.minMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXMMQSDIFF=$(get_config_value wes.varscan.fpfilter.maxMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFMAPQUAL=$(get_config_value wes.varscan.fpfilter.minRefMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARMAPQUAL=$(get_config_value wes.varscan.fpfilter.minVarMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXMAPQUALDIFF=$(get_config_value wes.varscan.fpfilter.maxMapQualDiff "${PARAM_DIR_PATIENT}")

# mutect2
readonly CFG_MUTECT_CALLABLEDEPTH=$(get_config_value wes.mutect.callableDepth "${PARAM_DIR_PATIENT}")

# Panel Parameter
# sametools mpileup
readonly CFG_PANEL_SAMTOOLS_MPILEUP_MINMQ=$(get_config_value panel.samtools.mpileup.minMQ "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_ADJUSTMQ=$(get_config_value panel.samtools.mpileup.adjustMQ "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_MAXDEPTH=$(get_config_value panel.samtools.mpileup.maxDepth "${PARAM_DIR_PATIENT}")

# mpileup flags
readonly CFG_PANEL_SAMTOOLS_MPILEUP_ILLUMINA=$(get_config_value panel.samtools.mpileup.illumina "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_COUNTORPHANS=$(get_config_value panel.samtools.mpileup.countOrphans "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_NOBAQ=$(get_config_value panel.samtools.mpileup.noBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_REDOBAQ=$(get_config_value panel.samtools.mpileup.redoBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_EXCLUDERG=$(get_config_value panel.samtools.mpileup.excludeRG "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_POSITION=$(get_config_value panel.samtools.mpileup.position "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_REGION=$(get_config_value panel.samtools.mpileup.region "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_IGNORERG=$(get_config_value panel.samtools.mpileup.ignoreRG "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_SAMTOOLS_MPILEUP_IGNOREOVERLAPS=$(get_config_value panel.samtools.mpileup.ignoreOverlaps "${PARAM_DIR_PATIENT}")

CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_ILLUMINA,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -6"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_COUNTORPHANS,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -A"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_NOBAQ,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -B"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_REDOBAQ,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -E"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_EXCLUDERG,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -G"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_POSITION,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -l"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_REGION,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -r"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_IGNORERG,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -R"
CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_PANEL_SAMTOOLS_MPILEUP_IGNOREOVERLAPS,,}" == "true" ]] && CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP="${CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP} -x"

# mpileup2snp
readonly CFG_PANEL_VARSCAN_MPILEUP2SNP_MINCOVERAGE=$(get_config_value panel.varscan.mpileup2snp.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2SNP_MINREADS2=$(get_config_value panel.varscan.mpileup2snp.minReads2 "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2SNP_MINFREQFORHOM=$(get_config_value panel.varscan.mpileup2snp.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2SNP_PVALUE=$(get_config_value panel.varscan.mpileup2snp.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2SNP_STRANDFILTER=$(get_config_value panel.varscan.mpileup2snp.strandFilter "${PARAM_DIR_PATIENT}")

# mpileup2indel
readonly CFG_PANEL_VARSCAN_MPILEUP2INDEL_MINCOVERAGE=$(get_config_value panel.varscan.mpileup2indel.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2INDEL_MINREADS2=$(get_config_value panel.varscan.mpileup2indel.minReads2 "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM=$(get_config_value panel.varscan.mpileup2indel.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2INDEL_PVALUE=$(get_config_value panel.varscan.mpileup2indel.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_MPILEUP2INDEL_STRANDFILTER=$(get_config_value panel.varscan.mpileup2indel.strandFilter "${PARAM_DIR_PATIENT}")

# varscan fpfilter
readonly CFG_PANEL_VARSCAN_FPFILTER_MINVARCOUNT=$(get_config_value panel.varscan.fpfilter.minVarCount "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINVARCOUNTLC=$(get_config_value panel.varscan.fpfilter.minVarCountLC "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXSOMATICP=$(get_config_value panel.varscan.fpfilter.maxSomaticP "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXSOMATICPDEPTH=$(get_config_value panel.varscan.fpfilter.maxSomaticPDepth "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINREFREADPOS=$(get_config_value panel.varscan.fpfilter.minRefReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINVARREADPOS=$(get_config_value panel.varscan.fpfilter.minVarReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINREFDIST3=$(get_config_value panel.varscan.fpfilter.minRefDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINVARDIST3=$(get_config_value panel.varscan.fpfilter.minVarDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINSTRANDEDNESS=$(get_config_value panel.varscan.fpfilter.minStrandedness "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINSTRANDREADS=$(get_config_value panel.varscan.fpfilter.minStrandReads "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXBASEQUALDIFF=$(get_config_value panel.varscan.fpfilter.maxBasequalDiff "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINREFAVGRL=$(get_config_value panel.varscan.fpfilter.minRefAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINVARAVGRL=$(get_config_value panel.varscan.fpfilter.minVarAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXRLDIFF=$(get_config_value panel.varscan.fpfilter.maxRlDiff "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXREFMMQS=$(get_config_value panel.varscan.fpfilter.maxRefMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXVARMMQS=$(get_config_value panel.varscan.fpfilter.maxVarMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINMMQSDIFF=$(get_config_value panel.varscan.fpfilter.minMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXMMQSDIFF=$(get_config_value panel.varscan.fpfilter.maxMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINREFMAPQUAL=$(get_config_value panel.varscan.fpfilter.minRefMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MINVARMAPQUAL=$(get_config_value panel.varscan.fpfilter.minVarMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_VARSCAN_FPFILTER_MAXMAPQUALDIFF=$(get_config_value panel.varscan.fpfilter.maxMapQualDiff "${PARAM_DIR_PATIENT}")

# mutect2
readonly CFG_PANEL_MUTECT_CALLABLEDEPTH=$(get_config_value panel.mutect.callableDepth "${PARAM_DIR_PATIENT}")

# Fusions
readonly CFG_FUSION_GENES="${DIR_SEQUENCING}/$(get_config_value panel.fusions.fusionGenes "${PARAM_DIR_PATIENT}")"

# Amplifications
readonly CFG_AMPLIFICATION_GENES="${DIR_SEQUENCING}/$(get_config_value panel.amplification.amplificationGenes "${PARAM_DIR_PATIENT}")"

## tumorOnly Parameter
# samtools mpileup
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_MINMQ=$(get_config_value tumorOnly.samtools.mpileup.minMQ "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_ADJUSTMQ=$(get_config_value tumorOnly.samtools.mpileup.adjustMQ "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_MAXDEPTH=$(get_config_value tumorOnly.samtools.mpileup.maxDepth "${PARAM_DIR_PATIENT}")

# samtools mpileup flags
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_ILLUMINA=$(get_config_value tumorOnly.samtools.mpileup.illumina "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_COUNTORPHANS=$(get_config_value tumorOnly.samtools.mpileup.countOrphans "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_NOBAQ=$(get_config_value tumorOnly.samtools.mpileup.noBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_REDOBAQ=$(get_config_value tumorOnly.samtools.mpileup.redoBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_EXCLUDERG=$(get_config_value tumorOnly.samtools.mpileup.excludeRG "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_POSITION=$(get_config_value tumorOnly.samtools.mpileup.position "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_REGION=$(get_config_value tumorOnly.samtools.mpileup.region "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_IGNORERG=$(get_config_value tumorOnly.samtools.mpileup.ignoreRG "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_SAMTOOLS_MPILEUP_IGNOREOVERLAPS=$(get_config_value tumorOnly.samtools.mpileup.ignoreOverlaps "${PARAM_DIR_PATIENT}")

CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_ILLUMINA,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -6"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_COUNTORPHANS,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -A"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_NOBAQ,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -B"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_REDOBAQ,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -E"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_EXCLUDERG,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -G"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_POSITION,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -l"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_REGION,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -r"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_IGNORERG,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -R"
CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_IGNOREOVERLAPS,,}" == "true" ]] && CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP="${CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP} -x"

# varscan mpileup2snp
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINCOVERAGE=$(get_config_value tumorOnly.varscan.mpileup2snp.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINREADS2=$(get_config_value tumorOnly.varscan.mpileup2snp.minReads2 "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINFREQFORHOM=$(get_config_value tumorOnly.varscan.mpileup2snp.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_PVALUE=$(get_config_value tumorOnly.varscan.mpileup2snp.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_STRANDFILTER=$(get_config_value tumorOnly.varscan.mpileup2snp.strandFilter "${PARAM_DIR_PATIENT}")

# varscan mpileup2indel
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINCOVERAGE=$(get_config_value tumorOnly.varscan.mpileup2indel.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINREADS2=$(get_config_value tumorOnly.varscan.mpileup2indel.minReads2 "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM=$(get_config_value tumorOnly.varscan.mpileup2indel.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_PVALUE=$(get_config_value tumorOnly.varscan.mpileup2indel.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_STRANDFILTER=$(get_config_value tumorOnly.varscan.mpileup2indel.strandFilter "${PARAM_DIR_PATIENT}")

# varscan fpfilter
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARCOUNT=$(get_config_value tumorOnly.varscan.fpfilter.minVarCount "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARCOUNTLC=$(get_config_value tumorOnly.varscan.fpfilter.minVarCountLC "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXSOMATICP=$(get_config_value tumorOnly.varscan.fpfilter.maxSomaticP "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXSOMATICPDEPTH=$(get_config_value tumorOnly.varscan.fpfilter.maxSomaticPDepth "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFREADPOS=$(get_config_value tumorOnly.varscan.fpfilter.minRefReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARREADPOS=$(get_config_value tumorOnly.varscan.fpfilter.minVarReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFDIST3=$(get_config_value tumorOnly.varscan.fpfilter.minRefDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARDIST3=$(get_config_value tumorOnly.varscan.fpfilter.minVarDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINSTRANDEDNESS=$(get_config_value tumorOnly.varscan.fpfilter.minStrandedness "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINSTRANDREADS=$(get_config_value tumorOnly.varscan.fpfilter.minStrandReads "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXBASEQUALDIFF=$(get_config_value tumorOnly.varscan.fpfilter.maxBasequalDiff "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFAVGRL=$(get_config_value tumorOnly.varscan.fpfilter.minRefAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARAVGRL=$(get_config_value tumorOnly.varscan.fpfilter.minVarAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXRLDIFF=$(get_config_value tumorOnly.varscan.fpfilter.maxRlDiff "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXREFMMQS=$(get_config_value tumorOnly.varscan.fpfilter.maxRefMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXVARMMQS=$(get_config_value tumorOnly.varscan.fpfilter.maxVarMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINMMQSDIFF=$(get_config_value tumorOnly.varscan.fpfilter.minMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXMMQSDIFF=$(get_config_value tumorOnly.varscan.fpfilter.maxMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFMAPQUAL=$(get_config_value tumorOnly.varscan.fpfilter.minRefMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARMAPQUAL=$(get_config_value tumorOnly.varscan.fpfilter.minVarMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_TUMORONLY_VARSCAN_FPFILTER_MAXMAPQUALDIFF=$(get_config_value tumorOnly.varscan.fpfilter.maxMapQualDiff "${PARAM_DIR_PATIENT}")

# mutect2
readonly CFG_TUMORONLY_MUTECT_CALLABLEDEPTH=$(get_config_value tumorOnly.mutect.callableDepth "${PARAM_DIR_PATIENT}")

# ANNOVAR Databases
readonly CFG_ANNOVAR_PROTOCOL=$(get_config_value annovar.protocol "${PARAM_DIR_PATIENT}")
readonly CFG_ANNOVAR_ARGOP=$(get_config_value annovar.argop "${PARAM_DIR_PATIENT}")

## Tools and paths
# Paths
readonly BIN_JAVA="java -Djava.io.tmpdir=${DIR_TMP} " # path to java

# Pre-Processing
readonly BIN_FASTQC="${DIR_TOOLS}/FastQC/bin/fastqc -t ${CFG_COMMON_CPUCORES} --extract "

readonly BIN_TRIM="${BIN_JAVA} -Xcompressedrefs -jar ${DIR_TOOLS}/Trimmomatic/trimmomatic.jar PE -threads ${CFG_COMMON_CPUCORES} -phred33 "
readonly DIR_TRIMMOMATIC_ADAPTER="${DIR_TOOLS}/Trimmomatic/adapters"
readonly BIN_CUT="cut -f1,2,3"

# Alignment
readonly BIN_BWAMEM="${DIR_TOOLS}/bwa/bwa mem -M "

# BAM-Readcount
readonly BIN_BAM_READCOUNT="${DIR_TOOLS}/bam-readcount/bin/bam-readcount -q ${CFG_SAMTOOLS_MPILEUP_MINMQ} -b ${CFG_GENERAL_MINBASEQUAL} -w 1 -f ${FILE_GENOME} "

# SAMTOOLS
readonly BIN_SAMTOOLS="${DIR_TOOLS}/samtools/samtools" # path to samtools
readonly BIN_SAMVIEW="${BIN_SAMTOOLS} view -@ ${CFG_COMMON_CPUCORES} "
readonly BIN_SAMSORT="${BIN_SAMTOOLS} sort -@ ${CFG_COMMON_CPUCORES} "
readonly BIN_SAMRMDUP="${BIN_SAMTOOLS} rmdup "
readonly BIN_SAMINDEX="${BIN_SAMTOOLS} index "
readonly BIN_MPILEUP="${BIN_SAMTOOLS} mpileup ${CFG_SAMTOOLS_FLAGS_MPILEUP} "
readonly BIN_STATS="${BIN_SAMTOOLS} stats "

# GATK
readonly BIN_GATK="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/gatk/GenomeAnalysisTK.jar"
readonly BIN_REALIGNER_TARGER_CREATOR="${BIN_GATK} -T RealignerTargetCreator -R ${FILE_GENOME} -nt ${CFG_COMMON_CPUCORES} "
readonly BIN_INDEL_REALIGNER="${BIN_GATK} -R ${FILE_GENOME} -T IndelRealigner "
readonly BIN_BASE_RECALIBRATOR="${BIN_GATK} -T BaseRecalibrator -l INFO -R ${FILE_GENOME} -knownSites ${CFG_REFERENCE_DBSNP} -nct ${CFG_COMMON_CPUCORES} "
readonly BIN_PRINT_READS="${BIN_GATK} -T PrintReads -R ${FILE_GENOME} -nct ${CFG_COMMON_CPUCORES} "

# GATK4
readonly BIN_GATK4="${DIR_TOOLS}/gatk4/gatk" # --java-options '-Xmx${CFG_COMMON_MEMORY}'"

# PICARD
readonly BIN_FIX_MATE="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -Dpicard.useLegacyParser=false -jar ${DIR_TOOLS}/picard/picard.jar FixMateInformation "

# VARSCAN
readonly BIN_VAR_SCAN="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/varscan/VarScan.jar"
readonly BIN_SOMATIC="${BIN_VAR_SCAN} somatic"
readonly BIN_PROCESSSOMATIC="${BIN_VAR_SCAN} processSomatic"

# ANNOVAR
readonly DIR_ANNOVAR="${DIR_TOOLS}/annovar/"
readonly DIR_ANNOVAR_DATA="${DIR_ANNOVAR}/humandb"
readonly CONVERT2ANNOVAR2="${DIR_ANNOVAR}/convert2annovar.pl --format vcf4old --outfile "
readonly CONVERT2ANNOVAR3="${DIR_ANNOVAR}/convert2annovar.pl --format vcf4old --includeinfo --comment --outfile "
readonly CONVERT2ANNOVAR="${DIR_ANNOVAR}/convert2annovar.pl --format vcf4 --includeinfo --comment --withzyg --outfile "
readonly TABLEANNOVAR="${DIR_ANNOVAR}/table_annovar.pl"

# COVERAGE
readonly BIN_COVERAGE="${DIR_TOOLS}/bedtools2/bin/bedtools coverage -hist -g ${FILE_GENOME}.fai -sorted "

# SNPEFF
readonly BIN_SNPEFF="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/snpEff/snpEff.jar GRCh37.75 -c ${DIR_TOOLS}/snpEff/snpEff.config -canon -v -noLog -noStats"

# FREEC
readonly BIN_FREEC="${DIR_TOOLS}/FREEC/bin/freec "
readonly FILE_REFERENCE_MAPPABILITY="${DIR_REF}/mappability/$(get_config_value reference.mappability "${PARAM_DIR_PATIENT}")"

# R
readonly BIN_RSCRIPT=$(command -v Rscript)

# fusioncatcher
readonly FUSIONCATCHER_DB="${DIR_TOOLS}/fusioncatcher/data/current"
readonly BIN_FUSIONCATCHER="${DIR_TOOLS}/fusioncatcher/bin/fusioncatcher.py -p ${CFG_COMMON_CPUCORES} -d ${FUSIONCATCHER_DB} "

# msisnesor2
readonly MSISENSOR2="${DIR_TOOLS}/msisensor2/msisensor2 msi -b ${CFG_COMMON_CPUCORES} -M ${DIR_TOOLS}/msisensor2/models_hg19_GRCh37"

# msisensor-pro
readonly MSISENSOR_PRO="${DIR_TOOLS}/msisensor-pro/binary/msisensor-pro msi -b ${CFG_COMMON_CPUCORES}"
readonly MSISENSOR_PRO_SCAN="${DIR_TOOLS}/msisensor-pro/binary/msisensor-pro scan"
readonly MICROSATELLITE_SITES="${DIR_DATABASE}/$(get_config_value reference.microsatelliteSites "${PARAM_DIR_PATIENT}")"

# HRD
readonly SEQUENZA_UTILS=$(command -v sequenza-utils)
readonly SEQUENZA_WINDOW=$(get_config_value sequenza.window "${PARAM_DIR_PATIENT}")
readonly SEQUENZA_NON_MATCHING_NORMAL="${DIR_REF}/sequenza/$(get_config_value sequenza.nonMatchingNormal "${PARAM_DIR_PATIENT}")"
readonly SEQUENZA_CHROMOSOMES=$(get_config_value sequenza.chromosomes "${PARAM_DIR_PATIENT}")

# export parameters
export CFG_AUTHOR
export CFG_CENTER
export CFG_PROTOCOL
export CFG_ENTITY

export CFG_FILE_TUMOR_R1
export CFG_FILE_TUMOR_R2
export CFG_FILE_GERMLINE_R1
export CFG_FILE_GERMLINE_R2

export CFG_PANEL_FILE_TUMOR
export CFG_PANEL_FILE_NUMBER

export CFG_FOLDER_RNA

export CFG_UCSC_SERVER

export CFG_CNV_ANNOTATION

export DIR_TARGET
export DIR_WES
export DIR_ANALYSES
export DIR_RNA
export DIR_FUSIONS

export FILE_GENOME
export CFG_REFERENCE_LENGTH
export HRD_REF_WIG

export CFG_REFERENCE_CAPTUREREGIONS
export CFG_REFERENCE_CAPTUREGENES
export CFG_REFERENCE_COVEREDREGION
export CFG_REFERENCE_CAPTUREREGIONNAME
export CFG_REFERENCE_CAPTURECORFACTORS
export CFG_REFERENCE_COVERED_EXONS
export CFG_REFERENCE_ACTIONABLEGENES

export CFG_REFERENCE_DBSNP

export CFG_COMMON_CPUCORES
export CFG_COMMON_MEMORY

export CFG_SAMTOOLS_MPILEUP_MINMQ
export CFG_SAMTOOLS_MPILEUP_ADJUSTMQ
export CFG_SAMTOOLS_MPILEUP_MAXDEPTH
export CFG_SAMTOOLS_MPILEUP_ILLUMINA
export CFG_SAMTOOLS_MPILEUP_COUNTORPHANS
export CFG_SAMTOOLS_MPILEUP_NOBAQ
export CFG_SAMTOOLS_MPILEUP_REDOBAQ
export CFG_SAMTOOLS_MPILEUP_EXCLUDERG
export CFG_SAMTOOLS_MPILEUP_POSITION
export CFG_SAMTOOLS_MPILEUP_REGION
export CFG_SAMTOOLS_MPILEUP_IGNORERG
export CFG_SAMTOOLS_MPILEUP_IGNOREOVERLAPS
export CFG_SAMTOOLS_FLAGS_MPILEUP

export CFG_GENERAL_MINBASEQUAL
export CFG_GENERAL_MAFCUTOFF
export CFG_GENERAL_MINVAF

export CFG_VARSCAN_SOMATIC_MINCOVERAGE
export CFG_VARSCAN_SOMATIC_TUMORPURITY
export CFG_VARSCAN_SOMATIC_MINFREQFORHOM
export CFG_VARSCAN_SOMATIC_MINCOVERAGENORMAL
export CFG_VARSCAN_SOMATIC_MINCOVERAGETUMOR
export CFG_VARSCAN_SOMATIC_NORMALPURITY
export CFG_VARSCAN_SOMATIC_PVALUE
export CFG_VARSCAN_SOMATIC_SOMATICPVALUE
export CFG_VARSCAN_SOMATIC_VALIDATION

export CFG_VARSCAN_PROCESSSOMATIC_MAXNORMALFREQ
export CFG_VARSCAN_PROCESSSOMATIC_PVALUE

export CFG_VARSCAN_FPFILTER_MINVARCOUNT
export CFG_VARSCAN_FPFILTER_MINVARCOUNTLC
export CFG_VARSCAN_FPFILTER_MAXSOMATICP
export CFG_VARSCAN_FPFILTER_MAXSOMATICPDEPTH
export CFG_VARSCAN_FPFILTER_MINREFREADPOS
export CFG_VARSCAN_FPFILTER_MINVARREADPOS
export CFG_VARSCAN_FPFILTER_MINREFDIST3
export CFG_VARSCAN_FPFILTER_MINVARDIST3
export CFG_VARSCAN_FPFILTER_MINSTRANDEDNESS
export CFG_VARSCAN_FPFILTER_MINSTRANDREADS
export CFG_VARSCAN_FPFILTER_MAXBASEQUALDIFF
export CFG_VARSCAN_FPFILTER_MINREFAVGRL
export CFG_VARSCAN_FPFILTER_MINVARAVGRL
export CFG_VARSCAN_FPFILTER_MAXRLDIFF
export CFG_VARSCAN_FPFILTER_MAXREFMMQS
export CFG_VARSCAN_FPFILTER_MAXVARMMQS
export CFG_VARSCAN_FPFILTER_MINMMQSDIFF
export CFG_VARSCAN_FPFILTER_MAXMMQSDIFF
export CFG_VARSCAN_FPFILTER_MINREFMAPQUAL
export CFG_VARSCAN_FPFILTER_MINVARMAPQUAL
export CFG_VARSCAN_FPFILTER_MAXMAPQUALDIFF

export CFG_PANEL_SAMTOOLS_MPILEUP_MINMQ
export CFG_PANEL_SAMTOOLS_MPILEUP_ADJUSTMQ
export CFG_PANEL_SAMTOOLS_MPILEUP_MAXDEPTH
export CFG_PANEL_SAMTOOLS_MPILEUP_ILLUMINA
export CFG_PANEL_SAMTOOLS_MPILEUP_COUNTORPHANS
export CFG_PANEL_SAMTOOLS_MPILEUP_NOBAQ
export CFG_PANEL_SAMTOOLS_MPILEUP_REDOBAQ
export CFG_PANEL_SAMTOOLS_MPILEUP_EXCLUDERG
export CFG_PANEL_SAMTOOLS_MPILEUP_POSITION
export CFG_PANEL_SAMTOOLS_MPILEUP_REGION
export CFG_PANEL_SAMTOOLS_MPILEUP_IGNORERG
export CFG_PANEL_SAMTOOLS_MPILEUP_IGNOREOVERLAPS
export CFG_PANEL_SAMTOOLS_FLAGS_MPILEUP

export CFG_PANEL_VARSCAN_MPILEUP2SNP_MINCOVERAGE
export CFG_PANEL_VARSCAN_MPILEUP2SNP_MINREADS2
export CFG_PANEL_VARSCAN_MPILEUP2SNP_MINFREQFORHOM
export CFG_PANEL_VARSCAN_MPILEUP2SNP_PVALUE
export CFG_PANEL_VARSCAN_MPILEUP2SNP_STRANDFILTER

export CFG_PANEL_VARSCAN_MPILEUP2INDEL_MINCOVERAGE
export CFG_PANEL_VARSCAN_MPILEUP2INDEL_MINREADS2
export CFG_PANEL_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM
export CFG_PANEL_VARSCAN_MPILEUP2INDEL_PVALUE
export CFG_PANEL_VARSCAN_MPILEUP2INDEL_STRANDFILTER

export CFG_PANEL_VARSCAN_FPFILTER_MINVARCOUNT
export CFG_PANEL_VARSCAN_FPFILTER_MINVARCOUNTLC
export CFG_PANEL_VARSCAN_FPFILTER_MAXSOMATICP
export CFG_PANEL_VARSCAN_FPFILTER_MAXSOMATICPDEPTH
export CFG_PANEL_VARSCAN_FPFILTER_MINREFREADPOS
export CFG_PANEL_VARSCAN_FPFILTER_MINVARREADPOS
export CFG_PANEL_VARSCAN_FPFILTER_MINREFDIST3
export CFG_PANEL_VARSCAN_FPFILTER_MINVARDIST3
export CFG_PANEL_VARSCAN_FPFILTER_MINSTRANDEDNESS
export CFG_PANEL_VARSCAN_FPFILTER_MINSTRANDREADS
export CFG_PANEL_VARSCAN_FPFILTER_MAXBASEQUALDIFF
export CFG_PANEL_VARSCAN_FPFILTER_MINREFAVGRL
export CFG_PANEL_VARSCAN_FPFILTER_MINVARAVGRL
export CFG_PANEL_VARSCAN_FPFILTER_MAXRLDIFF
export CFG_PANEL_VARSCAN_FPFILTER_MAXREFMMQS
export CFG_PANEL_VARSCAN_FPFILTER_MAXVARMMQS
export CFG_PANEL_VARSCAN_FPFILTER_MINMMQSDIFF
export CFG_PANEL_VARSCAN_FPFILTER_MAXMMQSDIFF
export CFG_PANEL_VARSCAN_FPFILTER_MINREFMAPQUAL
export CFG_PANEL_VARSCAN_FPFILTER_MINVARMAPQUAL
export CFG_PANEL_VARSCAN_FPFILTER_MAXMAPQUALDIFF

export CFG_TUMORONLY_SAMTOOLS_MPILEUP_MINMQ
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_ADJUSTMQ
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_MAXDEPTH
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_ILLUMINA
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_COUNTORPHANS
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_NOBAQ
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_REDOBAQ
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_EXCLUDERG
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_POSITION
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_REGION
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_IGNORERG
export CFG_TUMORONLY_SAMTOOLS_MPILEUP_IGNOREOVERLAPS
export CFG_TUMORONLY_SAMTOOLS_FLAGS_MPILEUP

export CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINCOVERAGE
export CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINREADS2
export CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINFREQFORHOM
export CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_PVALUE
export CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_STRANDFILTER

export CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINCOVERAGE
export CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINREADS2
export CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM
export CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_PVALUE
export CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_STRANDFILTER

export CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARCOUNT
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARCOUNTLC
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXSOMATICP
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXSOMATICPDEPTH
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFREADPOS
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARREADPOS
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFDIST3
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARDIST3
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINSTRANDEDNESS
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINSTRANDREADS
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXBASEQUALDIFF
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFAVGRL
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARAVGRL
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXRLDIFF
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXREFMMQS
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXVARMMQS
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINMMQSDIFF
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXMMQSDIFF
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFMAPQUAL
export CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARMAPQUAL
export CFG_TUMORONLY_VARSCAN_FPFILTER_MAXMAPQUALDIFF

export CFG_PANEL_MUTECT_CALLABLEDEPTH
export CFG_TUMORONLY_MUTECT_CALLABLEDEPTH
export CFG_MUTECT_CALLABLEDEPTH

export CFG_FUSION_GENES
export CFG_AMPLIFICATION_GENES

export CFG_ANNOVAR_PROTOCOL
export CFG_ANNOVAR_ARGOP

export BIN_JAVA

export BIN_FASTQC
export BIN_TRIM
export DIR_TRIMMOMATIC_ADAPTER
export BIN_CUT

export BIN_BWAMEM

export BIN_BAM_READCOUNT

export BIN_SAMTOOLS
export BIN_SAMVIEW
export BIN_SAMSORT
export BIN_SAMRMDUP
export BIN_SAMINDEX
export BIN_MPILEUP
export BIN_STATS

export BIN_GATK
export BIN_REALIGNER_TARGER_CREATOR
export BIN_INDEL_REALIGNER
export BIN_BASE_RECALIBRATOR
export BIN_PRINT_READS

export BIN_GATK4

export BIN_FIX_MATE

export BIN_VAR_SCAN
export BIN_SOMATIC
export BIN_PROCESSSOMATIC

export DIR_ANNOVAR
export DIR_ANNOVAR_DATA
export CONVERT2ANNOVAR2
export CONVERT2ANNOVAR3
export CONVERT2ANNOVAR
export TABLEANNOVAR

export BIN_COVERAGE

export BIN_SNPEFF

export BIN_FREEC
export FILE_REFERENCE_MAPPABILITY

export BIN_RSCRIPT

export FUSIONCATCHER_DB
export BIN_FUSIONCATCHER

export MSISENSOR2
export MSISENSOR_PRO
export MSISENSOR_PRO_SCAN
export MICROSATELLITE_SITES

export SEQUENZA_UTILS
export SEQUENZA_WINDOW
export SEQUENZA_NON_MATCHING_NORMAL
export SEQUENZA_CHROMOSOMES
