#!/usr/bin/env bash

# common settings
readonly CFG_AUTHOR=$(get_config_value common.author "${PARAM_DIR_PATIENT}")
readonly CFG_CENTER=$(get_config_value common.center "${PARAM_DIR_PATIENT}")
readonly CFG_PROTOCOL=$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")

# temporary folder
readonly DIR_TMP="$(get_config_value common.dirTmp "${PARAM_DIR_PATIENT}")/${PARAM_DIR_PATIENT}"

readonly CFG_FILE_TUMOR_R1=$(get_config_value common.files.tumor_R1 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_TUMOR_R2=$(get_config_value common.files.tumor_R2 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE_R1=$(get_config_value common.files.germline_R1 "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE_R2=$(get_config_value common.files.germline_R2 "${PARAM_DIR_PATIENT}")

# folder containing patient output
readonly DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
readonly DIR_WES="${DIR_TARGET}/WES"
readonly DIR_ANALYSES="${DIR_TARGET}/Analyses"

# end paths

## Genome
readonly FILE_GENOME="${DIR_REF}/genome/$(get_config_value reference.genome "${PARAM_DIR_PATIENT}")"

readonly CFG_REFERENCE_LENGTH="${DIR_CHROMOSOMES}/$(get_config_value reference.length "${PARAM_DIR_PATIENT}")"

# depending on measurement machine
## SureSelect (Capture Kit)
readonly CFG_REFERENCE_CAPTUREREGIONS="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.captureRegions "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTUREGENES="${DIR_SEQUENCING}/$(get_config_value reference.sequencing.captureGenes "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_COVEREDREGION="$(get_config_value reference.sequencing.coveredRegion "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTUREREGIONNAME="$(get_config_value reference.sequencing.captureRegionName "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTURECORFACTORS="${DIR_DATABASE}/$(get_config_value reference.sequencing.captureCorFactors "${PARAM_DIR_PATIENT}")"


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

# general tool parameters
readonly CFG_GENERAL_MINBASEQUAL=$(get_config_value tools.general.minBaseQual "${PARAM_DIR_PATIENT}")
readonly CFG_GENERAL_MAFCUTOFF=$(get_config_value tools.general.maf_cutoff "${PARAM_DIR_PATIENT}")

# sametools mpileup
readonly CFG_SAMTOOLS_MPILEUP_MINMQ=$(get_config_value tools.samtools.mpileup.minMQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_ADJUSTMQ=$(get_config_value tools.samtools.mpileup.adjustMQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_MAXDEPTH=$(get_config_value tools.samtools.mpileup.maxDepth "${PARAM_DIR_PATIENT}")

# mpileup flags
readonly CFG_SAMTOOLS_MPILEUP_ILLUMINA=$(get_config_value tools.samtools.mpileup.illumina "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_COUNTORPHANS=$(get_config_value tools.samtools.mpileup.countOrphans "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_NOBAQ=$(get_config_value tools.samtools.mpileup.noBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_REDOBAQ=$(get_config_value tools.samtools.mpileup.redoBAQ "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_EXCLUDERG=$(get_config_value tools.samtools.mpileup.excludeRG "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_POSITION=$(get_config_value tools.samtools.mpileup.position "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_REGION=$(get_config_value tools.samtools.mpileup.region "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_IGNORERG=$(get_config_value tools.samtools.mpileup.ignoreRG "${PARAM_DIR_PATIENT}")
readonly CFG_SAMTOOLS_MPILEUP_IGNOREOVERLAPS=$(get_config_value tools.samtools.mpileup.ignoreOverlaps "${PARAM_DIR_PATIENT}")

CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_ILLUMINA,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -6"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_COUNTORPHANS,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -A"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_NOBAQ,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -B"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_REDOBAQ,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -E"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_EXCLUDERG,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -G"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_POSITION,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -l"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_REGION,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -r"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_IGNORERG,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -R"
CFG_SAMTOOLS_FLAGS_MPILEUP=[[ "${CFG_SAMTOOLS_MPILEUP_IGNOREOVERLAPS,,}" == "true" ]] && CFG_SAMTOOLS_FLAGS_MPILEUP="${CFG_SAMTOOLS_FLAGS_MPILEUP} -x"


# global Varscan parameters
readonly CFG_VARSCAN_MINVAF=$(get_config_value tools.varscan.minVAF "${PARAM_DIR_PATIENT}")

# VarScan somatic
readonly CFG_VARSCAN_SOMATIC_MINCOVERAGE=$(get_config_value tools.varscan.somatic.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_TUMORPURITY=$(get_config_value tools.varscan.somatic.tumorPurity "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_MINFREQFORHOM=$(get_config_value tools.varscan.somatic.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_MINCOVERAGENORMAL=$(get_config_value tools.varscan.somatic.minCoverageNormal "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_MINCOVERAGETUMOR=$(get_config_value tools.varscan.somatic.minCoverageTumor "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_NORMALPURITY=$(get_config_value tools.varscan.somatic.normalPurity "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_PVALUE=$(get_config_value tools.varscan.somatic.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_SOMATICPVALUE=$(get_config_value tools.varscan.somatic.somaticPValue "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_STRANDFILTER=$(get_config_value tools.varscan.somatic.strandFilter "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_VALIDATION=$(get_config_value tools.varscan.somatic.validation "${PARAM_DIR_PATIENT}")

# Varscan processSomatic
readonly CFG_VARSCAN_PROCESSSOMATIC_MAXNORMALFREQ=$(get_config_value tools.varscan.processSomatic.maxNormalFreq "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PROCESSSOMATIC_PVALUE=$(get_config_value tools.varscan.processSomatic.pValue "${PARAM_DIR_PATIENT}")

# VarScan fpfilter
readonly CFG_VARSCAN_FPFILTER_MINVARCOUNT=$(get_config_value tools.varscan.fpfilter.minVarCount "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARCOUNTLC=$(get_config_value tools.varscan.fpfilter.minVarCountLC "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXSOMATICP=$(get_config_value tools.varscan.fpfilter.maxSomaticP "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXSOMATICPDEPTH=$(get_config_value tools.varscan.fpfilter.maxSomaticPDepth "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFREADPOS=$(get_config_value tools.varscan.fpfilter.minRefReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARREADPOS=$(get_config_value tools.varscan.fpfilter.minVarReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFDIST3=$(get_config_value tools.varscan.fpfilter.minRefDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARDIST3=$(get_config_value tools.varscan.fpfilter.minVarDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINSTRANDEDNESS=$(get_config_value tools.varscan.fpfilter.minStrandedness "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINSTRANDREADS=$(get_config_value tools.varscan.fpfilter.minStrandReads "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXBASEQUALDIFF=$(get_config_value tools.varscan.fpfilter.maxBasequalDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFAVGRL=$(get_config_value tools.varscan.fpfilter.minRefAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARAVGRL=$(get_config_value tools.varscan.fpfilter.minVarAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXRLDIFF=$(get_config_value tools.varscan.fpfilter.maxRlDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXREFMMQS=$(get_config_value tools.varscan.fpfilter.maxRefMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXVARMMQS=$(get_config_value tools.varscan.fpfilter.maxVarMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINMMQSDIFF=$(get_config_value tools.varscan.fpfilter.minMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXMMQSDIFF=$(get_config_value tools.varscan.fpfilter.maxMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINREFMAPQUAL=$(get_config_value tools.varscan.fpfilter.minRefMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MINVARMAPQUAL=$(get_config_value tools.varscan.fpfilter.minVarMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_FPFILTER_MAXMAPQUALDIFF=$(get_config_value tools.varscan.fpfilter.maxMapQualDiff "${PARAM_DIR_PATIENT}")

# Panel Parameter for Varscan
readonly CFG_PANEL_MINVAF=$(get_config_value tools.varscan.panel.minVAF "${PARAM_DIR_PATIENT}")
readonly CFG_PANEL_MINBASEQUAL=$(get_config_value tools.varscan.panel.minBaseQual "${PARAM_DIR_PATIENT}")

# mpileup2snp
readonly CFG_VARSCAN_MPILEUP2SNP_MINCOVERAGE=$(get_config_value tools.varscan.panel.mpileup2snp.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2SNP_MINREADS2=$(get_config_value tools.varscan.panel.mpileup2snp.minReads2 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2SNP_MINFREQFORHOM=$(get_config_value tools.varscan.panel.mpileup2snp.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2SNP_PVALUE=$(get_config_value tools.varscan.panel.mpileup2snp.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2SNP_STRANDFILTER=$(get_config_value tools.varscan.panel.mpileup2snp.strandFilter "${PARAM_DIR_PATIENT}")

# mpileup2indel
readonly CFG_VARSCAN_MPILEUP2INDEL_MINCOVERAGE=$(get_config_value tools.varscan.panel.mpileup2indel.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2INDEL_MINREADS2=$(get_config_value tools.varscan.panel.mpileup2indel.minReads2 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM=$(get_config_value tools.varscan.panel.mpileup2indel.minFreqForHom "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2INDEL_PVALUE=$(get_config_value tools.varscan.panel.mpileup2indel.pValue "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MPILEUP2INDEL_STRANDFILTER=$(get_config_value tools.varscan.panel.mpileup2indel.strandFilter "${PARAM_DIR_PATIENT}")

readonly CFG_VARSCAN_PANEL_FPFILTER_MINVARCOUNT=$(get_config_value tools.varscan.panel.fpfilter.minVarCount "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINVARCOUNTLC=$(get_config_value tools.varscan.panel.fpfilter.minVarCountLC "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXSOMATICP=$(get_config_value tools.varscan.panel.fpfilter.maxSomaticP "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXSOMATICPDEPTH=$(get_config_value tools.varscan.panel.fpfilter.maxSomaticPDepth "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINREFREADPOS=$(get_config_value tools.varscan.panel.fpfilter.minRefReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINVARREADPOS=$(get_config_value tools.varscan.panel.fpfilter.minVarReadpos "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINREFDIST3=$(get_config_value tools.varscan.panel.fpfilter.minRefDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINVARDIST3=$(get_config_value tools.varscan.panel.fpfilter.minVarDist3 "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINSTRANDEDNESS=$(get_config_value tools.varscan.panel.fpfilter.minStrandedness "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINSTRANDREADS=$(get_config_value tools.varscan.panel.fpfilter.minStrandReads "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXBASEQUALDIFF=$(get_config_value tools.varscan.panel.fpfilter.maxBasequalDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINREFAVGRL=$(get_config_value tools.varscan.panel.fpfilter.minRefAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINVARAVGRL=$(get_config_value tools.varscan.panel.fpfilter.minVarAVGRL "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXRLDIFF=$(get_config_value tools.varscan.panel.fpfilter.maxRlDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXREFMMQS=$(get_config_value tools.varscan.panel.fpfilter.maxRefMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXVARMMQS=$(get_config_value tools.varscan.panel.fpfilter.maxVarMMQS "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINMMQSDIFF=$(get_config_value tools.varscan.panel.fpfilter.minMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXMMQSDIFF=$(get_config_value tools.varscan.panel.fpfilter.maxMMQSDiff "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINREFMAPQUAL=$(get_config_value tools.varscan.panel.fpfilter.minRefMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MINVARMAPQUAL=$(get_config_value tools.varscan.panel.fpfilter.minVarMapQual "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_PANEL_FPFILTER_MAXMAPQUALDIFF=$(get_config_value tools.varscan.panel.fpfilter.maxMapQualDiff "${PARAM_DIR_PATIENT}")

# ANNOVAR Databases
readonly CFG_ANNOVAR_PROTOCOL=$(get_config_value tools.annovar.protocol "${PARAM_DIR_PATIENT}")
readonly CFG_ANNOVAR_ARGOP=$(get_config_value tools.annovar.argop "${PARAM_DIR_PATIENT}")

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

# PICARD
readonly BIN_FIX_MATE="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/picard/picard.jar FixMateInformation "

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

# cnvkit
readonly FILE_FLAT_REFERENCE="${DIR_REF}/cnvkit/$(get_config_value tools.cnvkit.flatReference "${PARAM_DIR_PATIENT}")"
readonly BIN_CNVKIT=$(command -v cnvkit.py)

# R
readonly BIN_RSCRIPT=$(command -v Rscript)

# export parameters
export CFG_AUTHOR
export CFG_CENTER
export CFG_PROTOCOL

export CFG_FILE_TUMOR_R1
export CFG_FILE_TUMOR_R2
export CFG_FILE_GERMLINE_R1
export CFG_FILE_GERMLINE_R2

export DIR_TARGET
export DIR_WES
export DIR_ANALYSES

export FILE_GENOME
export CFG_REFERENCE_LENGTH

export CFG_REFERENCE_CAPTUREREGIONS
export CFG_REFERENCE_CAPTUREGENES
export CFG_REFERENCE_COVEREDREGION
export CFG_REFERENCE_CAPTUREREGIONNAME
export CFG_REFERENCE_CAPTURECORFACTORS

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

export CFG_VARSCAN_MINVAF

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

export CFG_PANEL_MINVAF
export CFG_PANEL_MINBASEQUAL

export CFG_VARSCAN_MPILEUP2SNP_MINCOVERAGE
export CFG_VARSCAN_MPILEUP2SNP_MINREADS2
export CFG_VARSCAN_MPILEUP2SNP_MINFREQFORHOM
export CFG_VARSCAN_MPILEUP2SNP_PVALUE
export CFG_VARSCAN_MPILEUP2SNP_STRANDFILTER

export CFG_VARSCAN_MPILEUP2INDEL_MINCOVERAGE
export CFG_VARSCAN_MPILEUP2INDEL_MINREADS2
export CFG_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM
export CFG_VARSCAN_MPILEUP2INDEL_PVALUE
export CFG_VARSCAN_MPILEUP2INDEL_STRANDFILTER

export CFG_VARSCAN_PANEL_FPFILTER_MINVARCOUNT
export CFG_VARSCAN_PANEL_FPFILTER_MINVARCOUNTLC
export CFG_VARSCAN_PANEL_FPFILTER_MAXSOMATICP
export CFG_VARSCAN_PANEL_FPFILTER_MAXSOMATICPDEPTH
export CFG_VARSCAN_PANEL_FPFILTER_MINREFREADPOS
export CFG_VARSCAN_PANEL_FPFILTER_MINVARREADPOS
export CFG_VARSCAN_PANEL_FPFILTER_MINREFDIST3
export CFG_VARSCAN_PANEL_FPFILTER_MINVARDIST3
export CFG_VARSCAN_PANEL_FPFILTER_MINSTRANDEDNESS
export CFG_VARSCAN_PANEL_FPFILTER_MINSTRANDREADS
export CFG_VARSCAN_PANEL_FPFILTER_MAXBASEQUALDIFF
export CFG_VARSCAN_PANEL_FPFILTER_MINREFAVGRL
export CFG_VARSCAN_PANEL_FPFILTER_MINVARAVGRL
export CFG_VARSCAN_PANEL_FPFILTER_MAXRLDIFF
export CFG_VARSCAN_PANEL_FPFILTER_MAXREFMMQS
export CFG_VARSCAN_PANEL_FPFILTER_MAXVARMMQS
export CFG_VARSCAN_PANEL_FPFILTER_MINMMQSDIFF
export CFG_VARSCAN_PANEL_FPFILTER_MAXMMQSDIFF
export CFG_VARSCAN_PANEL_FPFILTER_MINREFMAPQUAL
export CFG_VARSCAN_PANEL_FPFILTER_MINVARMAPQUAL
export CFG_VARSCAN_PANEL_FPFILTER_MAXMAPQUALDIFF

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

export FILE_FLAT_REFERENCE
export BIN_CNVKIT

export BIN_RSCRIPT
