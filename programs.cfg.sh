#!/usr/bin/env bash

 # temporary folder
readonly DIR_TMP="$(get_config_value common.dirTmp "${PARAM_DIR_PATIENT}")/${PARAM_DIR_PATIENT}"

readonly CFG_FILE_TUMOR=$(get_config_value common.files.tumor "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE=$(get_config_value common.files.germline "${PARAM_DIR_PATIENT}")

# folder containing patient output
readonly DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
readonly DIR_WES="${DIR_TARGET}/WES"
readonly DIR_ANALYSIS="${DIR_TARGET}/Analysis"
readonly DIR_RSCRIPT="${DIR_SCRIPT}/RScripts"
readonly DIR_DATABASE="${DIR_SCRIPT}/databases"

# end paths

## Genome
readonly FILE_GENOME="${DIR_REF}/Genome/$(get_config_value reference.genome "${PARAM_DIR_PATIENT}")"

readonly CFG_REFERENCE_LENGTH="${DIR_CHROMOSOMES}/$(get_config_value reference.length "${PARAM_DIR_PATIENT}")"

# depending on measurement machine
## SureSelect (Capture Kit)
readonly CFG_REFERENCE_CAPTUREREGIONS="${DIR_REF}/$(get_config_value reference.sequencing.captureRegions "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_CAPTUREGENES="${DIR_REF}/$(get_config_value reference.sequencing.captureGenes "${PARAM_DIR_PATIENT}")"
readonly CFG_REFERENCE_COVEREDREGION="$(get_config_value reference.sequencing.captureRegion)"


# database for known variants
## dbSNP vcf File
readonly CFG_REFERENCE_DBSNP="${DIR_DBSNP}/$(get_config_value reference.dbSNP "${PARAM_DIR_PATIENT}")"
# END variables

# if no parallel computation, set processes to 1
if [[ -z ${PARALLEL_PROCESSES} ]]; then
  readonly TMP_PROCESSES=1
else
  readonly TMP_PROCESSES="{PARALLEL_PROCESSES}"
fi

# take cpucores/2
readonly CFG_COMMON_CPUCORES=$(($(get_config_value common.cpucores "${PARAM_DIR_PATIENT}")/${TMP_PROCESSES}))

readonly tmp_memory=$(get_config_value common.memory "${PARAM_DIR_PATIENT}")
# take memory/2

readonly CFG_COMMON_MEMORY="$(("${tmp_memory//[^0-9.]/}"/${TMP_PROCESSES}))${tmp_memory//[^a-zA-Z]/}"
readonly CFG_VARSCAN_MINBASEQUAL=$(get_config_value varscan.minBaseQual "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_MINVAF=$(get_config_value varscan.minVAF "${PARAM_DIR_PATIENT}")

# VarScan somatic
readonly CFG_VARSCAN_SOMATIC_MINCOVERAGE=$(get_config_value varscan.somatic.minCoverage "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_TUMORPURITY=$(get_config_value varscan.somatic.tumorPurity "${PARAM_DIR_PATIENT}")
readonly CFG_VARSCAN_SOMATIC_MINFREQFORHOM=$(get_config_value varscan.somatic.minFreqForHom "${PARAM_DIR_PATIENT}")

# VarScan fpfilter
readonly CFG_VARSCAN_FPLFILTER_MINVARCOUNT=$(get_config_value varscan.fpfilter.minVarCount "${PARAM_DIR_PATIENT}")

# ANNOVAR Databases
readonly CFG_ANNOVAR_PROTOCOL=$(get_config_value annovar.protocol "${PARAM_DIR_PATIENT}")
readonly CFG_ANNOVAR_ARGOP=$(get_config_value annovar.argop "${PARAM_DIR_PATIENT}")

## Tools and paths
# Paths
readonly BIN_JAVA="java -Xcompressedrefs -Djava.io.tmpdir=${DIR_TMP} " # path to java

# Pre-Processing
readonly BIN_FASTQC="${DIR_TOOLS}/FastQC/bin/fastqc -t ${CFG_COMMON_CPUCORES} --extract "

readonly BIN_TRIM="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/Trimmomatic/trimmomatic.jar PE -threads ${CFG_COMMON_CPUCORES} -phred33 "
readonly DIR_TRIMMOMATIC_ADAPTER="${DIR_TOOLS}/Trimmomatic/adapters"
readonly BIN_CUT="cut -f1,2,3"

# Alignment
readonly BIN_BWAMEM="${DIR_TOOLS}/bwa/bwa mem -M "

# BAM-Readcount
readonly BIN_BAM_READCOUNT="${DIR_TOOLS}/bam-readcount/bin/bam-readcount -q 1 -b 20 -w 1 -f ${FILE_GENOME} "

# SAMTOOLS
readonly BIN_SAMTOOLS="${DIR_TOOLS}/samtools/samtools" # path to samtools
readonly BIN_SAMVIEW="${BIN_SAMTOOLS} view -@ ${CFG_COMMON_CPUCORES} "
readonly BIN_SAMSORT="${BIN_SAMTOOLS} sort -@ ${CFG_COMMON_CPUCORES} "
readonly BIN_SAMRMDUP="${BIN_SAMTOOLS} rmdup "
readonly BIN_SAMINDEX="${BIN_SAMTOOLS} index "
readonly BIN_MPILEUP="${BIN_SAMTOOLS} mpileup -B -C 50 -f ${FILE_GENOME} -q 1 --min-BQ ${CFG_VARSCAN_MINBASEQUAL}"
readonly BIN_STATS="${BIN_SAMTOOLS} stats "

# GATK
readonly BIN_GATK="java -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/gatk/GenomeAnalysisTK.jar"
readonly BIN_REALIGNER_TARGER_CREATOR="${BIN_GATK} -T RealignerTargetCreator -R ${FILE_GENOME} -nt ${CFG_COMMON_CPUCORES} "
readonly BIN_INDEL_REALIGNER="${BIN_GATK} -R ${FILE_GENOME} -T IndelRealigner "
readonly BIN_BASE_RECALIBRATOR="${BIN_GATK} -T BaseRecalibrator -l INFO -R ${FILE_GENOME} -knownSites ${CFG_REFERENCE_DBSNP} -nct ${CFG_COMMON_CPUCORES} "
readonly BIN_PRINT_READS="${BIN_GATK} -T PrintReads -R ${FILE_GENOME} -nct ${CFG_COMMON_CPUCORES} "

# PICARD
readonly BIN_FIX_MATE="java -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/picard/picard.jar FixMateInformation "

# VARSCAN
readonly BIN_VAR_SCAN="java -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/varscan/VarScan.jar"
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
readonly BIN_SNPEFF="${BIN_JAVA} -Xmx${CFG_COMMON_MEMORY} -jar ${DIR_TOOLS}/snpEff/snpEff.jar GRCh37.75 -c ${DIR_TOOLS}/snpEff/snpEff.config -canon -v"

# FREEC
readonly BIN_FREEC="${DIR_TOOLS}/FREEC/bin/freec "

readonly FILE_REFERENCE_MAPPABILITY="${DIR_REF}/mappability/$(get_config_value reference.mappability "${PARAM_DIR_PATIENT}")"

# R
readonly BIN_RSCRIPT=$(command -v Rscript)

export CFG_FILE_TUMOR
export CFG_FILE_GERMLINE

export DIR_TARGET
export DIR_WES
export DIR_ANALYSIS
export DIR_RSCRIPT
export DIR_DATABASE


export FILE_GENOME
export CFG_REFERENCE_LENGTH

export CFG_REFERENCE_CAPTUREREGIONS
export CFG_REFERENCE_COVEREDREGION

export CFG_REFERENCE_DBSNP

export CFG_COMMON_CPUCORES
export CFG_COMMON_MEMORY
export CFG_VARSCAN_MINBASEQUAL
export CFG_VARSCAN_MINVAF

export CFG_VARSCAN_SOMATIC_MINCOVERAGE
export CFG_VARSCAN_SOMATIC_TUMORPURITY
export CFG_VARSCAN_SOMATIC_MINFREQFORHOM

export CFG_VARSCAN_FPLFILTER_MINVARCOUNT

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

export BIN_RSCRIPT
