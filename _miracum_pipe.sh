#!/usr/bin/env bash

###########################################
## WES Pipeline for somatic and germline ##
###########################################
# script to run the actual analysis
# Version 31.07.2019

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)

## load settings
# shellcheck source=common.cfg.sh
. "${DIR_SCRIPT}"/common.cfg.sh


function usage() {
  echo "usage: miracum_pipe.sh -d dir [-h]"
  echo "  -d  dir             specify relative folder of patient"
  echo "  -t  task            specify task"
  echo "  -h                  show this help screen"
  exit 1
}

while getopts d:t:h option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
  t) readonly PARAM_TASK=$OPTARG ;;
  h) usage ;;
  \?)
    echo "Unknown option: -$OPTARG" >&2
    exit 1
    ;;
  :)
    echo "Missing option argument for -$OPTARG" >&2
    exit 1
    ;;
  *)
    echo "Unimplemented option: -$OPTARG" >&2
    exit 1
    ;;
  esac
done

# if no patient is defined
if [[ -z "${PARAM_DIR_PATIENT}" ]]; then
  echo "no patient defined."
  echo "--"
  usage
fi


# load patient yaml
readonly CFG_SEX=$(get_config_value sex "${PARAM_DIR_PATIENT}")
if [[ "$(get_config_value annotation.germline "${PARAM_DIR_PATIENT}")" = "True" ]]; then
  readonly CFG_CASE=somaticGermline
else
  readonly CFG_CASE=somatic
fi

# check inputs
readonly VALID_TASKS=("GD TD VC CNV Report")
readonly VALID_SEXES=("XX XY")

for value in "${VALID_TASKS[@]}"
do
  [[ "${PARAM_TASK}" = "${value}" ]] && \
    echo "unknown task: ${PARAM_TASK}" && \
    echo "use one of the following values: $(join_by ' ' ${VALID_TASKS})" && \
    exit 1
done

for value in "${VALID_SEXES[@]}"
do
  [[ "${CFG_SEX}" = "${value}" ]] && \
    echo "unknown sex: ${CFG_SEX}" && \
    echo "use one of the following values: $(join_by ' ' ${VALID_SEXES})" && \
    exit 1
done

##################################################################################################################
#### Parameters which have to be adjusted accoridng the the environment or the users needs

readonly CFG_FILE_TUMOR=$(get_config_value common.files.tumor "${PARAM_DIR_PATIENT}")
readonly CFG_FILE_GERMLINE=$(get_config_value common.files.germline "${PARAM_DIR_PATIENT}")

# folder containing patient output
readonly DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
readonly DIR_WES="${DIR_TARGET}/WES"
readonly DIR_ANALYSIS="${DIR_TARGET}/Analysis"
readonly DIR_RSCRIPT="${DIR_SCRIPT}/RScripts"
readonly DIR_DATABASE="${DIR_SCRIPT}/Databases"

# end paths

## Genome
readonly FILE_GENOME="${DIR_REF}/Genome/$(get_config_value reference.genome "${PARAM_DIR_PATIENT}")"

readonly CFG_REFERENCE_LENGTH="${DIR_CHROMOSOMES}/$(get_config_value reference.length "${PARAM_DIR_PATIENT}")"

# depending on measurement machine
## SureSelect (Capture Kit)
readonly CFG_REFERENCE_CAPTUREREGIONS="${DIR_REF}/$(get_config_value reference.sequencing.capture_regions "${PARAM_DIR_PATIENT}")"

# database for known variants
## dbSNP vcf File
readonly CFG_REFERENCE_DBSNP="${DIR_DBSNP}/$(get_config_value reference.dbSNP "${PARAM_DIR_PATIENT}")"
# END variables

# take cpucores/2
readonly CFG_COMMON_CPUCORES=$(($(get_config_value common.cpucores "${PARAM_DIR_PATIENT}")/2))
readonly CFG_COMMON_MEMORY=$(get_config_value common.memory "${PARAM_DIR_PATIENT}")
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

##################################################################################################################

##########
## MAIN ##
##########

# TODO: required in all 4 files
echo "STARTING ${CFG_CASE}_${PARAM_DIR_PATIENT}_${PARAM_TASK}";
touch "${DIR_TARGET}/.STARTING_MARKER_${PARAM_TASK}"

### program calls
case ${PARAM_TASK} in

## alignment -----------------------------------------------------------------------------------------------------
# TODO: make_alignment.sh
GD | TD)
  if [[ ! -d "${DIR_TMP}" ]]; then
    mkdir -p "${DIR_TMP}"
  fi

  # SAMPLE
  readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_${PARAM_TASK}

  if [[ ${PARAM_TASK} = GD ]]; then
    readonly FILE_FASTQ_1="${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_FILE_GERMLINE}1.fastq.gz"
    readonly FILE_FASTQ_2="${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_FILE_GERMLINE}2.fastq.gz"
  else
    readonly FILE_FASTQ_1="${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_FILE_TUMOR}1.fastq.gz"
    readonly FILE_FASTQ_2="${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_FILE_TUMOR}2.fastq.gz"
  fi

  # temp files
  readonly fastq_o1_p_t=${DIR_TMP}/${NameD}_output1_paired_trimmed.fastq.gz
  readonly fastq_o1_u_t=${DIR_TMP}/${NameD}_output1_unpaired_trimmed.fastq.gz
  readonly fastq_o2_p_t=${DIR_TMP}/${NameD}_output2_paired_trimmed.fastq.gz
  readonly fastq_o2_u_t=${DIR_TMP}/${NameD}_output2_unpaired_trimmed.fastq.gz
  readonly bam=${DIR_TMP}/${NameD}_output.bam
  readonly prefixsort=${DIR_TMP}/${NameD}_output.sort
  readonly sortbam=${DIR_TMP}/${NameD}_output.sort.bam
  readonly rmdupbam=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.bam
  readonly bai=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.bai
  readonly bamlist=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.bam.list
  readonly realignedbam=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.bam
  readonly realignedbai=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.bai
  readonly fixedbam=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.bam
  readonly fixedbai=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.bai
  readonly csv=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.recal_data.csv

  recalbam=${DIR_WES}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
  readonly statstxt=${DIR_WES}/${NameD}_stats.txt
  readonly coveragetxt=${DIR_WES}/${NameD}_coverage.all.txt

  # fastqc zip to WES
  ${BIN_FASTQC} "${FILE_FASTQ_1}" -o "${DIR_WES}"
  ${BIN_FASTQC} "${FILE_FASTQ_2}" -o "${DIR_WES}"

  # trim fastq
  ${BIN_TRIM} "${FILE_FASTQ_1}" "${FILE_FASTQ_2}" "${fastq_o1_p_t}" "${fastq_o1_u_t}" "${fastq_o2_p_t}" "${fastq_o2_u_t}" \
    ILLUMINACLIP:"${DIR_TRIMMOMATIC_ADAPTER}"/TruSeq3-PE-2.fa:2:30:10 HEADCROP:3 TRAILING:10 MINLEN:25

  ${BIN_FASTQC} "${fastq_o1_p_t}" -o "${DIR_WES}"
  ${BIN_FASTQC} "${fastq_o2_p_t}" -o "${DIR_WES}"

  # make bam
  ${BIN_BWAMEM} -R "@RG\tID:${NameD}\tSM:${NameD}\tPL:illumina\tLB:lib1\tPU:unit1" -t 12 "${FILE_GENOME}" \
    "${fastq_o1_p_t}" "${fastq_o2_p_t}" | ${BIN_SAMVIEW} -bS - >"${bam}"

  # stats
  ${BIN_STATS} "${bam}" >"${statstxt}"

  # sort bam
  ${BIN_SAMSORT} "${bam}" -T "${prefixsort}" -o "${sortbam}"

  # rmdup bam
  ${BIN_SAMVIEW} -b -f 0x2 -q1 "${sortbam}" | ${BIN_SAMRMDUP} - "${rmdupbam}"

  # make bai
  ${BIN_SAMINDEX} "${rmdupbam}" "${bai}"

  # make bam list
  ${BIN_REALIGNER_TARGER_CREATOR} -o "${bamlist}" -I "${rmdupbam}"

  # realign bam
  ${BIN_INDEL_REALIGNER} -I "${rmdupbam}" -targetIntervals "${bamlist}" -o "${realignedbam}"

  # fix bam
  ${BIN_FIX_MATE} INPUT="${realignedbam}" OUTPUT="${fixedbam}" SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

  # make csv
  ${BIN_BASE_RECALIBRATOR} -I "${fixedbam}" \
    -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o "${csv}"

  # recal bam
  ${BIN_PRINT_READS} -I "${fixedbam}" -BQSR "${csv}" -o "${recalbam}"

  # coverage
  ${BIN_COVERAGE} -b "${recalbam}" -a "${CFG_REFERENCE_CAPTUREREGIONS}" | grep '^all' >"${coveragetxt}"

  # zip
  ${BIN_FASTQC} "${recalbam}" -o "${DIR_WES}"
  ;;
# eo alignment

## variantCalling ------------------------------------------------------------------------------------------------
# TODO: make_vc.sh
VC)
  if [[ ! -d "${DIR_TMP}" ]]; then
    mkdir -p "${DIR_TMP}"
  fi

  readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_${PARAM_TASK}
  readonly NameGD=${CFG_CASE}_${PARAM_DIR_PATIENT}_GD
  readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_TD
  readonly recalbamGD=${DIR_WES}/${NameGD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
  readonly recalbamTD=${DIR_WES}/${NameTD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
  readonly snpvcf=${DIR_WES}/${NameD}.output.snp.vcf
  readonly indelvcf=${DIR_WES}/${NameD}.output.indel.vcf

  ${BIN_MPILEUP} "${recalbamGD}" "${recalbamTD}" | ${BIN_SOMATIC} --output-snp "${snpvcf}" --output-indel "${indelvcf}" \
    --min-coverage "${CFG_VARSCAN_SOMATIC_MINCOVERAGE}" --tumor-purity "${CFG_VARSCAN_SOMATIC_TUMORPURITY}" \
    --min-var-freq "${CFG_VARSCAN_MINVAF}" --min-freq-for-hom "${CFG_VARSCAN_SOMATIC_MINFREQFORHOM}" \
    --min-avg-qual "${CFG_VARSCAN_MINBASEQUAL}" --output-vcf 1 --mpileup 1

  # Processing of somatic mutations

  ${BIN_PROCESSSOMATIC} "${snpvcf}" --min-tumor-freq "${CFG_VARSCAN_MINVAF}"
  ${BIN_PROCESSSOMATIC} "${indelvcf}" --min-tumor-freq "${CFG_VARSCAN_MINVAF}"

  # FP Filter:  snp.Somatic.hc snp.LOH.hc snp.Germline.hc
  # FP Filter:  indel.Somatic.hc indel.LOH.hc indel.Germline.hc

  readonly names1="snp indel"
  for name1 in ${names1}; do

    if [[ "${CFG_CASE}" = somatic ]]; then
      names2="Somatic LOH"
    else
      names2="Somatic LOH Germline"
    fi

    for name2 in ${names2}; do
      hc_vcf=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.vcf
      hc_avi=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.avinput
      hc_rci=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.readcount.input
      hc_rcs=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.readcounts
      hc_fpf=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.fpfilter.vcf
      if [[ "${name2}" = Somatic ]]; then
        recalbam=${recalbamTD}
      else
        recalbam=${recalbamGD}
      fi
      ${CONVERT2ANNOVAR2} "${hc_avi}" "${hc_vcf}"
      ${BIN_CUT} "${hc_avi}" > "${hc_rci}"
      ${BIN_BAM_READCOUNT} -l "${hc_rci}" "${recalbam}" > "${hc_rcs}"
      ${BIN_VAR_SCAN} fpfilter "${hc_vcf}" "${hc_rcs}" --output-file "${hc_fpf}" --keep-failures 1 \
        --min-ref-basequal "${CFG_VARSCAN_MINBASEQUAL}" --min-var-basequal "${CFG_VARSCAN_MINBASEQUAL}" \
        --min-var-count "${CFG_VARSCAN_FPLFILTER_MINVARCOUNT}" --min-var-freq "${CFG_VARSCAN_MINVAF}"
    done
  done

  readonly data=${DIR_WES}
  for name1 in ${names1}; do
    # Annotation snp.Somatic.hc $data/NameD.output.snp.Somatic.hc.fpfilter.vcf
    # Annotation indel.Somatic.hc $data/NameD.output.indel.Somatic.hc.fpfilter.vcf
    hc_=${data}/${NameD}.output.${name1}.Somatic.hc
    hc_fpf=${data}/${NameD}.output.${name1}.Somatic.hc.fpfilter.vcf
    hc_T_avi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput
    hc_T_avi_multi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput.hg19_multianno.csv
    ${CONVERT2ANNOVAR} "${hc_}" "${hc_fpf}" -allsample
    ${TABLEANNOVAR} "${hc_T_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
        -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

    hc_snpeff="${data}/${NameD}.output.${name1}.Somatic.SnpEff.vcf"
    ${BIN_SNPEFF} "${hc_fpf}" > "${hc_snpeff}"

    if [[ "${CFG_CASE}" = somaticGermline ]]; then
      # Annotation snp.Germline.hc $data/NameD.output.snp.Germline.hc.fpfilter.vcf
      # Annotation indel.Germline.hc $data/NameD.output.indel.Germline.hc.fpfilter.vcf
      hc_=${data}/${NameD}.output.${name1}.Germline.hc
      hc_fpf=${data}/${NameD}.output.${name1}.Germline.hc.fpfilter.vcf
      hc_N_avi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput
      hc_N_avi_multi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput.hg19_multianno.csv
      ${CONVERT2ANNOVAR} "${hc_}" "${hc_fpf}" -allsample
      ${TABLEANNOVAR} "${hc_N_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
          -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

      hc_N_snpeff=${data}/${NameD}.output.${name1}.NORMAL.SnpEff.vcf
      ${BIN_SNPEFF} "${hc_fpf}" > "${hc_N_snpeff}"
    fi

    # Annotation snp.LOH.hc
    # Annotation indel.LOH.hc
    hc_vcf=${data}/${NameD}.output.${name1}.LOH.hc.vcf
    hc_fpf=${data}/${NameD}.output.${name1}.LOH.hc.fpfilter.vcf
    hc_avi=${data}/${NameD}.output.${name1}.LOH.hc.avinput
    hc_avi_multi=${data}/${NameD}.output.${name1}.LOH.hc.avinput.hg19_multianno.csv
    ${CONVERT2ANNOVAR3} "${hc_avi}" "${hc_fpf}"
    ${TABLEANNOVAR} "${hc_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
        -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

    hc_L_snpeff=${data}/${NameD}.output.${name1}.LOH.SnpEff.vcf
    ${BIN_SNPEFF} "${hc_fpf}" > "${hc_L_snpeff}"
  done

  rm -r "${DIR_TMP:?}/*"
  ;;
# eo VC

## CNV  ----------------------------------------------------------------------------------------------------------
# TODO: make_cnv.sh
CNV)
  readonly output="${DIR_WES}/CNV"

  if [[ ! -d "${output}" ]]; then
    mkdir -p "${output}"
  fi

  cat >"${DIR_WES}"/CNV_config.txt <<EOI
[general]

chrFiles = ${DIR_CHROMOSOMES}
chrLenFile = ${CFG_REFERENCE_LENGTH}
breakPointType = 4
breakPointThreshold = 1.2
forceGCcontentNormalization = 1
gemMappabilityFile = ${FILE_REFERENCE_MAPPABILITY}
intercept = 0
minCNAlength = 3
maxThreads = 12
noisyData = TRUE
outputDir = ${output}
ploidy = 2
printNA = FALSE
readCountThreshold = 50
samtools = ${BIN_SAMTOOLS}
sex = ${CFG_SEX}
step = 0
window = 0
uniqueMatch = TRUE
contaminationAdjustment = TRUE

[sample]

mateFile = ${DIR_WES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_TD_output.sort.filtered.rmdup.realigned.fixed.recal.bam
inputFormat = BAM
mateOrientation = FR

[control]

mateFile = ${DIR_WES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_GD_output.sort.filtered.rmdup.realigned.fixed.recal.bam
inputFormat = BAM
mateOrientation = FR

[target]

captureRegions = ${CFG_REFERENCE_CAPTUREREGIONS}
EOI

  export PATH=${PATH}:${BIN_SAMTOOLS}
  ${BIN_FREEC} -conf "${DIR_WES}"/CNV_config.txt
  ;;
  # eo CNV

## Report  -------------------------------------------------------------------------------------------------------
# TODO: make_report.sh
Report)
  cd "${DIR_ANALYSIS}" || exit 1

  # TODO: refactor R-Script to use different folders and also add DIR_REF (to link to Target file)
  ${BIN_RSCRIPT} "${DIR_ANALYSIS}"/Main.R "${CFG_CASE}" "${PARAM_DIR_PATIENT}" "${CFG_FILE_GERMLINE}" "${CFG_FILE_TUMOR}" \
    "${DIR_TARGET}" "${DIR_RSCRIPT}" "${DIR_DATABASE}"

  ${BIN_RSCRIPT} -e "library(knitr); knit('${DIR_ANALYSIS}/Report.Rnw')"
  pdflatex -interaction=nonstopmode Report.tex
  pdflatex -interaction=nonstopmode Report.tex
  ;;
  # eo Report
esac
# eo program calls

# TODO: required in all 4 files
echo "task ${PARAM_TASK} for ${PARAM_DIR_PATIENT} finished"
rm "${DIR_TARGET}/.STARTING_MARKER_${PARAM_TASK}"