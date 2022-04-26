#!/usr/bin/env bash

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
  echo "  -p                  computing as parallel process"
  echo "  -h                  show this help screen"
  exit 1
}

while getopts d:ph option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
  p) readonly PARALLEL_PROCESSES=2 ;;
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
if [[ "$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")" = "tumorOnly" ]]; then
  readonly CFG_CASE=tumorOnly
fi

# check inputs
readonly VALID_SEXES=("XX XY")

if [[ ! " ${VALID_SEXES[@]} " =~ " ${CFG_SEX} " ]]; then
  echo "unknown sex: ${CFG_SEX}"
  echo "use one of the following values: $(join_by ' ' ${VALID_SEXES})"
  exit 1
fi

##################################################################################################################

## load programs
# shellcheck source=programs.cfg.sh
. "${DIR_SCRIPT}/programs.cfg.sh"

##################################################################################################################

readonly DIR_CNV_OUTPUT="${DIR_WES}/CNV"

[[ -d "${DIR_CNV_OUTPUT}" ]] || mkdir -p "${DIR_CNV_OUTPUT}"

readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td
readonly recalbam=${DIR_WES}/${NameTD}_output.sort.rmdup.realigned.fixed.recal.bam
readonly HRD_OUTPUT=${DIR_WES}/${NameD}.seqz.gz
readonly HRD_OUTPUT_SMALL=${DIR_WES}/${NameD}.small.seqz.gz
readonly HRD_INFO=${DIR_WES}/${NameD}.small.seqz._info_seg.txt

# HRD
if [ ! -f "${HRD_REF_WIG}" ]; then
    echo "${HRD_REF_WIG} does not exist. Generating ..."
    ${SEQUENZA_UTILS} gc_wiggle --fasta "${FILE_GENOME}" -w "${SEQUENZA_WINDOW}" -o "${HRD_REF_WIG}"
fi

${SEQUENZA_UTILS} bam2seqz -S "${BIN_SAMTOOLS}" -gc "${HRD_REF_WIG}" --fasta "${FILE_GENOME}" -n "${recalbam}" --tumor "${recalbam}" --normal2 "${SEQUENZA_NON_MATCHING_NORMAL}" \
  -C ${SEQUENZA_CHROMOSOMES} -o "${HRD_OUTPUT}"
${SEQUENZA_UTILS} seqz_binning -s "${HRD_OUTPUT}" -w "${SEQUENZA_WINDOW}" -o "${HRD_OUTPUT_SMALL}"
${BIN_RSCRIPT} "${DIR_RSCRIPT}/HRD.R" "${NameD}" "${DIR_WES}"

# Extract ploidy and purity from sequenza
PLOIDY=$(${BIN_RSCRIPT} "--vanilla" "-e" "cat(round(read.table('${HRD_INFO}', header = T)\$ploidy))")
CONTAMINATION=$(${BIN_RSCRIPT} "--vanilla" "-e" "cat(1-read.table('${HRD_INFO}', header = T)\$cellularity)")

cat >"${DIR_WES}"/CNV_config.txt <<EOI
[general]

chrFiles = ${DIR_CHROMOSOMES}
chrLenFile = ${CFG_REFERENCE_LENGTH}
breakPointType = 4
breakPointThreshold = 1.2
forceGCcontentNormalization = 1
gemMappabilityFile = ${FILE_REFERENCE_MAPPABILITY}
intercept = 1
minCNAlength = 3
maxThreads = ${CFG_COMMON_CPUCORES}
noisyData = TRUE
outputDir = ${DIR_CNV_OUTPUT}
ploidy = ${PLOIDY}
contamination = ${CONTAMINATION}
contaminationAdjustment = TRUE
printNA = FALSE
readCountThreshold = 50
samtools = ${BIN_SAMTOOLS}
bedtools = ${DIR_TOOLS}/bedtools2/bin/bedtools
sex = ${CFG_SEX}
step = 0
window = 0
uniqueMatch = TRUE

[sample]

mateFile = ${DIR_WES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_td_output.sort.rmdup.realigned.fixed.recal.bam
inputFormat = BAM
mateOrientation = FR

[BAF]

makePileup = ${CFG_REFERENCE_DBSNP}
fastaFile = ${FILE_GENOME}
minimalCoveragePerPosition = ${CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINCOVERAGE}
minimalQualityPerPosition = ${CFG_GENERAL_MINBASEQUAL}
SNPfile = ${CFG_REFERENCE_DBSNP}

[target]

captureRegions = ${CFG_REFERENCE_CAPTUREREGIONS}
EOI

export PATH="${BIN_SAMTOOLS}:${PATH}"
${BIN_FREEC} -conf "${DIR_WES}"/CNV_config.txt