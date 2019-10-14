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

## load programs
# shellcheck source=programs.cfg.sh
. "${DIR_SCRIPT}"/programs.cfg.sh

##################################################################################################################

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