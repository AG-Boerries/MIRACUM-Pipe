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
if [[ "$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")" = "panel" ]]; then
  readonly CFG_CASE=panelTumor
fi

##################################################################################################################

## load programs
# shellcheck source=programs.cfg.sh
. "${DIR_SCRIPT}/programs.cfg.sh"

##################################################################################################################

[[ -d "${DIR_RNA}" ]] || mkdir -p "${DIR_RNA}"
[[ -d "${DIR_FUSIONS}" ]] || mkdir -p "${DIR_FUSIONS}"

readonly input="${DIR_INPUT}/${PARAM_DIR_PATIENT}/${CFG_FOLDER_RNA}"

# SAMPLE
# readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_${PARAM_TASK}

# FASTQC
for f in ${input}/*.fastq.gz; do
    ${BIN_FASTQC} "${f}" -o "${DIR_RNA}"
done

# Pseudoalignment
#	${KALLISTO} -o ${Out} --plaintext -t 12 ${InputPath}/${fastq1} ${InputPath}/${fastq2}

# Fusioncatcher
	${BIN_FUSIONCATCHER} -i "${input}" -o "${DIR_FUSIONS}"
	