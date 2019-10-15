#!/usr/bin/env bash

###########################################
## WES Pipeline for somatic and germline ##
###########################################
# script to run the actual analysis

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

while getopts d:h option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
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

cd "${DIR_ANALYSIS}" || exit 1

# TODO: refactor R-Script to use different folders and also add DIR_REF (to link to Target file)
${BIN_RSCRIPT} "${DIR_RSCRIPT}"/Main.R "${CFG_CASE}" "${PARAM_DIR_PATIENT}" "${CFG_FILE_GERMLINE}" "${CFG_FILE_TUMOR}" \
  "${DIR_TARGET}" "${DIR_RSCRIPT}" "${DIR_DATABASE}" "${CFG_REFERENCE_CAPTUREGENES}

${BIN_RSCRIPT} --vanilla -e "load('${DIR_ANALYSIS}/WES.Rdata'); library(knitr); knit('${DIR_RSCRIPT}/Report.Rnw')"

pdflatex -interaction=nonstopmode Report.tex
pdflatex -interaction=nonstopmode Report.tex
