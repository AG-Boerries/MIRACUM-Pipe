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

cd "${DIR_ANALYSES}" || exit 1

# execute Main.R, i.e. complete analyses in R
${BIN_RSCRIPT} "${DIR_RSCRIPT}/Main.R" "${CFG_CASE}" "${PARAM_DIR_PATIENT}" "${CFG_FILE_GERMLINE_R1}" "${CFG_FILE_TUMOR_R1}" \
  "${DIR_TARGET}" "${DIR_RSCRIPT}" "${DIR_DATABASE}" "${CFG_REFERENCE_CAPTUREGENES}" "${CFG_REFERENCE_COVEREDREGION}" \
  "${CFG_AUTHOR}" "${CFG_CENTER}" "${CFG_REFERENCE_CAPTUREREGIONS}" "${CFG_REFERENCE_CAPTUREREGIONNAME}" "${FILE_GENOME}" \
  "${CFG_REFERENCE_CAPTURECORFACTORS}" "${CFG_GENERAL_MINVAF}" "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARCOUNT}" "${CFG_GENERAL_MAFCUTOFF}" \
  "${CFG_REFERENCE_ACTIONABLEGENES}" "${CFG_REFERENCE_COVERED_EXONS}" "${CFG_ENTITY}" "${CFG_SEX}" "${CFG_FUSION_GENES}" \
  "${CFG_AMPLIFICATION_GENES}" "${CFG_UCSC_SERVER}" "${CFG_CNV_ANNOTATION}"

# translate to tex
${BIN_RSCRIPT} --vanilla -e "load('${DIR_ANALYSES}/Report.RData'); library(knitr); knit('${DIR_RSCRIPT}/Report_tumorOnly.Rnw');"

# PDF report
mv "${DIR_ANALYSES}/Report_tumorOnly.tex" "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.tex"
pdflatex -interaction=nonstopmode "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.tex" \
  --output-directory="${DIR_ANALYSES}"
pdflatex -interaction=nonstopmode "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.tex" \
  --output-directory="${DIR_ANALYSES}"

# clean up
# remove aux files which are created while pdflatex
rm -f "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.aux" \
      "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.toc" \
      "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.log" \
      "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.out"
