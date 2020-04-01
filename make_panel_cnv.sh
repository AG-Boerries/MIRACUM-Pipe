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
readonly CFG_PROTOCOL=$(get_config_value common.protocol# "${PARAM_DIR_PATIENT}")
if [[ "${CFG_PROTOCOL}" == "panel" ]]; then
  readonly CFG_CASE=panelTumor
fi
#if [[ "$(get_config_value annotation.germline "${PARAM_DIR_PATIENT}")" = "True" ]]; then
#  readonly CFG_CASE=somaticGermline
#else
#  readonly CFG_CASE=somatic
#fi

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

readonly DIR_CNV_OUTPUT="${DIR_WES}"

[[ -d "${DIR_CNV_OUTPUT}" ]] || mkdir -p "${DIR_CNV_OUTPUT}"

# TODO: install cnvkit
# TODO: install FlatReference
if [ ${sureselect} = "TruSight_Tumor" ]
then
readonly FlatReference="/home/mhess/Panel/FlatReference2.cnn"
fi
if [ ${sureselect} = "TruSight_Amplicon" ]
then
readonly FlatReference="/home/mhess/Panel/FlatReference.cnn"
fi
if [ ${sureselect} = "V6" ]
then
readonly FlatReference="/home/mhess/Panel/FlatReferenceV6.cnn"
fi

readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_cnv
readonly cnr=${DIR_WES}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.cnr
readonly cns=${DIR_WES}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.cns

${BIN_CNVKIT} batch -m amplicon -r ${FlatReference} $WES_DIR}/*.bam
${BIN_CNVKIT} segment ${cnr} -o ${cns} --rscript-path ${BIN_RSCRIPT}
