#!/usr/bin/env bash

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)

## load settings
# shellcheck source=common.cfg.sh
. "${DIR_SCRIPT}"/common.cfg.sh

possible_tasks=("gd td vc cnv report")

function usage() {
  echo "usage: run.sh [-d dir] [-t task] [-f] [-h]"
  echo "  -t  task             specify task ${possible_tasks}"
  echo "  -f                   specify forced run, i.e. neglecting if a patient already was computed"
  echo "  -d  dir              specify directory in which the patient is located"
  echo "  -h                   show this help screen"
  exit 1
}

# run_pipe relative_patient_dir dir_log
function run_pipe() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_log="${dir_target}/log"

  mkdir -p "${dir_log}"

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_alignment.sh -t TD -d "${dir_patient}" &> ${dir_log}/TD.log) &
  ("${DIR_SCRIPT}"/make_alignment.sh -t GD -d "${dir_patient}" &> ${dir_log}/GD.log) &
  wait

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_vc.sh -d "${dir_patient}"  &> ${dir_log}/VC.log) &
  ("${DIR_SCRIPT}"/make_cnv.sh -d "${dir_patient}" &> ${dir_log}/CNV.log) &
  wait

  # create report based on the results of the processes above
  ("${DIR_SCRIPT}"/make_report.sh -d "${dir_patient}" &> ${dir_log}/Report.log) &
}

# get_case relative_patient_dir
function get_case() {
  local dir_patient="${1}"

  if [[ "$(get_config_value annotation.germline "${dir_patient}")" = "True" ]]; then
    echo "somaticGermline"
  else
    echo "somatic"
  fi
}

while getopts d:t:fh option; do
  case "${option}" in
  t) readonly PARAM_TASK=$OPTARG;;
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
  f) readonly PARAM_FORCE=true;;
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

for value in "${possible_tasks[@]}"
do
  [[ "${PARAM_TASK}" = "${value}" ]] && \
    echo "unknown task: ${PARAM_TASK}" && \
    echo "use one of the following values: $(join_by ' ' ${possible_tasks})" && \
    exit 1
done

# create temporary folder if not existent
[[ -d "${DIR_TMP}" ]] || mkdir -p "${DIR_TMP}"

# run script
if [[ -z "${PARAM_DIR_PATIENT}" && -z "${PARAM_TASK}" ]]; then
  for dir in "${DIR_SCRIPT}"/assets/input/*; do
    if [[ ! -f "${PARAM_DIR_PATIENT}"/.processed || ${PARAM_FORCE} ]]; then
      DIR_PATIENT=${dir##*/}
      echo "computing ${DIR_PATIENT}"

      CFG_CASE=$(get_case ${DIR_PATIENT})
      DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${DIR_PATIENT}"

      run_pipe "${DIR_PATIENT}" "${DIR_TARGET}"

      # check if report was generated successfully
      if [[ -f "${DIR_ANALYSIS}/Report.pdf" ]]; then
        touch "${DIR_TARGET}/.processed"
      else
        echo "${DIR_PATIENT} failed"
      fi
    fi
  done

else
  CFG_CASE=$(get_case ${PARAM_DIR_PATIENT})
  DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
  DIR_LOG="${DIR_TARGET}/log"

  mkdir -p "${DIR_LOG}"

  # if already computed, i.e. file .processed exists, only compute again if forced
  if [[ ! -f "${PARAM_DIR_PATIENT}/.processed" || "${PARAM_FORCE}" ]]; then
    if [[ ! -z ${PARAM_TASK} ]]; then
      # possibility to comfortably run tasks separately
      case "${PARAM_TASK}" in
        td)
          "${DIR_SCRIPT}"/make_alignment.sh -t td -d "${dir}"
        ;;
        gd)
          "${DIR_SCRIPT}"/make_alignment.sh -t gd -d "${dir}"
        ;;

        cnv)
          "${DIR_SCRIPT}"/make_cnv.sh -d "${dir}"
        ;;

        vc)
          "${DIR_SCRIPT}"/make_vc.sh -d "${dir}"
        ;;
        report)
          "${DIR_SCRIPT}"/make_report.sh -d "${dir}"
        ;;
      esac
    else
      echo "computing ${PARAM_DIR_PATIENT}"

      CFG_CASE=$(get_case ${PARAM_DIR_PATIENT})
      DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"

      run_pipe "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"

      # check if report was generated successfully
      if [[ -f "${DIR_ANALYSIS}/Report.pdf" ]]; then
        touch "${DIR_TARGET}/.processed"
      else
        echo "${PARAM_DIR_PATIENT} failed"
      fi
    fi
  fi
fi
