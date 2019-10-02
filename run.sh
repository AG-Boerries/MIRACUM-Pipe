#!/usr/bin/env bash

# createSet.sh and make_alignment_VC_CNV.sh need to be in the base directory
# first germline (GD), then tumor (TD)

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)

## load settings
# shellcheck source=global.sh
. "${DIR_SCRIPT}"/global.sh

possible_tasks=("GD TD VC CNV Report")

function usage() {
  echo "usage: run.sh [-d dir] [-t task] [-f] [-h]"
  echo "  -t  task             specify task ${possible_tasks}"
  echo "  -f                   specify forced run, i.e. neglecting if a patient already was computed"
  echo "  -d  dir              specify directory in which the patient is located"
  echo "  -h                   show this help screen"
  exit 1
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
      "${DIR_SCRIPT}"/createSet.sh -d "${DIR_PATIENT}"

      if [[ "$(get_config_value annotation.germline "${DIR_PATIENT}")" = "True" ]]; then
        readonly CFG_CASE=somaticGermline
      else
        readonly CFG_CASE=somatic
      fi

      ("${DIR_OUTPUT}/${CFG_CASE}_${DIR_PATIENT}"/run_jobs.sh 
        &> "${DIR_OUTPUT}/${CFG_CASE}_${DIR_PATIENT}/log/run.log") &
    fi
  done

else

  if [[ "$(get_config_value annotation.germline "${PARAM_DIR_PATIENT}")" = "True" ]]; then
    readonly CFG_CASE=somaticGermline
  else
    readonly CFG_CASE=somatic
  fi

  if [[ ! -z ${PARAM_TASK} ]]; then
    # possibility to comfortably run tasks separately
    case "${PARAM_TASK}" in
      cnv)
        # TODO: ${DOCKER_COMMAND} ${DIR_MIRACUM}/make_cnv.sh -d "${dir}"
        echo "not yet implemented"
      ;;

      vc)
        # TODO: ${DOCKER_COMMAND} ${DIR_MIRACUM}/make_vc.sh -d "${dir}"
        echo "not yet implemented"
      ;;
      report)
        # TODO: ${DOCKER_COMMAND} ${DIR_MIRACUM}/make_report.sh -d "${dir}"
        echo "not yet implemented"      
      ;;
    esac
  else
    if [[ ! -f "${PARAM_DIR_PATIENT}/.processed" || "${PARAM_FORCE}" ]]; then
      echo "computing ${PARAM_DIR_PATIENT}"
      "${DIR_SCRIPT}"/createSet.sh -d "${PARAM_DIR_PATIENT}"
      ("${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"/run_jobs.sh \ 
        &> "${DIR_OUTPUT}/${CFG_CASE}_${DIR_PATIENT}/log/run.log" ) &
    fi
  fi
fi
