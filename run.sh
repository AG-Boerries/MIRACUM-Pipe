#!/usr/bin/env bash

# createSet.sh and make_alignment_VC_CNV.sh need to be in the base directory
# first germline (GD), then tumor (TD)

DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd -P
)

DIR_OUTPUT="${DIR_SCRIPT}/assets/output"

possible_tasks=("GD TD VC CNV Report")

function usage() {
  echo "usage: run.sh [-d dir] [-t task] [-f] [-h]"
  echo "  -t  task             specify task ${possible_tasks}"
  echo "  -f                   specify forced run, i.e. neglecting if a patient already was computed"
  echo "  -d  dir              specify directory in which the patient is located"
  echo "  -h                   show this help screen"
  exit 1
}

while getopts d:t:h option; do
  case "${option}" in
  t) task=$OPTARG;;
  d) DIR_PATIENT=$OPTARG ;;
  f) force=true;;
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
  [[ "${task}" = "${value}" ]] && \
    echo "unknown task: ${task}" && \
    echo "use one of the following values: $(join_by ' ' ${possible_tasks})" && \
    exit 1
done


if [[ -z "${dir}" && -z "${task}" ]]; then
  for DIR_PATIENT in ${DIR_SCRIPT}/input/*; do
    if [[ ! -f ${DIR_PATIENT}/.processed || ${force} ]]; then
      ${DIR_SCRIPT}/createSet.sh -d "${DIR_PATIENT##*/}"

      (${DIR_OUTPUT}/${case}_${DIR_PATIENT##*/}/run_jobs.sh > out && touch .processed) &
    fi
  done

#  ${DOCKER_COMMAND} ${DIR_MIRACUM}/run.sh
else
  # possibility to comfortably run tasks separately
  case ${task} in
    cnv)
      # ${DOCKER_COMMAND} ${DIR_MIRACUM}/make_cnv.sh -d "${dir}"
      echo "not yet implemented"
    ;;

    vc)
      # ${DOCKER_COMMAND} ${DIR_MIRACUM}/make_vc.sh -d "${dir}"
      echo "not yet implemented"
    ;;
    report)
      # ${DOCKER_COMMAND} ${DIR_MIRACUM}/make_report.sh -d "${dir}"
      echo "not yet implemented"
    ;;

    default)
      if [[ ! -f ${DIR_PATIENT}/.processed || ${force} ]]; then
        ${DIR_OUTPUT}/createSet.sh -d "${DIR_PATIENT##*/}"
        (${DIR_OUTPUT}/${case}_${DIR_PATIENT##*/}/run_jobs.sh > out && touch .processed) &
      fi
    ;;
  esac
fi
