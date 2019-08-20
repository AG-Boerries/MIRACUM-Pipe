#!/usr/bin/env bash

# createSet.sh and make_alignment_VC_CNV.sh need to be in the base directory
# first germline (GD), then tumor (TD)

SCRIPT_PATH=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd -P
)

possible_cases=("somatic somaticGermline")
possible_sex=("XX XY")

function usage() {
  echo "usage: createSet.sh -s sex -c case -n num -fg filename -ft filename -g folder -t folder [-h]"
  echo "  -d  dir              specify directory in which the patients are located"
  echo "  -h                   show this help screen"
  exit 1
}

#dir_tumor=$1
#dir_germline=$2
#sex=$3 # XX or XY
#filename_germline=$4
#filename_tumor=$5
#method=$6 # somatic or somaticGermline
#ID=$7

# in patient's config.yaml:
#  echo "  -s  sex              specify sex (${possible_sex})"
#  echo "  -c  case             specify case (${possible_cases})"

while getopts dh option; do
  case "${option}" in
  d) patient_dir=$OPTARG ;;
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

if [[ -z ${patient_dir} ]]; then
  echo "no patient dir given"
  exit 1
fi

#if [[ ! " ${possible_sex[@]} " =~ " ${sex} " ]]; then
#  echo "unknown sex: ${sex}"
#  echo "use one of the following values: ${possible_sex}"
#  exit 1
#fi

# TODO: finalize call with patient dir as param and nothing else :)
for patient_dir in ${SCRIPT_PATH}/input/*; do
  if [[ ! -f ${patient_dir}/.processed ]]; then
    ${SCRIPT_PATH}/createSet.sh -p ${patient_dir}

    (${SCRIPT_PATH}/assets/input/${patient_dir}/run_jobs.sh > out && touch .processed) &
  fi
done




