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
  echo "  -s  sex              specify sex (${possible_sex})"
  echo "  -c  case             specify case (${possible_cases})"
  echo "  -n  num              specify num"
  echo "  -fg filename         specify filename germline without file extension"
  echo "  -ft filename         specify filename tumor without file extension"
  echo "  -g folder            specify folder containing germline"
  echo "  -t folder            specify folder containing tumor"
  echo "  -h                   show this help screen"
  exit 1
}

dir_tumor=$1
dir_germline=$2
sex=$3 # XX or XY
filename_germline=$4
filename_tumor=$5
method=$6 # somatic or somaticGermline
ID=$7

while getopts h option; do
  case "${option}" in
  id) ID=$OPTARG ;;
  m) method=$OPTARG ;;
  n) num=$OPTARG ;;
  s) sex=$OPTARG ;;
  g) dir_germline=$OPTARG ;;
  t) dir_tumor=$OPTARG ;;
  fg) filename_germline=$OPTARG ;;
  ft) filename_tumor=$OPTARG ;;
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

if [[ ! " ${possible_sex[@]} " =~ " ${sex} " ]]; then
  # whatever you want to do when arr contains value
  echo "unknown sex: ${sex}"
  echo "use one of the following values: ${possible_sex}"
  exit 1
fi

${SCRIPT_PATH}/createSet.sh -c ${method} ${ID} ${dir_tumor} ${dir_germline} \
  ${filename_germline} ${filename_tumor} ${sex}

cd ${method}_${ID}

./run_jobs.sh >out &
