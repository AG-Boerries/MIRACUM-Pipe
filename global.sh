#!/usr/bin/env bash

# This function gets configuration parameters from yaml files
# first trying a look-up inside the patient specific config and
# falls back to the global config if requested config key not available
#
#   usage: get_config_value some.key relative_patient_path
#
function get_config_value()
{
  SCRIPT_PATH=$(
    cd "$(dirname "${BASH_SOURCE[0]}")" || exit
    pwd -P
  )

  # first look inside the patient specific config
  patient_conf="${SCRIPT_PATH}"/assets/input/"${2}"/config.yaml
  [[ -f "${patient_conf}" ]] && value=$(shyaml get-value "${1}" < "${patient_conf}" 2> /dev/null)

  # if value not available in patient specific config, take global config's value
  if [[ -z "${value}" ]]; then
    value=$(shyaml get-value "${1}" < "${SCRIPT_PATH}"/conf/default.yaml)
  fi

  echo "${value}"
}

function join_by { local IFS="$1"; shift; echo "$*"; }

## Directories
## General
DIR_ASSETS="${DIR_SCRIPT}/assets"

DIR_INPUT="${DIR_ASSETS}/input"             # folder contatining the raw data (.fastq files)
DIR_OUTPUT="${DIR_ASSETS}/output"
DIR_REF="${DIR_ASSETS}/references"        # reference genome
DIR_TMP="/tmp" # temporary folder
