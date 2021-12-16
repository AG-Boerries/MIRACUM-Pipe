#!/usr/bin/env bash

readonly DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
  pwd -P
)

## load settings
# shellcheck source=common.cfg.sh
. "${DIR_SCRIPT}"/common.cfg.sh

readonly VALID_TASKS=("gd td vc cnv report td_gd_parallel vc_cnv_parallel rna_fusions")
readonly VALID_PROTOCOLS=("wes panel tumorOnly")

function usage() {
  echo "usage: miracum_pipe.sh [-d dir] [-t task] [-p protocol] [-f] [-h]"
  echo "  -t  task             specify task ${VALID_TASKS}"
  echo "  -p  protocol         specify protocol ${VALID_PROTOCOLS}"
  echo "  -f                   specify forced run, i.e. neglecting if a patient already was computed"
  echo "  -d  dir              specify directory in which the patient is located"
  echo "  -s                   computing sequentially; might be required if ressources are limited"
  echo "  -h                   show this help screen"
  exit 1
}


# setup relative_patient_dir dir_target
function setup() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_tmp="$(get_config_value common.dirTmp "${dir_patient}")/${dir_patient}"
  mkdir -p ${dir_tmp}

  local dir_log="${dir_target}/log"
  mkdir -p "${dir_log}"
}

# cleanup relative_patient_dir
function cleanup() {
  local dir_patient="${1}"
  local dir_tmp="$(get_config_value common.dirTmp "${dir_patient}")/${dir_patient}"
  rm -rf "${dir_tmp:?}"
}

# WES protocol
# run_pipe relative_patient_dir dir_target
function run_pipe() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_log="${dir_target}/log"

  setup "${dir_patient}" "${dir_target}"

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_alignment.sh -p -t td -d "${dir_patient}" &> "${dir_log}/td.log") &
  ("${DIR_SCRIPT}"/make_alignment.sh -p -t gd -d "${dir_patient}" &> "${dir_log}/gd.log") &
  wait

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_vc.sh  -p -d "${dir_patient}" &> "${dir_log}/vc.log") &
  ("${DIR_SCRIPT}"/make_cnv.sh -p -d "${dir_patient}" &> "${dir_log}/cnv.log") &
  wait

  # create report based on the results of the processes above
  ("${DIR_SCRIPT}"/make_report.sh -d "${dir_patient}" &> "${dir_log}/report.log")

  cleanup "${dir_patient}"
}

# run_pipe_seq relative_patient_dir dir_target
function run_pipe_seq() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_log="${dir_target}/log"

  setup "${dir_patient}" "${dir_target}"

  ("${DIR_SCRIPT}"/make_alignment.sh -t td -d "${dir_patient}" &> "${dir_log}/td.log")
  ("${DIR_SCRIPT}"/make_alignment.sh -t gd -d "${dir_patient}" &> "${dir_log}/gd.log")
  ("${DIR_SCRIPT}"/make_vc.sh  -d "${dir_patient}" &> "${dir_log}/vc.log")
  ("${DIR_SCRIPT}"/make_cnv.sh -d "${dir_patient}" &> "${dir_log}/cnv.log")
  ("${DIR_SCRIPT}"/make_report.sh -d "${dir_patient}" &> "${dir_log}/report.log")

  cleanup "${dir_patient}"
}


# Panel protocol
# run_pipe_panel relative_patient_dir dir_target
function run_panel_pipe() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_log="${dir_target}/log"

  setup "${dir_patient}" "${dir_target}"

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_panel_alignment.sh -t td -d "${dir_patient}" &> "${dir_log}/td.log") &
  ("${DIR_SCRIPT}"/make_rna_fusions.sh -d "${dir_patient}" &> "${dir_log}/fusions.log") &
  wait

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_panel_vc.sh -d "${dir_patient}" &> "${dir_log}/vc.log") &
  ("${DIR_SCRIPT}"/make_panel_cnv.sh -p -d "${dir_patient}" &> "${dir_log}/cnv.log") &
  wait

  # create report based on the results of the processes above
  ("${DIR_SCRIPT}"/make_panel_report.sh -d "${dir_patient}" &> "${dir_log}/report.log")

  cleanup "${dir_patient}"
}

# run_pipe_seq relative_patient_dir dir_target
function run_panel_pipe_seq() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_log="${dir_target}/log"

  setup "${dir_patient}" "${dir_target}"

  ("${DIR_SCRIPT}"/make_panel_alignment.sh -t td -d "${dir_patient}" &> "${dir_log}/td.log")
  ("${DIR_SCRIPT}"/make_panel_vc.sh  -d "${dir_patient}" &> "${dir_log}/vc.log")
  ("${DIR_SCRIPT}"/make_panel_cnv.sh -d "${dir_patient}" &> "${dir_log}/cnv.log")
  ("${DIR_SCRIPT}"/make_rna_fusions.sh -d "${dir_patient}" &> "${dir_log}/fusions.log")
  ("${DIR_SCRIPT}"/make_panel_report.sh -d "${dir_patient}" &> "${dir_log}/report.log")

  cleanup "${dir_patient}"
}

# tumorOnly protocol
# run_pipe_panel relative_patient_dir dir_target
function run_tumorOnly_pipe() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_log="${dir_target}/log"

  setup "${dir_patient}" "${dir_target}"

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_tumorOnly_alignment.sh -t td -d "${dir_patient}" &> "${dir_log}/td.log")

  # use parallel shell scripting
  ("${DIR_SCRIPT}"/make_tumorOnly_vc.sh -d "${dir_patient}" &> "${dir_log}/vc.log") &
  ("${DIR_SCRIPT}"/make_tumorOnly_cnv.sh -p -d "${dir_patient}" &> "${dir_log}/cnv.log") &
  wait

  # create report based on the results of the processes above
  ("${DIR_SCRIPT}"/make_tumorOnly_report.sh -d "${dir_patient}" &> "${dir_log}/report.log")

  cleanup "${dir_patient}"
}

# run_pipe_seq relative_patient_dir dir_target
function run_tumorOnly_pipe_seq() {
  local dir_patient="${1}"
  local dir_target="${2}"

  local dir_log="${dir_target}/log"

  setup "${dir_patient}" "${dir_target}"

  ("${DIR_SCRIPT}"/make_tumorOnly_alignment.sh -t td -d "${dir_patient}" &> "${dir_log}/td.log")
  ("${DIR_SCRIPT}"/make_tumorOnly_vc.sh  -d "${dir_patient}" &> "${dir_log}/vc.log")
  ("${DIR_SCRIPT}"/make_tumorOnly_cnv.sh -d "${dir_patient}" &> "${dir_log}/cnv.log")
  ("${DIR_SCRIPT}"/make_tumorOnly_report.sh -d "${dir_patient}" &> "${dir_log}/report.log")

  cleanup "${dir_patient}"
}

# get_case relative_patient_dir
function get_case() {
  local dir_patient="${1}"

  if [[ "$(get_config_value common.protocol "${dir_patient}")" = "wes" ]]; then
    if [[ "$(get_config_value common.germline "${dir_patient}")" = "True" ]]; then
      echo "somaticGermline"
    else
      echo "somatic"
    fi
  fi

  if [[ "$(get_config_value common.protocol "${dir_patient}")" = "panel" ]]; then
    echo "panelTumor"
  fi

  if [[ "$(get_config_value common.protocol "${dir_patient}")" = "tumorOnly" ]]; then
    echo "tumorOnly"
  fi
}

while getopts d:t:p:fhs option; do
  case "${option}" in
  t) readonly PARAM_TASK=$OPTARG;;
  p) readonly PARAM_PROTOCOL=$OPTARG;;
  d) readonly PARAM_DIR_PATIENT=$OPTARG;;
  f) readonly PARAM_FORCE=true;;
  s) readonly PARAM_SEQ=true;;
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

# run script
if [[ ! -z "${PARAM_PROTOCOL}" ]]; then
  if [[ ! " ${VALID_PROTOCOLS[@]} " =~ " ${PARAM_PROTOCOL} " ]]; then
    echo "unknown protocol: ${PARAM_PROTOCOL}"
    echo "use one of the following values: $(join_by ' ' ${VALID_PROTOCOLS})"
    exit 1
  fi

  if [[ !("${PARAM_FORCE}") && -e "${DIR_TARGET}/.processed" ]]; then
    echo "Analyses already exist. If you want to re-analyze it, it has to be forced (-f)."
  fi

  if [[ -z "${PARAM_DIR_PATIENT}" && -z "${PARAM_TASK}" && -n "${PARAM_PROTOCOL}" ]]; then
    for dir in "${DIR_SCRIPT}"/assets/input/*; do
      # get relative dir patient
      DIR_PATIENT=${dir##*/}
    
      # estimate output dir
      CFG_CASE=$(get_case ${DIR_PATIENT})
      DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${DIR_PATIENT}"
      DIR_ANALYSES="${DIR_TARGET}/Analyses"

      if [[ "${PARAM_FORCE}" || ! -f "${DIR_TARGET}/.processed" ]]; then
        echo "computing ${DIR_PATIENT}"
        echo "${PARAM_PROTOCOL}"

        if [[ "${PARAM_PROTOCOL}" == "wes" ]]; then
          echo "WES Protocol"

          if [[ -z "${PARAM_SEQ}" ]]; then
            run_pipe "${DIR_PATIENT}" "${DIR_TARGET}"
          else
            run_pipe_seq "${DIR_PATIENT}" "${DIR_TARGET}"
          fi
          # check if report was generated successfully
          if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${DIR_PATIENT}_Report.pdf" ]]; then
            touch "${DIR_TARGET}/.processed"
            echo "${DIR_PATIENT} finished"
          else
            echo "${DIR_PATIENT} failed"
          fi
        fi

        if [[ "${PARAM_PROTOCOL}" == "panel" ]]; then
          echo "Panel Protocol"

          if [[ -z "${PARAM_SEQ}" ]]; then
            run_panel_pipe "${DIR_PATIENT}" "${DIR_TARGET}"
          else
            run_panel_pipe_seq "${DIR_PATIENT}" "${DIR_TARGET}"
          fi
          # check if report was generated successfully
          if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${DIR_PATIENT}_Report_Panel.pdf" ]]; then
            touch "${DIR_TARGET}/.processed"
            echo "${DIR_PATIENT} finished"
          else
            echo "${DIR_PATIENT} failed"
          fi
        fi

        if [[ "${PARAM_PROTOCOL}" == "tumorOnly" ]]; then
          echo "tumorOnly Protocol"

          if [[ -z "${PARAM_SEQ}" ]]; then
            run_tumorOnly_pipe "${DIR_PATIENT}" "${DIR_TARGET}"
          else
            run_tumorOnly_pipe_seq "${DIR_PATIENT}" "${DIR_TARGET}"
          fi
          # check if report was generated successfully
          if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${DIR_PATIENT}_Report_tumorOnly.pdf" ]]; then
            touch "${DIR_TARGET}/.processed"
            echo "${DIR_PATIENT} finished"
          else
            echo "${DIR_PATIENT} failed"
          fi
        fi
      fi
    done

  else
    CFG_CASE=$(get_case ${PARAM_DIR_PATIENT})
    DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
    DIR_LOG="${DIR_TARGET}/log"

    readonly DIR_TMP="$(get_config_value common.dirTmp "${PARAM_DIR_PATIENT}")/${PARAM_DIR_PATIENT}"

    # create temporary folder if not existent
    [[ -d "${DIR_TMP}" ]] || mkdir -p "${DIR_TMP}"
    [[ -d "${DIR_LOG}" ]] || mkdir -p "${DIR_LOG}"
    
    # estimate output dir
    CFG_CASE=$(get_case ${PARAM_DIR_PATIENT})
    DIR_TARGET="${DIR_OUTPUT}/${CFG_CASE}_${PARAM_DIR_PATIENT}"
    DIR_ANALYSES="${DIR_TARGET}/Analyses"
    
    # if already computed, i.e. file .processed exists, only compute again if forced
    if [[ "${PARAM_FORCE}" || ! -f "${DIR_TARGET}/.processed" ]]; then
      if [[ ! -z ${PARAM_TASK} ]]; then
        if [[ ! " ${VALID_TASKS[@]} " =~ " ${PARAM_TASK} " ]]; then
          echo "unknown task: ${PARAM_TASK}"
          echo "use one of the following values: $(join_by ' ' ${VALID_TASKS})"
          exit 1
        fi

        # setup processing
        setup "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"

        # possibility to comfortably run tasks separately
        case "${PARAM_PROTOCOL}" in
          wes)
            case "${PARAM_TASK}" in
              td)
                "${DIR_SCRIPT}"/make_alignment.sh -t td -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/td.log"
              ;;
              gd)
                "${DIR_SCRIPT}"/make_alignment.sh -t gd -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/gd.log"
              ;;

              cnv)
		"${DIR_SCRIPT}"/make_cnv.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/cnv.log"
              ;;

              vc)
                "${DIR_SCRIPT}"/make_vc.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/vc.log"
              ;;

              report)
                "${DIR_SCRIPT}"/make_report.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/report.log"
                if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report.pdf" ]]; then
                  touch "${DIR_TARGET}/.processed"
                  echo "${PARAM_DIR_PATIENT} Report finished"
                else
                  echo "${PARAM_DIR_PATIENT} failed"
                fi
              ;;

              td_gd_parallel)
                ("${DIR_SCRIPT}"/make_alignment.sh -p -t td -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/td.log") &
                ("${DIR_SCRIPT}"/make_alignment.sh -p -t gd -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/gd.log") &
                wait
              ;;

              vc_cnv_parallel)
                ("${DIR_SCRIPT}"/make_vc.sh  -p -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/vc.log") &
                ("${DIR_SCRIPT}"/make_cnv.sh -p -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/cnv.log") &
                wait
              ;;
            esac
          ;;

          panel)
            case "${PARAM_TASK}" in
              td)
                "${DIR_SCRIPT}"/make_panel_alignment.sh -t td -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/td.log"
              ;;

              cnv)
                "${DIR_SCRIPT}"/make_panel_cnv.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/cnv.log"
              ;;

              vc)
                "${DIR_SCRIPT}"/make_panel_vc.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/vc.log"
              ;;

              rna_fusions)
                "${DIR_SCRIPT}"/make_rna_fusions.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/fusions.log"
              ;;

              report)
                "${DIR_SCRIPT}"/make_panel_report.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/report.log"
                if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_Panel.pdf" ]]; then
                  touch "${DIR_TARGET}/.processed"
                  echo "${PARAM_DIR_PATIENT} Report finished"
                else
                  echo "${PARAM_DIR_PATIENT} failed"
                fi
              ;;

              td_rna_fusions_parallel)
                ("${DIR_SCRIPT}"/make_panel_alignment.sh -t td -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/td.log") &
                ("${DIR_SCRIPT}"/make_rna_fusions.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/fusions.log") &
                wait
              ;;

              vc_cnv_parallel)
                ("${DIR_SCRIPT}"/make_panel_vc.sh  -p -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/vc.log") &
                ("${DIR_SCRIPT}"/make_panel_cnv.sh -p -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/cnv.log") &
                wait
              ;;
            esac
          ;;

          tumorOnly)
            case "${PARAM_TASK}" in
              td)
                "${DIR_SCRIPT}"/make_tumorOnly_alignment.sh -t td -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/td.log"
              ;;

              cnv)
                "${DIR_SCRIPT}"/make_tumorOnly_cnv.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/cnv.log"
              ;;

              vc)
                "${DIR_SCRIPT}"/make_tumorOnly_vc.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/vc.log"
              ;;

              report)
                "${DIR_SCRIPT}"/make_tumorOnly_report.sh -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/report.log"
                if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.pdf" ]]; then
                  touch "${DIR_TARGET}/.processed"
                  echo "${PARAM_DIR_PATIENT} Report finished"
                else
                  echo "${PARAM_DIR_PATIENT} failed"
                fi
              ;;

              vc_cnv_parallel)
                ("${DIR_SCRIPT}"/make_tumorOnly_vc.sh  -p -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/vc.log") &
                ("${DIR_SCRIPT}"/make_tumorOnly_cnv.sh -p -d "${PARAM_DIR_PATIENT}" &> "${DIR_LOG}/cnv.log") &
                wait
              ;;
            esac
          ;;
        esac

        cleanup "${PARAM_DIR_PATIENT}"
      else
        echo "computing ${PARAM_DIR_PATIENT}"
        
        if [[ ${PARAM_PROTOCOL} == "wes" ]]; then
          echo "WES Protocol"
          if [[ -z "${PARAM_SEQ}" ]]; then
            run_pipe "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"
          else
            run_pipe_seq "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"
          fi

          # check if report was generated successfully
          if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report.pdf" ]]; then
            touch "${DIR_TARGET}/.processed"
            echo "${PARAM_DIR_PATIENT} finished"
          else
            echo "${PARAM_DIR_PATIENT} failed"
          fi
        fi
        if [[ ${PARAM_PROTOCOL} == "panel" ]]; then
          echo "Panel Protocol"
          if [[ -z "${PARAM_SEQ}" ]]; then
            run_panel_pipe "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"
          else
            run_panel_pipe_seq "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"
          fi

          # check if report was generated successfully
          if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_Panel.pdf" ]]; then
            touch "${DIR_TARGET}/.processed"
            echo "${PARAM_DIR_PATIENT} finished"
          else
            echo "${PARAM_DIR_PATIENT} failed"
          fi
        fi
        if [[ ${PARAM_PROTOCOL} == "tumorOnly" ]]; then
          echo "tumorOnly Protocol"
          if [[ -z "${PARAM_SEQ}" ]]; then
            run_tumorOnly_pipe "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"
          else
            run_tumorOnly_pipe_seq "${PARAM_DIR_PATIENT}" "${DIR_TARGET}"
          fi

          # check if report was generated successfully
          if [[ -f "${DIR_ANALYSES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_Report_tumorOnly.pdf" ]]; then
            touch "${DIR_TARGET}/.processed"
            echo "${PARAM_DIR_PATIENT} finished"
          else
            echo "${PARAM_DIR_PATIENT} failed"
          fi
        fi
      fi
    fi
  fi
fi
