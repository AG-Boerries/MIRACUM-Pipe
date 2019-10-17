#!/usr/bin/env bash

###########################################
## WES Pipeline for somatic and germline ##
###########################################
# script to run the actual analysis
# Version 31.07.2019

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
  echo "  -p                  computing as parallel process"
  echo "  -h                  show this help screen"
  exit 1
}

while getopts d:t:ph option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
  t) readonly PARAM_TASK=$OPTARG ;;
  h) usage ;;
  p) readonly PARALLEL_PROCESSES=2 ;;
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

[[ -d "${DIR_ANALYSIS}" ]] || mkdir -p "${DIR_ANALYSIS}"

readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_${PARAM_TASK}
readonly NameGD=${CFG_CASE}_${PARAM_DIR_PATIENT}_gd
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td
readonly recalbamGD=${DIR_WES}/${NameGD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
readonly recalbamTD=${DIR_WES}/${NameTD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
readonly snpvcf=${DIR_WES}/${NameD}.output.snp.vcf
readonly indelvcf=${DIR_WES}/${NameD}.output.indel.vcf

${BIN_MPILEUP} "${recalbamGD}" "${recalbamTD}" | ${BIN_SOMATIC} --output-snp "${snpvcf}" --output-indel "${indelvcf}" \
  --min-coverage "${CFG_VARSCAN_SOMATIC_MINCOVERAGE}" --tumor-purity "${CFG_VARSCAN_SOMATIC_TUMORPURITY}" \
  --min-var-freq "${CFG_VARSCAN_MINVAF}" --min-freq-for-hom "${CFG_VARSCAN_SOMATIC_MINFREQFORHOM}" \
  --min-avg-qual "${CFG_VARSCAN_MINBASEQUAL}" --output-vcf 1 --mpileup 1

# Processing of somatic mutations

${BIN_PROCESSSOMATIC} "${snpvcf}" --min-tumor-freq "${CFG_VARSCAN_MINVAF}"
${BIN_PROCESSSOMATIC} "${indelvcf}" --min-tumor-freq "${CFG_VARSCAN_MINVAF}"

# FP Filter:  snp.Somatic.hc snp.LOH.hc snp.Germline.hc
# FP Filter:  indel.Somatic.hc indel.LOH.hc indel.Germline.hc

readonly names1="snp indel"
for name1 in ${names1}; do

  if [[ "${CFG_CASE}" = somatic ]]; then
    names2="Somatic LOH"
  else
    names2="Somatic LOH Germline"
  fi

  for name2 in ${names2}; do
    hc_vcf=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.vcf
    hc_avi=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.avinput
    hc_rci=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.readcount.input
    hc_rcs=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.readcounts
    hc_fpf=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.fpfilter.vcf
    if [[ "${name2}" = Somatic ]]; then
      recalbam=${recalbamTD}
    else
      recalbam=${recalbamGD}
    fi
    ${CONVERT2ANNOVAR2} "${hc_avi}" "${hc_vcf}"
    ${BIN_CUT} "${hc_avi}" > "${hc_rci}"
    ${BIN_BAM_READCOUNT} -l "${hc_rci}" "${recalbam}" > "${hc_rcs}"
    ${BIN_VAR_SCAN} fpfilter "${hc_vcf}" "${hc_rcs}" --output-file "${hc_fpf}" --keep-failures 1 \
      --min-ref-basequal "${CFG_VARSCAN_MINBASEQUAL}" --min-var-basequal "${CFG_VARSCAN_MINBASEQUAL}" \
      --min-var-count "${CFG_VARSCAN_FPLFILTER_MINVARCOUNT}" --min-var-freq "${CFG_VARSCAN_MINVAF}"
  done
done

readonly data=${DIR_WES}
for name1 in ${names1}; do
  # Annotation snp.Somatic.hc $data/NameD.output.snp.Somatic.hc.fpfilter.vcf
  # Annotation indel.Somatic.hc $data/NameD.output.indel.Somatic.hc.fpfilter.vcf
  hc_=${data}/${NameD}.output.${name1}.Somatic.hc
  hc_fpf=${data}/${NameD}.output.${name1}.Somatic.hc.fpfilter.vcf
  hc_T_avi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput
  hc_T_avi_multi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput.hg19_multianno.csv
  ${CONVERT2ANNOVAR} "${hc_}" "${hc_fpf}" -allsample
  ${TABLEANNOVAR} "${hc_T_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
      -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

  hc_snpeff="${data}/${NameD}.output.${name1}.Somatic.SnpEff.vcf"
  ${BIN_SNPEFF} "${hc_fpf}" > "${hc_snpeff}"

  if [[ "${CFG_CASE}" = somaticGermline ]]; then
    # Annotation snp.Germline.hc $data/NameD.output.snp.Germline.hc.fpfilter.vcf
    # Annotation indel.Germline.hc $data/NameD.output.indel.Germline.hc.fpfilter.vcf
    hc_=${data}/${NameD}.output.${name1}.Germline.hc
    hc_fpf=${data}/${NameD}.output.${name1}.Germline.hc.fpfilter.vcf
    hc_N_avi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput
    hc_N_avi_multi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput.hg19_multianno.csv
    ${CONVERT2ANNOVAR} "${hc_}" "${hc_fpf}" -allsample
    ${TABLEANNOVAR} "${hc_N_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
        -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

    hc_N_snpeff=${data}/${NameD}.output.${name1}.NORMAL.SnpEff.vcf
    ${BIN_SNPEFF} "${hc_fpf}" > "${hc_N_snpeff}"
  fi

  # Annotation snp.LOH.hc
  # Annotation indel.LOH.hc
  hc_vcf=${data}/${NameD}.output.${name1}.LOH.hc.vcf
  hc_fpf=${data}/${NameD}.output.${name1}.LOH.hc.fpfilter.vcf
  hc_avi=${data}/${NameD}.output.${name1}.LOH.hc.avinput
  hc_avi_multi=${data}/${NameD}.output.${name1}.LOH.hc.avinput.hg19_multianno.csv
  ${CONVERT2ANNOVAR3} "${hc_avi}" "${hc_fpf}"
  ${TABLEANNOVAR} "${hc_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
      -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

  hc_L_snpeff=${data}/${NameD}.output.${name1}.LOH.SnpEff.vcf
  ${BIN_SNPEFF} "${hc_fpf}" > "${hc_L_snpeff}"
done