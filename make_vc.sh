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

while getopts d:t:ph option; do
  case "${option}" in
  d) readonly PARAM_DIR_PATIENT=$OPTARG ;;
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
if [[ "$(get_config_value common.germline "${PARAM_DIR_PATIENT}")" = "True" ]]; then
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

[[ -d "${DIR_ANALYSES}" ]] || mkdir -p "${DIR_ANALYSES}"

# names
readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_vc
readonly NameGD=${CFG_CASE}_${PARAM_DIR_PATIENT}_gd
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td

# input
readonly recalbamGD=${DIR_WES}/${NameGD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
readonly recalbamTD=${DIR_WES}/${NameTD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam

# keep
readonly snpvcf=${DIR_WES}/${NameD}.output.snp.vcf
readonly indelvcf=${DIR_WES}/${NameD}.output.indel.vcf
readonly MSI_OUTPUT=${DIR_WES}/${NameD}_MSI
readonly BAMMATCHER_OUTPUT=${DIR_WES}/${CFG_CASE}_${PARAM_DIR_PATIENT}_bam-matcher.txt

${BIN_MPILEUP} --adjust-MQ "${CFG_SAMTOOLS_MPILEUP_ADJUSTMQ}" --min-MQ "${CFG_SAMTOOLS_MPILEUP_MINMQ}" --min-BQ "${CFG_GENERAL_MINBASEQUAL}" --max-depth "${CFG_SAMTOOLS_MPILEUP_MAXDEPTH}" -f "${FILE_GENOME}" -l "${CFG_REFERENCE_CAPTUREREGIONS}" "${recalbamGD}" "${recalbamTD}" | ${BIN_SOMATIC} --output-snp "${snpvcf}" --output-indel "${indelvcf}" \
  --min-coverage "${CFG_VARSCAN_SOMATIC_MINCOVERAGE}" --tumor-purity "${CFG_VARSCAN_SOMATIC_TUMORPURITY}" \
  --min-var-freq "${CFG_GENERAL_MINVAF}" --min-freq-for-hom "${CFG_VARSCAN_SOMATIC_MINFREQFORHOM}" \
  --min-avg-qual "${CFG_GENERAL_MINBASEQUAL}" --normal-purity "${CFG_VARSCAN_SOMATIC_NORMALPURITY}" \
  --min-coverage-normal "${CFG_VARSCAN_SOMATIC_MINCOVERAGENORMAL}" --min-coverage-tumor "${CFG_VARSCAN_SOMATIC_MINCOVERAGETUMOR}" \
  --p-value "${CFG_VARSCAN_SOMATIC_PVALUE}" --somatic-p-value "${CFG_VARSCAN_SOMATIC_SOMATICPVALUE}" \
  --strand-filter "${CFG_VARSCAN_SOMATIC_STRANDFILTER}" --validation "${CFG_VARSCAN_SOMATIC_VALIDATION}" --output-vcf 1 --mpileup 1

# Processing of somatic mutations

${BIN_PROCESSSOMATIC} "${snpvcf}" --min-tumor-freq "${CFG_GENERAL_MINVAF}" \
 --max-normal-freq "${CFG_VARSCAN_PROCESSSOMATIC_MAXNORMALFREQ}" --p-value "${CFG_VARSCAN_PROCESSSOMATIC_PVALUE}"
${BIN_PROCESSSOMATIC} "${indelvcf}" --min-tumor-freq "${CFG_GENERAL_MINVAF}" \
 --max-normal-freq "${CFG_VARSCAN_PROCESSSOMATIC_MAXNORMALFREQ}" --p-value "${CFG_VARSCAN_PROCESSSOMATIC_PVALUE}"

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
    # temp
    hc_avi=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.avinput
    hc_rci=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.readcount.input
    hc_rcs=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.readcounts

    # keep
    hc_vcf=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.vcf
    hc_fpf=${DIR_WES}/${NameD}.output.${name1}.${name2}.hc.fpfilter.vcf
    
    if [[ "${name2}" = "Somatic" ]]; then
      recalbam=${recalbamTD}
    else
      recalbam=${recalbamGD}
    fi
    
    ${CONVERT2ANNOVAR2} "${hc_avi}" "${hc_vcf}"
    ${BIN_CUT} "${hc_avi}" > "${hc_rci}"
    ${BIN_BAM_READCOUNT} -l "${hc_rci}" "${recalbam}" > "${hc_rcs}"
    ${BIN_VAR_SCAN} fpfilter "${hc_vcf}" "${hc_rcs}" --output-file "${hc_fpf}" --keep-failures 1 \
      --min-ref-basequal "${CFG_GENERAL_MINBASEQUAL}" --min-var-basequal "${CFG_GENERAL_MINBASEQUAL}" \
      --min-var-count "${CFG_VARSCAN_FPFILTER_MINVARCOUNT}" --min-var-freq "${CFG_GENERAL_MINVAF}" \
      --min-var-count-lc "${CFG_VARSCAN_FPFILTER_MINVARCOUNTLC}" --max-somatic-p "${CFG_VARSCAN_FPFILTER_MAXSOMATICP}" \
      --max-somatic-p-depth "${CFG_VARSCAN_FPFILTER_MAXSOMATICPDEPTH}" --min-ref-readpos "${CFG_VARSCAN_FPFILTER_MINREFREADPOS}" \
      --min-var-readpos "${CFG_VARSCAN_FPFILTER_MINVARREADPOS}" --min-ref-dist3 "${CFG_VARSCAN_FPFILTER_MINREFDIST3}" \
      --min-var-dist3 "${CFG_VARSCAN_FPFILTER_MINVARDIST3}" --min-strandedness "${CFG_VARSCAN_FPFILTER_MINSTRANDEDNESS}" \
      --min-strand-reads "${CFG_VARSCAN_FPFILTER_MINSTRANDREADS}" --max-basequal-diff "${CFG_VARSCAN_FPFILTER_MAXBASEQUALDIFF}" \
      --min-ref-avgrl "${CFG_VARSCAN_FPFILTER_MINREFAVGRL}" --min-var-avgrl "${CFG_VARSCAN_FPFILTER_MINVARAVGRL}" \
      --max-rl-diff "${CFG_VARSCAN_FPFILTER_MAXRLDIFF}" --max-ref-mmqs "${CFG_VARSCAN_FPFILTER_MAXREFMMQS}" \
      --max-var-mmqs "${CFG_VARSCAN_FPFILTER_MAXVARMMQS}" --min-mmqs-diff "${CFG_VARSCAN_FPFILTER_MINMMQSDIFF}" \
      --max-mmqs-diff "${CFG_VARSCAN_FPFILTER_MAXMMQSDIFF}" --min-ref-mapqual "${CFG_VARSCAN_FPFILTER_MINREFMAPQUAL}" \
      --min-var-mapqual "${CFG_VARSCAN_FPFILTER_MINVARMAPQUAL}" --max-mapqual-diff "${CFG_VARSCAN_FPFILTER_MAXMAPQUALDIFF}"
  done
done

readonly data=${DIR_WES}
for name1 in ${names1}; do
  # Annotation snp.Somatic.hc $data/NameD.output.snp.Somatic.hc.fpfilter.vcf
  # Annotation indel.Somatic.hc $data/NameD.output.indel.Somatic.hc.fpfilter.vcf
  # temp
  hc=${DIR_WES}/${NameD}.output.${name1}.Somatic.hc
  hc_T_avi=${DIR_WES}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput

  # keep
  hc_fpf=${DIR_WES}/${NameD}.output.${name1}.Somatic.hc.fpfilter.vcf
  hc_T_avi_multi=${DIR_WES}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput.hg19_multianno.csv
  hc_snpeff="${DIR_WES}/${NameD}.output.${name1}.Somatic.SnpEff.vcf"
  
  # annovar annotation
  ${CONVERT2ANNOVAR} "${hc}" "${hc_fpf}" -allsample
  ${TABLEANNOVAR} "${hc_T_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
      -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

  # snpEff; identify canonical transcript 
  ${BIN_SNPEFF} "${hc_fpf}" > "${hc_snpeff}"

  if [[ "${CFG_CASE}" = somaticGermline ]]; then
    # Annotation snp.Germline.hc $data/NameD.output.snp.Germline.hc.fpfilter.vcf
    # Annotation indel.Germline.hc $data/NameD.output.indel.Germline.hc.fpfilter.vcf
    # temp
    hc=${DIR_WES}/${NameD}.output.${name1}.Germline.hc
    hc_N_avi=${DIR_WES}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput

    # keep
    hc_fpf=${DIR_WES}/${NameD}.output.${name1}.Germline.hc.fpfilter.vcf
    hc_N_avi_multi=${DIR_WES}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput.hg19_multianno.csv
    hc_N_snpeff=${DIR_WES}/${NameD}.output.${name1}.NORMAL.SnpEff.vcf

    # annovar annotation
    ${CONVERT2ANNOVAR} "${hc}" "${hc_fpf}" -allsample
    ${TABLEANNOVAR} "${hc_N_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
        -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

    # snpEff; identify canonical transcripts
    ${BIN_SNPEFF} "${hc_fpf}" > "${hc_N_snpeff}"
  fi

  # Annotation snp.LOH.hc
  # Annotation indel.LOH.hc
  # temp
  hc_avi=${DIR_WES}/${NameD}.output.${name1}.LOH.hc.avinput

  # keep
  hc_fpf=${DIR_WES}/${NameD}.output.${name1}.LOH.hc.fpfilter.vcf
  hc_avi_multi=${DIR_WES}/${NameD}.output.${name1}.LOH.hc.avinput.hg19_multianno.csv
  hc_L_snpeff=${DIR_WES}/${NameD}.output.${name1}.LOH.SnpEff.vcf

  # annovar annotation
  ${CONVERT2ANNOVAR3} "${hc_avi}" "${hc_fpf}"
  ${TABLEANNOVAR} "${hc_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
      -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

  # snpEff; identify canonical transcript
  ${BIN_SNPEFF} "${hc_fpf}" > "${hc_L_snpeff}"
done

# MSI
if [ ! -f "${MICROSATELLITE_SITES}" ]; then
    echo "${MICROSATELLITE_SITES} does not exist. Generating ..."
    ${MSISENSOR_PRO_SCAN} -d "${FILE_GENOME}" -o "${MICROSATELLITE_SITES}"
fi

${MSISENSOR_PRO} -d ${MICROSATELLITE_SITES} -n "${recalbamGD}" -t "${recalbamTD}" -o "${MSI_OUTPUT}"

${BAM_MATCHER} --bam1 ${recalbamGD} --bam2 ${recalbamTD} --output "${BAMMATCHER_OUTPUT}"

#eo VC