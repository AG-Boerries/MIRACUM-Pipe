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
# readonly CFG_PROTOCOL=$(get_config_value common.protocol# "${PARAM_DIR_PATIENT}")
if [[ "$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")" = "panel" ]]; then
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

[[ -d "${DIR_ANALYSES}" ]] || mkdir -p "${DIR_ANALYSES}"

# names
readonly NameD=${CFG_CASE}_${PARAM_DIR_PATIENT}_vc
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td

# temp
readonly mpileup=${DIR_WES}/${NameD}_mpileup # TODO change back to DIR_TMP!

# keep
readonly recalbam=${DIR_WES}/${NameTD}_output.sort.realigned.fixed.recal.bam
readonly snpvcf=${DIR_WES}/${NameD}.output.snp.vcf
readonly indelvcf=${DIR_WES}/${NameD}.output.indel.vcf

${BIN_MPILEUP} --adjust-MQ "${CFG_SAMTOOLS_MPILEUP_ADJUSTMQ}" --min-MQ "${CFG_SAMTOOLS_MPILEUP_MINMQ}" --min-BQ "${CFG_PANEL_MINBASEQUAL}" --max-depth "${CFG_SAMTOOLS_MPILEUP_MAXDEPTH}" -f "${FILE_GENOME}" -l "${CFG_REFERENCE_CAPTUREREGIONS}" "${recalbam}" > "${mpileup}"
${BIN_VAR_SCAN} mpileup2snp "${mpileup}" --min-coverage "${CFG_VARSCAN_MPILEUP2SNP_MINCOVERAGE}" --min-reads2 "${CFG_VARSCAN_MPILEUP2SNP_MINREADS2}" \
    --min-freq-for-hom "${CFG_VARSCAN_MPILEUP2SNP_MINFREQFORHOM}" --p-value "${CFG_VARSCAN_MPILEUP2SNP_PVALUE}" --min-avg-qual "${CFG_PANEL_MINBASEQUAL}" \
    --strand-filter "${CFG_VARSCAN_MPILEUP2SNP_STRANDFILTER}" --min-var-freq "${CFG_PANEL_MINVAF}" --output-vcf 1 > "${snpvcf}"
${BIN_VAR_SCAN} mpileup2indel "${mpileup}" --min-coverage "${CFG_VARSCAN_MPILEUP2INDEL_MINCOVERAGE}" --min-reads2 "${CFG_VARSCAN_MPILEUP2INDEL_MINREADS2}" \
    --min-freq-for-hom "${CFG_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM}" --p-value "${CFG_VARSCAN_MPILEUP2INDEL_PVALUE}" --min-avg-qual "${CFG_PANEL_MINBASEQUAL}" \
    --strand-filter "${CFG_VARSCAN_MPILEUP2INDEL_STRANDFILTER}" --min-var-freq "${CFG_PANEL_MINVAF}" --output-vcf 1 > "${indelvcf}"


readonly names1="snp indel"
for name1 in ${names1}; do
  # temp
  hc_avi=${DIR_WES}/${NameD}.output.${name1}.avinput
  hc_rci=${DIR_WES}/${NameD}.output.${name1}.readcount.input
  hc_rcs=${DIR_WES}/${NameD}.output.${name1}.readcounts

  # keep
  hc_vcf=${DIR_WES}/${NameD}.output.${name1}.vcf
  hc_fpf=${DIR_WES}/${NameD}.output.${name1}.fpfilter.vcf
    
    
  ${CONVERT2ANNOVAR3} "${hc_avi}" "${hc_vcf}"
  ${BIN_CUT} "${hc_avi}" > "${hc_rci}"
  ${BIN_BAM_READCOUNT} -l "${hc_rci}" "${recalbam}" > "${hc_rcs}"
  ${BIN_VAR_SCAN} fpfilter "${hc_vcf}" "${hc_rcs}" --output-file "${hc_fpf}" --keep-failures 1 \
      --min-ref-basequal "${CFG_PANEL_MINBASEQUAL}" --min-var-basequal "${CFG_PANEL_MINBASEQUAL}" \
      --min-var-count "${CFG_VARSCAN_PANEL_FPFILTER_MINVARCOUNT}" --min-var-freq "${CFG_PANEL_MINVAF}" \
      --min-var-count-lc "${CFG_VARSCAN_PANEL_FPFILTER_MINVARCOUNTLC}" --max-somatic-p "${CFG_VARSCAN_PANEL_FPFILTER_MAXSOMATICP}" \
      --max-somatic-p-depth "${CFG_VARSCAN_PANEL_FPFILTER_MAXSOMATICPDEPTH}" --min-ref-readpos "${CFG_VARSCAN_PANEL_FPFILTER_MINREFREADPOS}" \
      --min-var-readpos "${CFG_VARSCAN_PANEL_FPFILTER_MINVARREADPOS}" --min-ref-dist3 "${CFG_VARSCAN_PANEL_FPFILTER_MINREFDIST3}" \
      --min-var-dist3 "${CFG_VARSCAN_PANEL_FPFILTER_MINVARDIST3}" --min-strandedness "${CFG_VARSCAN_PANEL_FPFILTER_MINSTRANDEDNESS}" \
      --min-strand-reads "${CFG_VARSCAN_PANEL_FPFILTER_MINSTRANDREADS}" --max-basequal-diff "${CFG_VARSCAN_PANEL_FPFILTER_MAXBASEQUALDIFF}" \
      --min-ref-avgrl "${CFG_VARSCAN_PANEL_FPFILTER_MINREFAVGRL}" --min-var-avgrl "${CFG_VARSCAN_PANEL_FPFILTER_MINVARAVGRL}" \
      --max-rl-diff "${CFG_VARSCAN_PANEL_FPFILTER_MAXRLDIFF}" --max-ref-mmqs "${CFG_VARSCAN_PANEL_FPFILTER_MAXREFMMQS}" \
      --max-var-mmqs "${CFG_VARSCAN_PANEL_FPFILTER_MAXVARMMQS}" --min-mmqs-diff "${CFG_VARSCAN_PANEL_FPFILTER_MINMMQSDIFF}" \
      --max-mmqs-diff "${CFG_VARSCAN_PANEL_FPFILTER_MAXMMQSDIFF}" --min-ref-mapqual "${CFG_VARSCAN_PANEL_FPFILTER_MINREFMAPQUAL}" \
      --min-var-mapqual "${CFG_VARSCAN_PANEL_FPFILTER_MINVARMAPQUAL}" --max-mapqual-diff "${CFG_VARSCAN_PANEL_FPFILTER_MAXMAPQUALDIFF}"
done

for name1 in ${names1}; do
  # Annotation
  # temp
  hc=${DIR_WES}/${NameD}.output.${name1}
  hc_T_avi=${DIR_WES}/${NameD}.output.${name1}.Sample1.avinput

  # keep
  hc_fpf=${DIR_WES}/${NameD}.output.${name1}.fpfilter.vcf
  #hc_T_avi_multi=${DIR_WES}/${NameD}.output.${name1}.Sample1.avinput.hg19_multianno.csv
  hc_snpeff="${DIR_WES}/${NameD}.output.${name1}.SnpEff.vcf"
  
  # annovar annotation
  ${CONVERT2ANNOVAR} "${hc}" "${hc_fpf}" -allsample
  ${TABLEANNOVAR} "${hc_T_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
      -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

  # snpEff; identify canonical transcript
  ${BIN_SNPEFF} "${hc_fpf}" > "${hc_snpeff}"
done