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
if [[ "$(get_config_value common.protocol "${PARAM_DIR_PATIENT}")" = "tumorOnly" ]]; then
  readonly CFG_CASE=tumorOnly
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
readonly NameTD=${CFG_CASE}_${PARAM_DIR_PATIENT}_td

# temp
readonly mpileup=${DIR_TMP}/${NameD}_mpileup # DIR_WES!

# keep
readonly recalbam=${DIR_WES}/${NameTD}_output.sort.rmdup.realigned.fixed.recal.bam
readonly snpvcf=${DIR_WES}/${NameD}.output.snp.vcf
readonly indelvcf=${DIR_WES}/${NameD}.output.indel.vcf
readonly OUTPUT_GZ=${DIR_WES}/${NameTD}_gatk4_mutect2.vcf.gz
readonly OUTPUT_FILTERED_GZ=${DIR_WES}/${NameTD}_gatk4_mutect2_filtered.vcf.gz
readonly OUTPUT=${DIR_WES}/${NameTD}_gatk4_mutect2_filtered
readonly ANNOVAR_OUTPUT=${DIR_WES}/${NameTD}.hg19_multianno.vcf
readonly MSI_OUTPUT=${DIR_WES}/${NameTD}_MSI

${BIN_MPILEUP} --adjust-MQ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_ADJUSTMQ}" --min-MQ "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_MINMQ}" --min-BQ "${CFG_GENERAL_MINBASEQUAL}" --max-depth "${CFG_TUMORONLY_SAMTOOLS_MPILEUP_MAXDEPTH}" -f "${FILE_GENOME}" -l "${CFG_REFERENCE_CAPTUREREGIONS}" "${recalbam}" > "${mpileup}"
${BIN_VAR_SCAN} mpileup2snp "${mpileup}" --min-coverage "${CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINCOVERAGE}" --min-reads2 "${CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINREADS2}" \
    --min-freq-for-hom "${CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_MINFREQFORHOM}" --p-value "${CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_PVALUE}" --min-avg-qual "${CFG_GENERAL_MINBASEQUAL}" \
    --strand-filter "${CFG_TUMORONLY_VARSCAN_MPILEUP2SNP_STRANDFILTER}" --min-var-freq "${CFG_GENERAL_MINVAF}" --output-vcf 1 > "${snpvcf}"
${BIN_VAR_SCAN} mpileup2indel "${mpileup}" --min-coverage "${CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINCOVERAGE}" --min-reads2 "${CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINREADS2}" \
    --min-freq-for-hom "${CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_MINFREQFORHOM}" --p-value "${CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_PVALUE}" --min-avg-qual "${CFG_GENERAL_MINBASEQUAL}" \
    --strand-filter "${CFG_TUMORONLY_VARSCAN_MPILEUP2INDEL_STRANDFILTER}" --min-var-freq "${CFG_GENERAL_MINVAF}" --output-vcf 1 > "${indelvcf}"

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
      --min-ref-basequal "${CFG_GENERAL_MINBASEQUAL}" --min-var-basequal "${CFG_GENERAL_MINBASEQUAL}" \
      --min-var-count "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARCOUNT}" --min-var-freq "${CFG_GENERAL_MINVAF}" \
      --min-var-count-lc "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARCOUNTLC}" --max-somatic-p "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXSOMATICP}" \
      --max-somatic-p-depth "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXSOMATICPDEPTH}" --min-ref-readpos "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFREADPOS}" \
      --min-var-readpos "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARREADPOS}" --min-ref-dist3 "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFDIST3}" \
      --min-var-dist3 "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARDIST3}" --min-strandedness "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINSTRANDEDNESS}" \
      --min-strand-reads "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINSTRANDREADS}" --max-basequal-diff "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXBASEQUALDIFF}" \
      --min-ref-avgrl "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFAVGRL}" --min-var-avgrl "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARAVGRL}" \
      --max-rl-diff "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXRLDIFF}" --max-ref-mmqs "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXREFMMQS}" \
      --max-var-mmqs "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXVARMMQS}" --min-mmqs-diff "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINMMQSDIFF}" \
      --max-mmqs-diff "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXMMQSDIFF}" --min-ref-mapqual "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINREFMAPQUAL}" \
      --min-var-mapqual "${CFG_TUMORONLY_VARSCAN_FPFILTER_MINVARMAPQUAL}" --max-mapqual-diff "${CFG_TUMORONLY_VARSCAN_FPFILTER_MAXMAPQUALDIFF}"
done

for name1 in ${names1}; do
  # Annotation
  # temp
  hc=${DIR_WES}/${NameD}.output.${name1}
  hc_T_avi=${DIR_WES}/${NameD}.output.${name1}.Sample1.avinput

  # keep
  hc_fpf=${DIR_WES}/${NameD}.output.${name1}.fpfilter.vcf
  hc_snpeff="${DIR_WES}/${NameD}.output.${name1}.SnpEff.vcf"

  # annovar annotation
  ${CONVERT2ANNOVAR} "${hc}" "${hc_fpf}" -allsample
  ${TABLEANNOVAR} "${hc_T_avi}" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" -buildver hg19 \
      -operation "${CFG_ANNOVAR_ARGOP}" -csvout -otherinfo -remove -nastring NA

  # snpEff; identify canonical transcript
  ${BIN_SNPEFF} "${hc_fpf}" > "${hc_snpeff}"
done

# GATK4 Mutect2
# ANNOVAR settings
CODINGARG="--includesnp --onlyAltering --mrnaseq --tolerate"
CONVERTARG="--includeinfo"

# Mutect2
${BIN_GATK4} Mutect2 -R ${FILE_GENOME} -I ${recalbam} -O ${OUTPUT_GZ} \
 --callable-depth "${CFG_TUMORONLY_MUTECT_CALLABLEDEPTH}" --intervals "${CFG_REFERENCE_CAPTUREREGIONS}" --min-base-quality-score "${CFG_GENERAL_MINBASEQUAL}" --base-quality-score-threshold "${CFG_GENERAL_MINBASEQUAL}"

# Filter
${BIN_GATK4} FilterMutectCalls -V ${OUTPUT_GZ} -R ${FILE_GENOME} -O ${OUTPUT_FILTERED_GZ} --intervals "${CFG_REFERENCE_CAPTUREREGIONS}" --min-median-base-quality "${CFG_GENERAL_MINBASEQUAL}" --min-allele-fraction "${CFG_GENERAL_MINVAF}"
gunzip "${OUTPUT_FILTERED_GZ}"

# Annovar
${TABLEANNOVAR} "${OUTPUT}.vcf" "${DIR_ANNOVAR_DATA}" -protocol "${CFG_ANNOVAR_PROTOCOL}" \
 --buildver hg19 --outfile "${OUTPUT}" --operation "${CFG_ANNOVAR_ARGOP}" \
 --nastring . --vcfinput --thread "${CFG_COMMON_CPUCORES}" --maxgenethread "${CFG_COMMON_CPUCORES}" \
 --otherinfo --remove --verbose --polish \
 --convertarg "${CONVERTARG}"
rm "${OUTPUT}.avinput"

# snpEff
${BIN_SNPEFF} "${OUTPUT}.vcf" > "${OUTPUT}_SnpEff.vcf"

# MSI
#${MSISENSOR2} -t "${recalbam}" -o "${MSI_OUTPUT}"

#eo VC
