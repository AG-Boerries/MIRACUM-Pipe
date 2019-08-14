#!/usr/bin/env bash

###########################################
## WES Pipeline for somatic and germline ##
###########################################
# script to run the actual analysis
# Version 31.07.2019



DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")" || exit
  pwd -P
)

DIR_TOOLS="${DIR_SCRIPT}/tools"

## load settings
# shellcheck source=global.sh
source "${DIR_SCRIPT}"/global.sh




function usage() {
  echo "usage: make_alignment_VC_CNV.sh -d dir [-h]"
  echo "  -d  dir             specify file a"
  echo "  -h                  show this help screen"
  exit 1
}

while getopts d:h option; do
  case "${option}" in
  d) DIR_PATIENT=$OPTARG ;;
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

#  c) case=$OPTARG ;;
#  n) num=$OPTARG ;;
#  s) sex=$OPTARG ;;
#  fa) file_a=$OPTARG ;;
#  fb) file_b=$OPTARG ;;

# load patient yaml
sex=get_config_value sex "${DIR_PATIENT}"
if [[ "$(get_config_value annotation.germline "${DIR_PATIENT}")" = "True" ]]; then
  case=somaticGermline
else
  case=somatic
fi

# check inputs
possible_tasks=("GD TD VC CNV Report")
possible_sex=("XX XY")

for value in "${possible_tasks[@]}"
do
  [[ "${task}" = "${value}" ]] && \
    echo "unknown task: ${task}" && \
    echo "use one of the following values: $(join_by ' ' ${possible_tasks})" && \
    exit 1
done

for value in "${possible_sex[@]}"
do
  [[ "${sex}" = "${value}" ]] && \
    echo "unknown sex: ${sex}" && \
    echo "use one of the following values: $(join_by ' ' ${possible_sex})" && \
    exit 1
done

##################################################################################################################
#### Parameters which have to be adjusted accoridng the the environment or the users needs

## everything in assets is intended to be linked as a volume
## General
DIR_ASSETS="${DIR_SCRIPT}/assets"

DIR_DATA="${DIR_ASSETS}/data"             # folder contatining the raw data (.fastq files)
DIR_OUTPUT="${DIR_ASSETS}/output"
DIR_REF="${DIR_ASSETS}/references"        # reference genome

mtb="${DIR_OUTPUT}/${DIR_PATIENT}"        # folder containing output
wes="${mtb}/WES"
ana="${mtb}/Analysis"
RscriptPath="${DIR_SCRIPT}/RScripts"
DatabasePath="${DIR_SCRIPT}/Databases"
DIR_TMP="/tmp" # temporary folder
# end paths

## Genome
GENOME="${DIR_REF}/Genome/$(get_config_value reference.genome "${DIR_PATIENT}")"
Chromosomes="${DIR_REF}/chromosomes"
ChromoLength="${DIR_REF}/chromosomes/$(get_config_value reference.length "${DIR_PATIENT}")"

# depending on measurement machine
## SureSelect (Capture Kit)
CaptureRegions="${DIR_REF}/$(get_config_value reference.sequencing.genome "${DIR_PATIENT}")"

# database for known variants
## dbSNP vcf File
dbSNPvcf="${DIR_REF}/dbSNP/$(get_config_value reference.dbSNP "${DIR_PATIENT}")"
# END variables


nCore=$(get_config_value common.cpucores "${DIR_PATIENT}")
minBaseQual=$(get_config_value varscan.minBaseQual "${DIR_PATIENT}")
minVAF=$(get_config_value varscan.minVAF "${DIR_PATIENT}")

# VarScan somatic
minCoverage=$(get_config_value varscan.somatic.minCoverage "${DIR_PATIENT}")
TumorPurity=$(get_config_value varscan.somatic.tumorPurity "${DIR_PATIENT}")
minFreqForHom=$(get_config_value varscan.somatic.minFreqForHom "${DIR_PATIENT}")

# VarScan fpfilter
minVarCount=$(get_config_value vascan.fpfilter.minVarCount "${DIR_PATIENT}")

# ANNOVAR Databases
protocol=$(get_config_value annovar.protocol "${DIR_PATIENT}")
argop=$(get_config_value annovar.argop "${DIR_PATIENT}")

## Tools and paths
# Paths
soft=${DIR_TOOLS}                       # folder containing all used tools
java="java -Djava.io.tmpdir=${DIR_TMP} " # path to java

# Pre-Processing
FASTQC="${soft}/FastQC/fastqc -t ${nCore} --extract "

TRIM="${java} -Xmx150g -jar ${soft}/Trimmomatic/trimmomatic.jar PE -threads ${nCore} -phred33 "
TrimmomaticAdapter="${soft}/Trimmomatic/adapters"
CUT="cut -f1,2,3"

# Alignment
BWAMEM="${soft}/bwa/bwa mem -M "

# BAM-Readcount
BamReadcount="${soft}/bam-readcount/bin/bam-readcount -q 1 -b 20 -w 1 -f ${GENOME} "

# SAMTOOLS
SAMTOOLS="${soft}/samtools/samtools" # path to samtools
SAMVIEW="${SAMTOOLS} view -@ ${nCore} "
SAMSORT="${SAMTOOLS} sort -@ ${nCore} "
SAMRMDUP="${SAMTOOLS} rmdup "
SAMINDEX="${SAMTOOLS} index "
MPILEUP="${SAMTOOLS} mpileup -B -C 50 -f ${GENOME} -q 1 --min-BQ ${minBaseQual}"
STATS="${SAMTOOLS} stats "

# GATK
GATK="${soft}/gatk/gatk"
RealignerTargetCreator="${GATK} -T RealignerTargetCreator -R ${GENOME} -nt ${nCore} "
IndelRealigner="${GATK} -R ${GENOME} -T IndelRealigner "
BaseRecalibrator="${GATK} -T BaseRecalibrator -l INFO -R ${GENOME} -knownSites ${dbSNPvcf} -nct ${nCore} "
PrintReads="${GATK} -T PrintReads -R ${GENOME} -nct ${nCore} "

# PICARD
FixMate="${soft}/bin/picard FixMateInformation "

# VARSCAN
VarScan="${soft}/varscan/varscan.py"
SOMATIC="${VarScan} somatic"
PROCESSSOMATIC="${VarScan} processSomatic"

# ANNOVAR
ANNOVAR="${soft}/annovar/bin/annovar"
ANNOVARData="${ANNOVAR}/humandb"
CONVERT2ANNOVAR2="${ANNOVAR}/convert2annovar.pl --format vcf4old --outfile "
CONVERT2ANNOVAR3="${ANNOVAR}/convert2annovar.pl --format vcf4old --includeinfo --comment --outfile "
CONVERT2ANNOVAR="${ANNOVAR}/convert2annovar.pl --format vcf4 --includeinfo --comment --withzyg --outfile "
TABLEANNOVAR="${ANNOVAR}/table_annovar.pl"

# COVERAGE
COVERAGE="${soft}/bedtools2/bin/bedtools coverage -hist -g ${GENOME}.fai -sorted "

# SNPEFF
SNPEFF="${java} -Xmx150g -jar ${soft}/snpEff/snpEff.jar GRCh37.75 -c ${soft}/snpEff/snpEff.config -canon -v"

# ControlFREEC
freec="${soft}/bin/freec "

gemMappabilityFile="${DIR_REF}/mappability/$(get_config_value reference.mappability "${DIR_PATIENT}")"

# R
Rscript=$(command -v Rscript)

##################################################################################################################

##########
## MAIN ##
##########

### program calls
case ${task} in

## alignment -----------------------------------------------------------------------------------------------------
# TODO: make_alignment.sh
GD | TD)
  if [ ! -d ${DIR_TMP} ]; then
    mkdir ${DIR_TMP}
  fi

  # SAMPLE
  NameD=${case}_${DIR_PATIENT}_${task}
  xx=${sex}
  InputPath=${DIR_DATA}/${xx} ## change later !!!
  Input1File=${file_a}1     # filename without extension
  Input2File=${file_a}2     # filename without extension

  # temp files
  fastq1=${InputPath}/${Input1File}.fastq.gz
  fastq2=${InputPath}/${Input2File}.fastq.gz
  fastq_o1_p_t=${DIR_TMP}/${NameD}_output1_paired_trimmed.fastq.gz
  fastq_o1_u_t=${DIR_TMP}/${NameD}_output1_unpaired_trimmed.fastq.gz
  fastq_o2_p_t=${DIR_TMP}/${NameD}_output2_paired_trimmed.fastq.gz
  fastq_o2_u_t=${DIR_TMP}/${NameD}_output2_unpaired_trimmed.fastq.gz
  bam=${DIR_TMP}/${NameD}_output.bam
  prefixsort=${DIR_TMP}/${NameD}_output.sort
  sortbam=${DIR_TMP}/${NameD}_output.sort.bam
  rmdupbam=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.bam
  bai=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.bai
  bamlist=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.bam.list
  realignedbam=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.bam
  realignedbai=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.bai
  fixedbam=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.bam
  fixedbai=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.bai
  csv=${DIR_TMP}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.recal_data.csv

  recalbam=${wes}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
  statstxt=${wes}/${NameD}_stats.txt
  coveragetxt=${wes}/${NameD}_coverage.all.txt

  # fastqc zip to WES
  ${FASTQC} "${fastq1}" -o "${wes}"
  ${FASTQC} "${fastq2}" -o "${wes}"

  # trim fastq
  ${TRIM} "${fastq1}" "${fastq2}" "${fastq_o1_p_t}" "${fastq_o1_u_t}" "${fastq_o2_p_t}" "${fastq_o2_u_t}" \
    ILLUMINACLIP:"${TrimmomaticAdapter}"/TruSeq3-PE-2.fa:2:30:10 HEADCROP:3 TRAILING:10 MINLEN:25

  ${FASTQC} "${fastq_o1_p_t}" -o "${wes}"
  ${FASTQC} "${fastq_o2_p_t}" -o "${wes}"

  # make bam
  ${BWAMEM} -R "@RG\tID:${NameD}\tSM:${NameD}\tPL:illumina\tLB:lib1\tPU:unit1" -t 12 "${GENOME}" \
    "${fastq_o1_p_t}" "${fastq_o2_p_t}" | ${SAMVIEW} -bS - >"${bam}"

  # stats
  ${STATS} "${bam}" >"${statstxt}"

  # sort bam
  ${SAMSORT} "${bam}" -T "${prefixsort}" -o "${sortbam}"

  # rmdup bam
  ${SAMVIEW} -b -f 0x2 -q1 "${sortbam}" | ${SAMRMDUP} - "${rmdupbam}"

  # make bai
  ${SAMINDEX} "${rmdupbam}" "${bai}"

  # make bam list
  ${RealignerTargetCreator} -o "${bamlist}" -I "${rmdupbam}"

  # realign bam
  ${IndelRealigner} -I "${rmdupbam}" -targetIntervals "${bamlist}" -o "${realignedbam}"

  # fix bam
  ${FixMate} INPUT="${realignedbam}" OUTPUT="${fixedbam}" SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

  # make csv
  ${BaseRecalibrator} -I "${fixedbam}" \
    -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o "${csv}"

  # recal bam
  ${PrintReads} -I "${fixedbam}" -BQSR "${csv}" -o "${recalbam}"

  # coverage
  ${COVERAGE} -b "${recalbam}" -a "${CaptureRegions}" | grep '^all' >"${coveragetxt}"

  # zip
  ${FASTQC} "${recalbam}" -o "${wes}"
  ;;
# eo alignment

## variantCalling ------------------------------------------------------------------------------------------------
# TODO: make_vc.sh
VC)
  if [ ! -d "${DIR_TMP}" ]; then
    mkdir "${DIR_TMP}"
  fi

  NameD=${case}_${DIR_PATIENT}_${task}
  NameGD=${case}_${DIR_PATIENT}_GD
  NameTD=${case}_${DIR_PATIENT}_TD
  recalbamGD=${wes}/${NameGD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
  recalbamTD=${wes}/${NameTD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
  snpvcf=${wes}/${NameD}.output.snp.vcf
  indelvcf=${wes}/${NameD}.output.indel.vcf

  ${MPILEUP} "${recalbamGD}" "${recalbamTD}" | ${SOMATIC} --output-snp "${snpvcf}" --output-indel "${indelvcf}" \
    --min-coverage "${minCoverage}" --tumor-purity "${TumorPurity}" --min-var-freq "${minVAF}" \
    --min-freq-for-hom "${minFreqForHom}" --min-avg-qual "${minBaseQual}" --output-vcf 1 --mpileup 1

  # Processing of somatic mutations

  ${PROCESSSOMATIC} "${snpvcf}" --min-tumor-freq "${minVAF}"
  ${PROCESSSOMATIC} "${indelvcf}" --min-tumor-freq "${minVAF}"

  # FP Filter:  snp.Somatic.hc snp.LOH.hc snp.Germline.hc
  # FP Filter:  indel.Somatic.hc indel.LOH.hc indel.Germline.hc

  names1="snp indel"
  for name1 in ${names1}; do

    if [ "${case}" = somatic ]; then
      names2="Somatic LOH"
    else
      names2="Somatic LOH Germline"
    fi

    for name2 in ${names2}; do
      hc_vcf=${wes}/${NameD}.output.${name1}.${name2}.hc.vcf
      hc_avi=${wes}/${NameD}.output.${name1}.${name2}.hc.avinput
      hc_rci=${wes}/${NameD}.output.${name1}.${name2}.hc.readcount.input
      hc_rcs=${wes}/${NameD}.output.${name1}.${name2}.hc.readcounts
      hc_fpf=${wes}/${NameD}.output.${name1}.${name2}.hc.fpfilter.vcf
      if [ "${name2}" = Somatic ]; then
        recalbam=${recalbamTD}
      else
        recalbam=${recalbamGD}
      fi
      ${CONVERT2ANNOVAR2} "${hc_avi}" "${hc_vcf}"
      ${CUT} "${hc_avi}" > "${hc_rci}"
      ${BamReadcount} -l "${hc_rci}" "${recalbam}" > "${hc_rcs}"
      ${VarScan} fpfilter "${hc_vcf}" "${hc_rcs}" --output-file "${hc_fpf}" --keep-failures 1 \
        --min-ref-basequal "${minBaseQual}" --min-var-basequal "${minBaseQual}" --min-var-count "${minVarCount}" \
        --min-var-freq "${minVAF}"
    done
  done

  data=${wes}
  for name1 in ${names1}; do
    # Annotation snp.Somatic.hc $data/NameD.output.snp.Somatic.hc.fpfilter.vcf
    # Annotation indel.Somatic.hc $data/NameD.output.indel.Somatic.hc.fpfilter.vcf
    hc_=${data}/${NameD}.output.${name1}.Somatic.hc
    hc_fpf=${data}/${NameD}.output.${name1}.Somatic.hc.fpfilter.vcf
    hc_T_avi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput
    hc_T_avi_multi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput.hg19_multianno.csv
    ${CONVERT2ANNOVAR} "${hc_}" "${hc_fpf}" -allsample
    ${TABLEANNOVAR} "${hc_T_avi}" "${ANNOVARData}" -protocol "${protocol}" -buildver hg19 -operation "${argop}" -csvout \
      -otherinfo -remove -nastring NA

    hc_snpeff=$data/${NameD}.output.$name1.Somatic.SnpEff.vcf
    ${SNPEFF} "${hc_fpf}" > "${hc_snpeff}"

    if [ "${case}" = somaticGermline ]; then
      # Annotation snp.Germline.hc $data/NameD.output.snp.Germline.hc.fpfilter.vcf
      # Annotation indel.Germline.hc $data/NameD.output.indel.Germline.hc.fpfilter.vcf
      hc_=${data}/${NameD}.output.${name1}.Germline.hc
      hc_fpf=${data}/${NameD}.output.${name1}.Germline.hc.fpfilter.vcf
      hc_N_avi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput
      hc_N_avi_multi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput.hg19_multianno.csv
      ${CONVERT2ANNOVAR} "${hc_}" "${hc_fpf}" -allsample
      ${TABLEANNOVAR} "${hc_N_avi}" "${ANNOVARData}" -protocol "${protocol}" -buildver hg19 -operation "${argop}" -csvout \
        -otherinfo -remove -nastring NA

      hc_N_snpeff=${data}/${NameD}.output.${name1}.NORMAL.SnpEff.vcf
      ${SNPEFF} "${hc_fpf}" >"${hc_N_snpeff}"
    fi

    # Annotation snp.LOH.hc
    # Annotation indel.LOH.hc
    hc_vcf=${data}/${NameD}.output.${name1}.LOH.hc.vcf
    hc_fpf=${data}/${NameD}.output.${name1}.LOH.hc.fpfilter.vcf
    hc_avi=${data}/${NameD}.output.${name1}.LOH.hc.avinput
    hc_avi_multi=${data}/${NameD}.output.${name1}.LOH.hc.avinput.hg19_multianno.csv
    ${CONVERT2ANNOVAR3} "${hc_avi}" "${hc_fpf}"
    ${TABLEANNOVAR} "${hc_avi}" "${ANNOVARData}" -protocol "${protocol}" -buildver hg19 -operation "${argop}" -csvout \
      -otherinfo -remove -nastring NA

    hc_L_snpeff=${data}/${NameD}.output.${name1}.LOH.SnpEff.vcf
    ${SNPEFF} "${hc_fpf}" >"${hc_L_snpeff}"
  done

  rm -r ${DIR_TMP}
  ;;
# eo VC

## CNV  ----------------------------------------------------------------------------------------------------------
# TODO: make_cnv.sh
CNV)
  output="${wes}/CNV"

  if [ ! -d "${output}" ]; then
    mkdir "${output}"
  fi

  cat >"${wes}"/CNV_config.txt <<EOI
[general]

chrFiles = ${Chromosomes}
chrLenFile = ${ChromoLength}
breakPointType = 4
breakPointThreshold = 1.2
forceGCcontentNormalization = 1
gemMappabilityFile = ${gemMappabilityFile}
intercept = 0
minCNAlength = 3
maxThreads = 12
noisyData = TRUE
outputDir = ${output}
ploidy = 2
printNA = FALSE
readCountThreshold = 50
samtools = ${SAMTOOLS}
sex = ${sex}
step = 0
window = 0
uniqueMatch = TRUE
contaminationAdjustment = TRUE

[sample]

mateFile = ${wes}/${case}_${DIR_PATIENT}_TD_output.sort.filtered.rmdup.realigned.fixed.recal.bam
inputFormat = BAM
mateOrientation = FR

[control]

mateFile = ${wes}/${case}_${DIR_PATIENT}_GD_output.sort.filtered.rmdup.realigned.fixed.recal.bam
inputFormat = BAM
mateOrientation = FR

[target]

captureRegions = ${CaptureRegions}
EOI

  export PATH=${PATH}:${SAMTOOLS}
  "${freec}"-conf "${wes}"/CNV_config.txt
  ;;
  # eo CNV

## Report  -------------------------------------------------------------------------------------------------------
# TODO: make_report.sh
Report)
  cd "${ana}" || exit

  ${Rscript} "${ana}"/Main.R "${case}" "${DIR_PATIENT}" "${file_a}" "${file_b}" "${mtb}" "${RscriptPath}" "${DatabasePath}"

  ${Rscript} -e "library(knitr); knit('Report.Rnw')"
  pdflatex -interaction=nonstopmode Report.tex
  pdflatex -interaction=nonstopmode Report.tex
  ;;
  # eo Report
esac
# eo program calls

# TODO: required in all 4 files
echo "task ${task} for ${DIR_PATIENT} finished"
rm .STARTING_MARKER_"${task}"
exit
