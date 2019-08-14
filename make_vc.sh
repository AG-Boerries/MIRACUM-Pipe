#!/usr/bin/env bash

if [ ! -d ${DIR_TMP} ]; then
  mkdir ${DIR_TMP}
fi

NameD=${case}_${num}_${task}
NameGD=${case}_${num}_GD
NameTD=${case}_${num}_TD
recalbamGD=${wes}/${NameGD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
recalbamTD=${wes}/${NameTD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
snpvcf=${wes}/${NameD}.output.snp.vcf
indelvcf=${wes}/${NameD}.output.indel.vcf

${MPILEUP} ${recalbamGD} ${recalbamTD} | ${SOMATIC} --output-snp ${snpvcf} --output-indel ${indelvcf} --min-coverage ${minCoverage} --tumor-purity ${TumorPurity} --min-var-freq ${minVAF} --min-freq-for-hom ${minFreqForHom} --min-avg-qual ${minBaseQual} --output-vcf 1 --mpileup 1

# Processing of somatic mutations

${PROCESSSOMATIC} ${snpvcf} --min-tumor-freq ${minVAF}
${PROCESSSOMATIC} ${indelvcf} --min-tumor-freq ${minVAF}

# FP Filter:  snp.Somatic.hc snp.LOH.hc snp.Germline.hc
# FP Filter:  indel.Somatic.hc indel.LOH.hc indel.Germline.hc

names1="snp indel"
for name1 in ${names1}; do

  if [ ${case} = somatic ]; then
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
    if [ ${name2} = Somatic ]; then
      recalbam=${recalbamTD}
    else
      recalbam=${recalbamGD}
    fi
    ${CONVERT2ANNOVAR2} ${hc_avi} ${hc_vcf}
    ${CUT} ${hc_avi} >${hc_rci}
    ${BamReadcount} -l ${hc_rci} ${recalbam} >${hc_rcs}
    ${VarScan} fpfilter ${hc_vcf} ${hc_rcs} --output-file ${hc_fpf} --keep-failures 1 --min-ref-basequal ${minBaseQual} --min-var-basequal ${minBaseQual} --min-var-count ${minVarCount} --min-var-freq ${minVAF}
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
  ${CONVERT2ANNOVAR} ${hc_} ${hc_fpf} -allsample
  ${TABLEANNOVAR} ${hc_T_avi} ${ANNOVARData} -protocol ${protocol} -buildver hg19 -operation ${argop} -csvout -otherinfo -remove -nastring NA
  hc_snpeff=$data/${NameD}.output.$name1.Somatic.SnpEff.vcf
  ${SNPEFF} ${hc_fpf} >${hc_snpeff}

  if [ $case = somaticGermline ]; then
    # Annotation snp.Germline.hc $data/NameD.output.snp.Germline.hc.fpfilter.vcf
    # Annotation indel.Germline.hc $data/NameD.output.indel.Germline.hc.fpfilter.vcf
    hc_=${data}/${NameD}.output.${name1}.Germline.hc
    hc_fpf=${data}/${NameD}.output.${name1}.Germline.hc.fpfilter.vcf
    hc_N_avi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput
    hc_N_avi_multi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput.hg19_multianno.csv
    ${CONVERT2ANNOVAR} ${hc_} ${hc_fpf} -allsample
    ${TABLEANNOVAR} ${hc_N_avi} ${ANNOVARData} -protocol ${protocol} -buildver hg19 -operation ${argop} -csvout -otherinfo -remove -nastring NA
    hc_N_snpeff=${data}/${NameD}.output.${name1}.NORMAL.SnpEff.vcf
    ${SNPEFF} ${hc_fpf} >${hc_N_snpeff}
  fi

  # Annotation snp.LOH.hc
  # Annotation indel.LOH.hc
  hc_vcf=${data}/${NameD}.output.${name1}.LOH.hc.vcf
  hc_fpf=${data}/${NameD}.output.${name1}.LOH.hc.fpfilter.vcf
  hc_avi=${data}/${NameD}.output.${name1}.LOH.hc.avinput
  hc_avi_multi=${data}/${NameD}.output.${name1}.LOH.hc.avinput.hg19_multianno.csv
  ${CONVERT2ANNOVAR3} ${hc_avi} ${hc_fpf}
  ${TABLEANNOVAR} ${hc_avi} ${ANNOVARData} -protocol ${protocol} -buildver hg19 -operation ${argop} -csvout -otherinfo -remove -nastring NA
  hc_L_snpeff=${data}/${NameD}.output.${name1}.LOH.SnpEff.vcf
  ${SNPEFF} ${hc_fpf} >${hc_L_snpeff}
done

rm -r ${DIR_TMP}