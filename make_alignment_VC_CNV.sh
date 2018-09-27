#!/bin/bash
###########################################
## WES Pipeline for somatic and germline ##
###########################################
# script to run the actual analysis
# Version 05.09.2018


case=$1;  # somatic or somaticGermline
task=$2;  # GD or TD (alignment)
num=$3;   # 73
sex=$4;

##################################################################################################################
#### Parameters which have to be adjusted accoridng the the environment or the users needs

## General
homedata="/path/to/data" # folder contatining the raw data (.fastq files)
annot="/path/to/annotation" # fodler containing annotation files like captureRegions.bed
motherpath="/path/to/output"
mtb="${motherpath}/${case}_${num}" # folder containing output
wes="${mtb}/WES"
ana="${mtb}/Analysis"
RscriptPath="${motherpath}/RScripts"
DatabasePath="${motherpath}/Databases"
tempdir="/scratch" # temporary folder


## Genome
GENOME="/path/to/ref/genome/including/index/hg19.fa"
Chromosomes="/path/to/ref/genome/chromosomes"
ChromoLength="/path/to/ref/chromosomes/length/hg19_chr.len"

## SureSelect (Capture Kit)
CaptureRegions="${annot}/Agilent/SureSelectV5UTR/V5UTR.bed"

## dbSNP vcf File
dbSNPvcf="/path/to/dbSNP/dbSNP/snp150hg19.vcf.gz "

### Software
## Cores to use
nCore=12

soft="/path/to/tools/software" # folder containing all used tools
java="${soft}/bin/java" # path to java

FASTQC="${soft}/FastQC/fastqc -t ${nCore} --extract "
TRIM="${java} -Xmx150g -jar ${soft}/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads ${nCore} -phred33 "
CUT="cut -f1,2,3"

BWAMEM="${soft}/bin/bwa mem -M "
BamReadcount="${soft}/bin/bam-readcount -q 1 -b 20 -w 1 -f ${GENOME} "

# SAMTOOLS
SAMTOOLS="${soft}/bin/samtools" # path to samtools
SAMVIEW="${SAMTOOLS} view -@ ${nCore} "
SAMSORT="${SAMTOOLS} sort -@ ${nCore} "
SAMRMDUP="${SAMTOOLS} rmdup "
SAMINDEX="${SAMTOOLS} index "
MPILEUP="${SAMTOOLS} mpileup -B -C50 -f ${GENOME} -q 1 "
STATS="${SAMTOOLS} stats "

# GATK
GATK="${soft}/bin/gatk"
RealignerTargetCreator="${GATK} -T RealignerTargetCreator -R ${GENOME} -nt ${nCore} "
IndelRealigner="${GATK} -R ${GENOME} -T IndelRealigner "
BaseRecalibrator="${GATK} -T BaseRecalibrator -l INFO -R ${GENOME} -knownSites ${dbSNPvcf} -nct ${nCore} "
PrintReads="${GATK} -T PrintReads -R ${GENOME} -nct ${nCore} "

# PICARD
FixMate="${soft}/bin/picard FixMateInformation "

# VARSCAN
VarScan="${soft}/anaconda2/bin/varscan"
SOMATIC="${VarScan} somatic"
PROCESSSOMATIC="${VarScan} processSomatic"


# ANNOVAR
CONVERT2ANNOVAR2="${soft}/annovar/convert2annovar.pl --format vcf4old --outfile "
CONVERT2ANNOVAR3="${soft}/annovar/convert2annovar.pl --format vcf4old --includeinfo --comment --outfile "
CONVERT2ANNOVAR="${soft}/annovar/convert2annovar.pl --format vcf4 --includeinfo --comment --withzyg --outfile "
TABLEANNOVAR="${soft}/annovar/table_annovar.pl"

# COVERAGE
COVERAGE="${soft}/bin/bedtools coverage -hist "

# SNPEFF
SNPEFF="${java} -Xmx150g -jar ${soft}/snpEff/snpEff.jar GRCh37.75 -c ${soft}/snpEff/snpEff.config -canon -v"

# ControlFREEC
freec="${soft}/bin/freec "
gemMappabilityFile="${soft}/FREEC-11.0/mappability/out100m2_hg19.gem"

# R
Rscript="${soft}/bin/Rscript"

### Software Parameters
# VarScan somatic
minCoverageNormal="8" #default
minCoverageTumor="8" #default
TumorPurity="0.5" # -> I would rather adjust the min-var-freq; tumnor-purity effectively changes (reduces) min-var-freq, i.e. VAF <- tumor-purity * VAF (0.5*0.1=0.05)
minVarFreq="0.10" # VAF, default
minFreqForHom="0.75" # VAF to call homozygote

# VarScan processSomatic
minTumorFreq="0.10" # default

# VarScan fpfilter
minRefBasequal="28"
minVarBasequal="28"
minVarCount="4" #default

# ANNOVAR Databases
protocol='refGene,gnomad_exome,exac03,esp6500siv2_ea,1000g2015aug_eur,avsnp150,clinvar_20170905,cadd13,intervar_20170202,dbnsfp33a,cosmic84_coding,cosmic84_noncoding'
argop='g,f,f,f,f,f,f,f,f,f,f,f' 

 
##################################################################################################################

##########
## MAIN ##
##########

# SAMPLE
NameD=${case}_${num}_${task}
xx=${4};
InputPath=${homedata}/ngs/$xx/fastq ## change later !!!
Input1File=${5}1; # filename without extension
Input2File=${5}2; # filename without extension

# temp files
fastq1=${InputPath}/${Input1File}.fastq.gz
fastq2=${InputPath}/${Input2File}.fastq.gz
fastq_o1_p_t=${tempdir}/${NameD}_output1_paired_trimmed.fastq.gz 
fastq_o1_u_t=${tempdir}/${NameD}_output1_unpaired_trimmed.fastq.gz 
fastq_o2_p_t=${tempdir}/${NameD}_output2_paired_trimmed.fastq.gz 
fastq_o2_u_t=${tempdir}/${NameD}_output2_unpaired_trimmed.fastq.gz 
bam=${tempdir}/${NameD}_output.bam
prefixsort=${tempdir}/${NameD}_output.sort
sortbam=${tempdir}/${NameD}_output.sort.bam
rmdupbam=${tempdir}/${NameD}_output.sort.filtered.rmdup.bam
bai=${tempdir}/${NameD}_output.sort.filtered.rmdup.bai
bamlist=${tempdir}/${NameD}_output.sort.filtered.rmdup.bam.list
realignedbam=${tempdir}/${NameD}_output.sort.filtered.rmdup.realigned.bam
realignedbai=${tempdir}/${NameD}_output.sort.filtered.rmdup.realigned.bai
fixedbam=${tempdir}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.bam
fixedbai=${tempdir}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.bai
csv=${tempdir}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.recal_data.csv

recalbam=${wes}/${NameD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
statstxt=${wes}/${NameD}_stats.txt
coveragetxt=${wes}/${NameD}_coverage.all.txt


### program calls 


## alignment -----------------------------------------------------------------------------------------------------

if [ ${task} = GD ] || [ ${task} = TD ] 
then

# fastqc zip to WES
     ${FASTQC} ${fastq1} -o ${wes}
     ${FASTQC} ${fastq2} -o ${wes}

# trim fastq
     ${TRIM} ${fastq1} ${fastq2} ${fastq_o1_p_t} ${fastq_o1_u_t}  ${fastq_o2_p_t} ${fastq_o2_u_t} ILLUMINACLIP:${soft}/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 HEADCROP:3 TRAILING:10 MINLEN:25
     ${FASTQC}  ${fastq_o1_p_t} -o ${wes}
     ${FASTQC}  ${fastq_o2_p_t} -o ${wes}

# make bam
     ${BWAMEM} -R "@RG\tID:${NameD}\tSM:${NameD}\tPL:illumina\tLB:lib1\tPU:unit1" -t 12 ${GENOME}  ${fastq_o1_p_t}  ${fastq_o2_p_t} | ${SAMVIEW} -bS - > ${bam}
     rm ${fastq_o1_p_t} ${fastq_o1_u_t} 
     rm ${fastq_o2_p_t} ${fastq_o2_u_t} 

# stats
     ${STATS} ${bam} > ${statstxt}

# sort bam
     ${SAMSORT} ${bam} -T ${prefixsort} -o ${sortbam}
     rm ${bam}
# rmdup bam
     ${SAMVIEW} -b -f 0x2 -q1 ${sortbam} | ${SAMRMDUP} - ${rmdupbam}
     rm ${sortbam}

# make bai
     ${SAMINDEX} ${rmdupbam} ${bai}
# make bam list
     ${RealignerTargetCreator} -o ${bamlist} -I ${rmdupbam}

# realign bam
     ${IndelRealigner} -I ${rmdupbam} -targetIntervals ${bamlist} -o ${realignedbam}
     rm ${bamlist} ${bai} ${rmdupbam}
     rm -r ${tempdir}

# fix bam
     ${FixMate} INPUT=${realignedbam} OUTPUT=${fixedbam} SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
     rm ${realignedbam} ${realignedbai}
     rm -r ${tempdir}

# make csv
     ${BaseRecalibrator} -I ${fixedbam} -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${csv}

# recal bam
     ${PrintReads} -I ${fixedbam} -BQSR ${csv} -o ${recalbam}
     rm ${csv} ${fixedbam} ${fixedbai}
     rm -r ${tempdir}

# coverage
     ${COVERAGE} -b ${recalbam} -a ${CaptureRegions} | grep '^all' > ${coveragetxt}
# zip
     ${FASTQC} ${recalbam} -o ${wes}


fi
# eo alignment 



## variantCalling ------------------------------------------------------------------------------------------------
if [ $task = VC ]
then

NameD=${case}_${num}_${task} 
NameGD=${case}_${num}_GD 
NameTD=${case}_${num}_TD 
recalbamGD=${wes}/${NameGD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
recalbamTD=${wes}/${NameTD}_output.sort.filtered.rmdup.realigned.fixed.recal.bam
snpvcf=${wes}/${NameD}.output.snp.vcf 
indelvcf=${wes}/${NameD}.output.indel.vcf 

   ${MPILEUP} ${recalbamGD} ${recalbamTD} | ${SOMATIC} --output-snp ${snpvcf} --output-indel ${indelvcf} --min-coverage-normal ${minCoverageNormal} --min-coverage-tumor ${minCoverageTumor} --tumor-purity ${TumorPurity} --min-var-freq ${minVarFreq} --min-freq-for-hom ${minFreqForHom} --output-vcf 1 --mpileup 1 

# Processing of somatic mutations

   ${PROCESSSOMATIC} ${snpvcf} --min-tumor-freq ${minTumorFreq}
   ${PROCESSSOMATIC} ${indelvcf} --min-tumor-freq ${minTumorFreq}

# FP Filter:  snp.Somatic.hc snp.LOH.hc snp.Germline.hc
# FP Filter:  indel.Somatic.hc indel.LOH.hc indel.Germline.hc


names1="snp indel"
for name1 in ${names1}
do

   if [ ${case} = somatic ]; then
   names2="Somatic LOH" 
   else
   names2="Somatic LOH Germline"
   fi

   for name2 in ${names2}
   do
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
   ${CONVERT2ANNOVAR2} ${hc_avi}   ${hc_vcf}
   ${CUT}              ${hc_avi} > ${hc_rci}
   ${BamReadcount}  -l ${hc_rci}   ${recalbam} > ${hc_rcs}
   ${VarScan} fpfilter ${hc_vcf}   ${hc_rcs} --output-file ${hc_fpf} --keep-failures 1 --min-ref-basequal ${minRefBasequal} --min-var-basequal ${minVarBasequal} --min-var-count ${minVarCount}
   done
done


data=${wes}
for name1 in ${names1}
do
# Annotation snp.Somatic.hc $data/NameD.output.snp.Somatic.hc.fpfilter.vcf
# Annotation indel.Somatic.hc $data/NameD.output.indel.Somatic.hc.fpfilter.vcf
hc_=${data}/${NameD}.output.${name1}.Somatic.hc 
hc_fpf=${data}/${NameD}.output.${name1}.Somatic.hc.fpfilter.vcf 
hc_T_avi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput 
hc_T_avi_multi=${data}/${NameD}.output.${name1}.Somatic.hc.TUMOR.avinput.hg19_multianno.csv 
   ${CONVERT2ANNOVAR} ${hc_} ${hc_fpf} -allsample
   ${TABLEANNOVAR}    ${hc_T_avi} ${soft}/annovar/humandb/ -protocol ${protocol} -buildver hg19 -operation ${argop} -csvout -otherinfo -remove -nastring NA
hc_snpeff=$data/${NameD}.output.$name1.Somatic.SnpEff.vcf
   ${SNPEFF} ${hc_fpf} > ${hc_snpeff}

if [ $case = somaticGermline ]; then
# Annotation snp.Germline.hc $data/NameD.output.snp.Germline.hc.fpfilter.vcf
# Annotation indel.Germline.hc $data/NameD.output.indel.Germline.hc.fpfilter.vcf
hc_=${data}/${NameD}.output.${name1}.Germline.hc 
hc_fpf=${data}/${NameD}.output.${name1}.Germline.hc.fpfilter.vcf 
hc_N_avi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput 
hc_N_avi_multi=${data}/${NameD}.output.${name1}.Germline.hc.NORMAL.avinput.hg19_multianno.csv 
   ${CONVERT2ANNOVAR} ${hc_} ${hc_fpf} -allsample
   ${TABLEANNOVAR}    ${hc_N_avi} ${soft}/annovar/humandb/ -protocol ${protocol} -buildver hg19 -operation ${argop} -csvout -otherinfo -remove -nastring NA
hc_N_snpeff=${data}/${NameD}.output.${name1}.NORMAL.SnpEff.vcf
      ${SNPEFF} ${hc_fpf} > ${hc_N_snpeff}
fi

# Annotation snp.LOH.hc
# Annotation indel.LOH.hc
hc_vcf=${data}/${NameD}.output.${name1}.LOH.hc.vcf 
hc_fpf=${data}/${NameD}.output.${name1}.LOH.hc.fpfilter.vcf 
hc_avi=${data}/${NameD}.output.${name1}.LOH.hc.avinput 
hc_avi_multi=${data}/${NameD}.output.${name1}.LOH.hc.avinput.hg19_multianno.csv 
   ${CONVERT2ANNOVAR3} ${hc_avi} ${hc_fpf} 
   ${TABLEANNOVAR}     ${hc_avi} ${soft}/annovar/humandb/ -protocol ${protocol} -buildver hg19 -operation ${argop} -csvout -otherinfo -remove -nastring NA
hc_L_snpeff=${data}/${NameD}.output.${name1}.LOH.SnpEff.vcf
   ${SNPEFF} ${hc_fpf} > ${hc_L_snpeff}
done


fi
# eo VC 



## CNV  ----------------------------------------------------------------------------------------------------------
if [ ${task} = CNV ]; then

output="${wes}/CNV"

if [ ! -d ${output} ]; then
   mkdir ${output}
fi

cat >  ${wes}/CNV_config.txt <<EOI
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

mateFile = ${wes}/${case}_${num}_TD_output.sort.filtered.rmdup.realigned.fixed.recal.bam
inputFormat = BAM
mateOrientation = FR

[control]

mateFile = ${wes}/${case}_${num}_GD_output.sort.filtered.rmdup.realigned.fixed.recal.bam
inputFormat = BAM
mateOrientation = FR

[target]

captureRegions = ${CaptureRegions}
EOI

export PATH=${PATH}:${SAMTOOLS}
${freec}-conf ${wes}/CNV_config.txt

fi
# eo CNV

## Report  -------------------------------------------------------------------------------------------------------
if [ ${task} = Report ]; then

cd ${ana}

${Rscript} ${ana}/Main.R ${case} ${num} ${4} ${5} ${mtb} ${RscriptPath} ${DatabasePath}

${Rscript} -e "library(knitr); knit('Report.Rnw')"
${soft}/pdflatex -interaction=nonstopmode Report.tex
${soft}/pdflatex -interaction=nonstopmode Report.tex

fi
# eo Report

echo "task ${task} for ${num} finished"
rm .STARTING_MARKER_${task} 
exit
