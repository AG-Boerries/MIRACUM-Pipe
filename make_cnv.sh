#!/usr/bin/env bash

output="${wes}/CNV"

if [ ! -d ${output} ]; then
  mkdir ${output}
fi

cat >${wes}/CNV_config.txt <<EOI
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
