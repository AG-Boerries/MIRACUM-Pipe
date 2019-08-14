#!/usr/bin/env bash

# TODO: finalize this
# required params: task


if [ ! -d ${DIR_TMP} ]; then
  mkdir ${DIR_TMP}
fi

# SAMPLE
NameD=${case}_${num}_${task}
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
${FASTQC} ${fastq1} -o ${wes}
${FASTQC} ${fastq2} -o ${wes}

# trim fastq
${TRIM} ${fastq1} ${fastq2} ${fastq_o1_p_t} ${fastq_o1_u_t} ${fastq_o2_p_t} ${fastq_o2_u_t} ILLUMINACLIP:${TrimmomaticAdapter}/TruSeq3-PE-2.fa:2:30:10 HEADCROP:3 TRAILING:10 MINLEN:25
${FASTQC} ${fastq_o1_p_t} -o ${wes}
${FASTQC} ${fastq_o2_p_t} -o ${wes}

# make bam
${BWAMEM} -R "@RG\tID:${NameD}\tSM:${NameD}\tPL:illumina\tLB:lib1\tPU:unit1" -t 12 ${GENOME} ${fastq_o1_p_t} ${fastq_o2_p_t} | ${SAMVIEW} -bS - >${bam}

# stats
${STATS} ${bam} >${statstxt}

# sort bam
${SAMSORT} ${bam} -T ${prefixsort} -o ${sortbam}

# rmdup bam
${SAMVIEW} -b -f 0x2 -q1 ${sortbam} | ${SAMRMDUP} - ${rmdupbam}

# make bai
${SAMINDEX} ${rmdupbam} ${bai}

# make bam list
${RealignerTargetCreator} -o ${bamlist} -I ${rmdupbam}

# realign bam
${IndelRealigner} -I ${rmdupbam} -targetIntervals ${bamlist} -o ${realignedbam}

# fix bam
${FixMate} INPUT=${realignedbam} OUTPUT=${fixedbam} SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

# make csv
${BaseRecalibrator} -I ${fixedbam} -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -o ${csv}

# recal bam
${PrintReads} -I ${fixedbam} -BQSR ${csv} -o ${recalbam}

# coverage
${COVERAGE} -b ${recalbam} -a ${CaptureRegions} | grep '^all' >${coveragetxt}

# zip
${FASTQC} ${recalbam} -o ${wes}

# TODO: required in all 4 files
echo "task ${task} finished"
rm .STARTING_MARKER_${task}
exit