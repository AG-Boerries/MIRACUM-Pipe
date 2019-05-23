#!/usr/bin/env bash

# script to create working directory and qsub job submissions for the WES samples
# Version 05.02.2019

########################
## Example run script ##
########################
# example
#  ./createSet.sh somaticGermline ID 'folder/containing/germline' 'folder/containing/tumor' filename_germline_without_file_extension filename_tumor_without_file_extension XY
# examples for paired-end fastq files (AS-258824-LR-37679_R1.fastq.gz and AS-258824-LR-37679_R2.fastq.gz)
# filename_germline_without_file_extension: AS-258824-LR-37679_R

# somatic: do not do call germline varaints only use germline for identifying somatic variants
# somaticGermline: call somatic variants as well as germline variants

case=$1; # somatic or somaticGermline
num=$2; # Patient ID
xx1=$3; # folder/containing/germline
xx2=$4; # folder/containing/tumor
f1=$5; # filename_germline_without_file_extension          ## Inputf1=$homedata/ngs/$xx1/fastq/$f1 -> R1,R2
f2=$6; # filename_tumor_without_file_extension          ## Inputf2=$homedata/ngs/$xx2/fastq/$f2 -> R1,R2
sex=$7; # gender

##################################################################################################################
## Parameters which need to be adjusted to the local environment
 # TODO: same volumes as in make_alignment
mtb="/path/to/output/folder/${case}_${num}" # path to output folder, subfolder with type of analysis and ID is automatically created
wes="${mtb}/WES" # subfolder containing the alignemnt, coverage, copy number variation and variant calling results
ana="${mtb}/Analysis" # subfolder containing PDF Report, annotated copy number variations and annotated variants
 
##################################################################################################################
 

##########
## MAIN ##
##########

pathtothis="."
if [ ! -d ${mtb} ]; then
  mkdir ${mtb}
  cp ${pathtothis}/make_alignment_VC_CNV.sh ${mtb}
  mkdir ${wes}
  mkdir ${ana}
  cp ${pathtothis}/RScripts/Main.R ${ana}
  cp ${pathtothis}/RScripts/Report.Rnw ${ana}
fi

# cycle on tasks

for d in GD TD VC CNV Report
do

# create sh scripts
dname=${case}_${num}_${d}
jobname=${num}_${d}
runname=${mtb}/run_${dname}.sh

cat > ${runname} <<EOI
#!/bin/bash
#PBS -S /bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=12,walltime=48:00:00
#PBS -m a
#PBS -N $jobname
#PBS -e ${mtb}/err.${dname}.log
#PBS -o ${mtb}/out.${dname}.log
#PBS -d ${mtb}/

EOI

# create sh run scripts
    case ${d} in
    GD)
   echo "bash ${mtb}/make_alignment_VC_CNV.sh ${case} ${d} ${num} ${xx1} ${f1} " >> ${runname} 
   ;;
    TD)
   echo "bash ${mtb}/make_alignment_VC_CNV.sh ${case} ${d} ${num} ${xx2} ${f2} " >> $runname 
   ;;
    VC)
   echo "bash ${mtb}/make_alignment_VC_CNV.sh ${case} ${d} ${num} " >> ${runname} 
   ;;
    CNV)
   echo "bash ${mtb}/make_alignment_VC_CNV.sh ${case} ${d} ${num} ${sex}" >> ${runname} 
   ;;
	 Report)
   echo "bash ${mtb}/make_alignment_VC_CNV.sh ${case} ${d} ${num} ${f1} ${f2} ${mtb} " >> ${runname}
   ;;
    esac
chmod a+x ${runname}
done
# cycle on d

# create jobs submissions script
cat > ${mtb}/run_jobs.sh <<EOI2
#!/bin/bash

maxhours=48
# -------------------------------------------
    date
    echo "Submitting  GD TD"

for task in GD TD; do
   dname=${case}_${num}_\${task}
   cd $mtb
   if [ -f .STARTING_MARKER_\${task} ]; then
       exit "Previous job uncompleted. Aborting!"
   else
       echo "STARTING \$dname"; echo "\$dname" > .STARTING_MARKER_\${task}
       qsub run_\${dname}.sh
   fi
done

    echo "Waiting for GD TD"

# wait no more than maxhours
count=1
while [ \${count} -lt \${maxhours} ]; do
   if [ -e ${mtb}/.STARTING_MARKER_GD ] || [ -e ${mtb}/.STARTING_MARKER_TD ]; then
     sleep 1h
     date
     (( count++ ))
   else
     break
   fi
done
    if [ \${count} -ge \${maxhours} ]; then
    date
    echo "Waited to long for job results for GD TD; exiting" 
    exit
    fi

    echo "Finished GD TD"

# -------------------------------------------

    date
    echo "Submitting  VC CNV"

for task in VC CNV; do
   dname=${case}_${num}_\${task}
   cd ${mtb}
   if [ -f .STARTING_MARKER_\${task} ]; then
       exit "Previous job uncompleted. Aborting!"
   else
       echo "STARTING \$dname"; echo "\$dname" > .STARTING_MARKER_\${task}
       qsub run_\${dname}.sh
   fi
done

    echo "Waiting for VC CNV"
count=1
while [ \${count} -lt \${maxhours} ]; do
   if [ -e ${mtb}/.STARTING_MARKER_VC ] || [ -e ${mtb}/.STARTING_MARKER_CNV ]; then
     sleep 1h
     date
     (( count++ ))
   else
     break
   fi
done

    echo "Finished VC CNV"
	
# -------------------------------------------

	date
	echo "Submitting  Report"

for task in Report; do
	dname=${case}_${num}_\${task}
	cd $mtb
	if [ -f .STARTING_MARKER_\${task} ]; then
	    exit "Previous job uncompleted. Aborting!"
	else
	    echo "STARTING \${dname}"; echo "\${dname}" > .STARTING_MARKER_\${task}
	    qsub run_\${dname}.sh
		fi
done

	 echo "Waiting for Report"
count=1
while [ \${count} -lt \${maxhours} ]; do
	if [ -e ${mtb}/.STARTING_MARKER_Report]; then
     sleep 1h
     date
     (( count++ ))
   else
     break
   fi
done

    echo "Finished Report"
    date
    echo "Finished all jobs for $num "
# -------------------------------------------
exit
EOI2
chmod a+x ${mtb}/run_jobs.sh




exit
# --------------------------------------------------------------------------------------------------------
