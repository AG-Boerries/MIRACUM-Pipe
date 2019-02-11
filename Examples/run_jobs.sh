#!/bin/bash

maxhours=48
# -------------------------------------------
    date
    echo "Submitting  GD TD"

for task in GD TD; do
   dname=somaticGermline_TCRBOA6_VCRome_${task}
   cd /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome
   if [ -f .STARTING_MARKER_${task} ]; then
       exit "Previous job uncompleted. Aborting!"
   else
       echo "STARTING $dname"; echo "$dname" > .STARTING_MARKER_${task}
       qsub run_${dname}.sh
   fi
done

    echo "Waiting for GD TD"

# wait no more than maxhours
count=1
while [ ${count} -lt ${maxhours} ]; do
   if [ -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/.STARTING_MARKER_GD ] || [ -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/.STARTING_MARKER_TD ]; then
     sleep 1h
     date
     (( count++ ))
   else
     break
   fi
done
    if [ ${count} -ge ${maxhours} ]; then
    date
    echo "Waited to long for job results for GD TD; exiting" 
    exit
    fi

    echo "Finished GD TD"

# -------------------------------------------

    date
    echo "Submitting  VC CNV"

for task in VC CNV; do
   dname=somaticGermline_TCRBOA6_VCRome_${task}
   cd /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome
   if [ -f .STARTING_MARKER_${task} ]; then
       exit "Previous job uncompleted. Aborting!"
   else
       echo "STARTING $dname"; echo "$dname" > .STARTING_MARKER_${task}
       qsub run_${dname}.sh
   fi
done

    echo "Waiting for VC CNV"
count=1
while [ ${count} -lt ${maxhours} ]; do
   if [ -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/.STARTING_MARKER_VC ] || [ -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/.STARTING_MARKER_CNV ]; then
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
	dname=somaticGermline_TCRBOA6_VCRome_${task}
	cd /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome
	if [ -f .STARTING_MARKER_${task} ]; then
	    exit "Previous job uncompleted. Aborting!"
	else
	    echo "STARTING ${dname}"; echo "${dname}" > .STARTING_MARKER_${task}
	    qsub run_${dname}.sh
		fi
done

	 echo "Waiting for Report"
count=1
while [ ${count} -lt ${maxhours} ]; do
	if [ -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/.STARTING_MARKER_Report]; then
     sleep 1h
     date
     (( count++ ))
   else
     break
   fi
done

    echo "Finished Report"
    date
    echo "Finished all jobs for TCRBOA6_VCRome "
# -------------------------------------------
exit
