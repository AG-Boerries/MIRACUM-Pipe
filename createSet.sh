#!/usr/bin/env bash

# script to create working directory and sh job submissions for the WES samples
# Version 31.07.2019

DIR_SCRIPT=$(
  cd "$(dirname "${BASH_SOURCE[0]}")"
  pwd -P
)

## load settings
# shellcheck source=global.sh
source "${DIR_SCRIPT}"/global.sh

########################
## Example run script ##
########################
# example
#  ./createSet.sh somaticGermline ID 'folder/containing/germline' 'folder/containing/tumor' filename_germline_without_file_extension filename_tumor_without_file_extension XY
# examples for paired-end fastq files (AS-258824-LR-37679_R1.fastq.gz and AS-258824-LR-37679_R2.fastq.gz)
# filename_germline_without_file_extension: AS-258824-LR-37679_R

# somatic: do not do call germline varaints only use germline for identifying somatic variants
# somaticGermline: call somatic variants as well as germline variants

possible_cases=("somatic somaticGermline")
possible_sex=("XX XY")

function usage() {
  echo "usage: createSet.sh -s sex -c case -n num -fg filename -ft filename -g folder -t folder [-h]"
  echo "  -d  dir              specify relative folder of patient"
  echo "  -h                   show this help screen"
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


# load patient yaml
sex=$(get_config_value sex "${DIR_PATIENT}")
if [[ "$(get_config_value annotation.germline "${DIR_PATIENT}")" = "True" ]]; then
  case=somaticGermline
else
  case=somatic
fi

for value in "${possible_sex[@]}"
do
  [[ "${sex}" = "${value}" ]] && \
    echo "unknown sex: ${sex}" && \
    echo "use one of the following values: $(join_by ' ' ${possible_sex})" && \
    exit 1
done

#case=$1                                     # somatic or somaticGermline
#sex=$7                                      # gender

# gets killed with
#num=$2                                      # Patient ID
#xx1=$3                                      # folder/containing/germline
#xx2=$4                                      # folder/containing/tumor
#f1=$5                                       # filename_germline_without_file_extension       ## Inputf1=$homedata/ngs/$xx1/fastq/$f1 -> R1,R2
#f2=$6                                       # filename_tumor_without_file_extension          ## Inputf2=$homedata/ngs/$xx2/fastq/$f2 -> R1,R2


##################################################################################################################
## Parameters which need to be adjusted to the local environment

DIR_OUTPUT="${DIR_SCRIPT}/assets/output"

mtb="${DIR_OUTPUT}/${case}_${DIR_PATIENT}" # path to output folder, subfolder with type of analysis and ID is automatically created
wes="${mtb}/WES"                   # subfolder containing the alignemnt, coverage, copy number variation and variant calling results
ana="${mtb}/Analysis"              # subfolder containing PDF Report, annotated copy number variations and annotated variants

##################################################################################################################

##########
## MAIN ##
##########

# reate folder if not exists
if [ ! -d ${mtb} ]; then
  mkdir -p ${wes}
  mkdir ${ana}
fi

# cycle on tasks

for d in GD TD VC CNV Report; do
  # create sh scripts
  dname=${case}_${DIR_PATIENT}_${d}
  jobname=${DIR_PATIENT}_${d}
  runname=${mtb}/run_${dname}.sh

  cat >${runname} <<EOI
#!/usr/bin/env bash

EOI
  # TODO: split make alignment to 4 files
  # create sh run scripts
  case ${d} in
  GD)
    echo "bash ${DIR_SCRIPT}/make_alignment_VC_CNV.sh -t ${d} -d ${DIR_PATIENT}" >> "${runname}"
    ;;
  TD)
    echo "bash ${DIR_SCRIPT}/make_alignment_VC_CNV.sh -t ${d} -d ${DIR_PATIENT}" >> "${runname}"
    ;;
  VC)
    echo "bash ${DIR_SCRIPT}/make_alignment_VC_CNV.sh -t ${d} -d ${DIR_PATIENT}" >> "${runname}"
    ;;
  CNV)
    echo "bash ${DIR_SCRIPT}/make_alignment_VC_CNV.sh -t ${d} -d ${DIR_PATIENT}" >> "${runname}"
    ;;
  Report)
    echo "bash ${DIR_SCRIPT}/make_alignment_VC_CNV.sh -t ${d} -d ${DIR_PATIENT}" >> "${runname}"
    ;;
  esac
  chmod a+x ${runname}
done
# end of cycle on d

# create jobs submissions script
cat >${mtb}/run_jobs.sh <<EOI2
#!/usr/bin/env bash

maxhours=48
# -------------------------------------------
date
echo "Submitting  GD TD"

for task in GD TD; do
   dname=${case}_${DIR_PATIENT}_\${task}
   cd ${mtb}
   if [ -f .STARTING_MARKER_\${task} ]; then
       exit "Previous job uncompleted. Aborting!"
   else
       echo "STARTING \$dname"; echo "\$dname" > .STARTING_MARKER_\${task}
       sh run_\${dname}.sh
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
   dname=${case}_${DIR_PATIENT}_\${task}
   cd ${mtb}
   if [ -f .STARTING_MARKER_\${task} ]; then
       exit "Previous job uncompleted. Aborting!"
   else
       echo "STARTING \$dname"; echo "\$dname" > .STARTING_MARKER_\${task}
       sh run_\${dname}.sh
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
	dname=${case}_${DIR_PATIENT}_\${task}
	cd ${mtb}
	if [ -f .STARTING_MARKER_\${task} ]; then
	    exit "Previous job uncompleted. Aborting!"
	else
	    echo "STARTING \${dname}"; echo "\${dname}" > .STARTING_MARKER_\${task}
	    sh run_\${dname}.sh
		fi
done

echo "Waiting for Report"
count=1
while [ \${count} -lt \${maxhours} ]; do
	if [ -e ${mtb}/.STARTING_MARKER_Report ]; then
     sleep 1h
     date
     (( count++ ))
   else
     break
   fi
done

echo "Finished Report"
date
echo "Finished all jobs for ${DIR_PATIENT} "
# -------------------------------------------

cd ${mtb} && touch .processed

exit
EOI2
chmod a+x ${mtb}/run_jobs.sh

exit
# --------------------------------------------------------------------------------------------------------
