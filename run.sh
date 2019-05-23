#!/usr/bin/env bash

# createSet.sh and make_alignment_VC_CNV.sh need to be in the base directory
# first germline (GD), then tumor (TD)

dir_tumor=$2
dir_germline=$3
gender=$4 # XX or XY
filename_germline_without_file_extension=$5
filename_tumor_without_file_extension=$6
method=$7 # somatic or somaticGermline
ID=$8

./createSet.sh ${method} ${ID} ${dir_tumor} ${dir_germline} ${filename_germline_without_file_extension} ${filename_tumor_without_file_extension} ${gender}

cd ${method}_${ID}

./run_jobs.sh > out &