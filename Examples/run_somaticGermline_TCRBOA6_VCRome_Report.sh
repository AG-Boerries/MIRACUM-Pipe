#!/bin/bash
#PBS -S /bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=12,walltime=48:00:00
#PBS -m a
#PBS -N TCRBOA6_VCRome_Report
#PBS -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/err.somaticGermline_TCRBOA6_VCRome_Report.log
#PBS -o /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/out.somaticGermline_TCRBOA6_VCRome_Report.log
#PBS -d /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/

bash /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/make_alignment_VC_CNV.sh somaticGermline Report TCRBOA6_VCRome TCRBOA6-N-WEX.read TCRBOA6-T-WEX.read /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome 
