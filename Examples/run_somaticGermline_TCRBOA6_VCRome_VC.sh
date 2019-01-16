#!/bin/bash
#PBS -S /bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=12,walltime=12:00:00
#PBS -m a
#PBS -N TCRBOA6_VCRome_VC
#PBS -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/err.somaticGermline_TCRBOA6_VCRome_VC.log
#PBS -o /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/out.somaticGermline_TCRBOA6_VCRome_VC.log
#PBS -d /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/

bash /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/make_alignment_VC_CNV.sh somaticGermline VC TCRBOA6_VCRome 
