#!/bin/bash
#PBS -S /bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=12,walltime=48:00:00
#PBS -m a
#PBS -N TCRBOA6_VCRome_GD
#PBS -e /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/err.somaticGermline_TCRBOA6_VCRome_GD.log
#PBS -o /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/out.somaticGermline_TCRBOA6_VCRome_GD.log
#PBS -d /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/

bash /home/miracum/Analysis/somaticGermline_TCRBOA6_VCRome/make_alignment_VC_CNV.sh somaticGermline GD TCRBOA6_VCRome P110/TexasBiobank TCRBOA6-N-WEX.read 
