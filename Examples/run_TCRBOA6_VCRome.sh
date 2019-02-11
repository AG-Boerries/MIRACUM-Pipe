# createSet.sh and make_alignment_VC_CNV.sh need to be in the base directory
# first germline (GD), then tumor (TD)

./createSet.sh somaticGermline TCRBOA6_VCRome 'P110/TexasBiobank' 'P110/TexasBiobank' TCRBOA6-N-WEX.read TCRBOA6-T-WEX.read XX

cd somaticGermline_TCRBOA6_VCRome

./run_jobs.sh > out &
