# createSet.sh and make_alignment_VC_CNV.sh need to be in the base directory
# first germline (GD), then tumor (TD)

./createSet.sh somaticGermline CYDTY9 'P117/data/180713_ST-K00265_0163_BHWNKTBBXX/AS-253876-LR-37206' 'P117/data/180713_ST-K00265_0163_BHWNKTBBXX/AS-253877-LR-37206' AS-253876-LR-37206_R AS-253877-LR-37206_R XY

cd somaticGermline_CYDTY9

./run_jobs.sh > out &