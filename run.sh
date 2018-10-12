# createSet.sh and make_alignment_VC_CNV.sh need to be in the base directory
# first germline (GD), then tumor (TD)

./createSet.sh somaticGermline ID 'folder/containing/germline' 'folder/containing/tumor' filename_germline_without_file_extension filename_tumor_without_file_extension XY

cd somaticGermline_ID

./run_jobs.sh > out &