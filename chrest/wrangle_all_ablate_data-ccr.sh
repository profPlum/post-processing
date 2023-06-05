#!/bin/bash

# NOTE: debugging breadcrumb-- I think you were done debugging? It was a pretty simple change to make the srun calls independent...
# It seems the problem was a problem with the data, but to be safe I left the debugging temp changes in here incase you need to debug something else too...

data_dir='/user/dwyerdei/data/ablate_data_processing/raw_ablate_data/flowField'
#'/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'

# TODO: use general-compute once done debugging! 
# simple srun that just runs 1 task! (i.e. allows reuse of amap & map!)
slurmi_qos=debug #general-compute
alias run_1task='srun --qos=$slurmi_qos --partition=$slurmi_qos'
#alias run_1task='' # dummy for when we are doing this locally

################### Wrangle Ablate: ################### 

ablate2chrest() {
	run_1task python ablateData.py --file "$1" --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_yi:Yi monitor_zMix:zmix
}
# TODO: comment once done debugging!
ablate2chrest "$data_dir"/flowField.00619.hdf5
echo if this worked without errors then come back and uncomment remainder of the script!

## TODO: uncomment once done debugging!
#amap ablate2chrest "$data_dir"*.hdf5
#wait
#
#################### Wrangle Chrest: ###################
#
#chrest2csv() {
#	run_1task python wrangle_chrest_dataset_for_UQ.py "$1"
#}
#chrest2csv "$data_dir"/flowField.00619.hdf5 
#
## TODO: uncomment once done debugging!
#amap chrest2csv "${data_dir}/flowField_mixtureFraction."*'.chrest/flowField_mixtureFraction.'*'.chrest.'*'.hdf5'
#wait
#######################################################
