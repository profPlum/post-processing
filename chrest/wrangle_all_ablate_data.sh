#!/bin/bash

# Debugging Breadcrumbs: problems were in data & with amap using a subshell which made `wait` not work and job terminal early.
# Also submitting each step as a job was a bad idea, high failure rate & long submission wait times.

data_dir='/user/dwyerdei/data/ablate_data_processing/raw_ablate_data/flowField_mixtureFraction/'
#'/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'

# do dependency injection for srun wrapper (or empty if local)
if [ -z $(which srun) ]; then
    run_1task='' # dummy for when we are doing this locally
else
    run_1task="srun -N 1 -n 1"
fi

################### Wrangle Ablate: ################### 

amap() { # asynchronous version of map!
    cmd="$1"; shift; for x in "$@"; do eval "$cmd $x" &>> amap_out.log & done
} # NOTE: we used to have the entire for-loop inside a sub-shell (i.e. inside parenthesis (for ... done)), why was this? Is it still important?

ablate2chrest() {
	run_1task python ablateData.py --file "$1" --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_yi:Yi monitor_zMix:zmix
}
## TODO: comment once done debugging!
#ablate2chrest "$data_dir"/flowField.00019.hdf5
#echo if this worked without errors then come back and uncomment remainder of the script!

# TODO: uncomment once done debugging!
amap ablate2chrest "$data_dir"*.hdf5
wait

################### Wrangle Chrest: ###################

chrest2csv() {
	run_1task python wrangle_chrest_dataset_for_UQ.py "$1"
}
## TODO: comment once done debugging!
#chrest2csv "$data_dir"/flowField.00019.hdf5 

amap chrest2csv "${data_dir}/flowField_mixtureFraction."*'.chrest/flowField_mixtureFraction.'*'.chrest.'*'.hdf5'
wait
#####################################################a
