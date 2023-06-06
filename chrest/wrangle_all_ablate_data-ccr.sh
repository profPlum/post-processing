#!/bin/bash

# NOTE: debugging breadcrumb-- I think you were done debugging? It was a pretty simple change to make the srun calls independent...
# It seems the problem was a problem with the data, but to be safe I left the debugging temp changes in here incase you need to debug something else too...

source /user/dwyerdei/.bashrc
data_dir='/user/dwyerdei/data/ablate_data_processing/raw_ablate_data/flowField_mixtureFraction/'
#'/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'

echo removing all existing chrest folders!
rm -r "$data_dir"/*.chrest

# do dependency injection for srun wrapper (or empty if local)
if [ -z $(which srun) ]; then
	run_1task='' # dummy for when we are doing this locally
else
	# simple srun wrapper that just runs 1 task! (i.e. allows reuse of amap & map!)
	#slurmi_qos=general-compute
	run_1task="srun -N 1 -n 1" # --time=00:15:00 --qos=$slurmi_qos --partition=$slurmi_qos"
fi

#run_1task='srun -N 1 -n 1'
#run_1task='' # dummy for local execution

process_ablate() {
	#ablate2chrest="python ablateData.py --file '$1' --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_Yi:Yi monitor_zMix:zmix"
	#chrest_data="${1//.hdf5/}.00000.chrest.hdf5"
	#chrest2csv="python wrangle_chrest_dataset_for_UQ.py '$chrest_data'"
	$run_1task ./process_single_ablate_file.sh "$1" #bash -c "$ablate2chrest; $chrest2csv"
	echo starting 1 job
}

amap() { # asynchronous version of map!
    cmd="$1"; shift; for x in "$@"; do eval "$cmd $x" &>> amap_out.log & done
} # NOTE: we used to have the entire for-loop inside a sub-shell (i.e. inside parenthesis (for ... done)), why was this? Is it still important?

## Iterate over each configuration file in the directory
#for config_file in "$data_dir"*.hdf5; do
#	process_ablate "$config_file" &
#	#$run_1task ./process_single_ablate_file.sh "$config_file" &
#done
#wait  # Wait for all the background processes to finish

amap process_ablate "$data_dir"*.hdf5
wait
echo done waiting!!

#################### Wrangle Ablate: ################### 
#ablate2chrest() {
#	$run_1task python ablateData.py --file "$1" --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_Yi:Yi monitor_zMix:zmix
#}
### TODO: comment once done debugging!
##ablate2chrest "$data_dir"/flowField_mixtureFraction.00019.hdf5
##echo if this worked without errors then come back and uncomment remainder of the script!
#
#amap ablate2chrest "$data_dir"*.hdf5
#wait
#echo done!!
#ls "${data_dir}/flowField_mixtureFraction."*'.chrest/flowField_mixtureFraction.'*'.chrest.'*'.hdf5'
#
##sleep $((60*5))

##################### Wrangle Chrest: ###################
#
#chrest2csv() {
#	$run_1task python wrangle_chrest_dataset_for_UQ.py "$1"
#}
##path="$data_dir"/flowField_mixtureFraction.00019.chrest/*.hdf5
##echo path = $path
##[ -e "$path" ] && echo it exists!
##chrest2csv "$path" 
#
#amap chrest2csv "${data_dir}/flowField_mixtureFraction."*'.chrest/flowField_mixtureFraction.'*'.chrest.'*'.hdf5'
#wait
#########################################################
