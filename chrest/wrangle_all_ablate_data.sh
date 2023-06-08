#!/bin/bash

source /user/dwyerdei/.bashrc
data_dir='/user/dwyerdei/data/ablate_data_processing/raw_ablate_data/flowField_mixtureFraction/'
#'/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'

echo removing all existing chrest csv files!
rm "$data_dir"/*.chrest/chrest_data.csv "$data_dir"/*.chrest/*.hdf5 "$data_dir"/*.chrest/*.xdmf

# do dependency injection for srun wrapper (or empty if local)
if [ -z $(which srun) ]; then
	run_1task='' # dummy for when we are doing this locally
else
	run_1task="srun -N 1 -n 1"
fi

process_ablate="$run_1task ./process_single_ablate_file.sh"

amap() { # asynchronous version of map!
    cmd="$1"; shift; for x in "$@"; do eval "$cmd $x" &>> amap_out.log & done
} # NOTE: we used to have the entire for-loop inside a sub-shell (i.e. inside parenthesis (for ... done)), why was this? Is it still important?

amap "$process_ablate" "$data_dir"*.hdf5
wait
echo job done sub-sub!!
