#!/bin/bash

#data_dir='/Users/dwyerdeighan/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleDiffusionFlame/flowField_mixtureFraction'
#data_dir='/Users/dwyerdeighan/ClionProjects/ablate/ablateInputs/chemtabDiffusionFlame/_1DSampleChemTabDiffusionFlame/flowField_chemTab'
#data_dir='/user/dwyerdei/data/ablate_data_processing/raw_ablate_data/flowField_mixtureFraction/'
#data_dir='/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'

if ! [ $data_dir ]; then
    echo Error: env var "data_dir" is unspecified! Please provide it.
    return 1 || exit 1
fi

chemtab_output="$(echo $data_dir | grep -i -e '.*chemTab$')"
export chemtab_output # so we have access inside process_single_ablate_file.sh

# do dependency injection for srun wrapper (or empty if local)
if [ $(which srun) ]; then
	run_1task="srun -N 1 -n 1"
else
	run_1task='' # dummy for when we are doing this locally
fi

process_ablate="$run_1task ./process_single_ablate_file.sh"

# NOTE: much easier than a bash loop!! e.g.: map echo 1 2 3
map() { cmd="$1"; shift; for x in "$@"; do eval "$cmd $x"; done }
amap() { # asynchronous version of map!
    cmd="$1"; shift; for x in "$@"; do eval "$cmd $x" & done
} # NOTE: we used to have the entire for-loop inside a sub-shell (i.e. inside parenthesis (for ... done)), why was this? Is it still important?

files=("$data_dir/"*.hdf5)
echo all files: ${files[*]}
echo num files: ${#files[*]}
amap "$process_ablate" ${files[*]}
wait
echo job done sub-sub!!
