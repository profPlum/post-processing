#!/bin/bash

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

alias xargs="xargs -n1 -P$(getconf _NPROCESSORS_ONLN)"

files=("$data_dir/"*.hdf5)
echo all files: ${files[*]}
echo num files: ${#files[*]}
echo "${files[*]}" | xargs -n1 -P$(getconf _NPROCESSORS_ONLN) $process_ablate
wait
echo job done sub-sub!!
