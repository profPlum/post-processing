#!/bin/bash

echo first arg: $1
final_output="${1%.*}.csv.gz"
intermediate_csv="$final_output"

#final_output="$intermediate_csv"
if ((should_rebalance)); then
    final_output="${1%.*}-rebalanced.csv.gz"
fi

echo final_output: $final_output

if [[ -e "$final_output" ]]; then
    echo exiting, final output already exists
    return 0 &> /dev/null || exit 0
fi
echo proceeding! with input hdf5: $1

echo checking if we have chemtab output...
fields_list="monitor_energySource:souener monitor_Yi:Yi aux_pressure:pressure aux_temperature:temp"
if [[ $chemtab_output ]]; then
    echo answer: True
    fields_list="$fields_list solution_densityProgress:density_ aux_Progress:mass_ monitor_progressSource:source_"
    # results in: mass_CPV_PC* & source_CPV_PC*, verified 10/4/23
else
    echo answer: False
    fields_list="$fields_list solution_densityYi:density_ monitor_yiSource:souspec monitor_zMix:zmix"
fi

# processes 1 ablate file
base_1="$(basename $1)"
echo base_1: ${base_1#*.}
domain_file="$(dirname $1)/../domain/domain.${base_1#*.}"
if ((should_rebalance)); then
    should_rebalance=--rebalance-flame
else
    should_rebalance=
fi

python wrangle_ablate_dataset_for_ChemTab.py $should_rebalance --file "$1" "$domain_file" --fields $fields_list #--n-cubes-per-dim 25

echo input_fn: "$1"
echo aux_fn: "$domain_file"

#if [ $should_rebalance ]; then
#    Rscript Rebalance_chrest_data.R "$intermediate_csv" #"$domain_file"
#    rm "$intermediate_csv"
#fi
