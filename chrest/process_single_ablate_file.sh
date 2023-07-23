#!/bin/bash

#should_rebalance=T # set this flag to T or "" (empty)
#chemtab_output=T # should be passed as an argument!

echo first arg: $1
echo should_rebalance: $should_rebalance

final_output="${1%.*}.csv.gz"
intermediate_csv="$final_output"

#final_output="$intermediate_csv"
if [ $should_rebalance ]; then
    final_output="${final_output%%.*}-rebalanced.csv.gz"
fi

echo final_output: $final_output
echo chrest_data: $chrest_data

if [ -e "$final_output" ]; then 
    echo exiting, final output already exists
    return 0 || exit 0 &> /dev/null
fi
echo proceeding! with input hdf5: $1

echo checking if we have chemtab output...
fields_list="monitor_energySource:souener monitor_Yi:Yi aux_temperature:temp"
if [ "$chemtab_output" ]; then
    echo answer: True
    fields_list="$fields_list monitor_progressSource:soucpv"
else
    echo answer: False
    fields_list="$fields_list monitor_yiSource:souspec monitor_zMix:zmix"
fi

# processes 1 ablate file
base_1="$(basename $1)"
echo base_1: ${base_1#*.}
domain_file="$(dirname $1)/../domain/domain.${base_1#*.}"
if [ $should_rebalance ]; then
    should_rebalance=--rebalance-flame
fi

python wrangle_chrest_dataset_for_UQ.py $should_rebalance --file "$1" "$domain_file" --fields $fields_list #--n-cubes-per-dim 25

echo input_fn: "$1"
echo aux_fn: "$domain_file"

#if [ $should_rebalance ]; then
#    Rscript Rebalance_chrest_data.R "$intermediate_csv" #"$domain_file"
#    rm "$intermediate_csv"
#fi
