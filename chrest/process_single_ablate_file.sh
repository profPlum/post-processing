#!/bin/bash

#should_rebalance=T # set this flag to T or "" (empty)
#chemtab_output=T # should be passed as an argument!

echo first arg: $1
echo should_rebalance: $should_rebalance

internal_hdf5="$(basename ${1//.hdf5/}).chrest.00000.hdf5"
chrest_folder="${1//.hdf5/}.chrest"
chrest_data="$chrest_folder/$internal_hdf5"
intermediate_csv="$(dirname $chrest_data)/chrest_data.csv.gz"

final_output="$intermediate_csv"
if [ $should_rebalance ]; then
    final_output="$(dirname $chrest_data)/chrest_data-rebalanced.csv.gz"
fi

echo final_output: $final_output
echo chrest_data: $chrest_data

if [ -e "$final_output" ]; then 
    echo exiting, final output already exists
    return 0
fi
echo proceeding! with input hdf5: $1

echo checking if we have chemtab output...
if [ "$chemtab_output" ]; then
    echo answer: True
    fields_list="monitor_energySource:souener monitor_Yi:Yi monitor_progressSource:soucpv"
else
    echo answer: False
    fields_list="monitor_energySource:souener monitor_yiSource:souspec monitor_Yi:Yi monitor_zMix:zmix"
fi

# processes 1 ablate file
python ablateData.py --file "$1" --fields $fields_list 
python wrangle_chrest_dataset_for_UQ.py --file "$chrest_data" --fields $fields_list

[ $should_rebalance ] && Rscript Rebalance_chrest_data.R "$intermediate_csv"

# move final output out of the directory (temporarily)
mv "$final_output" $(dirname $final_output)/..
rm "$chrest_folder"/* # clear the directoy
mv $(dirname $final_output)/../$(basename $final_output) "$final_output" # move final output back to where it was
# NOTE: it is necessary to retain th folder so that it is identifiable among other chrest_data[-rebalanced].csv.gz's
