#!/bin/bash

internal_hdf5="$(basename ${1//.hdf5/}).chrest.00000.hdf5"
chrest_data="${1//.hdf5/}.chrest/$internal_hdf5"
intermediate_csv="$(dirname $chrest_data)/chrest_data.csv.gz"
final_output="$(dirname $chrest_data)/chrest_data-rebalanced.csv.gz"

if [ -e "$final_output" ]; then 
	echo exiting, final output already exists
	return 0
fi
echo proceeding! with input hdf5: $1

# processes 1 ablate file
python ablateData.py --file "$1" --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_Yi:Yi monitor_zMix:zmix
python wrangle_chrest_dataset_for_UQ.py "$chrest_data"
Rscript Rebalance_chrest_data.R "$intermediate_csv"
rm "$intermediate_csv" "$chrest_data"
#"$(dirname $chrest_data)/chrest_data.csv"
#/user/dwyerdei/data/ablate_data_processing/raw_ablate_data/flowField_mixtureFraction/flowField_mixtureFraction.00199.chrest/chrest_data.csv
