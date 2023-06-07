#!/bin/bash

internal_hdf5="$(basename ${1//.hdf5/}).chrest.00000.hdf5"
chrest_data="${1//.hdf5/}.chrest/$internal_hdf5"

[ -e $chrest_data ] && return 0
echo proceeding! with chrest_data: $chrest_data

# processes 1 ablate file
python ablateData.py --file "$1" --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_Yi:Yi monitor_zMix:zmix
python wrangle_chrest_dataset_for_UQ.py "$chrest_data"
