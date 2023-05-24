#!/bin/bash

data_dir='/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'

################### Wrangle Ablate: ################### 

ablate2chrest() {
	python ablateData.py --file "$1" --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_yi:Yi monitor_zMix:zmix
}
#ablate2chrest '/Volumes/UB CHREST/v2-560x80x80-768/flowField_mixtureFraction/flowField_mixtureFraction.00001.hdf5'

amap ablate2chrest "$data_dir"*.hdf5
wait

################### Wrangle Chrest: ###################

chrest2csv() {
	python wrangle_chrest_dataset_for_UQ.py "$1"
}
#chrest2csv '/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/flowField_mixtureFraction.00000.chrest/flowField_mixtureFraction.00000.chrest.00000.hdf5'

amap chrest2csv "${data_dir}/flowField_mixtureFraction."*'.chrest/flowField_mixtureFraction.'*'.chrest.'*'.hdf5'
wait
#######################################################
