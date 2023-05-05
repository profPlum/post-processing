#!/bin/bash

ablate2chrest() {
	python ablateData.py --file "$1" --fields monitor_densityEnergySource:souener monitor_densityYiSource:souspec monitor_yi:Yi monitor_zMix:zmix
}

#ablate2chrest '/Volumes/UB CHREST/v2-560x80x80-768/flowField_mixtureFraction/flowField_mixtureFraction.00001.hdf5'

amap ablate2chrest '/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'*.hdf5 
