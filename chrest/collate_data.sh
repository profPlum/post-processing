#!/bin/bash

data_dir='/Volumes/UB_CHREST/v2-560x80x80-768/flowField_mixtureFraction/'
collated_fn="$data_dir/chrest_collated.csv"

collate() {
    tail -n +2 "$1" >> $collated_fn # skip the csv header!
}

# put in csv header! NOTE: 00000 is fake data! All 0's so Matt asked me to skip it & I will
head -n 1 "${data_dir}/flowField_mixtureFraction.00001.chrest"/chrest_data.csv > $collated_fn
map collate "${data_dir}/flowField_mixtureFraction."*'.chrest'/chrest_data.csv
