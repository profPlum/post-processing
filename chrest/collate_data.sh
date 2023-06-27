#!/bin/bash

data_dir='/user/dwyerdei/data/ablate_data_processing/raw_ablate_data/flowField_mixtureFraction'
collated_fn="$data_dir/chrest_collated.csv"

for i in $data_dir/*.chrest/*.gz; do zcat "$i" | sed '1d'; done > $collated_fn
return 0

#collate() {
#    tail -n +2 "$1" >> $collated_fn # skip the csv header!
#}
#
## put in csv header! NOTE: 00000 is fake data! All 0's so Matt asked me to skip it & I will
#head -n 1 "${data_dir}/flowField_mixtureFraction.00001.chrest"/chrest_data-rebalanced.csv > $collated_fn
#map collate "${data_dir}/flowField_mixtureFraction."*'.chrest'/chrest_data-rebalanced.csv
