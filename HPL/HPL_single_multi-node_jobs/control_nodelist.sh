#!/bin/bash

node_list=()
while IFS= read -r line; do
    node_list+=("$line")
done < nodelist_test.txt

# Print array elements
for element in "${node_list[@]}"; do
        ./sample_script.sh $element

done
