#!/bin/bash

output_file="nodelist_quana.txt"

i_start=1
i_end=10

j_start=1
j_end=60

echo "i_start=$i_start, i_end=$i_end"
echo "j_start=$j_start, j_end=$j_end"

> "output_file"

for i in $(seq $i_start $i_end); do
	for j in $(seq $j_start $j_end); do 
		echo "cpu-$i-$j" >> "output_file"
	done
done

echo "File '$output_filename' contains list of all possible nodes in Quanah"
