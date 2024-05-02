#!/bin/bash

# Use this script to extract peak RAM used by SPAdes from any spades.log files contained
# in subdirectories (SPAdes output directories) of the directory provided to this script.
# Example usage: ./spades_peakRAM.sh ~/metagenomes_watl/all_output/metaspades_output

spadesdir=$1

for i in $spadesdir/*; do
	name=$(basename -s _assembly $i)
	max=$(grep -E '/ [0-9]+G' $i/spades.log | awk '{print $4}' | sed 's/G//' | sort -n | tail -n 1)
	if [[ "$max" == "" ]]; then
		max=$(($(grep -E '/ [0-9]+M' $i/spades.log | awk '{print $4}' | sed 's/M//' | sort -n | tail -n 1) / 1000))
	fi
	echo $name $max
done
