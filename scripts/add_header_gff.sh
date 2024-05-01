#!/bin/bash

# first check if header attribute already exists
if grep -Pq "header=\w+;" $1; then
	echo "Error: header attribute already exists"
	exit 1
fi

# append header attribute to gff files created by prodigal
while read line; do
	contig=$(echo $line | grep -Po '^[\w.]+')
	number=$(echo $line | grep -Po 'ID=\d+_\K\d+')
	sed 's/ID=/header='$contig'_'$number';ID=/' <<< $line >> $1.tmp
done < $1

wait; mv $1.tmp $1
