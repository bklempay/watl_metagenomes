#!/bin/bash

# parse options from command line
while getopts 'i:o:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		o) outdir=$OPTARG ;;
        esac
done

mkdir -p $outdir
rm -f $outdir/merged_bins.fasta

for bin in $indir/*.fa; do
	name=$(basename -s .fa $bin)
	sed s/'>'/'>'$name'_'/g $bin >> $outdir/merged_bins.fasta
done
