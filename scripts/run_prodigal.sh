#!/bin/bash

# parse options from command line
while getopts 'i:o:n:m:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		o) outdir=$OPTARG ;;
		n) name=$OPTARG ;;
		m) mode=$OPTARG ;;
        esac
done

mkdir -p $outdir

prodigal -p $mode \
	-i $indir/$name'.fa' \
	-o $outdir/$name'.gff' \
	-a $outdir/$name'_proteins.faa' \
	-d $outdir/$name'_CDS.fna' \
	-f gff \
	&> $outdir/$name'_prodigal.log'
