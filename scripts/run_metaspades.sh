#!/bin/bash

# parse options from command line
while getopts 'i:o:n:t:m:k:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		o) outdir=$OPTARG ;;
                n) name=$OPTARG ;;
		t) threads=$OPTARG ;;
		m) memory=$OPTARG ;;
		k) kmers=$(echo $OPTARG | sed 's/,/ /g') ;;
        esac
done

mkdir -p $outdir

metaspades.py -o $outdir/$name'_assembly'\
	-1 $indir/$name'_R1_qc.fastq.gz' \
	-2 $indir/$name'_R2_qc.fastq.gz' \
	-s $indir/$name'_orphaned.fastq.gz' \
	--threads $threads \
	--memory $memory \
	-k $kmers \
	1> /dev/null
