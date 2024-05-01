#!/bin/bash

# parse options from command line
while getopts 'i:o:t:l:c:C:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		o) outdir=$OPTARG ;;
                t) threads=$OPTARG ;;
                l) min_length=$OPTARG ;;
                c) min_completeness=$OPTARG ;;
                C) max_contamination=$OPTARG ;;
        esac
done

mkdir -p $outdir

dRep dereplicate $outdir \
	-p $threads \
	-g $indir/*.fa \
	-l $min_length \
	-comp $min_completeness \
	-con $max_contamination \
	&> $outdir/drep.log
