#!/bin/bash

# parse options from command line
while getopts 'i:b:o:n:m:t:s:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		b) bwadir=$OPTARG ;;
		o) outdir=$OPTARG ;;
		n) name=$OPTARG ;;
		m) min_contig=$OPTARG ;;
		t) threads=$OPTARG ;;
		s) seed=$OPTARG ;;
        esac
done

mkdir -p $outdir

jgi_summarize_bam_contig_depths \
	--outputDepth $outdir/$name'_contigs_depth.txt' \
	$bwadir/$name'_alignments'/*.bam \
	&> $outdir/$name'_metabat.log'
metabat2 -i $indir/$name'_contigs.fasta' \
	-a $outdir/$name'_contigs_depth.txt' \
	-o $outdir/$name'_bin' \
	-m $min_contig \
	-t $threads \
	--seed $seed \
	-v \
	&>> $outdir/$name'_metabat.log'
