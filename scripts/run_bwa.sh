#!/bin/bash

# parse options from command line
while getopts 'i:r:o:n:m:t:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		r) refdir=$OPTARG ;;
		o) outdir=$OPTARG ;;
		n) refname=$OPTARG ;;
		m) mapname=$OPTARG ;;
		t) threads=$OPTARG ;;
        esac
done

mkdir -p $outdir/$refname'_alignments'

if [ ! -f $outdir/$refname'_index.amb' ]; then
	bwa index -p $outdir/$refname'_index' \
		$refdir/$refname*.fasta \
		&> $outdir/$refname'_index.log'
fi

bwa mem -t $threads \
	$outdir/$refname'_index' \
	$indir/$mapname'_R1_qc.fastq.gz' \
	$indir/$mapname'_R2_qc.fastq.gz' \
	2> $outdir/$refname'_alignments'/$mapname'_bwa.log' | \
samtools sort \
	-o $outdir/$refname'_alignments'/$mapname'_map.bam' \
	-O BAM \
	-@ $threads \
	&>> $outdir/$refname'_alignments'/$mapname'_bwa.log'
