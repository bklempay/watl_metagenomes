#!/bin/bash

# parse options from command line
while getopts 'i:o:n:f:r:m:q:u:b:t:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		o) outdir=$OPTARG ;;
		n) name=$OPTARG ;;
                f) forward_pattern=$OPTARG ;;
                r) reverse_pattern=$OPTARG ;;
		m) trim_quality=$OPTARG ;;
		q) filt_quality=$OPTARG ;;
		u) unqual_limit=$OPTARG ;;
		b) n_base_limit=$OPTARG ;;
		t) threads=$OPTARG ;;
        esac
done

mkdir -p $outdir

fastp	-i $indir/$name*$forward_pattern* \
	-o $outdir/$name'_R1_qc.fastq.gz' \
	-I $indir/$name*$reverse_pattern* \
	-O $outdir/$name'_R2_qc.fastq.gz' \
	--unpaired1 $outdir/$name'_orphaned.fastq.gz' \
	--unpaired2 $outdir/$name'_orphaned.fastq.gz' \
	--cut_front --cut_tail -M $trim_quality \
	--qualified_quality_phred $filt_quality \
	--unqualified_percent_limit $unqual_limit \
	--n_base_limit $n_base_limit \
	--json $outdir/$name'_qc.json' \
	--html $outdir/$name'_qc.html' \
	--thread $threads \
	&> $outdir/$name'_fastp.log'
