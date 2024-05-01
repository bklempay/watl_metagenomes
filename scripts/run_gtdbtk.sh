#!/bin/bash

# Specify locations of gtdb_to_ncbi_majority_vote.py and GTDB metadata directory
gtdb2ncbi_script=$(which gtdb_to_ncbi_majority_vote.py)
gtdb_meta_dir=$(echo $gtdb2ncbi_script | sed 's/GTDB-Tk\/.*/GTDB-Tk\/metadata/')
# I did my best to make these portable, but it might require some editting by hand
# e.g. if gtdb_to_ncbi_majority_vote.py is not in PATH

if [ ! -f $gtdb2ncbi_script ]; then
	echo "Error: gtdb_to_ncbi_majority_vote.py not found"
	echo "Error: Edit run_gtdbtk.sh before continuing"
	exit 1
fi

if [ ! -e $gtdb_meta_dir ]; then
	echo "Error: GTDB metadata directory not found"
	echo "Error: Edit run_gtdbtk.sh before continuing"
	exit 1
fi

# parse options from command line
while getopts 'i:o:t:' opt; do
	case $opt in
		i) indir=$OPTARG ;;
		o) outdir=$OPTARG ;;
		t) threads=$OPTARG ;;
	esac
done

mkdir -p $outdir

gtdbtk classify_wf \
	--genome_dir $indir \
	--out_dir $outdir \
	--mash_db $outdir/mash_db.msh \
	-x fa \
	--cpus $threads \
	--pplacer_cpus $threads \
	--force \
	&> /dev/null

$gtdb2ncbi_script \
	--gtdbtk_output_dir $outdir \
	--output_file $outdir/gtdb_to_ncbi.tsv \
	--ar53_metadata_file $gtdb_meta_dir/ar53_metadata.tsv.gz \
	--bac120_metadata_file $gtdb_meta_dir/bac120_metadata.tsv.gz \
	&>> $outdir/gtdbtk.log
