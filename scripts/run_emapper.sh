#!/bin/bash

# parse options from command line
while getopts 'i:o:n:t:m:e:b:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		o) outdir=$OPTARG ;;
		n) name=$OPTARG ;;
		t) threads=$OPTARG ;;
		m) mode=$OPTARG ;;
		e) evalue=$OPTARG ;;
		b) bitscore=$OPTARG ;;
        esac
done

mkdir -p $outdir

emapper.py --cpu $threads \
	-i $indir/$name'_proteins.faa' \
	--itype proteins \
	-m $mode \
	--evalue $evalue \
	--score $bitscore \
	--dbmem \
	--seed_ortholog_evalue $evalue \
	--seed_ortholog_score $bitscore \
	-o $name \
	--output_dir $outdir \
	--decorate_gff $indir/$name'.gff' \
	--decorate_gff_ID_field header \
	&> $outdir/$name'_emapper.log'
