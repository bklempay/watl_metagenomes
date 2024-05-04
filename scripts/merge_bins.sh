#!/bin/bash

# parse options from command line
while getopts 'i:o:s:t:' opt; do
        case $opt in
		i) indir=$OPTARG ;;
		o) outfile=$OPTARG ;;
		s) suffix=$OPTARG ;;
		t) type=$OPTARG ;;
        esac
done

mkdir -p $(dirname $outfile)
rm -f $outfile

# get column headers (ignored if type = fasta)
header=$( awk '!/^##/ {print; exit}' $(ls $indir/*$suffix | head -n 1))
# assumes that every file has the same header and number of columns
# ignores commented lines (beginning with ##)

if [ $type = "fasta" ]; then
	for bin in $indir/*$suffix; do
		name=$(basename -s "$suffix" $bin)
		sed 's/>/>'$name'_/g' $bin >> $outfile
	done
elif [ $type = "csv" ]; then
	echo "bin,$header" > $outfile
	for bin in $indir/*$suffix; do
		name=$(basename -s "$suffix" $bin)
		grep -v "^##" $bin | awk -F, 'NR > 1 {print "'$name'," $0}' >> $outfile
	done
elif [ $type = "tsv" ]; then
	echo -e "bin\t$header" > $outfile
	for bin in $indir/*$suffix; do
		name=$(basename -s "$suffix" $bin)
		grep -v "^##" $bin | awk 'NR > 1 {print "'$name'\t" $0}' >> $outfile
		# assumes that the first line is a header
		# ignores commented lines (beginning with ##)
	done
else
	echo $type type not recognized
	exit 1
fi
