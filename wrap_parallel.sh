#!/bin/bash

#### user defined parameters ####
cpus=60
indir=raw_reads
forward_pattern=_R1
reverse_pattern=_R2

# fastp parameters
trim_quality=30
filt_quality=25
unqual_limit=16
n_base_limit=0

# metaSPAdes parameters
kmers=21,33,55

# MetaBAT2 parameters
min_contig=2500
seed=4444

# dRep parameters
min_length=50000
min_completeness=50
max_contamination=10

# eggNOG-mapper parameters
search_mode=mmseqs
evalue=0.001
bitscore=60
# Note: default search algorithm parameters (e.g. sensitivity) will be kept for
# the sake of consistency and code compatibility between diamond and mmseqs modes


# parse options from command line
while getopts 'n:b:m:' opt; do
	case $opt in
		# REQUIRED sample names must be separated by spaces, and batches by linebreaks
		n) worklist_names=$OPTARG ;;
		# REQUIRED reference sample name followed by list of samples for read mapping
		b) bwa_maplist=$OPTARG ;;
		# OPTIONAL estimated RAM requirement for each metaspades assembly (Gb)
		m) worklist_mem=$OPTARG ;;
	esac
done
names=$(sed -z 's/\n/ /g' $worklist_names) # remove linebreaks from worklist


echo -e '\n\n###############################'
echo '#### Trim and QC raw reads ####'
echo '###############################'
date
fastp --version
echo -e '\nfastp paramters:'
echo '--cut_mean_quality 		'$trim_quality
echo '--qualified_quality_phred	'$filt_quality
echo '--unqualified_percent_limit	'$unqual_limit
echo '--n_base_limit			'$n_base_limit
echo '--thread			5'
# 5 threads seems to be pretty optimal for fastp. I could decrease runtime by sorting the worklist
# by size in order to start the biggest jobs first, but it won't save that much time to be worth it

# run fastp multiple jobs in parallel
parallel -j $(($cpus/5)) scripts/run_fastp.sh \
	-i $indir -o all_output/fastp_output -n {1} -f $forward_pattern -r $reverse_pattern \
	-m $trim_quality -q $filt_quality -u $unqual_limit -b $n_base_limit -t 5 \
	::: $names

echo -e '\nDone!'
echo 'See all_output/fastp_output/*fastp.log for details'


echo -e '\n\n##############################'
echo '#### Assemble metagenomes ####'
echo '##############################'
date
metaspades.py --version
echo -e '\nmetaspades.py parameters:'
echo '--threads			'$cpus' / n assemblies in batch (see '$worklist_names')'
echo '-k				'$(echo $kmers | sed 's/,/ /g')

# run metaSPAdes in sequential batches of parallel jobs
for i in $(seq $(wc -l < $worklist_names)); do
	batch=$(awk 'NR=='$i $worklist_names)
	memlim=$(printf '250 %.0s' $(seq $(wc -l < $worklist_names))) # default memory limit 250 Gb
	# memlim=$(awk 'NR=='$i $worklist_mem) # uncomment to allocate memory limit to each assembly
	threads=$(($cpus/$(echo $batch | wc -w))) # allocate all available threads to each batch

	parallel -j $cpus --link scripts/run_metaspades.sh \
		-i all_output/fastp_output -o all_output/metaspades_output -n {1} \
		-t $threads -m {2} -k $kmers \
		::: $batch ::: $memlim
done

echo -e '\nDone!'
echo 'See all_output/metaspades_output/*_assembly/spades.log for details'

# make symbolic links to contigs in ./assemblies/
mkdir -p assemblies
for name in $names; do
	if [ -f all_output/metaspades_output/$name'_assembly/contigs.fasta' ]; then
		ln -fs ../all_output/metaspades_output/$name'_assembly/contigs.fasta' \
			assemblies/$name'_contigs.fasta'
	else
		echo -e '\n'$name' assembly failed'
		echo -e 'See all_output/metaspades_output/'$name'_assembly/spades.log for details'
	fi
done


echo -e '\n\n########################################'
echo '#### Align QCed reads to assemblies ####'
echo '########################################'
date
bwa |& head -n 3
samtools |& head -n 3
echo -e '\nbwa mem parameters:'
echo '-t				'$cpus

# index multiple assemblies in parallel (not required but saves a lot of time)
parallel -j $cpus bwa index \
	-p all_output/bwa_output/{1}_index assemblies/{1}_contigs.fasta \
	::: $names &> /dev/null

# run alignments in sequence using all available threads
while read line; do
	refname=$(echo $line | awk '{print $1}')
	maplist=$(echo $line | awk '{$1=""; print $0}')
	for name in $maplist; do
		scripts/run_bwa.sh \
			-i all_output/fastp_output -r assemblies -o all_output/bwa_output \
			-n $refname -m $name -t $cpus
	done
done < $bwa_maplist

echo -e '\nDone!'
echo 'See all_output/bwa_output/*_alignments/*bwa.log for details'


echo -e '\n\n#####################'
echo '#### Bin contigs ####'
echo '#####################'
date
metabat2 |& head -n 2
echo -e '\nmetabat2 parameters:'
echo '--minContig			'$min_contig
echo '--numThreads			10'
echo '--seed				'$seed

# run MetaBAT2 multiple jobs in parallel
parallel -j $(($cpus/10)) scripts/run_metabat.sh \
	-i assemblies -b all_output/bwa_output -o all_output/metabat_output \
	-n {1} -m $min_contig -t 10 -s $seed \
	::: $names
# jgi_summarize_bam_contig_depths can only use a single thread per bam file, so we can
# parallelize more efficiently by keeping metabat2 to the same number of threads as samples

echo -e '\nDone!'
echo 'See all_output/metabat_output/*metabat.log for details'

# make symbolic links to bins in ./bins/
mkdir -p bins
for name in $names; do
	if [ -f all_output/metabat_output/$name'_bin.1.fa' ]; then
		for bin in all_output/metabat_output/$name'_bin'*; do
			ln -fs ../$bin bins
		done
	else
		echo -e '\n'$name' formed zero bins'
		echo -e 'See all_output/metabat_output/'$name'_metabat.log for details'
	fi
done


echo -e '\n\n#######################'
echo '#### Classify bins ####'
echo '#######################'
date
gtdbtk |& head -n 2
echo -e '\ngtdbtk classify_wf parameters:'
echo '--cpus				'$cpus
echo '--pplacer_cpus			'$cpus

scripts/run_gtdbtk.sh -i bins -o all_output/gtdbtk_output -t $cpus

echo -e '\nDone!'
echo 'See all_output/gtdbtk_output/gtdbtk.log for details'


echo -e '\n\n#################################'
echo '#### QC and dereplicate bins ####'
echo '#################################'
date
checkm |& head -n 2
dRep |& head -n 2
echo -e '\ndRep dereplicate parameters:'
echo '-p				'$cpus
echo '-l				'$min_length
echo '-comp				'$min_completeness
echo '-con				'$max_contamination

scripts/run_drep.sh -i bins -o all_output/drep_output -t $cpus \
	-l $min_length -c $min_completeness -C $max_contamination

echo -e '\nDone!'
echo 'See all_output/drep_output/drep.log for details'


echo -e '\n\n###############################################'
echo '#### Align QCed reads to dereplicated bins ####'
echo '###############################################'
date
bwa |& head -n 3
samtools |& head -n 3
echo -e '\nbwa mem parameters:'
echo '-t				'$cpus

# merge dereplicated bins into a single reference fasta
scripts/merge_bins.sh \
	-i all_output/drep_output/dereplicated_genomes \
	-o all_output/relabund_calc/merged_bins.fasta \
	-s .fa -t fasta

# run alignments in sequence using all available threads
for name in $names; do
	scripts/run_bwa.sh \
		-i all_output/fastp_output -r all_output/relabund_calc -o all_output/relabund_calc \
		-n merged_bins -m $name -t $cpus
done

echo -e '\nDone!'
echo 'See all_output/relabund_calc/merged_bins_alignments/*bwa.log for details'


echo -e '\n\n######################################'
echo '#### Predict protein-coding genes ####'
echo '######################################'
date
prodigal -v |& head -n 2
echo -e '\nprodigal parameters:'
echo '-p				single'

# get list of dereplicated bins
MAGs=all_output/drep_output/dereplicated_genomes/*.fa
names_MAGs=$(basename -a -s .fa $MAGs)

# run Prodigal multiple jobs in parallel
parallel -j $cpus scripts/run_prodigal.sh \
	-i all_output/drep_output/dereplicated_genomes -o all_output/prodigal_output \
	-n {1} -m single \
	::: $names_MAGs

echo -e '\nDone!'
echo 'See all_output/prodigal_output/*prodigal.log for details'

# eggNOG mapper cannot decorate gff files as created by prodigal
# append a header attribute corresponding to FASTA file headers
parallel -j $cpus scripts/add_header_gff.sh {1} \
	::: all_output/prodigal_output/*.gff


echo -e '\n\n###########################'
echo '#### Annotate proteins ####'
echo '###########################'
date
emapper.py --version
echo -e '\nemapper.py parameters:'
echo '--cpu				'$cpus
echo '-m				'$search_mode
echo '--evalue			'$evalue
echo '--score				'$bitscore
echo '--seed_ortholog_evalue		'$evalue
echo '--seed_ortholog_score		'$bitscore

# run eggNOG mapper in sequence using all available threads
for name in $names_MAGs; do
	scripts/run_emapper.sh \
		-i all_output/prodigal_output -o all_output/emapper_output \
		-n $name -t $cpus -m $search_mode -e $evalue -b $bitscore
done

echo -e '\nDone!'
echo 'See all_output/emapper_output/*emapper.log for details'

# merge emapper annotations into a single tsv
scripts/merge_bins.sh \
	-i all_output/emapper_output -o annotations_merged.tsv \
	-s .emapper.annotations -t tsv


#### Summarize bins and calculate relative abundance ####

# summarize bam contig depths
jgi_summarize_bam_contig_depths \
	--outputDepth all_output/relabund_calc/merged_bins_contigs_depth.txt \
	all_output/relabund_calc/merged_bins_alignments/*.bam \
	&> all_output/relabund_calc/merged_bins_jgi_summarize.log

# summarize bin assembly statistics and classification
# and calculate MAG relative abundance in each sample
scripts/summarize_bins.R \
	all_output/drep_output/data_tables \
	all_output/gtdbtk_output/gtdb_to_ncbi.tsv \
	all_output/relabund_calc/merged_bins_contigs_depth.txt


#### Calculate protein stats ####
mkdir -p all_output/protein_stats
parallel -j $cpus scripts/protein_stats.py {1} \
	::: all_output/prodigal_output/*.faa
mv all_output/prodigal_output/*.csv all_output/protein_stats
# merge protein stats into a single csv
scripts/merge_bins.sh \
	-i all_output/protein_stats -o bins_protein_stats.csv \
	-s proteins_stats.csv -t csv
