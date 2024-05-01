#!/usr/bin/env Rscript

# Required input: [1] dRep data_tables directory,
# [2] gtdb_to_ncbi_majority_vote.py output file, and
# [3] jgi_summarize_contigs_depth 
args <- commandArgs(trailingOnly = TRUE)


#### Summarize bin assembly statistics and classification ####
bins <- read.csv(paste0(args[1], "/genomeInfo.csv"), row.names = 1)
clusts <- read.csv(paste0(args[1], "/Cdb.csv"), row.names = 1)
winners <- read.csv(paste0(args[1], "/Wdb.csv"), row.names = 2)
taxa <- read.delim(args[2], row.names = 1)

# Summarize CheckM and dRep results
bins$drep_cluster <- clusts[row.names(bins), "secondary_cluster"]
bins$drep_winner <- winners[bins$drep_cluster, "genome"]
rownames(bins) <- gsub("\\.(fasta|fas|fa|fna)$", "", rownames(bins),
                       ignore.case = TRUE)
bins$drep_winner <- gsub("\\.(fasta|fas|fa|fna)$", "", bins$drep_winner,
                       ignore.case = TRUE)

# Summarize GTDB-Tk results
bins$gtdb_classification <- taxa[rownames(bins), "GTDB.classification"]
bins$ncbi_classification <- taxa[rownames(bins), "Majority.vote.NCBI.classification"]

# Output summary file
write.csv(bins, "bins_summary.csv", quote = FALSE)


#### Calculate MAG relative abundance in each sample ####
depth <- read.table(args[3], header = TRUE, check.names = FALSE)
depth$binName <- sub("_NODE.*$", "", depth$contigName)
sumabd <- depth[, grepl("_map$", names(depth))] * depth$contigLen
sumabd <- aggregate(. ~ depth$binName, data = sumabd, FUN = sum)
sumlen <- aggregate(depth$contigLen ~ depth$binName, FUN = sum)
avgabd <- sumabd[, -1] / sumlen [, 2]
relabd <- t(avgabd) / rowSums(t(avgabd))
rownames(relabd) <- sub("_map$", "", rownames(relabd))
colnames(relabd) <- sumabd$`depth$binName`

# Output relative abundance matrix
write.csv(relabd, "bins_relabund.csv", quote = FALSE)