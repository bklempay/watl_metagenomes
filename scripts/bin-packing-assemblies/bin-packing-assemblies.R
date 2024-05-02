#!/usr/bin/env Rscript

#### Bin packing assemblies ####

# Run script with path to your sample list: must have a `name` column with
# sample names and a `nreads` column with the number of reads in each sample
sample.list <- commandArgs(trailingOnly = TRUE)
# or to use the WATL samples, just run this script without any arguments
if (length(sample.list) == 0) sample.list <- "watl-assemblies-stats"

# Our lab has two remote servers that we use for our bioinformatics analyses:
# fram and gjoa. Both of these servers have massive amounts of RAM, 500 Gb and
# 1 Tb, respectively. But even that's not enough memory to run 125 metagenomic
# assemblies in parallel. So I developed this algorithm to group these 125
# samples into batches of samples that are small enough to be assembled in
# parallel, and the batches will then be run sequentially using all available
# CPUs. This strategy dramatically reduced assembly time from months to weeks.

fram.RAM <- 500
gjoa.RAM <- 1000


#### Estimate peak RAM requirement for WATL assemblies based on past runs ####

# Peak RAM usage for metaSPAdes assemblies scales up with the number of reads.
# We can model this relationship in order to predict the peak RAM requirements
# for our samples as a function of N reads. As far as I know, it's difficult to
# predict how much RAM metaSPAdes will require from one system to the next, so
# it might be necessary to run a test on your own system first.

# read in past assembly stats... replace this with your own data (you can use
# spades_peakRAM.sh to extract peak RAM used by SPAdes from spades.log files)
sbsw.dat <- read.table("sbsw-assemblies-stats", header = TRUE)
sbsw.dat$site <- "sbsw"
watl.dat <- read.table("watl-assemblies-stats", header = TRUE)
watl.dat$site <- "watl"
dat <- rbind(sbsw.dat, watl.dat)

# remove outliers (these were bad WMGS libraries which took up more RAM than
# their small number of reads would normally require)
dat <- dat[dat$peak.RAM.Gb < dat$nreads * 6e-6 + 5, ]

# predict peak RAM as a function of N reads (log-log scale linear model)
linreg <- summary(lm(log10(peak.RAM.Gb) ~ log10(nreads), data = dat))
a <- linreg$coefficients["log10(nreads)", 1]
b <- linreg$coefficients["(Intercept)", 1]
e <- linreg$sigma

# a pretty plot to admire just how excellent our model is:
library(ggplot2)
ggplot(dat, aes(x = nreads, y = peak.RAM.Gb, group = 1)) +
  geom_point(aes(col = site)) +
  scale_x_log10(limits = c(8e4, 3e8)) +
  scale_y_log10(limits = c(0.8, 1e3)) +
  geom_abline(slope = a, intercept = b) +
  geom_abline(slope = a, intercept = b + 2 * e, lty = "dotted") +
  geom_abline(slope = a, intercept = b - 2 * e, lty = "dotted") +
  geom_hline(yintercept = fram.RAM, col = "red", lty = "dashed") +
  geom_hline(yintercept = gjoa.RAM, col = "blue", lty = "dashed") +
  annotate("text", 1e5, fram.RAM, 
           label = "fram", 
           hjust = 1, 
           vjust = 1, 
           col = "red") +
  annotate("text", 1e5, gjoa.RAM,
           label = "gjoa", 
           hjust = 1, 
           vjust = 1, 
           col = "blue") +
  theme_bw()


# predict peak RAM for samples as a function of N reads
samples <- read.table(sample.list, header = TRUE)
samples$peak.RAM.Gb <- 10^(a * log10(samples$nreads) + b + 2 * e)


#### Assign assemblies to batches to be executed on multiple servers ####
# jobs in each batch will be run in parallel and batches will be run sequentially
# such that total RAM requirements never exceed available memory on either server
# this is a 1-dimensional bin-packing problem with multiple bin sizes

# variant on first-fit decreasing algorithm
pack.bins <- function(items, bin.sizes, labels = NULL) {
  if (max(items) > max(bin.sizes)) {
    stop("The largest item is too big for the largest bin!")
  }
  if (!is.null(labels) & length(items) != length(labels)) {
    stop("The number of labels must be equal to the number of items.")
  }
  
  sorted.labels <- labels[order(items, decreasing = TRUE)]
  sorted.items <- sort(items, decreasing = TRUE)
  bin.sizes <- sort(bin.sizes)
  
  # initiate 1 empty bin of each size
  bins <- labs <- emptybins <- vector("list", length = length(bin.sizes))
  remaining_capacity <- bin.sizes
  
  # for each item in the list, sorted in decreasing order
  for (i in 1:length(items)) {
    item <- sorted.items[i]
    label <- sorted.labels[i]
    added <- FALSE
    # iterate through the existing bins
    for (bin in 1:length(bins)) {
      # if item is less than or equal to bin's remaining capacity
      if (item <= remaining_capacity[bin]) {
        # add item to that bin
        bins[[bin]] <- c(bins[[bin]], item)
        if (!is.null(labels)) labs[[bin]] <- c(labs[[bin]], label)
        remaining_capacity[bin] <- remaining_capacity[bin] - item
        added <- TRUE
      }
      if (added) break
    }
    # if item can't be added to any existing bin
    if (!added) {
      # create new set of empty bins
      bins <- c(bins, emptybins)
      if (!is.null(labels)) labs <- c(labs, emptybins)
      remaining_capacity <- c(remaining_capacity, bin.sizes)
      # and add item to smallest possible bin
      bin <- bin + min(which(bin.sizes > item))
      bins[[bin]] <- c(bins[[bin]], item)
      if (!is.null(labels)) labs[[bin]] <- c(labs[[bin]], label)
      remaining_capacity[bin] <- remaining_capacity[bin] - item
    }
    # now sort bins by least remaining capacity, rinse, and repeat
    bins <- bins[order(remaining_capacity)]
    if (!is.null(labels)) labs <- labs[order(remaining_capacity)]
    remaining_capacity <- sort(remaining_capacity)
  }
  # toss empty bins, sort bins by size, and return list of packed bins
  bins <- bins[!sapply(bins, is.null)]
  labs <- labs[!sapply(labs, is.null)]
  labs <- labs[order(sapply(bins, sum), decreasing = TRUE)]
  bins <- bins[order(sapply(bins, sum), decreasing = TRUE)]
  if (!is.null(labels)) bins <- list(item.labels = labs, item.sizes = bins)
  return(bins)
}


#### Use pack.bins function to assign assemblies to batches ####

bins <- pack.bins(samples$peak.RAM.Gb, c(fram.RAM, gjoa.RAM),
                  labels = samples$name)

# create worklists for gjoa and fram

nbatches <- ceiling(length(bins$item.labels) / 2)

for (i in 1:nbatches) {
  write(paste(bins$item.labels[[i]], collapse = " "),
        file = "gjoa.worklist.names",
        append = TRUE)
  write(paste(bins$item.sizes[[i]], collapse = " "),
        file = "gjoa.worklist.mem",
        append = TRUE)
  write(paste(bins$item.labels[[i+nbatches]], collapse = " "),
        file = "fram.worklist.names",
        append = TRUE)
  write(paste(bins$item.sizes[[i+nbatches]], collapse = " "),
        file = "fram.worklist.mem",
        append = TRUE)
}