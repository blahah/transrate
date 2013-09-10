#!/usr/bin/env Rscript

# script to load output from compareassemblies.rb
# and create a merged table of gene-level coverage for all assemblies
library(plyr)

# load sample metadata
# filename, label, usearchresults, alignment
samples <- read.csv('samples.csv', head=TRUE, as.is=TRUE)

# load reference data
data <- read.csv('metricstable.csv', as.is=TRUE)

# merge in sample coverage data from RSEM files
mergeSample <- function(row) {   
    # extract filename and label from row
    filename <- row[3]
    label <- row[2]
    print(paste('merging in data for ', label))
    # read in the usearch results
    s <- read.csv(filename, head=FALSE, sep='\t', as.is=TRUE)
    # subset to only reference transcript ID + search ID + target coverage
    s <- s[,2:4]
    # take the maximum coverage for each transcript ID
    s <- ddply(s, .(V2), numcolwise(max))
    # subset to only transcript ID + target coverage
    s <- s[,c(1,3)]
    # label the columns
    names(s) <- c('Transcript.ID', label)
    # merge
    return(s)
}

for(i in 1:nrow(samples)) {
    row <- as.character(samples[i,])
    print(row)
    s <- mergeSample(row)
    data <- merge(data, s, by='Transcript.ID', all.x=TRUE)
}

# write the output
write.csv(data, 'allsamplesdata.csv', row.names=FALSE, na='0')
