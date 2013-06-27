source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR")) ## UNCOMMENT TO INSTALL

## USAGE:
# edgeR.R --args <commasep_samfiles> <commasep_samplenames> <gfffile>


# get arguments
samfiles <- strsplit(ARGV[1], ",")
samplenames <- strsplit(ARGV[2], ",")
gfffile <- ARGV[3]

######## get counts from a bam ############
# shamelessly pilfered from https://stat.ethz.ch/pipermail/bioc-devel/2012-September/003608.html
# biocLite(c("Rsubread")) ## UNCOMMENT TO INSTALL
library(Rsubread)

gfftab <- read.table(gfffile, sep="\t", stringsAsFactors=FALSE)
txtab <- gfftab[gfftab[,3]=="mRNA",] # change this to exon if you want isoforms

AGIs <- gsub(".+(AT.G\\d{5}\\.\\d).+","\\1",txtab[,9],perl=TRUE)  #get the agi's. dont ask how it works.

annotation <- data.frame(entrezid=as.integer(factor(AGIs)), 
                   chromosome=txtab$V1, chr_start=txtab$V4, chr_stop=txtab$V5)

# This sorts AGIs by chr position
ID <- sort(unique(annotation$entrezid))
AGIs <- AGIs[match(ID, annotation$entrezid)]


c_rsub <- featureCounts(SAMfiles=samfiles, annot=annotation, type = "gene")
counts <- c_rsub$counts

rownames(counts) <-AGIs
colnames(counts) <-basename(samfiles)
hist(counts[,1],breaks=1000)

library(edgeR)
edgeRUsersGuide()