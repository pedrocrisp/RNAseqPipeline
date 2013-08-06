source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR")) ## UNCOMMENT TO INSTALL

## USAGE:
# edgeR.R --args <commasep_countfiles> <commasep_samplenames> <gfffile>

library(edgeR)

keyfile.path <- ARGV[1]
#keyfile.path <- "~/hons-pipeline/km_pipeline/pete-exp277.keyfile"
keyfile <- read.delim(keyfile.path)

samples <- keyfile$Sample
countFiles <- paste("count/", samples, "/", samples, ".counts", sep="")

countDFs <- lapply(countFiles, function (fle) read.delim(fle))

# get num genes, == number of rows in first DF
numGenes <- nrow(countDFs[[1]])

geneNames <- as.character(countDFs[[1]]$geneid)
geneLengths <- countDFs[[1]]$length

countList <- lapply(countDFs, function (df) df$nreads)
countMatrix <- matrix(unlist(countList), nrow=numGenes)
sampleGroups <- as.character(keyfile$Treatment)
dge <- DGEList(counts=countMatrix, group=sampleGroups)

# TODO: DESeq style +1 transform to allow statisical comparison of low abundance -> high abundance transcripts.

dgen <- calcNormFactors(dge, method="TMM")
dgec <- estimateCommonDisp(dgen)
dget <- estimateTagwiseDisp(dgec)

# tad hacky, I'll have to come up with a better way of doing this.
groups <- as.character(unique(keyfile[,3]))
nonControlGroups <- groups[groups!="Control"]
testPairs <- lapply(nonControlGroups, function (x) c("Control", x))
testNames <- lapply(testPairs, function (x) paste(x, collapse=".vs."))
tests <- lapply(testPairs, function (x) exactTest(dget, pair=x))
names(tests) <- testNames

# incomplete: Need to write results out in a coherent table, like the output of cuffdiff. Also need to do plots.
