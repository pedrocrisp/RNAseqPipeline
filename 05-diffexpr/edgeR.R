source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR")) ## UNCOMMENT TO INSTALL

## USAGE:
# edgeR.R --args <commasep_countfiles> <commasep_samplenames> <gfffile>

library(edgeR)

keyfile.path <- ARGV[1]
keyfile <- read.delim(keyfile.path)

samples <- keyfile$Sample
countFiles <- paste("count/", samples, "/", samples, ".counts", sep="")

geneNames <- as.character(countDFs[[1]]$geneid)
sampleGroups <- as.character(keyfile$Treatment)
#geneLengths <- countDFs[[1]]$length

dge <- readDGE(
	       countFiles,
	       columns=c(1,3),
	       groups=sampleGroups,
	       labels=as.character(sampleGroups)
	       )

# TODO: DESeq style +1 transform to allow statisical comparison of low abundance -> high abundance transcripts.

dge <- calcNormFactors(dge, method="TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

groups <- as.character(unique(keyfile[,3]))
nonControlGroups <- groups[groups!="Control"]
testPairs <- lapply(nonControlGroups, function (x) c("Control", x))
testNames <- lapply(testPairs, function (x) paste(x, collapse=".vs."))
tests <- lapply(testPairs, function (x) exactTest(dge, pair=x))
names(tests) <- testNames

for (tst in tests) {
	fn <- paste(tst$comparison, collapse=".VS.")
	row.names(table) <- geneNames
	# next line doesn't work as it should. Row names are not changed after sort.
	#table <- table[with(table, order(PValue)),]
	write.csv(table, paste("de/", fn , ".csv", sep=""))
}
