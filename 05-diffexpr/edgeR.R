source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR")) ## UNCOMMENT TO INSTALL

## USAGE:
# edgeR.R --args <keyfile>
library(edgeR)

################################################################################
######################           Keyfile            ############################
################################################################################

keyfile.path <- ARGV[1]
keyfile <- read.delim(keyfile.path)
samples <- keyfile$Sample
analysisName <- unlist(strsplit(rev(unlist(strsplit(keyfile.path, "/")))[1], "\\."))[1]


################################################################################
######################           Data Entry         ############################
################################################################################

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


################################################################################
######################      DGE Normalisation       ############################
################################################################################

# TODO: DESeq style +1 transform to allow statisical comparison of low abundance -> high abundance transcripts.

dge <- calcNormFactors(dge, method="TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)


################################################################################
######################       Diff Exp Testing       ############################
################################################################################

groups <- as.character(unique(keyfile[,3]))
nonControlGroups <- groups[groups!="Control"]
testPairs <- lapply(nonControlGroups, function (x) c("Control", x))
testNames <- lapply(testPairs, function (x) paste(x, collapse=".vs."))
tests <- lapply(testPairs, function (x) exactTest(dge, pair=x))
names(tests) <- testNames


################################################################################
######################         Write Results        ############################
################################################################################

# global plots
pdf(paste0(analysisName, ".pdf"))
plotBCV(dge, main=analysisName)
plotMDS(dge, main=analysisName)
dev.off()

# write exact test tables out
xf <-  4 # fold change at which lines are ruled
for (tst in tests) {
	testName <- paste(tst$comparison, collapse=".VS.")
	baseName <-paste0("de/", testName)

	decision <- decideTestsDGE(et, p=0.05)
	detags <- geneNames[as.logical(decision)]
	table <- tests$table
	row.names(table) <- geneNames
	# next line doesn't work as it should. Row names are not changed after sort.
	#table <- table[with(table, order(PValue)),]
	write.csv(table, paste0(baseName, ".csv"))

	# plots
	pdf(paste(baseName, "_smear.pdf"))
	plotSmear(
		  dge,
		  de.tags=detags,
		  main=testName,
		  sub=paste("lines indicate", xf, "fold change")
		  )
	abline(h=c(-log2(xf), log2(xf)), col=blue)
	dev.off()
}
