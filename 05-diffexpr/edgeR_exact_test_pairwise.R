source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR")) ## UNCOMMENT TO INSTALL

## USAGE:
# edgeR.R --args <keyfile>
library(edgeR)
ARGV <- commandArgs(trailingOnly=TRUE)
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

sampleGroups <- as.character(keyfile$Treatment)
#geneLengths <- countDFs[[1]]$length

dge <- readDGE(
	       countFiles,
	       columns=c(1,3),
	       group=sampleGroups,
	       labels=as.character(samples)
	       )
geneNames <- as.character(rownames(dge$counts))


################################################################################
######################      DGE Normalisation       ############################
################################################################################

# TODO: DESeq style +1 transform to allow statisical comparison of low abundance -> high abundance transcripts.
min.reads <- 1
min.samples.with.min.reads <- 1
keep <- rowSums(cpm(dge)>min.reads)>min.samples.with.min.reads

dge <- dge[keep,]
dge$samples$lib.size <- colSums(dge$counts)

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

	decision <- decideTestsDGE(tst, p=0.05)
	detags <- geneNames[as.logical(decision)]
	table <- tst$table
	table$geneid <- geneNames
	row.names(table) <- NULL
	# next line doesn't work as it should. Row names are not changed after sort.
	table <- table[with(table, order(PValue)),]
	write.csv(table, paste0(baseName, ".csv"))

	# plots
	pdf(paste0(baseName, "_smear.pdf"))
	plotSmear(
		  dge,
		  de.tags=detags,
		  main=testName,
		  sub=paste("lines indicate", xf, "fold change")
		  )
	abline(h=c(-log2(xf), log2(xf)), col="blue")
	dev.off()
}
