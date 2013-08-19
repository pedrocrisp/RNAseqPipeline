library(edgeR)
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


