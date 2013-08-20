## UNCOMMENT TO INSTALL
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR"))

library(edgeR)
library(multicore)
library(reshape2)

## USAGE:
# edgeR.R --args <keyfile>
ARGV <- commandArgs(trailingOnly=TRUE)
if (length(ARGV) < 1) {stop("USAGE: edgeR.R --args <keyfile>"); q(status=-1)}

allpairs <- FALSE

################################################################################
######################           Keyfile            ############################
################################################################################

keyfile.path <- ARGV[1]
keyfile <- read.delim(keyfile.path)
samples <- keyfile$Sample
analysis.name <- unlist(strsplit(rev(unlist(strsplit(keyfile.path, "/")))[1], "\\."))[1]
rm(keyfile.path)

################################################################################
######################           Data Entry         ############################
################################################################################

count.files <- paste("count/", samples, "/", samples, ".counts", sep="")

sample.groups <- as.character(paste(keyfile[,3:ncol(keyfile)]))
#geneLengths <- countDFs[[1]]$length

dge <- readDGE(
  count.files,
  columns=c(1,3),
  group=sample.groups,
  labels=as.character(samples)
)
gene.names <- as.character(rownames(dge$counts))

rm(count.files)

################################################################################
######################      DGE Normalisation       ############################
################################################################################

# TODO: DESeq style variance stabilisng transform to allow statisical comparison
# of low abundance -> high abundance transcripts.

min.reads <- 1
min.samples.with.min.reads <- 1
num.samples.with.enough.reads <- rowSums(cpm(dge)>min.reads)
loci.2.keep <- num.samples.with.enough.reads > min.samples.with.min.reads

old.dge <- dge
dge <- old.dge[loci.2.keep,]
dge$samples$lib.size <- colSums(dge$counts)

gene.names.keep <- gene.names[loci.2.keep]

dge <- calcNormFactors(dge, method="TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

rm(min.reads)
rm(min.samples.with.min.reads)
rm(num.samples.with.enough.reads)

################################################################################
######################       Diff Exp Testing       ############################
################################################################################


groups <- unique(sample.groups)
if(allpairs) {
  testPairs <- combn(groups, 2)
} else {
  non.control.groups <- groups[2:length(groups)]
  test.pairs <- lapply(non.control.groups, function (x) c("Control", x))
}
test.names <- lapply(test.pairs, function (x) paste(x, collapse=".vs."))
tests <- mclapply(test.pairs, function (x) exactTest(dge, pair=x))
names(tests) <- test.names
names(tests)

################################################################################
######################         Write Results        ############################
################################################################################

short.sample.names <- sapply(
  strsplit(as.character(samples),  "_"),
  function (l) paste(l[4:5], collapse="_")
  )

nice.sample.names <- paste(sample.groups, rep(1:3,length(groups)), sep=".")

# global plots
pdf(paste0(analysisName, ".pdf"))
plotBCV(dge, main=analysisName)
plotMDS(dge,
        main=analysisName,
        labels=short.sample.names,
        col=rainbow(length(groups))
        )
dev.off()

# write exact test tables out
xf <-  4 # fold change at which lines are ruled
for (tst in tests) {
	test.name <- paste(tst$comparison, collapse=".VS.")
	base.name <-paste0("de/", test.name)

	decision <- decideTestsDGE(tst, p=0.05)
	detags <- gene.names.keep[as.logical(decision)]
	table <- tst$table
	table$geneid <- gene.names.keep
	row.names(table) <- NULL
	# next line doesn't work as it should. Row names are not changed after sort.           
	table <- table[with(table, order(PValue)),]
	write.csv(table, paste0(base.name, ".csv"))

	# plots
	pdf(paste0(base.name, "_smear.pdf"))
	plotSmear(
		  dge,
		  de.tags=detags,
		  main=test.name,
		  sub=paste("lines indicate", xf, "fold change")
		  )
	abline(h=c(-log2(xf), log2(xf)), col="blue")
	dev.off()
}

fc.matrix <- sapply(tests, function (t) t$table$logFC, simplify = "array")
p.matrix <- sapply(tests, function (t) t$table$logFC, simplify = "array")
fdr.matrix <- sapply(tests, function (t) p.adjust(t$table$PValue), simplify="array")
rownames(fc.matrix) <- gene.names.keep
rownames(p.matrix) <- gene.names.keep
rownames(fdr.matrix) <- gene.names.keep

cpm.matrix <- cpm(dge, log=T)
colnames(cpm.matrix) <- nice.sample.names

write.csv(cpm.matrix, file=paste0(analysis.name, "_cpm.csv"))
write.csv(fc.matrix, file=paste0(analysis.name, "_fc.csv"))
write.csv(fdr.matrix, file=paste0(analysis.name, "_fdr.csv"))
write.csv(p.matrix, file=paste0(analysis.name, "_p.csv"))

factor.names <- colnames(keyfile)[3:length(colnames(keyfile))]

cpm.melt <- melt(cpm.matrix, varnames=c("geneID", factor.names), value.name="CPM")
fdr.melt <- melt(fdr.matrix, varnames=c("geneID", factor.names), value.name="FDR")
fc.melt <- melt(fc.matrix, varnames=c("geneID", factor.names), value.name="logFC")
p.melt <- melt(p.matrix, varnames=c("geneID", factor.names), value.name="P")

write.csv(cpm.melt, file=paste0(analysis.name, "_cpm_melt.csv"), row.names=F)
write.csv(fdr.melt, file=paste0(analysis.name, "_fdr_melt.csv"), row.names=F)
write.csv(fc.melt, file=paste0(analysis.name, "_fc_melt.csv"), row.names=F)
write.csv(p.melt, file=paste0(analysis.name, "_p_melt.csv"), row.names=F)