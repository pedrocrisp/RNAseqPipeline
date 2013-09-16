## UNCOMMENT TO INSTALL
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR"))
# install.packages(c('reshape2))

library(edgeR)
library(multicore)
library(reshape2)
library(gplots)

## USAGE:
# edgeR-pairwise.R --args <keyfile>

ARGV <- commandArgs(trailingOnly=TRUE)

if (length(ARGV) < 1) {
  stop("USAGE: edgeR.R --args <keyfile>")
  q(status=-1)
}

allpairs <- FALSE


################################################################################
######################           Keyfile            ############################
################################################################################

# Load keyfile etc:
keyfile.path <- ARGV[1]
keyfile <- read.delim(keyfile.path)
samples <- keyfile$Sample
analysis.name <- unlist(
  strsplit(rev(unlist(strsplit(keyfile.path, "/")))[1], "\\."))[1]

out.base <- paste0("./de/", analysis.name, "/")
dir.create(out.base, recursive=T)


################################################################################
######################           Data Entry         ############################
################################################################################

count.files <- paste0("count/", samples, "/", samples, ".counts")

sample.groups <- as.character(paste(keyfile[,3]))

tmp <- read.delim(count.files[[1]])
gene.lengths <- as.numeric(tmp[,2])
names(gene.lengths) <- tmp[,1]
rm(tmp)

dge <- readDGE(
  count.files,
  columns=c(1,3),
  group=sample.groups,
  labels=as.character(samples)
)

gene.names <- as.character(rownames(dge$counts))


################################################################################
######################      DGE Normalisation       ############################
################################################################################

# TODO: DESeq style variance stabilisng transform to allow statisical comparison
# of low abundance -> high abundance transcripts.

n.samples <- length(sample.groups)
n.reps <- n.samples / length(groups)
min.avg.reads.per.sample <- 0.5 # cpm, invariant
min.total.reads <- min.reads.per.sample * n.samples
min.samples.with.reads <- ((n.samples / n.reps) * 0.5) + n.reps

mr.keep <- rowSums(cpm(dge)) > min.total.reads
table(mr.keep)

ms.keep <- rowSums(cpm(dge) > 0) > min.samples.with.reads
table(ms.keep)

loci.2.keep <- ms.keep & mr.keep
table(loci.2.keep)
n.tags <- sum(loci.2.keep) # sum of true values, i.e. count all genes to keep

old.dge <- dge
dge <- old.dge[loci.2.keep,]
dge$samples$lib.size <- colSums(dge$counts)

old.gene.names <- gene.names
gene.names <- gene.names[loci.2.keep]

dge <- calcNormFactors(dge, method="TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)


################################################################################
######################       Diff Exp Testing       ############################
################################################################################

groups <- unique(sample.groups)

if(allpairs) {
  test.pairs <- combn(groups, 2) # bug here
} else {
  non.control.groups <- groups[2:length(groups)]
  test.pairs <- lapply(non.control.groups, function (x) c(groups[1], x))
}

test.names <- lapply(test.pairs, function (x) paste(x, collapse=".vs."))
test.names
tests <- mclapply(test.pairs, function (x) exactTest(dge, pair=x))
names(tests) <- test.names
names(tests)

################################################################################
######################         Write Results        ############################
################################################################################

short.names <- paste(
  substr(sample.groups, 0, 4),
  rep_len(1:n.reps, length(sample.groups)),
  sep="R"
)
nice.sample.names <- paste(
  sample.groups,
  rep_len(1:n.reps, length(sample.groups)),
  sep=".R"
)

# global plots
pdf(paste0(out.base, analysis.name, "_bcv.pdf"))
plotBCV(dge, main=analysis.name)
dev.off()

pdf(paste0(out.base, analysis.name, "_mds.pdf"))
plotMDS(dge,
  main=analysis.name,
  labels=nice.sample.names,
  col=rep(rainbow(length(groups)), each=n.reps)
)
dev.off()


xf <-  4 # fold change at which lines are ruled in test writer below
n.hm.genes <- 500
# a nice red -> black -> blue gradient.
# Visualise w/
# image(matrix(1:255), col=redblue)
redblue <- paste0(
  "#",
  as.hexmode(c(rep(0, 128),1:127 * 2)),
  format(as.hexmode(rep(0, 255)), width=2),
  as.hexmode(c(127:1 * 2, rep(0, 128)))
)
# write exact test tables out
for (tst in tests) {
  test.name <- paste(tst$comparison, collapse=".VS.")
  test.base.dir <- paste0(out.base, test.name, "/")
  dir.create(test.base.dir)

  decision <- decideTestsDGE(tst, p=0.05)
  print(test.name)
  print(table(decision))
  detags <- gene.names[as.logical(decision)]

  # tables for this test
  tt <- topTags(tst, n=n.tags)
  write.csv(tt, paste0(test.base.dir, test.name, "_toptags.csv"))


  hm.cols <- sample.groups %in% tst$comparison
  hm.rows <- match(rownames(tt)[1:n.hm.genes], rownames(dge$counts))
  table <- dge$counts[hm.rows, hm.cols]
  write.csv(tt, paste0(test.base.dir, test.name, "_sampletable.csv"))

  # plots
  pdf(paste0(test.base.dir, test.name, "_smear.pdf"))
  plotSmear(
      dge,
      de.tags=detags,
      main=test.name,
      sub=paste("lines indicate", xf, "fold change")
      )
  abline(h=c(-log2(xf), log2(xf)), col="blue")
  dev.off()

  pdf(paste0(test.base.dir, test.name, "_heatmap.pdf"))
  heatmap.2(
    log(table+1),
    col=redblue,
    trace="none",
    scale='row',
    main=test.name,
    labRow=NA,
    density.info="histogram",
    cexCol=0.8,
    labCol=short.names,
    xlab="Samples",
    ylab=paste("Tags", "(n=", n.hm.genes, ")"),
    lmat=rbind(c(0,3), c(2,1), c(4,4)), # WHAT IS THIS BLACK MAGIC
    lwid=c(0.3,4), # http://stackoverflow.com/questions/15351575/moving-color-
    lhei=c(2,6,3) # url cont: # key-in-r-heatmap-2-function-of-gplots-package
  )
  dev.off()
}

pfc.matrix <- predFC(dge)
fc.matrix <- sapply(tests, function (t) t$table$logFC, simplify = "array")
p.matrix <- sapply(tests, function (t) t$table$PValue, simplify = "array")
fdr.matrix <- sapply(tests, function (t) p.adjust(t$table$PValue), simplify="array")
rownames(fc.matrix) <- gene.names
rownames(p.matrix) <- gene.names
rownames(fdr.matrix) <- gene.names

cpm.matrix <- cpm(dge, log=T)
colnames(cpm.matrix) <- nice.sample.names

write.csv(cpm.matrix, file=paste0(out.base, analysis.name, "_cpm.csv"))
write.csv(fc.matrix, file=paste0(out.base, analysis.name, "_fc.csv"))
write.csv(pfc.matrix, file=paste0(out.base, analysis.name, "_pfc.csv"))
write.csv(fdr.matrix, file=paste0(out.base, analysis.name, "_fdr.csv"))
write.csv(p.matrix, file=paste0(out.base, analysis.name, "_p.csv"))

factor.names <- colnames(keyfile)[3:ncol(keyfile)]

cpm.melt <- melt(cpm.matrix, varnames=c("geneID", factor.names), value.name="CPM")
fdr.melt <- melt(fdr.matrix, varnames=c("geneID", factor.names), value.name="FDR")
fc.melt <- melt(fc.matrix, varnames=c("geneID", factor.names), value.name="logFC")
p.melt <- melt(p.matrix, varnames=c("geneID", factor.names), value.name="P")
pfc.melt <- melt(p.matrix, varnames=c("geneID", factor.names), value.name="P")

write.csv(cpm.melt, file=paste0(out.base, analysis.name, "_cpm_melt.csv"), row.names=F)
write.csv(fdr.melt, file=paste0(out.base, analysis.name, "_fdr_melt.csv"), row.names=F)
write.csv(fc.melt, file=paste0(out.base, analysis.name, "_fc_melt.csv"), row.names=F)
write.csv(pfc.melt, file=paste0(out.base, analysis.name, "_pfc_melt.csv"), row.names=F)
write.csv(p.melt, file=paste0(out.base, analysis.name, "_p_melt.csv"), row.names=F)
