## UNCOMMENT TO INSTALL
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("edgeR"))
# install.packages(c('reshape2))

library(edgeR)
library(multicore)
library(reshape2)
library(gplots)

## USAGE:
# edgeR_glm.R --args <keyfile> <paramfile>

# <paramfile> should contain a model matrix, a list of contrasts, a cutoff, etc:
# modelmatrix <-  c(...)
# contrasts <- c(...)
# FDR.cutoff <- 0.05 # or other

ARGV <- commandArgs(trailingOnly=TRUE)

if (length(ARGV) < 2) {
  stop("USAGE: edgeR_glm.R --args <keyfile> <paramfile>")
  q(status=-1)
}


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
dir.create(out.base, recursive=T, showWarnings=F)

sample.groups <- apply(
  as.matrix(keyfile[,3:length(keyfile)]),
  1,
  paste,
  collapse=" "
)

# Model settings etc. Must happen after `keyfile` has been defined
source(ARGV[2], echo=T)


################################################################################
######################           Data Entry         ############################
################################################################################

count.files <- paste0("count/", samples, "/", samples, ".counts")

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

### find and remove rRNAs from analysis. ###
rRNA <-  c(
  "AT2G01010",
  "AT2G01020",
  "AT3G41768",
  "AT3G41979",
  "ATCG00920",
  "ATCG00950",
  "ATCG00960",
  "ATCG00970",
  "ATCG01160",
  "ATCG01170",
  "ATCG01180",
  "ATCG01210",
  "ATMG00020",
  "ATMG01380",
  "ATMG01390"
  )

# find rRNAs
rRNA.tags <- match(rRNA, gene.names)

# count rRNAs
rRNAs <- dge$counts[rRNA.tags,]
rRNA.summary <- colSums(rRNAs)
rRNA.rates <- rRNA.summary / dge$samples$lib.size
print(rRNA.rates)

# remove rRNA from dge matrix
dge$counts <- dge$counts[-rRNA.tags,]

### Remove Sparse Tags ###
groups <- unique(sample.groups)
n.samples <- length(sample.groups)
n.reps <- n.samples / length(groups)

fancy.filter <- FALSE

if (fancy.filter) {
  # A smart-arsed filter which is supposed to tolerate different
  # experimental designs, but really just confuses the user.
  # Use with caution.

  min.avg.reads.per.sample <- 0.5 # cpm, invariant
  min.total.reads <- min.avg.reads.per.sample * n.samples
  min.samples.with.reads <- ((n.samples / n.reps) * 0.5) + n.reps

  mr.keep <- rowSums(cpm(dge)) > min.total.reads
  table(mr.keep)

  ms.keep <- rowSums(cpm(dge) > 0) > min.samples.with.reads
  table(ms.keep)

  loci.2.keep <- ms.keep & mr.keep
} else {
  # A simple filter:
  # Only loci with more than 10 counts per million in at least `n.reps` samples
  # are kept
  loci.2.keep <- rowSums(dge$counts > 10) > n.reps
}

table(loci.2.keep)
n.tags <- sum(loci.2.keep) # sum of true values, i.e. count all genes to keep
gene.names.keep <- as.character(rownames(dge$counts))[loci.2.keep]
length(gene.names.keep)

old.dge <- dge
dge <- old.dge[loci.2.keep,]
dge$samples$lib.size <- colSums(dge$counts)

old.gene.names <- gene.names
gene.names <- gene.names[loci.2.keep]

## edgeR normalisation and dispersion calculation ###
dge <- calcNormFactors(dge, method="TMM")
dge <- estimateGLMCommonDisp(dge, modelmatrix, verbose=TRUE)
dge <- estimateGLMTrendedDisp(dge, modelmatrix)
dge <- estimateGLMTagwiseDisp(dge, modelmatrix)


################################################################################
######################       Diff Exp Testing       ############################
################################################################################

groups <- unique(sample.groups)

test.names <- colnames(contrasts)
n.tests <- ncol(contrasts)
test.names

glm <- glmFit(dge, design=modelmatrix)

# run all tests
tests <- mclapply(
  1:n.tests, 
  function (n) 
    glmLRT(glm, contrast=contrasts[,n])
  )


################################################################################
######################         Write Results        ############################
################################################################################

# make some constants #
n.reps <- length(sample.groups) / length(groups)
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

agi.fd <- read.delim("ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_functional_descriptions")
agi.fd$AGI <- substr(agi.fd$Model_name, 1, 9)


### global plots ###
pdf(paste0(out.base, analysis.name, "_bcv.pdf"))
plotBCV(dge, main=analysis.name)
dev.off()

pdf(paste0(out.base, analysis.name, "_mds.pdf"))
plotMDS(dge,
  main=analysis.name,
  labels=nice.sample.names,
  col=rep(rainbow(length(groups)), each=n.reps),
  top=n.tags
)
dev.off()


###  write and plot contrast  ###
xf <-  4 # fold change at which lines are ruled
n.hm.genes <- 500 # number of genes to heatmap (from toptags)
n.hm.genes <- n.tags
for (tn in 1:n.tests) {
  test.name <- test.names[tn]
  tst <- tests[[tn]]

  test.base.dir <- paste0(out.base, test.name, "/")
  dir.create(test.base.dir, showWarnings=F)

  decision <- decideTestsDGE(tst, p=0.05)
  print(test.name)
  print(table(decision))
  detags <- rownames(tst$table)[as.logical(decision)]

  # tables for this test
  tt <- topTags(tst, n=n.tags)
  agi.match <- match(detags, agi.fd$AGI)
  tt$short_description <- agi.fd$Short_description[agi.match]
  tt$curator_summary <- agi.fd$Curator_summary[agi.match]
  tt$computational_description <- agi.fd$Computational_description[agi.match]
  write.csv(tt, paste0(test.base.dir, test.name, "_toptags.csv"))

  hm.cols <- sample.groups %in% names(which(contrasts[,tn] !=0))
  hm.rows <- match(rownames(tt)[1:n.hm.genes], rownames(dge$counts))
  table <- dge$counts[, hm.cols]
  write.csv(table, paste0(test.base.dir, test.name, "_sampletable.csv"))

  # plots
  pdf(paste0(test.base.dir, test.name, "_smear.pdf"))
  plotSmear(
    tst,
    de.tags=detags,
    main=test.name,
    sub=paste("lines indicate", xf, "fold change")
  )
  abline(h=c(-log2(xf), log2(xf)), col="blue")
  dev.off()

  # FDR diagonstic plots
  test.name <- test.names[tn]
  pdf(paste0(test.base.dir, test.name, ".fdr_dist.pdf"))
  plot(
    tt$table$FDR,
    type="l",
    ylab="FDR",
    main="Distribution of FDR values",
    ylim=c(0,1)
  )
  abline(h=FDR.cutoff, col="red")
  abline(v=sum(decision != 0), col="blue")
  dev.off()

  pdf(paste0(test.base.dir, test.name, ".logfdr_dist.pdf"))
  plot(
    log(tt$table$FDR),
    type="l",
    ylab="log(FDR)",
    main="Distribution of log(FDR) values"
  )
  abline(h=log(FDR.cutoff), col="red")
  abline(v=sum(decision != 0), col="blue")
  dev.off()

  # density plots
  pdf(paste0(test.base.dir, test.name, ".fdr_dens.pdf"))
  plot(
    density(tt$table$FDR),
    type="l",
    ylab="Density",
    xlab="FDR",
    main="Density of FDR values",
  )
  abline(h=FDR.cutoff, col="red")
  abline(v=sum(decision != 0), col="blue")
  dev.off()

  pdf(paste0(test.base.dir, test.name, ".logfdr_dens.pdf"))
  plot(
    density(log(tt$table$FDR)),
    type="l",
    ylab="Density",
    xlab="log(FDR)",
    main="Density of log(FDR) values"
  )
  abline(h=log(FDR.cutoff), col="red")
  abline(v=sum(decision != 0), col="blue")
  dev.off()
}


# test-wise matricies
fc.matrix <- sapply(tests, function (t) t$table$logFC, simplify="array")
p.matrix <- sapply(tests, function (t) t$table$PValue, simplify="array")
fdr.matrix <- sapply(tests, function (t) p.adjust(t$table$PValue), simplify="array")
rownames(fc.matrix) <- gene.names.keep
rownames(p.matrix) <- gene.names.keep
rownames(fdr.matrix) <- gene.names.keep
colnames(fc.matrix) <- test.names
colnames(p.matrix) <- test.names
colnames(fdr.matrix) <- test.names

write.csv(fc.matrix, file=paste0(out.base, analysis.name, "_fc.csv"))
write.csv(fdr.matrix, file=paste0(out.base, analysis.name, "_fdr.csv"))
write.csv(p.matrix, file=paste0(out.base, analysis.name, "_p.csv"))


# conglomerate the above test stat matricies into a big table (3D array)
all.array <- array(c(fc.matrix, p.matrix, fdr.matrix), c(3, dim(fc.matrix)))
dimnames(all.array) <- list(
  c("fc", "p", "fdr"),
  rownames(fc.matrix),
  colnames(fc.matrix)
  )
all.df<-as.data.frame.table(all.array)
names(all.df) <- c("stat", "geneid", "test", "value")
rm(all.array)
write.csv(all.df, file=paste0(out.base, analysis.name, "_allstat.csv"))


# sample-wise matrices
pfc.matrix <- predFC(dge)
cpm.matrix <- cpm(dge, log=T)
colnames(cpm.matrix) <- nice.sample.names

write.csv(cpm.matrix, file=paste0(out.base, analysis.name, "_cpm.csv"))
write.csv(pfc.matrix, file=paste0(out.base, analysis.name, "_pfc.csv"))
