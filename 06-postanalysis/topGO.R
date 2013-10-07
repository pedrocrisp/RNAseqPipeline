## UNCOMMENT TO INSTALL
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("topgo"))
# install.packages(c('reshape2', 'gplots', 'multicore'))

library(topGO, quietly=T)
library(multicore)
library(GO.db)
library(reshape2)

## USAGE:
# edgeR_glm.R --args <keyfile> <paramfile>

ARGV <- commandArgs(trailingOnly=TRUE)
ARGV <- c("kevin_hons/kevin-hons.key", '~/prog/bio/honsPipeline/pipelines/de_glm/edgeR_kmhons.R')
if (length(ARGV) < 2) {
  stop("USAGE: edgeR_glm.R --args <keyfile> <paramfile>")
  q(status=-1)
}

#######################################################




p.matrix
names(p.matrix) <- c(test.names)

summary(p.matrix)
hist(p.matrix)

sum(fdr.tab)
t <- 1
tg <- new(
  "topGOdata",
  description=test.names[t],
  allGenes=p.matrix[,t],
  geneSelectionFun=function (p) p<0.05,
  ontology=)