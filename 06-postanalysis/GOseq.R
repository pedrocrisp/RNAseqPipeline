library(goseq)
library(multicore)

ARGV <- commandArgs(trailingOnly=TRUE)

if (length(ARGV) < 1) {
  stop("USAGE: GOseq.R --args <test.table> [<goslim mapping> [<gff>]]")
  q(status=-1)
}

out.base <- "go/"
dir.create(out.base, showWarnings=F)

test.path <- ARGV[1]

tair <- "ftp://ftp.arabidopsis.org/home/tair/"
if (length(ARGV) == 2) {
  go.path <- ARGV[2]
} else {
  go.path <- paste0(tair, "Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt")
}

if (length(ARGV) == 3) {
  gff.path <- ARGV[3]
} else {
  gff.path <- paste0(
    tair,
    "Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"
    )
}

# grab go table from tair
go <- read.delim(go.path, header=F)
names(go) <- c(
  "locusName", "TAIRAccession", "AGI", "GOrelationship", "GOterm",
  "GOID", "TAIRKeywordID", "Aspect", "GOslimTerm", "EvidenceCode",
  "EvidenceDesc", "EvidenceWith", "Reference", "Annotator", "Date"
  )
go.goseq <- go[,c("AGI", "GOID")]

# gff preparation
gff <- read.delim(gff.path, header=F)
gff <- gff[gff[,3] == "gene",]

agis <- sapply(
  as.character(gff[,9]),
  function(s)
    rev(unlist(strsplit(s, "=")))[1]
  )
names(agis) <- NULL
gene.lengths <- gff[,5] - gff[,4] + 1
names(gene.lengths) <- agis


test <- read.csv(test.path)
rownames(test) <- as.character(test[,1])
test <- as.matrix(test[,2:length(test)])

summary(test)

all.genes <- rownames(test)
gene.lengths <- gene.lengths[match(all.genes, names(gene.lengths))]

goTestIndex <- function(tn) {
  tst <- test[,tn]
  test.name <- colnames(test)[tn]
  de.genes <- as.integer(tst < 0.05)
  names(de.genes) <- all.genes
  pwf <- nullp(de.genes, "tair10", id=all.genes, bias.data=gene.lengths)
  go.all <- goseq(pwf, "tair10", gene2cat=go.goseq)
  this.go <- go[match(go.all$category, as.character(go$GOID)),c(4:7)]
  go.all$GOrelationship <- this.go$GOrelationship
  go.all$GOterm <- this.go$GOterm
  go.all$sig <- go.all$over_represented_pvalue < 0.05
  row.names(go.all) <- NULL

  out.name <- paste0(out.base, test.name, ".goall.csv")
  print(out.name)
  write.csv(go.all, out.name, row.names=F)
  go.all
}

GA <- mclapply(1:ncol(test), goTestIndex)

terms <- sapply(
  GA,
  function(x) {
    b <- x[,"GOterm"]
    b[!x[,"sig"]] <- ""
    b
    },
  simplify="array"
  )
colnames(terms) <- colnames(test)
write.csv(terms, paste0(out.base, "terms.csv"), row.names=F)
