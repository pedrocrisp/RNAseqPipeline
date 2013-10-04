library(multicore)
ARGV <- commandArgs(trailingOnly=TRUE)

if (length(ARGV) < 1) {
  stop("USAGE: getGOTerms.R --args <outbase>")
  q(status=-1)
}

outname <- ARGV[1]

go <- read.delim("ftp://ftp.arabidopsis.org/home/tair/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt", header=F)

names(go) <- c("locusName", "TAIRAccession", "AGI", "GOrelationship", "GOterm",
               "GOID", "TAIRKeywordID", "Aspect", "GOslimTerm", "EvidenceCode",
               "EvidenceDesc", "EvidenceWith", "Reference", "Annotator", "Date")

go.involvedin <- go[which(go$GOrelationship == "involved in"),]

annotLevel <- matrix(unlist(strsplit(as.character(go$TAIRAccession), ":")), ncol=2, byrow=T)[,1]

genes <- annotLevel == "gene"
go.out <- go.involvedin[genes,c(3:6, 9)]

agis <- unique(as.character(go.out$AGI))

go.map <-mclapply(
  agis,
  function (agi)
    as.character(go[which(go$AGI == agi), "GOID"])
  )

names(go.map) <- agis

write.csv(
	go.out,
  paste0(outname, ".table.csv"),
	row.names=F
	)

dump(
  "go.map",
  file=paste0(outname, ".gomap")
  )
