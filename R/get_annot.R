summ.f <- read.table("output/BC_NC_summary.txt", header=T)

# Create columns for uniprot accession numbers
summ.f$sp.id <- gsub("\\|", "", stringr::str_extract(summ.f$sp.hit, "\\|.+\\|"))
summ.f$tr.id <- gsub("\\|", "", stringr::str_extract(summ.f$tr.hit, "\\|.+\\|"))

# Function to get uniprot data for each row of data frame with a uniprot 
#   accession given in a column called "uid"
uniprotInfo <- function(x, col="id") {
  if (class(x)!="data.frame") x <- setNames(data.frame(x), col)
  if (is.na(x[,col])) {
    return(NULL)
  } else {
    raw <- RCurl::getURL(paste0("http://www.uniprot.org/uniprot/?query=accession:", x[,col],
                                "&columns=id,entry%20name,reviewed,protein%20names,go,organism,length&format=tab"))
    if (!raw=="") {
      parsed <- read.delim(textConnection(raw), sep="\t", stringsAsFactors=FALSE)
      return(parsed)
    } else {
      return(NULL)
    }
  }
}

# Get GenBank info from nt hits
annota <- function(acc) {
  tryCatch({
    choosebank("genbank")
    annot <- getAnnot(query("raw", acc))[[1]]
    start <- grep("DEFINITION", annot)
    end <- grep("ACCESSION", annot)
    def <- gsub("  ", "", paste(annot[c(seq(start,end-1,1))], collapse=""))
    def <- gsub("DEFINITION", "", def)
    return(def)
  }, error = function(e) {
    tryCatch({
      choosebank("refseq")
      annot <- getAnnot(query("raw", acc))[[1]]
      start <- grep("DEFINITION", annot)
      end <- grep("ACCESSION", annot)
      def <- gsub("  ", "", paste(annot[c(seq(start,end-1,1))], collapse=""))
      def <- gsub("DEFINITION", "", def)
      return(def)
    }, error = function(e) {
      return(NA)
    })
  })
}



# Apply the uniprotInfo function to each row of summ.f (each uniprot accession will be queried)
sprot.out  <- plyr::adply(summ.f, 1, .fun=function(x) uniprotInfo(x, "sp.id"))
sprot.out <- sprot.out[, c("contig", "Protein.names", "Gene.ontology..GO.", "Organism")]
colnames(sprot.out) <- c("contig", "sp.Protein.names", "sp.GO", "sp.Organism")
trembl.out <- plyr::adply(summ.f, 1, .fun=function(x) uniprotInfo(x, "tr.id"))
trembl.out <- trembl.out[, c("contig", "Protein.names", "Gene.ontology..GO.", "Organism")]
colnames(trembl.out) <- c("contig", "tr.Protein.names", "tr.GO", "tr.Organism")

all <- merge(summ.f, sprot.out, all=T)
all <- merge(all, trembl.out, all=T)

# Get GenBank annotations for contigs with nt hits
choosebank("genbank")
acc <- paste0("AC=", na.omit(gsub("\\..+$", "", all$nt.hit)))
gbannot <- sapply(acc, annota)
all$nt.annot <- rep(NA, nrow(all))
all$nt.annot[which(!is.na(all$nt.hit))] <- gbannot

# Subset columns...
all <- all[,c("contig", "nsnp", "meancov", "nt.hit", "nt.evalue", "nt.annot",
              "co.hit", "co.evalue", "sp.hit", "sp.evalue", "sp.Protein.names",
              "sp.GO", "sp.Organism", "tr.hit", "tr.evalue", 
              "tr.Protein.names", "tr.GO", "tr.Organism")]
all <- all[with(all, order(tr.evalue, -meancov)), ]


save(all, file="output/BC_NC_summary_plus.RData")
write.table(all, file="output/BC_NC_summary_plus.txt", row.names=F, quote=F)
