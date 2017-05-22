require(Biostrings)

# Import SNP data
snps <- read.delim("data/SNPs/SNPS_info", header=F, stringsAsFactors=F)
colnames(snps) <- c("contig", "pos", "NC", "ND", "BC", "FST")

# Convert nucleotide frequencies for each group to a list of numeric vectors
makeNumList <- function(x) lapply(strsplit(x, split=","), as.numeric)
snps <- within(snps, {
  BC <- makeNumList(BC)
  NC <- makeNumList(NC)
  ND <- makeNumList(ND)
})

# Extract SNPs with fixed differences between two groups with coverage threshold
getFixedSNPs <- function(groups, mincov) {
  numZeros <- function(x) sapply(x, function(n) sum(n==0))
  whichN <- function(x) sapply(x, function(n) which(n!=0))
  minCov <- function(x) sapply(x, function(n) any(n>=mincov))
  #fixed <- which(numZeros(snps[,groups[1]])==3 & numZeros(snps[,groups[2]])==3)
  fixed <- snps[numZeros(snps[,groups[1]])==3 & numZeros(snps[,groups[2]])==3, ]
  diff <- fixed[whichN(fixed[,groups[1]])!=whichN(fixed[,groups[2]]), ]
  cov <- diff[minCov(diff[,groups[1]]) & minCov(diff[,groups[2]]), ]
  return(cov)
}

BCNC <- getFixedSNPs(groups=c("NC", "BC"), mincov=5)

# Filter to only include SNPs from contigs with > 92% column conservation
## was the column conservation metric assessed for all of these unfiltered SNPs?
CSR92 <- readDNAStringSet("data/contigs/filtered_CSR_.92.fasta")
BCNC.f <- BCNC[BCNC$contig %in% names(CSR92), ]

# Filter based on visual assessment -- real SNPs vs. possible alignment errors...
real <- read.csv("data/BC_NC_SNPS_visual_assessment.csv")
BCNC.f <- BCNC.f[BCNC.f$contig %in% real[real$visual_assess=="real", "contig"], ]

# Write filtered SNP information to file
x <- as.matrix(format(BCNC.f))
write.table(x, file="output/BC_NC_snps.txt", sep="\t", row.names=F, quote=F)
save(BCNC.f, file="output/BC_NC_snps.RData")

# Write contigs containing filtered SNPs to file
contigs.f <- subset(CSR92, names(CSR92) %in% BCNC.f$contig)
writeXStringSet(contigs.f, "output/BC_NC_contigs.fasta")
