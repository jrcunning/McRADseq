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
## contigs92.txt is a list of the names of all contigs with > 92% conservation (see Makefile)
contigs92 <- read.table("data/contigs/contigs92.txt", header=F)
BCNC.f <- BCNC[BCNC$contig %in% contigs92$V1, ]

### Extract .sam files for visual assessment
write(paste0("samFiles/", unique(BCNC.f$contig), ".sam"), file = "data/contigs/contigs_cov5_csr92.txt")
system('cd data/contigs && gunzip < samFiles.tar.gz | tar -x -v --files-from contigs_cov5_csr92.txt -f -')

# Summarize SNP information and write to file
## Get position of each SNP in each contig
pos <- aggregate(data.frame(pos=BCNC.f$pos), by=list(contig=BCNC.f$contig), FUN=c)
## Get minimum coverage in both NC and BC groups of each SNP
BCNC.f$cov <- unlist(lapply(apply(BCNC.f[, c("NC", "BC")], 1, FUN=c), 
                            function(x) min(unlist(x)[unlist(x)!=0])))
cov <- aggregate(data.frame(cov=BCNC.f$cov), by=list(contig=BCNC.f$contig), FUN=c)
## Merge position and coverage information
snpstats <- merge(pos, cov)
snpstats <- snpstats[order(snpstats$contig), ]

# Write SNP stats to table and RData file
x <- as.matrix(format(snpstats))
write.table(x, file="output/BC_NC_snpstats.txt", sep="\t", row.names=F, quote=F)
save(snpstats, file="output/BC_NC_snpstats.RData")

