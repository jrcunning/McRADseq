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

# Function to extract SNPs with fixed differences between two groups with coverage threshold
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

x <- as.matrix(format(BCNC))
write.table(x, file="output/BC_NC_snps.txt", sep="\t", row.names=F, quote=F)

save(BCNC, file="output/BC_NC_snps.RData")
