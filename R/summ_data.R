# Load SNP data from filtered contigs
load("output/BC_NC_snps.RData")

# Aggregate SNP data for data frame to merge with various BLAST outputs
## Get position of each SNP in each contig
pos <- aggregate(data.frame(pos=BCNC$pos), by=list(contig=BCNC$contig), FUN=c)
## Get minimum coverage in both NC and BC groups of each SNP
BCNC$cov <- unlist(lapply(apply(BCNC[, c("NC", "BC")], 1, FUN=c), 
                          function(x) min(unlist(x)[unlist(x)!=0])))
cov <- aggregate(data.frame(cov=BCNC$cov), by=list(contig=BCNC$contig), FUN=c)
## Merge position and coverage information
snpstats <- merge(pos, cov)
## Get total number of SNPs in each contig
snpstats$nsnp <- sapply(snpstats$pos, length)
## Get the mean minimum coverage of all SNPs in each contig
snpstats$meancov <- sapply(snpstats$cov, mean)

# Load contig-to-nt blastn results
nt.annot <- read.table("output/BC_NC_contigs2nt_blastout_tabular.txt")
colnames(nt.annot) <- c("contig", "nt.hit", "nt.pident", "nt.length", "nt.mismatch", "nt.gapopen", "nt.qstart", "nt.qend", "nt.sstart", "nt.send", "nt.evalue", "nt.bitscore")
nt.annot.f <- nt.annot[!duplicated(nt.annot$contig), ]
snps.hits1 <- merge(snpstats, nt.annot.f, all=T)

# Load contig-to-coral blastx results
co.annot <- read.table("output/BC_NC_contigs2co_blastout_tabular.txt")
colnames(co.annot) <- c("contig", "co.hit", "co.pident", "co.length", "co.mismatch", "co.gapopen", "co.qstart", "co.qend", "co.sstart", "co.send", "co.evalue", "co.bitscore")
co.annot.f <- co.annot[!duplicated(co.annot$contig), ]
snps.hits2 <- merge(snps.hits1, co.annot.f, all=T)

# Load coral-to-uniprot-sprot blastp results
sp.annot <- read.table("output/BC_NC_cohits2sp_blastout_tabular.txt")
colnames(sp.annot) <- c("co.hit", "sp.hit", "sp.pident", "sp.length", "sp.mismatch", "sp.gapopen", "sp.qstart", "sp.qend", "sp.sstart", "sp.send", "sp.evalue", "sp.bitscore")
sp.annot.f <- sp.annot[!duplicated(sp.annot$sp.hit), ]
snps.hits3 <- merge(snps.hits2, sp.annot.f, all=T)

# Load coral-to-uniprot-TrEMBL blastp results
tr.annot <- read.table("output/BC_NC_cohits2tr_blastout_tabular.txt")
colnames(tr.annot) <- c("co.hit", "tr.hit", "tr.pident", "tr.length", "tr.mismatch", "tr.gapopen", "tr.qstart", "tr.qend", "tr.sstart", "tr.send", "tr.evalue", "tr.bitscore")
tr.annot.f <- tr.annot[!duplicated(tr.annot$tr.hit), ]
snps.hits4 <- merge(snps.hits3, tr.annot.f, all=T)

# Filter columns to summarize...
summ.f <- snps.hits4[,c("contig", "nsnp", "meancov", "nt.hit", "nt.evalue", "co.hit", "co.evalue", "sp.hit", "sp.evalue", "tr.hit", "tr.evalue")]
summ.f <- summ.f[with(summ.f, order(tr.evalue, -meancov)), ]

write.table(summ.f, file="output/BC_NC_summary.txt", row.names=F, quote=F)

