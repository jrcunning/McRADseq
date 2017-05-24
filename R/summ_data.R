require(Biostrings)

# Load filtered SNP data
load("output/BC_NC_snpstats.RData")
## Get total number of SNPs in each contig
snpstats$nsnp <- sapply(snpstats$pos, length)
## Get the mean minimum coverage of all SNPs in each contig
snpstats$meancov <- sapply(snpstats$cov, mean)


# Set column names for tabular blast outputs
blastcolnames <- function(y, x) {
  c(y, paste(x, c("hit", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"), sep="."))
}


# Load contig-to-nt blastn results
nt.blastout <- file.info("output/BC_NC_contigs2nt_blastout_tabular.txt")
if (nt.blastout$size!=0) {
  nt.annot <- read.table("output/BC_NC_contigs2nt_blastout_tabular.txt")
  colnames(nt.annot) <- blastcolnames("contig", "nt")
  nt.annot.f <- nt.annot[!duplicated(nt.annot$contig), ]
  annot <- plyr::join(snpstats, nt.annot.f)
}


# Load contig-to-bacteria blastn results
bac.blastout <- file.info("output/BC_NC_contigs2bac_blastout_tabular.txt")
if (bac.blastout$size!=0) {
  bac.annot <- read.table("output/BC_NC_contigs2bac_blastout_tabular.txt")
  colnames(bac.annot) <- blastcolnames("contig", "bac")
  bac.annot.f <- bac.annot[!duplicated(bac.annot$contig), ]
  annot <- plyr::join(annot, bac.annot.f)
}


# Load contig-to-virus blastn results
vir.blastout <- file.info("output/BC_NC_contigs2vir_blastout_tabular.txt")
if (vir.blastout$size!=0) {
  vir.annot <- read.table("output/BC_NC_contigs2vir_blastout_tabular.txt")
  colnames(vir.annot) <- blastcolnames("contig", "vir")
  vir.annot.f <- vir.annot[!duplicated(vir.annot$contig), ]
  annot <- plyr::join(annot, vir.annot.f)
}


# Load contig-to-coral blastx results
co.blastout <- file.info("output/BC_NC_contigs2co_blastout_tabular.txt")
if (co.blastout$size!=0) {
  co.annot <- read.table("output/BC_NC_contigs2co_blastout_tabular.txt")
  colnames(co.annot) <- blastcolnames("contig", "co")
  co.annot.f <- co.annot[!duplicated(co.annot$contig), ]
  annot <- plyr::join(annot, co.annot.f)
}


# Load coral-to-swissprot blastp results
sp.blastout <- file.info("output/BC_NC_cohits2sp_blastout_tabular.txt")
if (sp.blastout$size!=0) {
  sp.annot <- read.table("output/BC_NC_cohits2sp_blastout_tabular.txt")
  colnames(sp.annot) <- blastcolnames("co.hit", "sp")
  sp.annot.f <- sp.annot[!duplicated(sp.annot$sp.hit), ]
  annot <- plyr::join(annot, sp.annot.f)
}


# Load coral-to-TrEMBL blastp results
tr.blastout <- file.info("output/BC_NC_cohits2tr_blastout_tabular.txt")
if (tr.blastout$size!=0) {
  tr.annot <- read.table("output/BC_NC_cohits2tr_blastout_tabular.txt")
  colnames(tr.annot) <- blastcolnames("co.hit", "tr")
  tr.annot.f <- tr.annot[!duplicated(tr.annot$tr.hit), ]
  annot <- plyr::join(annot, tr.annot.f)
}


# Filter columns to summarize...
keep <- which(unlist(lapply(strsplit(colnames(annot), split = ".", fixed = T), "[", 2)) %in% c("hit", "evalue"))
summ.f <- annot[, c(1, 4, 5, keep)]
summ.f <- summ.f[with(summ.f, order(contig)), ]

# Write output to file
write.table(summ.f, file="output/BC_NC_summary.txt", sep="\t", row.names=F, quote=F)

