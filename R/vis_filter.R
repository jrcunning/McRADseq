load("output/BC_NC_snpstats.RData")

# Filter based on visual assessment -- real SNPs vs. possible alignment errors...
real <- read.csv("data/BC_NC_SNPS_visual_assessment.csv")
realsnps <- snpstats[snpstats$contig %in% real[real$visual=="real", "contig"], ]
save(realsnps, file="output/BC_NC_snpstats_real.RData")
write.table(as.matrix(format(realsnps)), file="output/BC_NC_snpstats_real.tsv",
            sep="\t", row.names=F, quote=F)
