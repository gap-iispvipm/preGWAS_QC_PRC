print("R script started")


# read data into R:
het <- read.table("5_pruned.data.het", head=TRUE)


#1st Check heterozygosity rate:
pdf("5_heterozygosity.pdf")
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity rate", ylab="Frequency", main= "Heterozygosity rate (+-3 SD)")
abline(v=c((mean(het$HET_RATE)-3*sd(het$HET_RATE)), (mean(het$HET_RATE)+3*sd(het$HET_RATE))), col="red")
invisible(dev.off())


#2nd Heterozygosity outliers list:
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, "5_fail-het-qc.txt", row.names=FALSE)

print("R script finished")