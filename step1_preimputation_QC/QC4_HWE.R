print("R script started")

hwe<-read.table (file="3B_plink.hwe", header=TRUE)
pdf("4_histhwe.pdf")
hist(hwe[,9],main="Histogram HWE", xlab = "P-VALUE")
invisible(dev.off())

hwe_below<-hwe[hwe$P<0.00001,]
pdf("4_histhwe_below_0.00001.pdf")
hist(hwe_below[,9],main="Histogram HWE: strongly deviating SNPs only (P<0.00001)", xlab = "P-VALUE")
invisible(dev.off())

hwe_above<-hwe[hwe$P>0.00001,]
pdf("4_histhwe_above_0.00001.pdf")
hist(hwe_above[,9],main="Histogram HWE: NOT strongly deviating SNPs only (P>0.00001)", xlab = "P-VALUE")
invisible(dev.off())

print("R script finished")