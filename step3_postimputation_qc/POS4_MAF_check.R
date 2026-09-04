print("R script started")

# read data into R:
fileslist <- list.files()
maf_freq <- fileslist[grep("\\MAF_.frq", fileslist)]
maf_freq <- tail(maf_freq, n=1)
maf_freq <- read.table(file=maf_freq, header=TRUE, as.is=TRUE)

pdf("MAF_distribution.pdf")
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
invisible(dev.off())

print("R script finished")