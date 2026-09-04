print("R script started")

# read data into R:
fileslist <- list.files()
gender <- fileslist[grep("\\plink.sexcheck", fileslist)]
gender <- tail(gender, n=1)
gender <- read.table(file=gender, header=TRUE, as.is=TRUE)

pdf("Gender_check.pdf")
hist(gender[,6],main="Gender", xlab="F estimates")
invisible(dev.off())

pdf("Men_check.pdf")
male=subset(gender, gender$PEDSEX==1)
hist(male[,6],main="Men",xlab="F estimates")
invisible(dev.off())

pdf("Women_check.pdf")
female=subset(gender, gender$PEDSEX==2)
hist(female[,6],main="Women",xlab="F estimates")
invisible(dev.off())

print("R script finished")