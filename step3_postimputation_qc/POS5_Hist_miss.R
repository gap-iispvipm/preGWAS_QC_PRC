print("R script started")

# read data into R:
fileslist <- list.files()
  #imiss:
  fileimiss <- fileslist[grep("\\.imiss", fileslist)]
  fileimiss <- tail(fileimiss, n=1)
  indmiss<-read.table(file=fileimiss, header=TRUE)

  #lmiss:
  filelmiss <- fileslist[grep("\\.lmiss", fileslist)]
  filelmiss <- tail(filelmiss, n=1)
  snpmiss <- read.table(file=filelmiss, header=TRUE)


# plot:

pdf("histimiss.pdf") #indicates pdf format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness", xlab="Missing rate") #selects column 6, names header of file

pdf("histlmiss.pdf") 
hist(snpmiss[,5],main="Histogram SNP missingness", xlab="Missing rate")  
invisible(dev.off()) # shuts down the current device

print("R script finished")