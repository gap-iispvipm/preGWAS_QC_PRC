print("R script started")
# SCRIPT PURPOUSE:
# Create the conversion file for FID and IID of individuals.
# Add phenotypes to target data using the .txt file.


# IMP: to free memory usage (when working with data requiring big amounts of memory)
rm(list=ls(all.names=TRUE))
invisible(gc())


#Libraries
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringi))


#Upload raw files and rearrange data (this may not be necessary in other databases):
  #1: all_data.txt
  files <- list.files(path=getwd(), pattern="\\.txt")
  filename_uploaded <- files[which(apply(as.array(sapply("Target_data_", grepl, files)), 1, all))]
  if (length(filename_uploaded)>1) {
    filename_uploaded <- filename_uploaded[which(apply(as.array(sapply("_ARRANGED_HEADING", grepl, filename_uploaded)), 1, all))]
  }

  alldataA <- read.delim(filename_uploaded)
  alldataB <- read.delim(filename_uploaded, sep=" ")
  
    #Rearrange data in a proper dataframe.
    all_data0 <- data.frame(matrix(ncol=length(alldataB), nrow=NROW(alldataA)))
    names(all_data0) <- colnames(alldataB)
    for (i in 1:NROW(alldataA)) {
      a <- unlist(stri_split_fixed(alldataA[i, ], pattern=" "))
      if (length(a)!=length(all_data0)) {
        a1 <- a[1:length(all_data0)-1]
        a2 <- paste0(a[length(all_data0):length(a)], collapse=" ")
        a <- c(a1, a2)
      }
      all_data0[i, ] <- a
    }
    all_data0[all_data0=="control"] <- "Control"
    
  
  #2: Target data upload
  files <- list.files(path=getwd(), pattern="\\.fam")
  filename_uploaded <- files[which(apply(as.array(sapply("\\Target_data.fam$", grepl, files)), 1, all))]
  
    #Rearrange data in a proper dataframe.
    fam <- read.table(filename_uploaded, h=F)
    colnames(fam)[1:2] <- c("FID","IID")
    fam1 <- fam[,1:5]
    

#Create a dataframe relating Individuals ID (IID) to Phenotypes (PHENOTYPE).
  all_data1 <- unique(all_data0[c("MUESTRA.PGC", "IID", "PHENOTYPE")])
  names(all_data1)[1] <- "FID"

  all_data1[which(all_data1$PHENOTYPE=="FALSE"),]$PHENOTYPE <- "NA"
  all_data1[which(all_data1$PHENOTYPE=="NA"),]$PHENOTYPE <- -9
  all_data1[which(all_data1$PHENOTYPE=="Control"),]$PHENOTYPE <- 1
  all_data1[which(all_data1$PHENOTYPE=="Paciente"),]$PHENOTYPE <- 2
  

#merge both files by IID.
dat <- left_join(fam1, all_data1, by=c("FID", "IID"))
if (length(dat[is.na(dat$PHENOTYPE),]$PHENOTYPE)>0) {
  dat[is.na(dat$PHENOTYPE),]$PHENOTYPE <- -9
}
dat <- dat[order(dat$FID),]
dat[,5] <- as.numeric(dat[,5])
dat[,6] <- as.numeric(dat[,6])
dat1 <- dat[,c("FID", "IID", "PHENOTYPE")]


#Save data in files
write.table(dat1, "0A_phenotypes.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
write_tsv(all_data0,  "Target_data_all_data_ARRANGED_ALLCONTENT.txt")

print("R script finished")