print("R script started")
# SCRIPT PRUPOUSE:
# Create the conversion file for FID and IID of individuals.


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
    filename_uploaded <- filename_uploaded[which(apply(as.array(sapply("_ARRANGED_ALLCONTENT", grepl, filename_uploaded)), 1, all))]
  }
  all_data0 <- read.delim(filename_uploaded)
    
  
  #2: Target data upload
  files <- list.files(path=getwd(), pattern="\\.fam")
  filename_uploaded <- files[which(apply(as.array(sapply("1H_IMPtargetdata.fam", grepl, files)), 1, all))]
  
    #Rearrange data in a proper dataframe.
    fam0 <- read.table(filename_uploaded, h=F)
    colnames(fam0)[c(1:2,6)] <- c("FID_wrong","IID_wrong", "PHENOTYPE")
    fam1 <- fam0
      # Arrange FID
      fam1$FID <- gsub("\\_.*", "", fam0$IID_wrong)
      # Arrange IID
      fam1$IID <- gsub("^[^_]+_", "", fam0$IID_wrong)
    fid_iid <- fam1[ , c(1:2,7:8)]

#Save data in files
write.table(fid_iid, "1H_FID.IID.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)

print("R script finished")