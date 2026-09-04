print("R script started")

# SCRIPT PRUPOUSE:
# Identify ids frequencies and save the lowest to be removed in later plink


# IMP: to free memory usage (when working with data requiring big amounts of memory)
rm(list=ls(all.names=TRUE))
invisible(gc())


# Upload files.
  #Target data:
  tar_dup <- read.table("0C_targetID.duplicated.dupvar", header=TRUE, sep="\t")
  tar_frq <- read.table("0C_targetID.duplicated.frq", header=TRUE)
  
#QUITAR ->   #Reference data:
#QUITAR ->   #ref_dup <- read.table("0C_referenceID.duplicated.dupvar", header=TRUE, sep="\t")
#QUITAR ->   #ref_frq <- read.table("0C_referenceID.duplicated.frq", header=TRUE)
  

# Create a function to identify ids frequencies and save the lowest:
fun_lowest_freq <- function(df_dup, df_frq) {
  ids_exclude <- character(0)
  if (nrow(df_dup) > 0 ) {
        for (i in 1:nrow(df_dup)) {
      ids <- unlist(strsplit(df_dup[i, "IDS"], " "))
      df_ids.frq <- df_frq[df_frq$SNP %in% ids, ]
      max_id <- df_ids.frq[df_ids.frq$NCHROBS == max(df_ids.frq$NCHROBS), "SNP"]
      if (length(max_id)>1) {
        max_id <- max_id[1]
      }
      ids_exclude <- c(ids_exclude, ids[!(ids %in% max_id)])
    }
  }
  ids_exclude <<- ids_exclude
}


# Apply the function to each cohort
  #Target data:
  fun_lowest_freq(tar_dup, tar_frq)
  tar_exclude <- ids_exclude

#QUITAR ->   #Reference data:
#QUITAR ->   #fun_lowest_freq(ref_dup, ref_frq)
#QUITAR ->   #ref_exclude <- ids_exclude
  
  
#Save data
write.table(tar_exclude, "0C_targetID.duplicated.exclude.txt", col.names=FALSE, quote = FALSE, row.names = FALSE)
#QUITAR -> write.table(ref_exclude, "0C_referenceID.duplicated.exclude.txt", col.names=FALSE, quote = FALSE, row.names = FALSE)
 

print("R script finished")