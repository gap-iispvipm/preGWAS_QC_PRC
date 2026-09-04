print("R script started")

# SCRIPT PRUPOUSE:
# Create a PCA for Target and Reference dataset.


# IMP: to free memory usage (when working with data requiring big amounts of memory)
rm(list=ls(all.names=TRUE))
invisible(gc())

#Install packages required
if (!require("inspectdf")) {
  install.packages("inspectdf")
}


# Libraries:
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(inspectdf))
suppressMessages(require("RColorBrewer"))

# Read the .eigenvec file, the target and reference IDs, and reshape data:
eigenvec <- read.table("3_PCA_results.eigenvec", header=FALSE)
  names(eigenvec) <- c("FID", "IID", paste0("PC", c(1:20)))
ind.ref <- read.table("1C_QCreferencedata_individualsID.txt", header=FALSE)
  names(ind.ref) <- c("FID", "IID")
ind.target <- read.table("1C_QCtargetdata_individualsID.txt", header=FALSE)
  names(ind.target) <- c("FID", "IID")

  # Identify reference and target individuals:
  eigenvec$DATASET <- NA
  eigenvec[(eigenvec$FID %in% ind.ref$FID) & (eigenvec$IID %in% ind.ref$IID), "DATASET"] <- "REFERENCE"
  eigenvec[(eigenvec$FID %in% ind.target$FID) & (eigenvec$IID %in% ind.target$IID), "DATASET"] <- "TARGET"

  # If there is any IID duplicated, create unique IID:
  eigenvec$FID_IID <- paste(eigenvec$FID, eigenvec$IID, sep="_")
  eigenvec$IID_unique <- eigenvec$IID
  duplicated_IID <- eigenvec[duplicated(eigenvec$IID), "IID"]
  if (length(duplicated_IID)>0) {
    for (i in seq_along(duplicated_IID)) {
      eigenvec[(eigenvec$IID==duplicated_IID[i]) & (eigenvec$DATASET=="REFERENCE"), "IID_unique"] <- eigenvec[(eigenvec$IID==duplicated_IID[i]) & (eigenvec$DATASET=="REFERENCE"), "FID_IID"]
    }
  }


# read in the PED data:
PED <- read.table('PEDreferencedata_1000G.ped', header = TRUE, skip = 0, sep = '\t')
  PED <- PED[which(PED$Individual.ID %in% ind.ref$IID), ]
  

# Set Population, group, and colors:
  # FOR POPULATION:
  #from: http://www.internationalgenome.org/category/population/
  PED$Population <- factor(PED$Population, levels=c(
      "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
      "CLM","MXL","PEL","PUR",
      "CDX","CHB","CHS","JPT","KHV",
      "CEU","FIN","GBR","IBS","TSI",
      "BEB","GIH","ITU","PJL","STU"))
  
  eigenvec$COLOR <- NA  
  eigenvec[eigenvec$DATASET=="REFERENCE", "COLOR"] <- colorRampPalette(c(
      "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
      "forestgreen","forestgreen","forestgreen","forestgreen",
      "grey","grey","grey","grey","grey",
      "royalblue","royalblue","royalblue","royalblue","royalblue",
      "purple","purple","purple","purple","purple"))(length(unique(PED$Population)))[factor(PED$Population)]
  eigenvec[eigenvec$DATASET=="TARGET", "COLOR"] <- "#FF0000"
  
  eigenvec$POPULATION <- NA
  eigenvec[eigenvec$DATASET=="REFERENCE", "POPULATION"] <-as.character(factor(PED$Population))
  eigenvec[eigenvec$DATASET=="TARGET", "POPULATION"] <- "TARGET"
  
  
  # FOR GROUPS: 
  eigenvec$GROUP <- NA
  eigenvec[eigenvec$DATASET=="TARGET", "GROUP"] <- "TARGET"
  eigenvec[grepl(paste(c("ACB","ASW","ESN","GWD","LWK","MSL","YRI"), collapse="|"), eigenvec$POPULATION), "GROUP"] <- "AFR" 
  eigenvec[grepl(paste(c("CLM","MXL","PEL","PUR"), collapse="|"), eigenvec$POPULATION), "GROUP"] <- "AMR"
  eigenvec[grepl(paste(c("CDX","CHB","CHS","JPT","KHV"), collapse="|"), eigenvec$POPULATION), "GROUP"] <- "EAS"
  eigenvec[grepl(paste(c("CEU","FIN","GBR","IBS","TSI"), collapse="|"), eigenvec$POPULATION), "GROUP"] <- "EUR"
  eigenvec[grepl(paste(c("BEB","GIH","ITU","PJL","STU"), collapse="|"), eigenvec$POPULATION), "GROUP"] <- "SAS"
  
  myColRamp1 <- colorRampPalette(c("lemonchiffon", "yellow2"))
  myColRamp2 <- colorRampPalette(c("palegreen", "forestgreen"))
  myColRamp3 <- colorRampPalette(c("whitesmoke", "ivory4"))
  myColRamp4 <- colorRampPalette(c("lightcyan", "darkblue"))
  myColRamp5 <- colorRampPalette(c("plum1", "purple4"))
   
  col.legend.pop <- c(myColRamp1(length(c("ACB","ASW","ESN","GWD","LWK","MSL","YRI"))), myColRamp2(length(c("CLM","MXL","PEL","PUR"))), myColRamp3(length(c("CDX","CHB","CHS","JPT","KHV"))), myColRamp4(length(c("CEU","FIN","GBR","IBS","TSI"))), myColRamp5(length(c("BEB","GIH","ITU","PJL","STU"))), "#FF0000")
  names(col.legend.pop) <- c("ACB","ASW","ESN","GWD","LWK","MSL","YRI",
    "CLM","MXL","PEL","PUR",
    "CDX","CHB","CHS","JPT","KHV",
    "CEU","FIN","GBR","IBS","TSI",
    "BEB","GIH","ITU","PJL","STU", "TARGET")
  
  col.legend.group <- c("yellow","forestgreen","grey","royalblue","purple","red")
  names(col.legend.group) <- c("AFR", "AMR", "EAS", "EUR", "SAS", "TARGET")


# Plot the first two principal components

  # A) TARGET VS REFERENCE:
  pdf("3_PCA.datasets.pdf")
  
  ggplot(data=eigenvec, aes(x=PC1, y=PC2, color=DATASET)) +
    geom_point() +
    theme_linedraw() +
    labs(title="PCA: Target vs Reference Dataset") +
    theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black")) +
    scale_color_manual(values=c("black", "red"))

  # B) FOR GROUPS:
  ggplot(data=eigenvec, aes(x=PC1, y=PC2, color=GROUP)) +
    geom_point() +
    theme_linedraw() +
    labs(title="PCA: Target vs Reference Dataset") +
    theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black")) +
    scale_color_manual(values=col.legend.group)
  
  # C) FOR POPULATIONS:
  ggplot(data=eigenvec, aes(x=PC1, y=PC2, color=POPULATION)) +
    geom_point() +
    theme_linedraw() +
    scale_color_manual(values=col.legend.pop) +
    labs(title="PCA: Target vs Reference Dataset") +
    theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black"))
  
  invisible(dev.off())

  
# CLUSTER REFERENENCE DATASET:
  
  # Convert data columns in factors and subset for reference group:
  eigenvec_ref <- eigenvec[eigenvec$DATASET=="REFERENCE" & eigenvec$POPULATION!="TARGET",]
  eigenvec_ref[, c("DATASET", "COLOR", "POPULATION", "GROUP")] <- lapply(eigenvec_ref[, c("DATASET", "COLOR", "POPULATION", "GROUP")], as.factor)
  
  # Check summary of reference data:
  pdf("3_Reference.summary.pdf")
  
    df_ref_group<- summary(eigenvec_ref$GROUP)
    df_ref_group %>% barplot(ylim=c(0, max(summary(df_ref_group)) + 100)) %>% text(x=barplot(df_ref_group, ylim=c(0, max(summary(df_ref_group)) + 1)), y=df_ref_group + 70, labels=df_ref_group)
    
    df_ref_percentage <- as.matrix(table(eigenvec_ref$GROUP) %>% prop.table() * 100)
    df_ref_percentage %>% barplot(legend=TRUE, ylim=c(0, max(df_ref_percentage) + 100), args.legend=list(bty="n", x="top", ncol=5))
    
    # feature distribution
    eigenvec_ref[, c("PC1", "PC2")] %>% inspect_num() %>% show_plot()
    
    # class distribution
    eigenvec_ref[, c("POPULATION", "GROUP", "DATASET")] %>%  inspect_cat() %>% show_plot()
    
    # group plotting
    ggplot(data=eigenvec[eigenvec$DATASET=="REFERENCE",], aes(x=PC1, y=PC2, color=POPULATION)) +
      geom_point() +
      theme_linedraw() +
      scale_color_manual(values=col.legend.pop) +
      labs(title="PCA: Reference Dataset") +
      theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black"))

  invisible(dev.off())
  
  # Search means and 3SD per each group of groups:
  groups.mean = groups.3sd <- data.frame(PC1=rep(NA, length(levels(eigenvec_ref$GROUP))), PC2=rep(NA, length(levels(eigenvec_ref$GROUP))), row.names=levels(eigenvec_ref$GROUP))
  for (i in seq_along(levels(eigenvec_ref$GROUP)))  {
    group <- levels(eigenvec_ref$GROUP)[i]
    b1 <- eigenvec_ref[eigenvec_ref$GROUP==group, c("PC1","PC2")]
    groups.mean[rownames(groups.mean)==group,] <- sapply(b1, function(x) mean(x))
    groups.3sd[rownames(groups.3sd)==group,] <- sapply(b1, function(x) 3*sd(x))
  }
  groups.max <- groups.mean + groups.3sd
  groups.min <- groups.mean - groups.3sd
  
  #plot centroids and 3SD:
    #prepare dataframe:
    eigenvec_ref_centroid <- eigenvec_ref[0,]
    eigenvec_ref_centroid[nrow(eigenvec_ref_centroid)+length(levels(eigenvec_ref$GROUP)), ] <- NA
    eigenvec_ref_centroid$PC1 <- groups.mean$PC1
    eigenvec_ref_centroid$PC2 <- groups.mean$PC2
    eigenvec_ref_centroid$GROUP <- "Group centroid"
    eigenvec_ref_centroid$COLOR <- "#000000"
  
    eigenvec_ref_centroid <- rbind(eigenvec_ref, eigenvec_ref_centroid)
    
    #prepare colors:
    col.legend.centroid <- "black"
    names(col.legend.centroid) <- "Group centroid"
    col.legend.centroid <- append(col.legend.group, col.legend.centroid)

    #plot:  
    pdf("3_PCA.centroids.pdf")
    
      #Reference Data Set Centroids
      ggplot(data=eigenvec_ref_centroid, aes(x=PC1, y=PC2, color=GROUP)) +
        geom_point() +
        theme_linedraw() +
        labs(title="PCA: Reference Dataset Centroids") +
        theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black")) +
        scale_color_manual(values=col.legend.centroid)+
        annotate("rect", xmin=groups.min$PC1, xmax=groups.max$PC1, ymin=groups.min$PC2, ymax=groups.max$PC2, alpha=0.2, color=col.legend.centroid[rownames(groups.min)], fill="grey")
      
      #European Reference Data Set Centroids
      xmin <- groups.min[rownames(groups.min)=="EUR",]$PC1
      xmax <- groups.max[rownames(groups.max)=="EUR",]$PC1
      ymin <- groups.min[rownames(groups.min)=="EUR",]$PC2
      ymax <- groups.max[rownames(groups.max)=="EUR",]$PC2
      
      eigenvec_ref_eur <- eigenvec_ref_centroid[eigenvec_ref_centroid$GROUP=="EUR", ]
      eigenvec_ref_eur$AREA <- "EUR ref. include"
      eur_ref_out <- c(eigenvec_ref_eur[eigenvec_ref_eur$PC1>xmax, "IID_unique"], eigenvec_ref_eur[eigenvec_ref_eur$PC1<xmin, "IID_unique"], eigenvec_ref_eur[eigenvec_ref_eur$PC2>ymax, "IID_unique"], eigenvec_ref_eur[eigenvec_ref_eur$PC2<ymin, "IID_unique"])
      if (length(eur_ref_out)>0) {
        eigenvec_ref_eur[eigenvec_ref_eur$IID_unique %in% eur_ref_out, ]$AREA <- "EUR ref. exclude"
      }

      col.legend.area1 <- c(col.legend.group[names(col.legend.group)=="EUR"], "cyan")
      names(col.legend.area1) <- c("EUR ref. include", "EUR ref. exclude")
      
      ggplot(data=eigenvec_ref_eur, aes(x=PC1, y=PC2, color=AREA)) +
        geom_point() +
        theme_linedraw() +
        labs(title="PCA: Reference Dataset - European") +
        theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black")) +
        scale_color_manual(values=col.legend.area1)+
        annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=0.2, color=col.legend.centroid[rownames(groups.min)=="EUR"], fill="grey")
    
      #European Target Data Set:
      eigenvec_target <- eigenvec[eigenvec$DATASET=="TARGET", ]
      eigenvec_target$AREA <- "EUR target include"
      eur_target_out <- c(eigenvec_target[eigenvec_target$PC1>xmax, "IID_unique"], eigenvec_target[eigenvec_target$PC1<xmin, "IID_unique"], eigenvec_target[eigenvec_target$PC2>ymax, "IID_unique"], eigenvec_target[eigenvec_target$PC2<ymin, "IID_unique"])
      if (length(eur_target_out)>0) {
        eigenvec_target[eigenvec_target$IID_unique %in% eur_target_out, ]$AREA <- "EUR target exclude"  
      }
      eur_target_in <- eigenvec_target[eigenvec_target$AREA=="EUR target include", "IID_unique"]
       
      col.legend.area2 <- c(col.legend.group[names(col.legend.group)=="TARGET"], "lightpink")
      names(col.legend.area2) <- c("EUR target include", "EUR target exclude")
        
      ggplot(data=eigenvec_target, aes(x=PC1, y=PC2, color=AREA)) +
        geom_point() +
        theme_linedraw() +
        labs(title="PCA: Target Dataset - European") +
        theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black")) +
        scale_color_manual(values=col.legend.area2)+
        annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=0.2, color=col.legend.centroid[rownames(groups.min)=="EUR"], fill="grey")
      
      #European Target & Reference Data Set:
      eigenvec_ref_target_eur <- rbind(eigenvec_ref_eur, eigenvec_target)
      col.legend.area <- append(col.legend.area1, col.legend.area2)
      
      ggplot(data=eigenvec_ref_target_eur, aes(x=PC1, y=PC2, color=AREA)) +
        geom_point() +
        theme_linedraw() +
        labs(title="PCA: Reference & Target Dataset - European") +
        theme(plot.title=element_text(face="bold", size=18, hjust=0.5), axis.title=element_text(size=14,face="bold"), axis.text=element_text(size=12), legend.title=element_blank(),  legend.position="right", legend.background=element_rect(linewidth=0.2, linetype="solid", colour ="black")) +
        scale_color_manual(values=col.legend.area)+
        annotate("rect", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=0.2, color=col.legend.centroid[rownames(groups.min)=="EUR"], fill="grey")
    
    invisible(dev.off())

     
  #Keep target individuals:
  target.individuals.excluded <- eigenvec_target[eigenvec_target$IID_unique %in% eur_target_out, c("FID","IID")]
  target.individuals.included <- eigenvec_target[eigenvec_target$IID_unique %in% eur_target_in, c("FID","IID")]  
  
  a <- target.individuals.excluded[NA,]
  a$IID <- paste(target.individuals.excluded$FID, target.individuals.excluded$IID, sep="_")
  a$FID <- 0
  
  write.table(target.individuals.excluded, paste0("3_target.individuals.excluded", ".txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
  write.table(target.individuals.included, paste0("3_target.individuals.included", ".txt"), sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
 
  
print("R script finished")