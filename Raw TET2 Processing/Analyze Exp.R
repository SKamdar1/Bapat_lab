simplemasscalc <- function(mytable, GOI, wilcox_ids, wilcox_names, StageCol, ScoreCol, GroupCol) {
  #Remove columns without data
  
  GeneTable <- mytable
  
  GeneTable <- GeneTable[, colSums(is.na(GeneTable)) !=nrow(GeneTable)]
  
  #Remove rows without data
  GeneTable <- GeneTable[rowSums(is.na(GeneTable)) !=ncol(GeneTable),]
  
  #Transform non-grouping variables to numeric
  for (i in ((IncrementCounter+1):ncol(GeneTable))) {
    GeneTable[,i] <- as.numeric(as.character(GeneTable[,i]))
  }
  
  GeneOutput <<- GeneTable
  GeneTable <<- GeneOutput
  
  WilcoxVals <<- as.numeric(vector())
  KruskalVals <<- as.numeric(vector())
  NameStorageQuant <<- as.character(vector())
  NameStorageYN <<- as.character(vector())
  DunnTable <<- list()
  
  if (length(wilcox_ids) != length(wilcox_names)) {
    stop("Your Wilcox input variables are not equal in length.")
  }
  
  for (i in 1:length(wilcox_ids)) {
    wilcox_print_PMR(wilcox_ids[i],wilcox_names[i])
  }

  
  #Always performs regular/path stage tests.
  SimpleStageKruskal(StageCol,"PathStageSimpleTest")
  SimpleStageKruskal(StageCol,"PathStageCollapsed",collapse=TRUE)
  GeneTable[,StageCol] <- as.character(GeneTable[,StageCol])
  GeneTable[,StageCol][GeneTable[,StageCol]=='T4'|GeneTable[,StageCol]=='T3b'] <- 'T3b+'
  GeneTable[,StageCol] <- as.factor(GeneTable[,StageCol])
  KruskalCalc(StageCol,"PathStageAll")
  
  #Gleason Score calculated below (blank values become Normals if present, with another series of tests run without the normal cases afterward.)
  
  #Tidy the GS data
  GeneTable[,ScoreCol]<- as.character(GeneTable[,ScoreCol])
  GeneTable[,ScoreCol][GeneTable[,ScoreCol]=='0'] <- "NULL"
  GeneTable[,ScoreCol][GeneTable[,ScoreCol] == '' & GeneTable[,GroupCol] == "Solid Tissue Normal"] <- 'Normal'
  GeneTable[,ScoreCol] <- as.factor(GeneTable[,ScoreCol])
  
  GleasonCalc(groupid=ScoreCol)
  
  #Remove the normal cases(if applicable) and re-run
  if (any(GeneTable[,ScoreCol]=='Normal')) {
    GeneTable[,ScoreCol] <- as.character(GeneTable[,ScoreCol])
    GeneTable[,ScoreCol][GeneTable[,ScoreCol] == 'Normal'] <- ''
    GeneTable[,ScoreCol] <- as.factor(GeneTable[,ScoreCol])
    GleasonCalc(prefix="TumorOnly",groupid=ScoreCol)
  }
  
  write.table(GeneOutput, file=sprintf("%sProbeStats.txt", GOI), sep="\t", col.names=TRUE)
}


TCGA_Exp <- read.table("TCGA_Exp.txt",sep="\t",header=TRUE)
TCGA_Exp_t <- as.data.frame(t(TCGA_Exp))
colnames(TCGA_Exp_t) <- as.character(TCGA_Exp$Hybridization.REF)
TETStorage <- TCGA_Exp_t$`TET2|54790`
TCGA_Exp_t$`TET2|54790` <- as.numeric(as.character(TCGA_Exp_t$`TET2|54790`))
TCGA_Exp_TET_sorted <- TCGA_Exp_t[order(TCGA_Exp_t$`TET2|54790`),]
TCGA_Exp_TET_sorted <- TCGA_Exp_TET_sorted[,-c(2,6)]

TCGA_Normals <- TCGA_Exp_TET_sorted[which(TCGA_Exp_TET_sorted$`Sample Type`=="Solid Tissue Normal"),]
TCGA_Cancers <- TCGA_Exp_TET_sorted[which(TCGA_Exp_TET_sorted$`Sample Type`=="Primary solid Tumor"),]

TET10_90 <- quantile(TCGA_Cancers$`TET2|54790`,probs=c(0.1,0.9),na.rm=TRUE)
TET_lowest10 <- TCGA_Cancers[which(TCGA_Cancers$`TET2|54790`<=283.0601),]
TET_lowest10 <- rbind(TET_lowest10,TCGA_Normals)
TET_highest10 <- TCGA_Cancers[which(TCGA_Cancers$`TET2|54790`>=830.1193),]
TET_highest10 <- rbind(TET_highest10,TCGA_Normals)

TET20_80 <- quantile(TCGA_Cancers$`TET2|54790`,probs=c(0.2,0.8),na.rm=TRUE)
TET_lowest20 <- TCGA_Cancers[which(TCGA_Cancers$`TET2|54790`<=374.3045),]
TET_lowest20 <- rbind(TET_lowest20,TCGA_Normals)
TET_highest20 <- TCGA_Cancers[which(TCGA_Cancers$`TET2|54790`>=701.8048),]
TET_highest20 <- rbind(TET_highest20,TCGA_Normals)

#Use functions from TET and TCGA Processing.R

#Generating heatmap (from start, for genes with lowered expression)
tetexp <- read.table("RNASeqCRISPRData.txt",sep="\t",header=TRUE)
tetexp_sigdown <- tetexp[which((tetexp$FCB<0.75 & tetexp$Pvalue_1vRW<0.05)|(tetexp$FCD<0.75 & tetexp$Pvalue_2vRW<0.05)),]

TCGA_Exp <- read.table("TCGA_Exp.txt",sep="\t",header=TRUE)
row.names(TCGA_Exp) <- TCGA_Exp$Hybridization.REF
TCGA_Exp$Hybridization.REF <- as.character(TCGA_Exp$Hybridization.REF)
TCGA_Exp$Hybridization.REF <- gsub("\\|.*","",TCGA_Exp$Hybridization.REF)

sig_genes_down <- c("Sample Type","Gleason Score", "pathologic_T", "Biochemical recurrence", "TET2",as.character(tetexp_sigdown$GeneName))
tcga_exp_down <- TCGA_Exp %>% rownames_to_column('bloop') %>% filter(Hybridization.REF %in% sig_genes_down) %>% column_to_rownames('bloop')
tcga_exp_down <- as.data.frame(t(tcga_exp_down))
tcga_exp_down <- tcga_exp_down[-1,]


TCGA_Exp_t <- tcga_exp_down
TCGA_Exp_t$`TET2|54790` <- as.numeric(as.character(TCGA_Exp_t$`TET2|54790`))
TCGA_Exp_TET_sorted <- TCGA_Exp_t
TCGA_Normals <- TCGA_Exp_TET_sorted[which(TCGA_Exp_TET_sorted$`Sample Type`=="Solid Tissue Normal"),]
TCGA_Cancers <- TCGA_Exp_TET_sorted[which(TCGA_Exp_TET_sorted$`Sample Type`=="Primary solid Tumor"),]
TET_lowest10 <- TCGA_Cancers[which(TCGA_Cancers$`TET2|54790`<=283.0601),]
TET_lowest10 <- rbind(TET_lowest10,TCGA_Normals)

#Run TET and TCGA processing mass pvalue calcs

library(dplyr)
library(tibble)

simplemasscalc(mytable=TET_lowest10,GOI="TET_TCGA_Exp_Down", wilcox_ids=c(1,4),wilcox_names=c("TvN","Recurrence"),StageCol=3,ScoreCol=2,GroupCol=1)

#Filter GeneOutput to find genes with significantly increased expression in TvN
#Significant expression
sig_exp_mass_down_indices <- which(as.numeric(GeneOutput[(nrow(tcga_exp_down)+1),])<0.05)
sig_exp_mass_down <- GeneOutput[,c(1:4,sig_exp_mass_down_indices)]

#Significantly increased expression
#Convert datacols to numeric
sig_exp_mass_down[,c(5:ncol(sig_exp_mass_down))] <- lapply(sig_exp_mass_down[,c(5:ncol(sig_exp_mass_down))],as.numeric)
#Summarize mean across groups
directiontest <- aggregate(sig_exp_mass_down[,c(5:ncol(sig_exp_mass_down))],list(sig_exp_mass_down$`Sample Type`),mean)
direction_indices <- which(directiontest[1,]<directiontest[2,])
sig_exp_mass_down_concordant <- sig_exp_mass_down[,c(1:4,direction_indices+3)]

#Output the concordant genes into a text file.
write.table(sig_exp_mass_down_concordant, file="BF Corrected Sig Down Exp Matched TCGA Samples 118Genes.txt",sep="\t",col.names=TRUE)

###########################################

#Gathering data for heatmap

library(tidyr)
library(RColorBrewer)
library(viridis)

TET_lowest10[,c(5:ncol(TET_lowest10))] <- lapply(TET_lowest10[,c(5:ncol(TET_lowest10))],function(x) as.numeric(as.character(x)))
TET_lowest10$Sample_Name <- row.names(TET_lowest10)

gathered_TET_data <- gather(data=TET_lowest10,key=Gene, value=Expression, -c(1:4,3551))

significant_genes <- colnames(sig_exp_mass_down_concordant)
TET_lowest10_copy <- TET_lowest10[,c(3551,which(colnames(TET_lowest10) %in% significant_genes))]

TET_samplematrix <- as.matrix(TET_lowest10_copy[,-c(1:5)])
TET_samplematrix2 <- t(TET_samplematrix)

classfactor <- as.factor(TET_lowest10_copy$`Sample Type`)
my_col <- brewer.pal(2,"RdBu")[classfactor]

#y_col <- brewer.pal(9,"Set1")[classfactor]

coul <- colorRampPalette(brewer.pal(8,"Purples"))(8)
coul <- colorRampPalette(c("black","red"))

#heatmap for stratification of TET samples alone vs. normals
heatmap(TET_samplematrix2,scale="row",labRow="",labCol=c(1:78),Rowv=NA,ColSideColors=my_col,col=coul)
heatmap(TET_samplematrix2,scale="row",labRow="",labCol="",Rowv=NA,ColSideColors=my_col,col=viridis(7))

#heatmap for stratification of all cancers vs. normals based on these ~700 genes
colnames(TCGA_Exp_t) <- gsub("\\|",".",colnames(TCGA_Exp_t))
significant_genes[c(1:4)] <- c("Sample Type", "Gleason Score","pathologic_T","Biochemical recurrence")
TCGA_Exp_t_copy <- TCGA_Exp_t[,which(colnames(TCGA_Exp_t) %in% significant_genes)]

TCGA_Exp_t_copy2 <- TCGA_Exp_t_copy[,-c(1:4)]
TCGA_Exp_t_copy2 <- data.frame(sapply(TCGA_Exp_t_copy2, function(x) as.numeric(as.character(x))))
TCGA_samplematrix <- as.matrix(TCGA_Exp_t_copy2)
TCGA_samplematrix2 <- t(TCGA_samplematrix)

classfactor <- as.factor(TCGA_Exp_t_copy$`Sample Type`)
my_col <- brewer.pal(2,"RdBu")[classfactor]

heatmap(TCGA_samplematrix2,scale="row",labRow="",Rowv=NA,ColSideColors=my_col,col=viridis(7))

#heatmap for stratification of TET samples alone for the trimmed list of ~120 genes

#heatmap for all TCGA samples on those ~120 genes
colnames(TCGA_Exp_t) <- gsub("\\|.*","",colnames(TCGA_Exp_t))
genes_120 <- read.table("60_genes.txt",sep="\t",header=TRUE)
TCGA_Exp_t_copy <- TCGA_Exp_t[,c(1:4,which(colnames(TCGA_Exp_t) %in% genes_120[,1]))]
classfactor <- as.factor(TCGA_Exp_t_copy$`Sample Type`)
my_col <- brewer.pal(2,"RdBu")[classfactor]

TCGA_Exp_t_copy2 <- data.frame(sapply(TCGA_Exp_t_copy[,-c(1:4)], function(x) as.numeric(as.character(x))))
TCGA_samplematrix <- as.matrix(TCGA_Exp_t_copy2)
TCGA_samplematrix2 <- t(TCGA_samplematrix)

TCGA_samplematrix2 <- TCGA_samplematrix2[,-1]
heatmap(TCGA_samplematrix2,scale="row",labCol="",Rowv=NA,ColSideColors=my_col,col=viridis(7))

#Function for methylation

TET_lowest10_Cancers <- TET_lowest10[-which(TET_lowest10$`Sample Type`=="Solid Tissue Normal"),]

AutoMeth_TCGA <- function(clonetableid) {
  #Initialize storage Data Frames
  file1 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file2 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file3 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file4 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file5 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file6 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file7 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file8 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file9 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  file10 <- as.data.frame(matrix(data=rep("NA"), nrow=1, ncol=ncol(clonetableid)),stringsAsFactors=FALSE)
  
  filestoragelist <- c("file1","file2","file3","file4","file5","file6","file7","file8","file9","file10")
  
  SigMethResultList <- list()
  
  print("Files initialized.")
  
  #Sort all methylation peaks from file into a list based on chromosomal location
  for (i in 1:nrow(clonetableid)) {
    if (clonetableid$chr[i]=="chr1") {
      file1[(nrow(file1))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr12") {
      file3[(nrow(file3))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr13") {
      file3[(nrow(file3))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr15") {
      file4[(nrow(file4))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr16") {
      file4[(nrow(file4))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr18") {
      file5[(nrow(file5))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr19") {
      file5[(nrow(file5))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr20") {
      file6[(nrow(file6))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr21") {
      file6[(nrow(file6))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr3") {
      file7[(nrow(file7))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr5") {
      file8[(nrow(file8))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr7") {
      file9[(nrow(file9))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr9") {
      file10[(nrow(file10))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chrX") {
      file10[(nrow(file10))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chrY") {
      file10[(nrow(file10))+1,] <- as.character(unlist(clonetableid[i,]))
    } else if (clonetableid$chr[i] == "chr10") {
      if (as.numeric(as.character(clonetableid$end[i])) < 12110680) {
        file1[(nrow(file1))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file2[(nrow(file2))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr11") {
      if (as.numeric(as.character(clonetableid$end[i])) < 134628435) {
        file2[(nrow(file2))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file3[(nrow(file3))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr14") {
      if (as.numeric(as.character(clonetableid$end[i])) < 103839089) {
        file3[(nrow(file3))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file4[(nrow(file4))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr17") {
      if (as.numeric(as.character(clonetableid$end[i])) < 36646404) {
        file4[(nrow(file4))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file5[(nrow(file5))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr2") {
      if (as.numeric(as.character(clonetableid$end[i])) < 3606332) {
        file5[(nrow(file5))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file6[(nrow(file6))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr22") {
      if (as.numeric(as.character(clonetableid$end[i])) < 24671138) {
        file6[(nrow(file6))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file7[(nrow(file7))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr4") {
      if (as.numeric(as.character(clonetableid$end[i])) < 175135938) {
        file7[(nrow(file7))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file8[(nrow(file8))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr6") {
      if (as.numeric(as.character(clonetableid$end[i])) < 68599197) {
        file8[(nrow(file8))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file9[(nrow(file9))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else if (clonetableid$chr[i] == "chr8") {
      if (as.numeric(as.character(clonetableid$end[i])) < 37553078) {
        file9[(nrow(file9))+1,] <- as.character(unlist(clonetableid[i,]))
      } else {
        file10[(nrow(file10))+1,] <- as.character(unlist(clonetableid[i,]))
      }
    } else {
      print("There were no matches found for this entry.")
    }
  }
  print("Files are sorted.")
  #For each file of the ten TCGA files:
  for (i in 1:length(filestoragelist)) {
    print("The file being processed is")
    print(i)
    #Import each of the files and then rm them, so each master file is only kept in memory while it is needed
    TCGA_file <- read.table(sprintf('%s_TCGA.txt',filestoragelist[i]),sep="\t",header=TRUE)
    row.names(TCGA_file) <- TCGA_file$Composite.Element.REF
    TCGA_file$Gene_Symbol <- as.character(TCGA_file$Gene_Symbol)
    TCGA_file$Gene_Symbol[c(1:6)] <- as.character(TCGA_file$Composite.Element.REF[c(1:6)])
    #Pull the significant methylation regions from the files for CR1 or CR2
    print("Pulling sig meth regions")
    sig_meth <- c("Sample Type","Gleason Score","pathologic_T","Biochemical recurrence", as.character(get(filestoragelist[i])$V4))
    
    #Perform the grepping and remove the original file from memory
    print("Performing the grep")
    sig_meth_TCGA <- TCGA_file %>% rownames_to_column('bloop') %>% filter(Gene_Symbol %in% sig_meth) %>% column_to_rownames('bloop')
    rm(TCGA_file)
    
    print("The grep is complete")
    #Transpose to acceptable format for mass pvalue calculations
    
    sig_meth_TCGA <- as.data.frame(t(sig_meth_TCGA))
    
    #Identify lowest 10% of TET2 expression tumors
    TCGA_sig_Normals <- sig_meth_TCGA[which(sig_meth_TCGA$`Sample Type`=="Solid Tissue Normal"),]

    IDs <- substr(row.names(TET_lowest10_Cancers),start=1,stop=16)
    my_IDs <- substr(row.names(sig_meth_TCGA),start=12,stop=27)
    TCGA_sig_Cancers <- sig_meth_TCGA[c(1:4,which(my_IDs %in% IDs)),]
    
    sig_meth_TCGA <- rbind(TCGA_sig_Cancers,TCGA_sig_Normals)
    
    #Add the CpG IDs to gene symbols
    
    CombiGene <- sprintf('%s_%s',as.character(unlist(sig_meth_TCGA[2,])),as.character(unlist(sig_meth_TCGA[1,])))
    
    #Store co-ordinates for each CpG in a separate object and replace the first row with CombiGene
    CoordStorage <- t(sig_meth_TCGA[c(1:4),])
    
    colnames(sig_meth_TCGA) <- CombiGene
    colnames(sig_meth_TCGA)[c(1:4)] <- as.character(unlist(sig_meth_TCGA[1,c(1:4)]))
    
    sig_meth_TCGA <- sig_meth_TCGA[-c(1:4),]
    
    print("Pre-processing complete; mass pvalue calcs beginning")
    
    #Process with mass pvalue calculations
    
    simplemasscalc(mytable=sig_meth_TCGA,GOI=sprintf('%s_TET_TCGA_Meth',filestoragelist[i]), wilcox_ids=c(1,4),wilcox_names=c("TvN","Recurrence"),StageCol=3,ScoreCol=2,GroupCol=1)
    
    #Filter to find genes with significantly increased methylation in TvN
    
    print("Filtration beginning")
    
    sig_meth_up_indices <- which(as.numeric(GeneOutput[(nrow(sig_meth_TCGA)+1),])<0.05)
    sig_meth_up <- GeneOutput[,c(1:4,sig_meth_up_indices)]
    
    sig_meth_up[,c(5:ncol(sig_meth_up))] <- lapply(sig_meth_up[,c(5:ncol(sig_meth_up))],as.numeric)
    
    #Significantly increased expression
    directiontest_meth <- aggregate(sig_meth_up[,c(5:ncol(sig_meth_up))],list(sig_meth_up$`Sample Type`),mean)
    directiontest_meth_indices <- which(directiontest_meth[1,]>directiontest_meth[2,])
    sig_meth_up_concordant <- sig_meth_up[,c(1,directiontest_meth_indices+3)]
    
    #Add the co-ordinates back to the sigmeth table, then store the resultant file in the same growing list every time.
    mysig_genelist <- colnames(sig_meth_up_concordant)
    mysig_genelist <- c(mysig_genelist[c(1:4)],gsub(".*\\_","",mysig_genelist[c(5:length(mysig_genelist))]))  
    
    filtered_coords <- as.data.frame(CoordStorage) %>% filter(Composite.Element.REF %in% mysig_genelist)
    
    sig_meth_with_coords <- as.data.frame(rbind(sig_meth_up_concordant,as.character(filtered_coords$Genomic_Coordinate)))
    
    SigMethResultList <- c(SigMethResultList,list(sig_meth_with_coords))
  }
  return(SigMethResultList)
}

gene_name_storage <- gsub("\\_.*","",as.character(colnames(CR1_up_meth_df)))
storagevector <- rep(NA,length(gene_name_storage))

for(i in 2:(length(unique(gene_name_storage)))) {
  print("The gene name being analysed is")
  print(i)
  TCGA_points <- which(gene_name_storage==unique(gene_name_storage)[i])
  TCGA_coords <- as.numeric(CR1_up_meth_df[nrow(CR1_up_meth_df),TCGA_points])
  storage_of_storagevector <- which(gene_name_storage==unique(gene_name_storage)[i])
  CR1_points <- which(CR1_final$Gene==unique(gene_name_storage)[i])
  for (k in 1:length(TCGA_coords)) {
    if (is.na(TCGA_coords[k])) {
      print("This coordinate is NA.")
    } else {
      for (j in 1:length(CR1_points)) {
        #Check to make sure that there is a match for the name
        if (length(CR1_points[j])==0) {
          print("This value was not found.")
        } else if (is.na(CR1_points[j])) {
          print("This value was NA.")
        } else {
          #Check if the gene is oriented backwards or not
          mystart <- as.numeric(as.character((CR1_final[CR1_points[j],2])))
          myend <- as.numeric(as.character((CR1_final[CR1_points[j],3])))
          if (mystart < myend) {
            #The gene is oriented forwards///////////////
            if ((mystart-1500)<(TCGA_coords[k])) {
              if ((TCGA_coords[k])<(myend+1500)) {
                #The coordinate is in-between
                print("TRUE FORWARD")
                storagevector[TCGA_points[k]] <- "Between"
              }
            }
          } else {
            #The gene is oriented backwards
            if ((myend-1500)<(TCGA_coords[k])) {
              if ((TCGA_coords[k])<(myend+1500)) {
                #The co-ordinate is in between
                print("TRUE BACKWARD")
                storagevector[TCGA_points[k]] <- "Between"
              }
            }
          }
        }
      }
    }
  }
  
}

#Place the points located within 500 bp of a TET2 DMR into a separate table and print it.

bp_500_CR1_meth_points <- CR1_meth_sig_TCGA_df[,c(1,which(storagevector=="Between"))]