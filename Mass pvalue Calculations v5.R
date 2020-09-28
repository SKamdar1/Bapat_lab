#Mass pvalue Calculations.R
#Enabling easy statistical analysis of data (methylation, expression, TCGA, non-TCGA) for GOI
#This script distinguishes between sub-T stages (ie: T2a) and GS7(3+4)/GS7(4+3).
#
#Purpose: Perform the following statistical tests for n probes for a given GOI:
#         Note that the recurrence and stage tests apply to tumor samples ONLY.
#         Mann-Whitney U Test (2 variables):
#         Tumor vs normal samples
#         Biochemical recurrence (Y/N)
#         Surgical margin status (R0 vs R1)
#         Kruskal-Wallis Test (multiple groups):
#         Tumor samples classified by stage
#           1.T2a
#             T2b
#             T2c
#             T3a
#             T3b+
#           2.T2, T3
#           3.T2,T3a,T3b+ 
#             
#         Tumor samples classified by Gleason Score vs normal and each other
#           Gleason scores tested:
#           1.6-, 7a, 7b, 8+, Normal
#           2.6-, 7a, 7b+, Normal
#           3.6-, 7, 8-10, Normal  
#           4.7a-, 7b+, Normal
#             If normal is present, a second round of GS testing with tumors
#             only will also be performed.
#
#Version: v5.0
#Date:    2017-10-19
#Author:  Shivani Kamdar
#
#Input:   .txt file with the following format:
#           Rows are samples (ie: each row is a tumor or normal sample)
#             No need to add rownames, R will automatically do this for you.
#           Columns are CONTINUOUS VALUES (ie probes, PMR values, anything - as long as they are CONTINUOUS)
#           IMPORTANT: Place all grouping columns (non-value columns, such as score, group, stage, recurrence status, etc.)
#                      as the FIRST columns. Otherwise you will run into errors with the way the script handles group
#                      stratification and will probably get false values.
#           The rest of the columns are your probes/PMRs/other values and can be titled anything you like.
#
#           I provide automatic correction for some errors. All blanks and NA
#           values are automatically corrected for. However, to make your life
#           easier, organize the input file appropriately before you run this
#           script.
#Output:  The same .txt file that you put in, with all the stats added to the 
#         bottom for each probe.
#
#How to Use:
#         1. Generate your input .txt file and name it.
#         2. Define "GOI" below as the name of your .txt file (ie: BAP1 is the name of BAP1.txt).
#         3. Specify your stage column(s), Gleason score column, and group column, if applicable.
#         4. Specify vectors of column numbers and output names for Wilcox test.
#         5. Run the script!
#ToDo:    Fix hardcoded GeneTable$Score and GeneTable$Stage specifications.
#Dependencies: dunnTest package
#Notes:   Requires dunnTest to work properly. Should install dependency
#         "dunnTest" automatically if you don't have it. Let me know if there are
#         any issues.
#Changelog: See "Functions" file
# ====================================================================

#Set GOI variable
GOI <- "CandGeneFull"

#Identify your stage column(s), group column, and your Gleason score column

#ClinStageCol <- 1
StageCol <- 3
ScoreCol <- 2
GroupCol <- 1

#Identify the column numbers and desired output names that you want to perform Wilcox tests on.
#Remember that all grouping values in these columns MUST BE BINARY.
#The names don't have to match the column names - they can be set to whatever you like!

#wilcox_ids <- c(1,4,5)
#wilcox_names <- c("TumorvNormal","Recurrence","Margins")
wilcox_ids <- c(4,5,6,7,8,9,10,11,12,13,14)
wilcox_names <- c("PSARecurrence", "Salvage", "BCRRecurrence", "Metastasis","Margins","AdjRT","FailedRP","RecurredBelow1.5Years", "RecurredBelow3Years","RecurredBelow5Years","RecurredBelow7Years")

# ====================================================================
#Source all files and install dependencies

source("Mass pvalue Calc Functions.R")
# ====================================================================
#Tidying the data

#Import dataset
GeneTable <- read.table(file=sprintf("%s.txt", GOI), header = TRUE, sep="\t")

#Remove columns without data
GeneTable <- GeneTable[, colSums(is.na(GeneTable)) !=nrow(GeneTable)]

#Identify columns with binary variables
IncrementCounter <- CallIncrement()

#Remove rows without data
GeneTable <- GeneTable[rowSums(is.na(GeneTable)) !=ncol(GeneTable),]

#Transform non-grouping variables to numeric
for (i in ((IncrementCounter+1):ncol(GeneTable))) {
  GeneTable[,i] <- as.numeric(as.character(GeneTable[,i]))
}

# ====================================================================
#Initializing variables

WilcoxVals <- as.numeric(vector())
KruskalVals <- as.numeric(vector())
NameStorageQuant <- as.character(vector())
NameStorageYN <- as.character(vector())
DunnTable <- list()

# ====================================================================
#Calculate yes/no methylation and quantiles and add to data.

#Quantiles
covariates <- colnames(GeneTable)[c((IncrementCounter+1):ncol(GeneTable))]
Quantiles <- sapply(GeneTable[,c((IncrementCounter+1):ncol(GeneTable))],function(x){quantile(x,na.rm=TRUE)})

for (i in 1:length(covariates)) {
  assign(paste(covariates[i],"Quantile",sep=""),ifelse((GeneTable[,(IncrementCounter+i)])>Quantiles[4,i],1,0))
  GeneTable <- cbind(GeneTable,get(paste(covariates[i],"Quantile",sep="")))
  NameStorageQuant[i] <- paste(covariates[i],"Quantile",sep="")
}

colnames(GeneTable) <- c((colnames(GeneTable[,c(1:IncrementCounter)])),covariates,NameStorageQuant)

#Yes/no
for (i in 1:length(covariates)) {
  assign(paste(covariates[i],"YN",sep=""),ifelse(GeneTable[,i+IncrementCounter]>0,1,0))
  GeneTable <- cbind(GeneTable,get(paste(covariates[i],"YN",sep="")))
  NameStorageYN[i] <- paste(covariates[i],"YN",sep="")
}

colnames(GeneTable) <- c((colnames(GeneTable[,c(1:IncrementCounter)])),covariates,NameStorageQuant,NameStorageYN)

#Initialize GeneOutput
GeneOutput <- GeneTable

# ====================================================================
#Perform Mann-Whitney U Tests

#Check to make sure Wilcox parameters are equal in length

if (length(wilcox_ids) != length(wilcox_names)) {
  stop("Your Wilcox input variables are not equal in length.")
}

for (i in 1:length(wilcox_ids)) {
  wilcox_print_PMR(wilcox_ids[i],wilcox_names[i])
}

# ====================================================================
#Perform Kruskal Tests

#Perform StageTests on StageCol

#Perform clinical stage tests if applicable
if (exists(ClinStageCol)) {
  SimpleStageKruskal(ClinStageCol,"ClinStageSimpleTest")
  SimpleStageKruskal(ClinStageCol,"ClinStageCollapsed",collapse=TRUE)
  GeneTable[,ClinStageCol] <- as.character(GeneTable[,ClinStageCol])
  GeneTable[,ClinStageCol][GeneTable[,ClinStageCol]=='T4'|GeneTable[,ClinStageCol]=='T3b'] <- 'T3b+'
  GeneTable[,ClinStageCol] <- as.factor(GeneTable[,ClinStageCol])
  KruskalCalc(ClinStageCol,"ClinStageAll")
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

#Write the stats out and you're done!
write.table(GeneOutput, file=sprintf("%sProbeStats.txt", GOI), sep="\t", col.names=TRUE)
