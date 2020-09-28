#CoxLasso.R
##Script for Cox-corrected ridge and lasso regression, as well as univariate/multivariate Cox for continuous and dichotomized variables
#
#Purpose: Perform the following functions for a given file:
#         1. Generate ridge and lasso coefficients based on the chosen selection variable for provided continuous variables (ie: PMR values)
#         2. Perform univariate and multivariate Cox regression on the chosen selection variable for provided continuous variables
#         3. Select all genes with non-zero lasso coefficients and analyze every possible gene model on:
#             a. The selection variable (training cohort)
#             b. All categorical variables available (validation cohort)
#            The output includes: PPV, NPV, chi-squared pvalue, sensitivity/specificity, AIC, BIC, and AUC.
#         4. Take the top gene models and combine them with provided clinicopathological variables to assess additive value performance
#
#VARIATIONS: This file is specifically can be used for either a Cox regression-based or binomial lasso model, depending on the number of arguments provided to GroupNames.
#
#Version: v5.5
#Date:    2019-04-11
#Author:  Shivani Kamdar
#
#Input:   .txt file with the following format:
#           Columns contain both your categorical variables (including selection Status variable), continous Time variable, and all continuous variables to analyze
#Output:  Multiple .csv files:
#         1. Ridge and lasso coefficient file
#         2. Gene selection training file, listing every possible model as described under Purpose.3, with coefficients and odds ratios provided
#         3. Gene selection validation file, listing every possible model as described under Purpose.3
#         4. Similar files to 2. and 3. above, with clinicopathological variables added for additive value
#         Multiple .pdf files:
#         1. ROC files (pdf format) for every possible gene model
#         2. Graphical representations of ridge and lasso coefficients and deviance vs log(lambda) for each respective model (pdf format)
#
#ToDo:  
#Dependencies: caTools, tidyverse, glmnet, ROCR, survival, survminer
#Notes:   
#Changelog: See LassoFuncs
# ====================================================================
#LOAD FUNCTIONS

source("LassoFuncsv5.R")

# ====================================================================
#SCRIPT

#############Initialize name variables####################

#General Note: Make sure any variables specifying column names are the same in both training and validation or you'll get an error.
#DataNames: The names of your variables - for example, gene methylation values. 
#GroupNames: The name(s) of your outcome variable(s). If performing Cox, time first, then yes/no status. Otherwise, this should be just one value. 
#VariableVector: All outcome variables that you wish to test the model performance on (including your binary outcome variable).
#ClinicalNames: Clinical variables that you wish to test in combination with chosen gene models
#ModelIDs: When first run, the script will automatically choose the top 10 models based on the highest combined AUC in train and validation for the binary outcome variable.
#          If you want to specify specific models, replace NULL with the desired model number(s) (ID number[s]) output from the All Training file. (ie: c(65, 29))
#MaxCutoffValue: By default, maximizes both sensitivity and specificity when choosing a cutoff.
#               For max sensitivity: replace "none" with "sens"
#               For max specificity: replace "none" with "spec"

DataNames <- c("APC.PMR","HOXD3.PMR","TGFB2.PMR","GSTP1.PMR","TBX15.PMR","GENX.PMR","HOXD8.PMR","KLF8","GBX2","HIST.1","C1ORF114","HAPLN3")
GroupNames <- c("days.to.PSA.recurrence","PSARecurrence")
VariableVector <- c("BCR","PSARecurrence","Salvage","ADJDEP","ADJ_Salvage","Metastasis","Margins","Rec1_5years","Rec3years","Rec5years","Rec7years","ADJDEP1_5years","ADJDEP3years","ADJDEP5years","ADJDEP7years","IDC_LC")
#ClinicalNames <- c("PSA.at.diagnosis","Age","PSA_density","GleasonGradeGroups","PathStageCont")
ClinicalNames <- c("CAPRA_S","Age","PSA_density","PathStageCont")
ModelIDs <- NULL
MaxCutoffValue <- "sens"

#############Data processing####################

#Load data table(s)
#If you have pregenerated train and test
traindata <- LoadData("CAN.txt")
testdata <- LoadData("UHN.txt")
#If you need to generate your own from some file (example)
#headdata <- LoadData("CAN.txt")
#mylist <- GetSplitDataFrame(headdata,headdata$BCR,0.66)
#traindata <- as.data.frame(mylist[[1]])
#testdata <- as.data.frame(mylist[[2]])

#Removing the adjRT patients from the analysis
#traindata <- traindata[-which(traindata$hadAdjRT==1),]
#testdata <- testdata[-which(testdata$hadAdjRT==1),]
#traindata$BCR_RECURRENCE_DAYS_FROM_RP <- as.numeric(traindata$BCR_RECURRENCE_DAYS_FROM_RP)

#Grab list of stratified continuous data (FinalDat2) and grouping variable (FinalGroup2) and create these variables
HoldingList <- StratifyData(traindata,DataNames,GroupNames)
FinalDat <- HoldingList[[1]]
FinalGroup <- HoldingList[[2]]

#Generate correlation analysis
corr <- cor(FinalDat, use="complete.obs")

#############Ridge and lasso analysis####################

#Generate ridge and lasso coefficients and graphs
if (length(GroupNames)==2) {
  LassoCoefs <- LassoRidge(id="lasso",family="cox")
  RidgeCoefs <- LassoRidge(id="ridge",family="cox")
} else {
  LassoCoefs <- LassoRidge(id="lasso",family="binomial")
  RidgeCoefs <- LassoRidge(id="ridge",family="binomial")
}


########################Cox Analyses###############################################

#############Continuous variable for PMR####################

UVTable <- UnivariateCoxAnalysis(DataNames,traindata,GroupNames)
MVTable <- MultivariateCoxAnalysis(traindata,UVTable,GroupNames)

#############Third quartile yes/no methylation for PMR####################

#Generate third quartile values
QuantileData <- GenerateQuantiles(traindata,DataNames,4)

#Perform UV and MV analysis
UVTableQuant <- UnivariateCoxAnalysis(c(paste0(DataNames,"Quantile")),QuantileData,GroupNames)
MVTableQuant <- MultivariateCoxAnalysis(QuantileData,UVTableQuant,GroupNames)


###########################Fitting the gene models################################

#Automatically picks non-zero coefficient genes for testing

if (length(GroupNames)==2) {
  CoefGenes <- IdentifyCoefs(LassoCoefs,DataNames)
} else {
  CoefGenes <- IdentifyCoefs(LassoCoefs,DataNames,family="binomial")
}

if(length(CoefGenes)==1) {
  StorageList <- list()
  StorageList[[1]] <- CoefGenes
} else {
  StorageList <- StorePermutations(CoefGenes,DataNames)
}

#Run for loop to get pvalues, ORs, AIC/BIC/AUC, chi-sq/PPV/NPV/sens/spec/TruePos/FalsePos, and cutoffs (if training) for each model in StorageList

#Generates single list based on lasso selection variables. Also outputs the ROC curves.
#Check if we are using Cox:
if (length(GroupNames)==2) {
  TrainLoop <- AnalyzeLassoModels(traindata,GroupNames[2])
} else {
  TrainLoop <- AnalyzeLassoModels(traindata,GroupNames)
}

#Generates lists analyzing coefficients derived from this training loop for all possible clinical variables.
for (k in 1:length(VariableVector)) {
assign((paste(VariableVector[k],"TrainLoop",sep="")),AnalyzeLassoModels(traindata,groupid=VariableVector[k],cutofflist=TrainLoop))
}

#########################################Validation cohort analysis#######################################################################

for (k in 1:length(VariableVector)) {
  assign((paste(VariableVector[k],"ValLoop",sep="")),AnalyzeLassoModels(testdata,groupid=VariableVector[k],cutofflist=TrainLoop))
}

#Can analyze as many dichotomized groups as we want as long as they are specified in VariableVector.

#########################################Choosing best gene models#######################################################################

#Choose top ten with highest combined AUCs in training and validation IF modelids=NULL
if (is.null(ModelIDs)) {
  #Add a check for whether this is Cox-based or not to identify the grouping variable
  if (length(GroupNames)==2) {
    GroupingVariable <- GroupNames[2]
  } else {
    GroupingVariable <- GroupNames
  }
  TopTenModels <- ChooseBestModel(10)
  #IdentifiedTrain <- as.data.frame(TrainLoop[1])[as.data.frame(TrainLoop[1])$NameStorage %in% TopTenModels$V6,]
  ModelIDs <- TopTenModels$V7
}

#########################################Adding clinicopathological variables#######################################################################

#Training cohort
#Use stratifydata
ClinDataTables <- ClinicalSums(traindata,GroupNames,DataNames,modelids=ModelIDs)
#If you are doing non-Cox:
#ClinDataTables <- ClinicalSums(traindata,GroupNames,DataNames)
SplitClinList <- split(ClinDataTables,1:2)

ClinicalTrainingLoop <- AnalyzeClinicalLassoModels()
#Loop[1] is the table
#Loop[2] is all the coefficients
#Loop[3] is ORs and confidence intervals
#Loop[4] is sens, spec, and cutoff values
#Loop[5] is the StorageList

#Output all of the outcomes for the training cohort as well
#Generates lists analyzing coefficients derived from this training loop for all possible clinical variables.
for (k in 1:length(VariableVector)) {
  #ClinTrainDataTables <- ClinicalSums(traindata,VariableVector[k],DataNames,coefficienttable=SplitClinList)
  ClinTrainDataTables <- ClinicalSums(traindata,VariableVector[k],DataNames,coefficienttable=SplitClinList,modelids=ModelIDs)
  assign((paste(VariableVector[k],"TrainClinLoop",sep="")),AnalyzeClinicalLassoModels(clinicaltables=ClinTrainDataTables,cutofflist=ClinicalTrainingLoop,groupid=VariableVector[k]))
}

#Validation cohort

for (k in 1:length(VariableVector)) {
  ClinValDataTables <- ClinicalSums(testdata,VariableVector[k],DataNames,coefficienttable=SplitClinList)
  assign((paste(VariableVector[k],"ValClinLoop",sep="")),AnalyzeClinicalLassoModels(clinicaltables=ClinValDataTables,cutofflist=ClinicalTrainingLoop,groupid=VariableVector[k]))
}

#########################################Print data as CSV files#######################################################################

#############Ridge and lasso analysis####################

#Output ridge and lasso results as tables in an Excel file
sink('Cox-Corrected Ridge and Lasso.csv')

cat('Ridge Regression Coefficients')
cat('\n')
write.csv(RidgeCoefs)
cat('_________________________________________________________')
cat('\n')
cat('\n')

cat('Lasso Regression Coefficients')
cat('\n')
write.csv(LassoCoefs)
cat('_________________________________________________________')
cat('\n')
cat('\n')

sink()

########################Cox Analyses###############################################

#Output all of uni/multivariate logistic regression analyses to one csv
sink('Cox Analysis Continuous and Quartile Based.csv')

cat('Univariate Continuous PMR Cox Analysis')
cat('\n')
write.csv(UVTable)
cat('_________________________________________________________')
cat('\n')
cat('\n')

cat('Multivariate Continuous PMR Cox Analyis')
cat('\n')
write.csv(MVTable)
cat('_________________________________________________________')
cat('\n')
cat('\n')

cat('Univariate Third Quartile Based PMR Cox Analysis')
cat('\n')
write.csv(UVTableQuant)
cat('_________________________________________________________')
cat('\n')
cat('\n')

cat('Multivariate Third Quartile Based PMR Cox Analysis')
cat('\n')
write.csv(MVTableQuant)
cat('_________________________________________________________')
cat('\n')
cat('\n')

sink()

########################Training Cohort Analyses###############################################

#Output this data
sink('Gene Model Analysis Training.csv')

cat('Prediction Table for All Models')
cat('\n')
write.csv(as.data.frame(TrainLoop[1]))
cat('_____________________________________')
cat('\n')
cat('\n')

cat('Formula Coefficients for All Models')
cat('\n')
for (i in 1:length(TrainLoop[2][[1]])){
  write.csv(as.data.frame(TrainLoop[2][[1]][[i]]))
  cat('_____________________________________')
  cat('\n')
  cat('\n')
}
cat('_____________________________________')
cat('\n')
cat('\n')

cat('Coefficients and OR for All Models')
cat('\n')
for (i in 1:length(TrainLoop[3][[1]])){
  cat(as.data.frame(TrainLoop[2][[1]][[i]])[nrow(TrainLoop[2][[1]][[i]]),])
  cat('\n')
  write.csv(TrainLoop[3][[1]][[i]])
  cat('_____________________________________')
  cat('\n')
  cat('\n')
}
cat('_____________________________________')
cat('\n')
cat('\n')

sink()

#Training cohort: output analysis with all loops

#Output this data
sink('Gene Model Analysis All Variables Tested Training.csv')

cat('Analysis Table for All Models')
cat('\n')
for (k in 1:length(VariableVector)) {
  cat("Prediction Table -",VariableVector[k],"\n")
  mytable <- get(paste(VariableVector[k],"TrainLoop",sep=""))[1]
  write.csv(mytable)
  cat('_____________________________________')
  cat('\n')
  cat('\n')
}

cat('Coefficients and OR for All Models')
cat('\n')
for (k in 1:length(VariableVector)) {
  mytable2 <- get(paste(VariableVector[k],"TrainLoop",sep=""))[2]
  cat("OR Table -",VariableVector[k],"\n")
  cat('_____________________________________')
  cat('\n')
  for (j in 1:length(mytable2[[1]])) {
    cat(as.data.frame(TrainLoop[2][[1]][[j]])[nrow(TrainLoop[2][[1]][[j]]),])
    cat('\n')
    write.csv(as.data.frame(mytable2[[1]][[j]]))
    cat('\n')
  }
}

sink()

#########################################Validation cohort analysis#######################################################################

#Output this data
sink('Gene Model Analysis Validation.csv')

cat('Analysis Table for All Models')
cat('\n')
for (k in 1:length(VariableVector)) {
  cat("Prediction Table -",VariableVector[k],"\n")
  mytable <- get(paste(VariableVector[k],"ValLoop",sep=""))[1]
  write.csv(mytable)
  cat('_____________________________________')
  cat('\n')
  cat('\n')
}

cat('Coefficients and OR for All Models')
cat('\n')
for (k in 1:length(VariableVector)) {
  mytable2 <- get(paste(VariableVector[k],"ValLoop",sep=""))[2]
  cat("OR Table -",VariableVector[k],"\n")
  cat('_____________________________________')
  cat('\n')
  for (j in 1:length(mytable2[[1]])) {
    cat(as.data.frame(TrainLoop[2][[1]][[j]])[nrow(TrainLoop[2][[1]][[j]]),])
    cat('\n')
    write.csv(as.data.frame(mytable2[[1]][[j]]))
    cat('\n')
  }
}

sink()

########################Training Cohort Analyses (Clinical Variables Added)###############################################

for (i in 1:length(ClinDataTables)) {
  #If TopTen Models is being used, output names from there; if custom ModelIDs, extract names from Storage List
  if (exists("TopTenModels")) {
    sink(sprintf("%s Clin Integration Analysis Training.csv",TopTenModels$V6[i]))
  } else {
    sink(sprintf("%s Clin Integration Analysis Training.csv",paste(DataNames[unlist(StorageList[as.numeric(as.character(ModelIDs[i]))])],collapse=",")))
  }
  
  
  cat('Prediction Table for All Models')
  cat('\n')
  write.csv(as.data.frame(ClinicalTrainingLoop[1+(5*(i-1))]))
  cat('_____________________________________')
  cat('\n')
  cat('\n')
  
  cat('Formula Coefficients for All Models')
  cat('\n')
  for (j in 1:length(ClinicalTrainingLoop[2+(5*(i-1))][[1]])) {
    write.csv(as.data.frame(ClinicalTrainingLoop[2+(5*(i-1))][[1]][[j]]))
    cat('_____________________________________')
    cat('\n')
    cat('\n')
  }
  cat('_____________________________________')
  cat('\n')
  cat('\n')
  
  cat('Coefficients and OR for All Models')
  cat('\n')
  for (j in 1:length(ClinicalTrainingLoop[3+(5*(i-1))][[1]])){
    cat(as.data.frame(ClinicalTrainingLoop[2+(5*(i-1))][[1]][[j]])[nrow(ClinicalTrainingLoop[2+(5*(i-1))][[1]][[j]]),])
    cat('\n')
    write.csv(ClinicalTrainingLoop[3+(5*(i-1))][[1]][[j]])
    cat('_____________________________________')
    cat('\n')
    cat('\n')
  }
  cat('_____________________________________')
  cat('\n')
  cat('\n')
  
  sink()
}

#Output the training cohort data for all clinical variables
sink('Clin Integration Analysis Training All.csv')

for (k in 1:length(VariableVector)) {
  cat("Prediction Table -",VariableVector[k],"\n")
  mylist <- get(paste(VariableVector[k],"TrainClinLoop",sep=""))
  for (y in 1:length(ModelIDs)) {
    if (exists("TopTenModels")) {
      cat(paste(as.character(TopTenModels$V6[j],sep=","),"\n"))
    } else {
      cat(paste(DataNames[unlist(StorageList[as.numeric(as.character(ModelIDs[y]))])],collapse=","),"\n")
    }
    mytable <- mylist[1+(2*(y-1))]
    write.csv(mytable)
    cat('\n')
  }
  cat('_____________________________________')
  cat('\n')
  cat('\n')
}

sink()

########################Validation Cohort Analyses (Clinical Variables Added)###############################################

#Output this data
sink('Clinical Integration Analysis Validation.csv')

for (k in 1:length(VariableVector)) {
  cat("Prediction Table -",VariableVector[k],"\n")
  mylist <- get(paste(VariableVector[k],"ValClinLoop",sep=""))
  for (y in 1:length(ModelIDs)) {
    if (exists("TopTenModels")) {
      cat(paste(as.character(TopTenModels$V6[j],sep=","),"\n"))
    } else {
      cat(paste(DataNames[unlist(StorageList[as.numeric(as.character(ModelIDs[y]))])],collapse=","),"\n")
    }
    mytable <- mylist[1+(2*(y-1))]
    write.csv(mytable)
    cat('\n')
  }
  cat('_____________________________________')
  cat('\n')
  cat('\n')
}

cat('Coefficients and OR for All Models')
cat('\n')
for (k in 1:length(VariableVector)) {
  mylist2 <- get(paste(VariableVector[k],"ValClinLoop",sep=""))
  cat("OR Table -",VariableVector[k],"\n")
  cat('_____________________________________')
  cat('\n')
  for (j in 1:length(ModelIDs)) {
    if (exists("TopTenModels")) {
      cat(paste(as.character(TopTenModels$V6[j],sep=","),"\n"))
    } else {
      cat(paste(DataNames[unlist(StorageList[as.numeric(as.character(ModelIDs[y]))])],collapse=","),"\n")
    }
    mytable2 <- mylist[2+(2*(j-1))]
    for (y in 1:length(mytable2[[1]])) {
      cat(as.data.frame(ClinicalTrainingLoop[2+(5*(j-1))][[1]][[y]])[nrow(ClinicalTrainingLoop[2+(5*(j-1))][[1]][[y]]),])
      cat('\n')
      write.csv(as.data.frame(mytable2[[1]][[y]]))
      cat('\n')
    }
  }
}

sink()