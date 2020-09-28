#Mass pvalue Calc Functions.R
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
#Changelog: v1.1 - Now performs correction for stage and GS Score automatically
#           so you don't have to!
#           v2.0 - Kruskal tests now provide p-values for all scores.
#                  Full support for various sub-pathological stages and primary GS scores added.
#           v3.0 - Functions for IncrementCounter automation and vast simplification of all calculations added.
#           v3.1 - Binary variables undergo chi-squared (unadj) test for Wilcoxon and chi-squared plus Bonferroni-
#                  adjusted comparisons for Kruskal-Wallis/Dunn tests.
#           v3.2 - Quantiles and Y/N values are now added automatically.
#           v3.3 - T3b and T4 are now grouped together.
#                  Gleason score comparisons are now performed once with normals included and once with only
#                  tumors (if applicable).
#                  An additional stage comparison (passing collapse=TRUE to SimpleStageKruskal()) where
#                  T2 is compared to T3a and T3b+ has now been added.
#           v4.0 - Now you can simply set all of your desired test parameters at the start and run! Ease of use version.
#           v5.0 - Fifer package no longer hosted on CRAN. Adjusting script so that it installs and downloads binary plus all necessary dependencies.
#                  Switched all functions to a different file for ease of use.
#                  Script now correctly identifies and removes all N/A variables.
#                  Fixed an error causing chi-squared post hoc tests to not run in certain cases.
# ====================================================================
#Install dependencies

#Fifer is no longer hosted on CRAN so this is a giant pain and requires an installation function

#Function to install multiple packages at once
instpacks <- function(packagenames) {
  for (i in 1:length(packagenames)) {
    if(require(packagenames[i], character.only=TRUE)==FALSE) {
      install.packages(packagenames[i])
      library(packagenames[i])
    }
  }
}

#Define all dependencies required for this script and for Fifer from binary
mypackages <- c("dunn.test","party","Hmisc","plotrix","randomForest","xtable","randomForestSRC","fields","knitr")


#Install Fifer from archived binary if required
if(require("fifer")==FALSE) {
  packageurl <- "http://cran.r-project.org/src/contrib/Archive/fifer/fifer_1.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
  instpacks=mypackages
}
library("fifer")

# ====================================================================
#FUNCTIONS

#The IncrementCounter will track the number of metadata columns preceding the actual data.
#It automatically finds all columns with binary or categorical variables and adds these to IncrementCounter.
#Categorical variables WITH 15 OR LESS VARIABLES are also added automatically.
#WARNING: IF YOUR CATEGORICAL VARIABLE HAS MORE THAN 15 CATEGORIES FOR SOME REASON, YOU MUST ADD IT MANUALLY TO INCREMENTCOUNTER.
CallIncrement <- function() {
  IncrementBool <- sapply((apply(GeneTable,2,unique)),function(x){length(x)<15})
  return(length(which(IncrementBool==TRUE)))
}

#Function: perform Wilcox test for binary groups. If the tested column of GeneTable has binary values (ie quantiles or YN), perform chi-squared test instead.
wilcox_print_PMR <- function(colid,colname) {
  if (nrow(GeneTable[-which(GeneTable[,colid]==""|is.na(GeneTable[,colid])|GeneTable[,colid]=="no info"|GeneTable[,colid]=="N/A"|GeneTable[,colid]=="NULL"|GeneTable[,colid]=="Unknown"),])==0) {
    tablename <- GeneTable
  } else {
    tablename <- GeneTable[-which(GeneTable[,colid]==""|is.na(GeneTable[,colid])|GeneTable[,colid]=="no info"|GeneTable[,colid]=="N/A"|GeneTable[,colid]=="NULL"|GeneTable[,colid]=="Unknown"),]
  }
  for (i in (IncrementCounter+1):ncol(tablename)) {
    if (length(unique((na.omit(GeneTable[,i]))))==2){
      WilcoxVals[i-IncrementCounter] <- chisq.test(GeneTable[,i],GeneTable[,colid],correct=FALSE)$p.value
    } else {
      WilcoxVals[i-IncrementCounter] <- wilcox.test(tablename[,i] ~ tablename[,colid])$p.value
    }
  }
  #Add this data to what will eventually be our output file.
  GeneOutput[(nrow(GeneOutput))+1,] <<- c((rep("NA", IncrementCounter)), WilcoxVals)
  rownames(GeneOutput) <<- c((rownames(GeneOutput[1:((nrow(GeneOutput))-1),])), paste(colname,"pvalue",sep=" "))
}

#Function: perform Kruskal and Dunn tests for nominal groups. If the tested column of GeneTable has binary values (ie quantiles, YN), perform chi-squared tests instead.
KruskalCalc <- function(kruskalid,identifiername) {
  if (nrow(GeneTable[-which(GeneTable[,kruskalid]==""|is.na(GeneTable[,kruskalid])|GeneTable[,kruskalid]=="no info"|GeneTable[,kruskalid]=="N/A"|GeneTable[,kruskalid]=="NULL"|GeneTable[,kruskalid]=="Unknown"),])==0) {
    print("No invalid values found")
  } else {
    GeneTable <- GeneTable[-which(GeneTable[,kruskalid]==""|is.na(GeneTable[,kruskalid])|GeneTable[,kruskalid]=="no info"|GeneTable[,kruskalid]=="N/A"|GeneTable[,kruskalid]=="NULL"|GeneTable[,kruskalid]=="Unknown"),]
    #Remove the empty factors which will otherwise cause errors when doing chi-squared tests
    GeneTable[,kruskalid] <- as.character(GeneTable[,kruskalid])
    GeneTable[,kruskalid] <- as.factor(GeneTable[,kruskalid])
    }
  DunnHold <- list()
  DunnTable <- list()
  for (i in (IncrementCounter+1):ncol(GeneTable)) {
    if (length(na.omit(unique((GeneTable[,i]))))==2){
      KruskalVals[i-IncrementCounter] <- chisq.test(GeneTable[,i],GeneTable[,kruskalid],correct=FALSE)$p.value
      DunnHold <- dunn.test(GeneTable[,i],GeneTable[,kruskalid],method="bonferroni",list=TRUE)
      #DunnTable[[i-IncrementCounter]] <- DunnHold$P.adjusted
      ChiStorage <- table(GeneTable[,i],GeneTable[,kruskalid])
      ChiHold <- chisq.post.hoc(ChiStorage,test=c("chisq.test"),popsInRows=FALSE,control="bonferroni")
      DunnTable[[i-IncrementCounter]] <- ChiHold$adj.p
    } else if (length(na.omit(unique(GeneTable[,i])))==1){
      print("Test cannot be conducted: Only one variable available")
      KruskalVals[i-IncrementCounter] <- "NA"
      DunnTable[[i-IncrementCounter]] <- rep("NA",length(DunnHold))
    } else {
      KruskalVals[i-IncrementCounter] <- kruskal.test(GeneTable[,i],GeneTable[,kruskalid])$p.value
      DunnHold <- dunn.test(GeneTable[,i],GeneTable[,kruskalid],method="bonferroni",list=TRUE)
      DunnTable[[i-IncrementCounter]] <- DunnHold$P.adjusted
    }
  }
  DunnIdentifiers <- DunnHold$comparisons
  GeneOutput[(nrow(GeneOutput))+1,] <<- c((rep("NA",IncrementCounter)), KruskalVals)
  DunnTransferTable <- data.frame(matrix(data=rep("NA"), nrow=length(DunnHold$P.adjusted), ncol=ncol(GeneOutput)))
  for (i in 1:length(KruskalVals)) {
    if (length(unlist(DunnTable[i]))==nrow(DunnTransferTable)) {
      DunnTransferTable[,i+IncrementCounter] <- unlist(DunnTable[i])
    } else {
      print(paste("For column",i+IncrementCounter,"the number of comparisons did not match."),sep=" ")
    }
  }
  colnames(DunnTransferTable) <- colnames(GeneOutput)
  GeneOutput <<- rbind(GeneOutput,DunnTransferTable)
  rownames(GeneOutput) <<- c(rownames(GeneOutput)[1:(nrow(GeneOutput)-(nrow(DunnTransferTable)+1))], paste(identifiername,"pvalue",sep=" "), paste(identifiername, DunnIdentifiers, sep=" "))
}

#Test for number of comparisons and perform the appropriate test.
comparisontest <- function(groupid,identifiername) {
  comparisonnumber <- length(unique(GeneTable[,groupid]))
  if (comparisonnumber>2) {
    KruskalCalc(groupid,identifiername)
  } else if (comparisonnumber==2) {
    wilcox_print_PMR(groupid,identifiername)
  } else {
    print("Error: This comparison has only one variable!")
  }
}

#Kruskal for adjusted Stage
SimpleStageKruskal <- function(groupid,identifiername,collapse=FALSE) {
  GeneHold <- GeneTable
  if (nrow(GeneTable[-which(GeneTable[,groupid]==""|is.na(GeneTable[,groupid])|GeneTable[,groupid]=="no info"|GeneTable[,groupid]=="N/A"|GeneTable[,groupid]=="NULL"),])==0) {
    print("No invalid values found")
  } else {
    GeneTable <<- GeneTable[-which(GeneTable[,groupid]==""|is.na(GeneTable[,groupid])|GeneTable[,groupid]=="no info"|GeneTable[,groupid]=="N/A"|GeneTable[,groupid]=="NULL"),]
  }
  GeneTable[,groupid] <<- as.character(GeneTable[,groupid])
  GeneTable[,groupid][(grep('T1',GeneTable[,groupid]))] <<- 'T1'
  GeneTable[,groupid][(grep('T2',GeneTable[,groupid]))] <<- 'T2'
  if (collapse==FALSE) {
    GeneTable[,groupid][(grep('T3',GeneTable[,groupid]))] <<- 'T3+'
    GeneTable[,groupid][(grep('T4',GeneTable[,groupid]))] <<- 'T3+'
  } else {
    GeneTable[,groupid][(grep('T3b',GeneTable[,groupid]))] <<- 'T3b+'
    GeneTable[,groupid][(grep('T4',GeneTable[,groupid]))] <<- 'T3b+'
  }
  GeneTable[,groupid] <<- as.factor(GeneTable[,groupid])
  comparisonnumber <- length(unique(GeneTable[,groupid]))
  #Do Wilcox if only two groups left and Kruskal if there are more
  if (comparisonnumber>2) {
    KruskalCalc(groupid,identifiername)
  } else if (comparisonnumber==2) {
    wilcox_print_PMR(groupid,identifiername)
  } else {
    print("Error: This comparison has only one variable!")
  }
  GeneTable <<- GeneHold
}

#Gleason score-specific calculations

GleasonCalc <- function(prefix='', groupid) {
  if (nrow(GeneTable[-which(GeneTable[,groupid]==""|is.na(GeneTable[,groupid])|GeneTable[,groupid]=="no info"|GeneTable[,groupid]=="NULL"),])==0) {
    print("No invalid values found")
  } else {
    GeneTable <<- GeneTable[-which(GeneTable[,groupid]==""|is.na(GeneTable[,groupid])|GeneTable[,groupid]=="no info"|GeneTable[,groupid]=="NULL"),]
  }
  
  #Set 1: <=6, 7a, 7b, 8-10, Normal
  GeneHold <- GeneTable
  GeneTable[,groupid] <<- as.character(GeneTable[,groupid])
  GeneTable[,groupid][GeneTable[,groupid] == '8'|GeneTable[,groupid] == '9'|GeneTable[,groupid] == '10'] <<- '8+'
  GeneTable[,groupid][GeneTable[,groupid]=='3'|GeneTable[,groupid] == '4'|GeneTable[,groupid] == '5'|GeneTable[,groupid] == '6'] <<- '<6'
  GeneTable[,groupid] <<- as.factor(GeneTable[,groupid])
  comparisontest(groupid, paste(prefix,"6 v7a v7b v8plus"))
  GeneTable <<- GeneHold
  
  #Set 2: <=6, 7a, 7b-10, Normal
  GeneHold <- GeneTable
  GeneTable[,groupid] <<- as.character(GeneTable[,groupid])
  GeneTable[,groupid][GeneTable[,groupid]=='3'|GeneTable[,groupid] == '4'|GeneTable[,groupid] == '5'|GeneTable[,groupid] == '6'] <<- '<6'
  GeneTable[,groupid][GeneTable[,groupid] == '7 (4+3)'|GeneTable[,groupid] == '7b'|GeneTable[,groupid] == '8'|GeneTable[,groupid] == '9'|GeneTable[,groupid] == '10'] <<- '7b+'
  GeneTable[,groupid] <<- as.factor(GeneTable[,groupid])
  comparisontest(groupid,paste(prefix,"lessthan6 v 7a v greaterthan7b"))
  GeneTable <<- GeneHold
  
  #Set 3: <=6, 7, 8-10, Normal
  GeneHold <- GeneTable
  GeneTable[,groupid] <<- as.character(GeneTable[,groupid])
  GeneTable[,groupid][GeneTable[,groupid]=='3'|GeneTable[,groupid] == '4'|GeneTable[,groupid] == '5'|GeneTable[,groupid] == '6'] <<- '<6'
  GeneTable[,groupid][GeneTable[,groupid] == '7 (3+4)'|GeneTable[,groupid] == '7a'|GeneTable[,groupid] == '7 (4+3)'|GeneTable[,groupid] == '7b'] <<- '7'
  GeneTable[,groupid][GeneTable[,groupid] == '8'|GeneTable[,groupid]=='9'|GeneTable[,groupid]=='10'] <<- '8+'
  GeneTable[,groupid] <<- as.factor(GeneTable[,groupid])
  comparisontest(groupid,paste(prefix,"low v intermediate v high"))
  GeneTable <<- GeneHold
  
  #Set 4: <=7a, >=7b, Normal
  GeneHold <- GeneTable
  GeneTable[,groupid] <<- as.character(GeneTable[,groupid])
  GeneTable[,groupid][GeneTable[,groupid]=='3'|GeneTable[,groupid] == '4'|GeneTable[,groupid] == '5'|GeneTable[,groupid] == '6'|GeneTable[,groupid] == '7 (3+4)'] <<- '<7a'
  GeneTable[,groupid][GeneTable[,groupid] == '7 (4+3)'|GeneTable[,groupid] == '7b'|GeneTable[,groupid] == '8'|GeneTable[,groupid]=='9'|GeneTable[,groupid]=='10'] <<- '>7b'
  GeneTable[,groupid] <<- as.factor(GeneTable[,groupid])
  comparisontest(groupid,paste(prefix,"lessthan7a v greaterthan7b"))
  GeneTable <<- GeneHold
}
