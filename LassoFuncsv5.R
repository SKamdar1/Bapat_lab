#LassoFuncs.R
##All functions and dependencies for lasso analysis, with specific classic, threefold cross-additive, and Cox-based scripts currently available
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
#How to Use:
#         Detailed instructions to be generated
#ToDo:    
#Dependencies: caTools, tidyverse, glmnet, ROCR, survival, survminer, pROC
#Notes:   
#Changelog:
#           v5.5 - Clinical functions now work with binomial and Cox-based lasso (binomial was previously not supported)
#           v5.4 - Fixed bug with 0<y<1 error due to R converting binary numerics to 1 and 2 instead of 0 and 1.
#           v5.3 - Added DeLong test function support
#           v5.2 - Added support for using binomial lasso or Cox based on GroupNames
#           v5.1 - Added support for manual input of ModelIDs; fixed hardcoding of TopTenModels usage
#           v5.0 - Added opt.cut.senspec function to generate cutoff points with sensitivity or specificity greater than the given value while 
#                  maximizing the other characteristic.
#           v4.0 - Now fully automated for Cox-based lasso. Need to add automation functions for regular lasso.
# ====================================================================

#Install and load dependencies

if (require("caTools")==FALSE) {
  install.packages("caTools")
}
library(caTools)

if(require("tidyverse")==FALSE) {
  install.packages("tidyverse")
}
library(tidyverse)
library(broom)

if(require("glmnet")==FALSE) {
  install.packages("glmnet")
}
library(glmnet)

if(require("ROCR")==FALSE) {
  install.packages("ROCR")
}
library(ROCR)

if(require("survival")==FALSE) {
  install.packages("survival")
}
library(survival)

if(require("survminer")==FALSE) {
  install.packages("survminer")
}
library(survminer)

if(require("pROC")==FALSE) {
  install.packages("pROC")
}

# ====================================================================
#FUNCTIONS

#Load data and place into data frame
LoadData <- function(filename) {
  GeneTable <- read.table(filename,sep="\t",header=TRUE)
  return(GeneTable)
}

#Split data, if necessary
GetSplitDataFrame <- function(mydata,splitvar,splitratio) {
  set.seed(1)
  sample=sample.split(splitvar,SplitRatio=splitratio)
  train <- subset(mydata,sample==TRUE)
  test <- subset(mydata,sample==FALSE)
  samplelist <- list(train,test)
  return(samplelist)
  #Training can be extracted by selecting samplelist[[1]] and test by selecting samplelist[[2]]
}

#Generate the continuous variable dataset and the grouping variable dataset
StratifyData <- function(mydata,contnames,groupnames) {
  StratData <- data.frame(mydata[,contnames])
  if (length(groupnames)==2) {
      if (length(unique(na.omit(mydata[,groupnames[2]])))==2) {
        group <- data.frame(mydata[,groupnames])
        colnames(group) <- c("time","status")
        #Remove the NA values from the group stratification, and delete the PMRs corresponding to those values in the StratData table.
        if (length(which(is.na(group$time)|is.na(group$status)))>0) {
          grouptrans <- group[-which(is.na(group$time)|is.na(group$status)),]
          StratDataTrans <- StratData[-which(is.na(group$time)|is.na(group$status)),]
        } else {
          StratDataTrans <- StratData
          grouptrans <- group
        }
        #Now we need to pare down the StratData table to only the PMR values, then delete all the rows where there is an NA PMR. We do the same to the group selection column.
        row.x2.na <- apply(StratDataTrans, 1, function(x){any(is.na(x))})
        FinalDat <- StratDataTrans[!row.x2.na,]
        FinalGroup <- grouptrans[!row.x2.na,]
        #Before doing the regression, we need to first make sure everything is numeric.
        FinalDat2 = as.matrix(as.data.frame(lapply(FinalDat, as.numeric)))
        FinalGroup2 <- as.matrix(as.data.frame(lapply(FinalGroup, as.numeric)))
      } else {
        stop("Please make sure that the second variable in groupnames, representing 'status', contains only two factor levels for grouping (ie: 0/1, yes/no, etc)")
      }
    } else {
    group <- (mydata[,groupnames])
    if (length(which(is.na(group)))>0) {
      grouptrans <- group[-which(is.na(group))]
      StratDataTrans <- StratData[-which(is.na(group)),]
    } else {
      StratDataTrans <- StratData
      grouptrans <- group
    }
    #Now we need to pare down the StratData table to only the PMR values, then delete all the rows where there is an NA PMR. We do the same to the group selection column.
    row.x2.na <- apply(StratDataTrans, 1, function(x){any(is.na(x))})
    FinalDat <- StratDataTrans[!row.x2.na,]
    FinalGroup <- grouptrans[!row.x2.na]
    #Before doing the regression, we need to first make sure everything is numeric.
    FinalDat2 = as.matrix(as.data.frame(lapply(FinalDat, as.numeric)))
    FinalGroup2 <- as.numeric(FinalGroup)
  }
  FinalList <- list(FinalDat2,FinalGroup2)
  return(FinalList)
}

#Perform ridge or lasso analysis
LassoRidge <- function(data=FinalDat,group=FinalGroup,id,family,output=TRUE) {
  #Optimal lambda determination using the glmnet package, setting a seed to make sure that results are the same each time
  if (id=="lasso") {
    set.seed(1)
    #Alpha is set to 1 for lasso and 0 for ridge
    cv.fit.lasso <- cv.glmnet(data, group, alpha=1, family=family, type.measure='deviance',  nlambda=100, standardize=FALSE)
    opt_lambda <- cv.fit.lasso$lambda.min
    if (output==TRUE) {
      #Plot and output deviance vs. log(lambda) as a png file
      plot(cv.fit.lasso)
      dev.copy(pdf, 'Lasso - Deviance vs log(lambda).pdf')
      dev.off()
      #Plot and output coefficients vs. log(lambda) as a png file
      plot(cv.fit.lasso$glmnet.fit, xvar="lambda", label=TRUE)
      dev.copy(pdf, 'Lasso - Coefficients vs log(lambda).pdf')
      dev.off()
    }
    Coefs <- as.matrix(coef(cv.fit.lasso, s= opt_lambda))
    return(Coefs)
  } else if (id=="ridge") {
    set.seed(1)
    cv.fit1 <- cv.glmnet(data, group, alpha=0, family=family, type.measure='deviance', nlambda=100, standardize=FALSE)
    if (output==TRUE) {
      #Plot and output deviance vs. log(lambda) as a png file
      plot(cv.fit1)
      dev.copy(pdf, 'Ridge - Deviance vs log(lambda).pdf')
      dev.off()
      opt_lambda <- cv.fit1$lambda.min
      plot(cv.fit1$glmnet.fit, xvar="lambda", label=TRUE)
      dev.copy(pdf, 'Ridge - Coefficients vs log(lambda).pdf')
      dev.off()
    }
    #This output is the most important. Here are all of the coefficients for each gene. It appears as if selection is based on which one has the highest absolute value.
    #We will output this into a table called "Ridge Regression Coefficients".
    #Later, we will combine this with lasso regression results and output as one CSV.
    Coefs <- as.matrix(coef(cv.fit1, s= opt_lambda))
    return(Coefs)
  } else {
    stop("ID must be either 'lasso' or 'ridge'")
  }
}

#Perform Univariate analysis for each gene
UnivariateCoxAnalysis <- function(covariates,data,groupnames) {
  univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(',groupnames,')~',x,sep="")))
  univ_models <- lapply(univ_formulas,function(x){coxph(x,data=data)})
  #Extract data
  univ_results <- lapply(univ_models, function(x){ 
    x <- summary(x)
    p.value<-signif(x$wald["pvalue"], digits=5)
    wald.test<-signif(x$wald["test"], digits=5)
    beta<-signif(x$coef[1], digits=5);#coeficient beta
    HR <-signif(x$coef[2], digits=5);#exp(beta)
    HR.confint.lower <- signif(x$conf.int[,"lower .95"], 5)
    HR.confint.upper <- signif(x$conf.int[,"upper .95"],5)
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-c(beta, HR, wald.test, p.value)
    names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                  "p.value")
    return(res)
    #return(exp(cbind(coef(x),confint(x))))
  })
  univariates <- t(as.data.frame(univ_results,check.names=FALSE))
  unicox1 <- as.data.frame(univariates)
  return(unicox1)
}

#Perform multivariate analysis for each gene
MultivariateCoxAnalysis <- function(data,UVData,groupnames) {
  #Extract those with p<0.05 and put into the multivariate analysis
  sigindex <- which(as.numeric(as.character(UVData$p.value))<0.05)
  signames <- rownames(UVData)[sigindex]
  #Check if there are no significant values from univariate
  if (signames=="") {
    cat("There are no significant genes from univariate.","\n")
  } else {
    MVID <- noquote(paste(signames,collapse="+"))
    #Perform BCR analysis on this
    mv.cox <- coxph(as.formula(paste('Surv(',groupnames,')~', MVID,sep="")), data=data)
    #Put pvalues, confidence intervals, and hazard ratios into a data table
    pval.mvcox <- as.data.frame(coef(summary(mv.cox))[,5])
    confint.mvcox <- as.data.frame(exp(confint(mv.cox)))
    HR.mvcox <- as.data.frame(exp(coef(summary(mv.cox))[,1]))
    mvcoxtable <- cbind(pval.mvcox,HR.mvcox,confint.mvcox)
    colnames(mvcoxtable) <- c("p value", "HR", "2.5%", "97.5%")
    return(mvcoxtable)
  }
}

#Generate desired quantile cutoff
GenerateQuantiles <- function(data,colids,desiredquantileplusone) {
  #Find third quartile values
  Quantiles <- sapply(data[,colids],function(x){quantile(x,na.rm=TRUE)})
  QuantileFrame <- data.frame(matrix(data=rep("NA"),nrow=nrow(data[,colids]),ncol=ncol(Quantiles)))
  colnames(QuantileFrame) <- paste0(colnames(Quantiles),"Quantile")
  for (i in 1:length(colids)) {
    QuantileFrame[,i] <- ifelse(data[,colids][,i]>Quantiles[desiredquantileplusone,i],1,0) 
  }
  QuantileFrame <- cbind(data,QuantileFrame)
  return(QuantileFrame)
}

#Use previously generated coefficients (ie LassoCoefs or RidgeCoefs) to find the genes with non-zero coefficients. Note that for ridge, this will be everything.
IdentifyCoefs <- function(coefficients,ids,family=NULL) {
  if (!is.null(family)) {
    coefficients <- coefficients[-1]
  }
  Coefs <- cbind(as.data.frame(coefficients),ids)
  colnames(Coefs) <- c("Coefs","ID")
  IDStorage <- which(!Coefs$Coefs==0)
  return(IDStorage)
}

#Generate all possible permutations of the provided genes
StorePermutations <- function(genes,namelist) {
  StorageList <- list()
  for (i in 1:length(genes)) {
    ComboList <- combn(genes,i)
    for (j in 1:ncol(ComboList)) {
      StorageList <- c(StorageList,list(assign(paste(namelist[ComboList[,j]],collapse=" "),ComboList[,j])))
    }
  }
  return(StorageList)
}

#Find the optimal cutoff point maximizing sensitivity and specificity, based on a performance and prediction object
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (1-x)+y
    ind = which(d == max(d))[1]
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

#Alternate function to find the optimal cutoff point with sensitivity or specificity set to some minimum value (default=90) while maximizing the other.
opt.cut.senspec = function(perf,pred,sens=TRUE,minval=0.90){
  cutoffs=data.frame(cut=perf@alpha.values[[1]],sens=perf@y.values[[1]],spec=1-(perf@x.values[[1]]))
  if (sens==TRUE) {
    cutoffs <- cutoffs[order(cutoffs$spec,decreasing=TRUE),]
    cutoffchoice <- subset(cutoffs,sens>minval)
  } else {
    cutoffs <- cutoffs[order(cutoffs$sens,decreasing=TRUE),]
    cutoffchoice <- subset(cutoffs,spec>minval)
  }
  return(cutoffchoice[1,])
}

#Perform chi-squared testing and extract PPV/NPV/true positive/false positive/sensitivity/specificity/pvalue/chi-statistic.
#If you run into bugs with needing more than 1 variable: check when you call the function that both of your parameters have more than 1 level.
#You can get a check for this criterion from the AnalyzeClinicalLassoModels function.
chisq.processing <- function(data){
  #FinalLasso <- FinalLasso[-which(is.na(FinalLasso[,columnid])),]
  ObsTable <- chisq.test(data[,ncol(data)],data[,1],correct=FALSE)$observed
  PPV <- ObsTable[4]/(ObsTable[2]+ObsTable[4])
  NPV <- ObsTable[1]/(ObsTable[1]+ObsTable[3])
  TruePos <- ObsTable[4]/(ObsTable[4]+ObsTable[3])
  FalsePos <- ObsTable[2]/(ObsTable[1]+ObsTable[2])
  sens <- TruePos
  specificity <- ObsTable[1]/(ObsTable[1]+ObsTable[2])
  pval.chisq <- chisq.test(data[,ncol(data)],data[,1],correct=FALSE)$p.value
  chistat <- chisq.test(data[,ncol(data)],data[,1],correct=FALSE)$statistic
  chilist <- list(cbind(pval.chisq,PPV,NPV,TruePos,FalsePos,sens,specificity,chistat))
  return(chilist)
}

#Generate a data table containing grouping variable, transformed PMR values multiplied by fit model coefficients, and final values
LassoModeling <- function(genedata,genenames=DataNames,groupnames,genelist=StorageList,id,validation=FALSE,coefvalues=NULL) {
  #Split up the data based on whether it belongs to data variables or your group. Remove all NA values from the group.
  StratData <- genedata[,genenames]
  StratGroup <- genedata[,groupnames]
  StratLasso <- StratData[,id]
  StratLasso[StratLasso==""] <- "NA"
  GroupTrans <- StratGroup[!is.na(StratGroup)]
  print(id)
  if (length(id)>1) {
    #If you have more than 1 gene/clinical variable from the StorageList permutation passed to this function
    #Modifies the data so that there are no NA values in either the group or the data. Makes sure all the data is numeric.
    StratLassoTrans1 <- StratLasso[!is.na(StratGroup),]
    row.x2.na <- apply(StratLassoTrans1,1,function(x) {any(is.na(x))})
    FinalDat <- StratLassoTrans1[!row.x2.na,]
    FinalGroup <- GroupTrans[!row.x2.na]
    FinalDat2 <- as.matrix(as.data.frame(lapply(FinalDat,as.numeric)))
    FinalGroup2 <- as.numeric(FinalGroup)
    if (length(which(FinalGroup2=="2"))>0) {
      FinalGroup2[which(FinalGroup2=="1")] <- "0"
      FinalGroup2[which(FinalGroup2=="2")] <- "1"
      FinalGroup2 <- as.numeric(FinalGroup2)
    }
  } else {
    #If you have only 1 gene/clinical variable from the StorageList permutation passed to this function
    #Modifies the data so that there are no NA values in either the group or the data. Makes sure all the data is numeric.
    StratLassoTrans1 <- StratLasso[!is.na(StratGroup)]
    FinalDat <- StratLassoTrans1[!is.na(StratLassoTrans1)]
    FinalGroup <- GroupTrans[!is.na(StratLassoTrans1)]
    FinalDat2 <- as.numeric(FinalDat)
    FinalGroup2 <- as.numeric(FinalGroup)
    if (length(which(FinalGroup2=="2"))>0) {
      FinalGroup2[which(FinalGroup2=="1")] <- "0"
      FinalGroup2[which(FinalGroup2=="2")] <- "1"
      FinalGroup2 <- as.numeric(FinalGroup2)
    }
  }
  #Perform logistic regression to get the desired model
  #To take time into account: coef.fit.model <- coxph(Surv(FinalGroup1)~FinalDat2)(theoretical)
  if (validation==FALSE) {
    #This is for the training cohort, so make sure to generate coefficients to use.
    coef.fit.model <- glm(FinalGroup2~FinalDat2,family=binomial())
    ModelCoefs <- coef.fit.model$coefficients
    #Mutliply the lasso models by the coefficients, generating new tables or vectors.
    if (length(id)>1) {
      CoefData <- sweep(FinalDat2,2,ModelCoefs[-1],FUN="*")
      LassoSums <- rowSums(CoefData) + ModelCoefs[1]
    } else {
      LassoSums <- FinalDat2*ModelCoefs[-1] + ModelCoefs[1]
      CoefData <- LassoSums
    }
    FinalLasso <- as.data.frame(cbind(FinalGroup2,CoefData,LassoSums))
    LassoList <- list(FinalLasso,ModelCoefs)
  } else {
    #This is for the validation cohort, so we import the previously output coefficient data from the first time we ran this and use it to sweep the validation dataset.
    if (length(id)>1) {
      CoefData <- sweep(FinalDat2,2,(coefvalues[!is.na(coefvalues)][-1]),FUN="*")
      LassoSums <- rowSums(CoefData) + coefvalues[1]
    } else {
      LassoSums <- FinalDat2*coefvalues[!is.na(coefvalues)][-1] + coefvalues[1]
      CoefData <- LassoSums
    }
    FinalLasso <- as.data.frame(cbind(FinalGroup2,CoefData,LassoSums))
    LassoList <- list(FinalLasso,LassoSums)
  }
  #LassoList has two elements: the first being FinalLasso table (including the group, coefficient-transformed values, and the sum totals (LassoSums).
  #                            The second is either the coefficients (for training cohort) or the sum totals (for the validation cohort).
  return(LassoList)
}

#Generate pvalue/OR for fit, AUC/AIC/BIC, and cutoff values and dichotomized data
FitModel <- function(data,datanames=DataNames,modellist=StorageList,id,validation=FALSE,cutoffvalues=NULL,senspec=NULL,senspecvalue=NULL) {
  #Generate fit model via glm
  fit.model <- glm(data[,1]~data$LassoSums,family=binomial())
  #Generate the pvalues and odds ratios of the glm model
  pval.model <- as.data.frame(coef(summary(fit.model))[,4])
  OR.model <- exp(cbind(Odds_Ratio=coef(fit.model),confint(fit.model)))
  #Generate probability, prediction, and performance objects
  prob2 <- predict(fit.model)
  pred2 <- prediction(prob2,data[,1])
  perf2 <- performance(pred2,measure="tpr",x.measure="fpr")
  #Generate AUC curve based on the prediction object, as well as AIC and BIC of the fit model
  auc <- performance(pred2,"auc")
  Model.AUC <- as.numeric(auc@y.values)
  AICModel <- AIC(fit.model)
  BICModel <- BIC(fit.model)
  if (validation==FALSE) {
    #In the training cohort:
    if(senspec=="none") {
      #Generate cutoff values maximizing sensitivity and specificity. Then check if the total sums are above or below the cutoff value
      CutoffValues <- print(opt.cut(perf2,pred2))
      CutoffTraining <- ifelse(data$LassoSums>CutoffValues[3],1,0)
    } else if (senspec=="sens") {
      #Maximizes sensitivity
      #Note that the default senspec value is 0.90
      CutoffValues <- print(opt.cut.senspec(perf2,pred2,minval=senspecvalue))
      CutoffTraining <- ifelse(data$LassoSums>as.numeric(CutoffValues[1]),1,0)
    } else if (senspec=="spec") {
      #Maximizes specificity
      CutoffValues <- print(opt.cut.senspec(perf2,pred2,sens=FALSE,minval=senspecvalue))
      CutoffTraining <- ifelse(data$LassoSums>as.numeric(CutoffValues[1]),1,0)
    } else {
      stop("ERROR: An incorrect value was entered for senspec. Please make sure that it is either none (default) for equal weight to sensitivity and specificity or either
            sens or spec for maximizing sensitivity or specificity above the cutpoint senspecvalue.")
    }
    #Plot AUC curve for this model
    plot(perf2,col=rainbow(7),main="ROC curve",xlab="Specificity",ylab="Sensitivity")
    abline(0,1)
    dev.copy(pdf, file=sprintf("ROC Curve for %s.pdf",paste(datanames[id],collapse=" ")))
    dev.off()
    ModelCutoffList <- list(cbind(pval.model,OR.model),cbind(AICModel,BICModel,Model.AUC),CutoffValues,CutoffTraining)
  } else {
    #In the validation cohort
    #Check if the total sums are above or below the imported cutoff value
    CutoffTraining <- ifelse(data$LassoSums>cutoffvalues,1,0)
    ModelCutoffList <- list(cbind(pval.model,OR.model),cbind(AICModel,BICModel,Model.AUC),cutoffvalues,CutoffTraining)
  }
  #ModelCutofflist has 4 elements:
  #1. the pvalue and OR of the fitted model as a table
  #2. The AIC, BIC, and AUC values as a table
  #3. The cutoffvalues (if training, stores the ones that were generated in this function. If validation, it's the imported values, kept there to ensure that the list has
  #   the same number of elements whether it's training or validation.)
  #4. A 0/1 vector telling you whether each element in LassoSums is greater than (1) or less than (0) the cutoff.
  return(ModelCutoffList)
}

#Generate a nicely formatted prediction table for printing
GeneratePredTable <- function(aucs,chis,cutoffs=NULL){
  #Bind all the AUC and chi-squared result values together into separate data tables
  PredictionTable <- data.frame(do.call(rbind,aucs))
  ChiTable <- data.frame(do.call(rbind,chis))
  if (!is.null(cutoffs)) {
    #This is the training cohort. Every third value from the cutoffs (aka: the actual cutoffs; values 1 and 2 are sensitivity and specificity from opt.cut) is added to the
    #output table.
    CutoffTable <- data.frame(do.call(rbind,cutoffs))
    trip_indices <- seq(3,nrow(CutoffTable),3)
    CutoffTable <- CutoffTable[trip_indices,]
    PredictionTable <- cbind(ChiTable,CutoffTable,PredictionTable)
  } else {
    PredictionTable <- cbind(ChiTable,PredictionTable)
  }
  #PredictionTables are ordered in descending order based on AUC (highest AUC will always be at the top.)
  PredictionTable <- PredictionTable[order(PredictionTable$Model.AUC,decreasing=TRUE),]
  return(PredictionTable)
}

#Generate a nicely formatted prediction table for printing, using either a sens or spec max (no trip_indices)
GenerateMaxPredTable <- function(aucs,chis,cutoffs){
  #Bind all the AUC and chi-squared result values together into separate data tables
  PredictionTable <- data.frame(do.call(rbind,aucs))
  ChiTable <- data.frame(do.call(rbind,chis))
  #No need for cutoff indices here since the cutoffs will be output as a regular data frame
  CutoffTable <- data.frame(do.call(rbind,cutoffs))
  PredictionTable <- cbind(ChiTable,CutoffTable,PredictionTable)
  #PredictionTables are ordered in descending order based on AUC (highest AUC will always be at the top.)
  PredictionTable <- PredictionTable[order(PredictionTable$Model.AUC,decreasing=TRUE),]
  return(PredictionTable)
}

#Analyze lasso models for cutoffs, OR, chi-sq tests, AIC/BIC/AUC
AnalyzeLassoModels <- function(data,groupid,datanames=DataNames,modellist=StorageList,cutofflist=NULL,maxsenspec=MaxCutoffValue,minsenspec=0.90) {
  #Initialize lists to hold output values
  ORList <- list()
  FitAUCList <- list()
  CutoffList <- list()
  CoefficientList <- list()
  ChiResultList <- list()
  #Run for loop to iterate through the StorageList
  if (is.null(cutofflist)) {
    #Since there's no cutofflist, this is the training cohort.
    for (i in 1:length(modellist)) {
      print("This model number is")
      print(i)
      #Perform modeling on n factor combinations (derived from every possible)
      #LassoModeling processes the data so it is all numeric and split into datapoints or grouping points, gets the coefficients of the glm (if training), and returns:
      ##A table of the group, coefficient-transformed values, and the total LassoSums (FinalLasso[1])
      ##The coefficients themselves (training) or LassoSums again (validation) (FinalLasso[2])
      FinalLasso <- LassoModeling(data,groupnames=groupid,id=as.numeric(as.character(unlist(modellist[i]))))
      #FitModel analyzes and returns attributes of the glm for the model:
      ##pvalue and odds ratio (ModelCutoffList[1]); AIC/BIC/AUC (ModelCutoffList[2]); cutoff values (ModelCutoffList[3]), binary vector for threshold cutoff (ModelCutoffList[4])
      #Maxsenspec, if none, will attempt to maximize both sensitivity and specificity when generating a cutoff. Sens will maximize sensitivity; spec specificity.
      ModelCutoffList <- FitModel(as.data.frame(FinalLasso[1]),id=as.numeric(as.character(unlist(modellist[i]))),senspec=maxsenspec,senspecvalue=minsenspec)
      #Binds the group table from FinalLasso to that binary vector of whether values exceed threshold
      CutoffLasso <- cbind(as.data.frame(FinalLasso[1]),unlist(ModelCutoffList[4]))
      #Store the names of the genes that comprise a given model (from the id/ith element of StorageList)
      NameStorage <- paste(datanames[unlist(modellist[i])],collapse=",")
      #For all of the storage lists, append the appropriate values from ModelCutoffList. This will allow them to be stored as we iterate through the loop.
      ORList <- c(ORList,list(ModelCutoffList[[1]]))
      FitAUCList <- c(FitAUCList,list(cbind(as.data.frame(ModelCutoffList[2]),NameStorage)))
      CutoffList <- c(CutoffList,list(ModelCutoffList[[3]]))
      CoefficientList <- c(CoefficientList,list(rbind(as.data.frame(unlist(FinalLasso[2])),NameStorage)))
      #A check to make sure there is are at least two binary values - if this is skipped and you have all 1s or all 0s, the chisq.test will throw an error.
      if (length(unique(ModelCutoffList[[4]]))==1) {
        chilist <- list(c(rep("NA",8)))
        print("Error: This entry is skipped due to all values being either above or below cutoff.")
      } else {
        #Generate PPV, NPV, sens, spec, chi-statistic, p-value, etc. via chi-squared test
        chilist <- chisq.processing(CutoffLasso)
      }
      ChiResultList <- c(ChiResultList,chilist)
    }
    #Since the CutoffList is formatted a little differently if we maximize sens or spec, we use a different function to format the final table containing all of our performance measures
    if (maxsenspec=="none") {
      TrainPredictionTable <- GeneratePredTable(FitAUCList,ChiResultList,CutoffList)
    } else {
      TrainPredictionTable <- GenerateMaxPredTable(FitAUCList,ChiResultList,CutoffList)
    }
    #Make that table part of a list, along with the odds ratios and actual cutoff values (which we may need to analyze models further down the line)
      LoopList <- list(TrainPredictionTable,CoefficientList,ORList,CutoffList)
    } else {
    #We are on the validation cohort
    for (i in 1:length(modellist)) {
      #The output is as above. We set validation to true in the call, and unlist the provided cutofflist exported from LoopList from the training cohort to get cutoffs
      FinalLasso <- LassoModeling(genedata=data,genenames=datanames,groupnames=groupid,id=as.numeric(as.character(unlist(modellist[i]))),validation=TRUE,coefvalues=as.numeric(unlist(cutofflist[2][[1]][i])))
      #Add a check to see if the cutofflist is formatted with maximizing sensitivity or specificity or not
      if (is.null(names(unlist(cutofflist[4][[1]][i]))[3])) {
        #Maximizing both sensitivity and specificity
        ModelCutoffList <- FitModel(as.data.frame(FinalLasso[1]),validation=TRUE,cutoffvalues=as.numeric(unlist(cutofflist[4][[1]][i]))[3])
      } else {
        #Maximizing either sens or spec
        ModelCutoffList <- FitModel(as.data.frame(FinalLasso[1]),validation=TRUE,cutoffvalues=as.numeric(unlist(cutofflist[4][[1]][i]))[1])
      }
      #Same as above - bind the binary vector of whether cutoff is passed to the group data
      CutoffLasso <- cbind(as.data.frame(FinalLasso[1]),unlist(ModelCutoffList[4]))
      #Get the gene names in models
      NameStorage <- paste(datanames[unlist(modellist[i])],collapse=",")
      ORList <- c(ORList,list(ModelCutoffList[[1]]))
      FitAUCList <- c(FitAUCList,list(cbind(as.data.frame(ModelCutoffList[2]),NameStorage)))
      if (length(unique(ModelCutoffList[[4]]))==1) {
        chilist <- list(c(rep("NA",8)))
        print("Error: This entry is skipped due to all values being either above or below cutoff.")
      } else {
      chilist <- chisq.processing(CutoffLasso)
      }
      ChiResultList <- c(ChiResultList,chilist)
      CutoffList <- c(CutoffList,list(ModelCutoffList[[3]]))
    }
      #The TrainPredictionTable and LoopLists are different (no need to output the cutoffs again; it's redundant)
      TrainPredictionTable <- GeneratePredTable(FitAUCList,ChiResultList)
      LoopList <- list(TrainPredictionTable,ORList)
    }
  return(LoopList)
}

#Find the top n best models (max AUC sum in both training and validation cohorts for the selection variable)
ChooseBestModel <- function(n,traintable=as.data.frame(TrainLoop[1]),testtable=as.data.frame(get(paste(GroupingVariable,"ValLoop",sep=""))[1])) {
  #NEED TO MATCH AND THEN ADD!!!!!!!!!!!!!
  traintable <- traintable[order(traintable$NameStorage),]
  testtable <- testtable[order(testtable$NameStorage),]
  sumAUCs <- traintable$Model.AUC+testtable$Model.AUC
  testingdata <- as.data.frame(cbind(as.numeric(as.character(traintable$pval.chisq)),as.numeric(as.character(testtable$pval.chisq)),traintable$Model.AUC,testtable$Model.AUC,sumAUCs,as.character(traintable$NameStorage),rownames(traintable)))
  testingdata <- testingdata[-which(as.numeric(as.character(testingdata[,2]))>0.05),]
  sortedtestingdata <- testingdata[order(-as.numeric(as.character(testingdata$sumAUCs))),]
  butwhymalemodels <- sortedtestingdata[c(1:n),]
  return(butwhymalemodels)
}

#Generate tables with clinicopathological variables and ModelSums added for usage as FinalDats
ClinicalSums <- function(data,groupnames,datanames,modelids=ModelIDs,modellist=StorageList,clinnames=ClinicalNames,coefficienttable=NULL) {
  ClinModelList <- list()
  #If you have any NA values in your model ID list (this is the number of the models, not their names), they will be removed
  if (any(is.na(as.numeric(as.character(modelids))))) {
    modelids <- modelids[-which(is.na(modelids))]
  }
  if (is.null(coefficienttable)) {
    #This is the training cohort
    for (i in 1:length(modelids)) {
      print(i)
      #Overall: Derive the model coefficients from the model IDs (actually easier than just importing from the previous cutoff table, sadly...)
      #Stratify and process the data to remove NA/missing values and separate data from grouping variables
      StratLasso <- StratifyData(data,datanames[unlist(modellist[as.numeric(as.character(modelids[i]))])],groupnames)
      LassoGroup <- as.matrix(as.data.frame(StratLasso[2]))
      LassoDat <- as.matrix(as.data.frame(StratLasso[1]))
      #Derive coefficients from the glm of the given gene model alone
      if (ncol(LassoGroup)==2) {
        if (length(which(LassoGroup[,2]=="2"))>0) {
          LassoGroup[,2][which(LassoGroup[,2]=="1")] <- "0"
          LassoGroup[,2][which(LassoGroup[,2]=="2")] <- "1"
          LassoGroup[,2] <- as.numeric(LassoGroup[,2])
        }
        coef.fit.model <- glm(LassoGroup[,2]~LassoDat,family=binomial())
      } else {
        if (length(which(LassoGroup=="2"))>0) {
          LassoGroup[which(LassoGroup=="1")] <- "0"
          LassoGroup[which(LassoGroup=="2")] <- "1"
          LassoGroup <- as.numeric(LassoGroup)
        }
        coef.fit.model <- glm(LassoGroup~LassoDat,family=binomial())
      }
      ModelCoefs <- coef.fit.model$coefficients
      #Now tidy and process the data again, but now with the clinical variables taken into account
      ClinicalList <- StratifyData(data,c(datanames[unlist(modellist[as.numeric(as.character(modelids[i]))])],clinnames),groupnames)
      ClinDat <- as.matrix(as.data.frame(ClinicalList[[1]]))
      ClinGroup <- as.matrix(as.data.frame(ClinicalList[[2]]))
      #Multiply the gene PMR values (we're only multiplying those!) by their coefficients
      StratClinTrans <- sweep(ClinDat[,c(1:length(datanames[unlist(modellist[as.numeric(as.character(modelids[i]))])]))],2,ModelCoefs[-1],FUN="*")
      #Generate LassoSums for the model, so now those gene values are in one variable
      ModelSums <- rowSums(StratClinTrans) + ModelCoefs[1]
      #Put that in a table with the clinical data variables and groups
      UpdatedTable <- cbind(ModelSums,ClinDat[,ClinicalNames],ClinGroup)
      #Create a list (appended to with each model from TopTenModels) outputting all data with the coefficients as well
      ClinModelList <- c(ClinModelList,list(UpdatedTable),list(ModelCoefs))
    }
  } else {
    #This is the validation cohort
    for (i in 1:length(modelids)) {
      #Process and organize data as above
      ClinicalList <- StratifyData(data,c(datanames[unlist(modellist[as.numeric(as.character(modelids[i]))])],clinnames),groupnames)
      ClinDat <- as.matrix(as.data.frame(ClinicalList[[1]]))
      ClinGroup <- as.matrix(as.data.frame(ClinicalList[[2]]))
      #Coefficients for all components are drawn from training
      ModelCoefs <- as.numeric(coefficienttable[[2]][[i]])
      StratClinTrans <- sweep(ClinDat[,c(1:length(datanames[unlist(modellist[as.numeric(as.character(modelids[i]))])]))],2,ModelCoefs[-1],FUN="*")
      ModelSums <- rowSums(StratClinTrans) + ModelCoefs[1]
      UpdatedTable <- cbind(ModelSums,ClinDat[,ClinicalNames],ClinGroup)
      ClinModelList <- c(ClinModelList,list(UpdatedTable))
    }
  }
  return(ClinModelList)
}

#Clinical version of AnalyzeLassoModels
AnalyzeClinicalLassoModels <- function(clinicaltables=SplitClinList[[1]],datanames=ClinicalNames,cutofflist=NULL,groupid=NULL,maxsenspec=MaxCutoffValue,minsenspec=0.90) {
  #Initialize all holding lists as with AnalyzeLassoModels
  ClinLassoCoefs <- list()
  GiantClinTrainList <- list()
  ClinORList <- list()
  ClinFitAUCList <- list()
  ClinicalCutoffList <- list()
  ClinCoefficientList <- list()
  ClinChiResultList <- list()
  StorageOfClinStorage <- list()
  if (is.null(cutofflist)) {
    #We are in the training cohort
    for (i in 1:length(clinicaltables)) {
      print("i is")
      print(i)
      #We are adding the ModelSums to our ClinicalNames, so we have mylength to extend our "data" by 1
      mylength <- length(datanames) + 1 
      #Grab the table for the i'th model (originally from TopTenModels)
      MakeTable <- as.data.frame(clinicaltables[i])
      #Separate the data from the grouping variables
      ClinData <- as.matrix(as.data.frame(lapply(MakeTable[,c(1:mylength)],as.numeric)))
      ClinGroup <- as.matrix(as.data.frame(lapply(MakeTable[,c((mylength+1):ncol(MakeTable))],as.numeric)))
      #Perform lasso to select variables for the model and store these coefficients in the ClinStorage list
      if (ncol(ClinGroup)==2) {
        if (length(which(ClinGroup[,2]=="2"))>0) {
          ClinGroup[,2][which(ClinGroup[,2]=="1")] <- "0"
          ClinGroup[,2][which(ClinGroup[,2]=="2")] <- "1"
          ClinGroup[,2] <- as.numeric(ClinGroup[,2])
        }
        ClinStorage <- LassoRidge(data=ClinData,group=ClinGroup,id="lasso",family="cox",output=FALSE)
        ClinLassoCoefs <- c(ClinLassoCoefs,list(ClinStorage))
        #Pick out the variables with non-zero coefficients and store all permutations of these variables in ClinStorageList
        ClinCoefGenes <- IdentifyCoefs(ClinLassoCoefs[i],colnames(as.data.frame(ClinData[i])))
      } else {
        if (length(which(ClinGroup=="2"))>0) {
          ClinGroup[which(ClinGroup=="1")] <- "0"
          ClinGroup[which(ClinGroup=="2")] <- "1"
          ClinGroup <- as.numeric(ClinGroup)
        }
        ClinStorage <- LassoRidge(data=ClinData,group=ClinGroup,id="lasso",family="binomial",output=FALSE)
        ClinLassoCoefs <- c(ClinLassoCoefs,list(ClinStorage))
        #Pick out the variables with non-zero coefficients and store all permutations of these variables in ClinStorageList
        coefficients <- ClinLassoCoefs[[i]][-1]
        Coefs <- cbind(as.data.frame(coefficients),colnames(ClinData))
        colnames(Coefs) <- c("Coefs","ID")
        ClinCoefGenes <- which(!Coefs$Coefs==0)
      }
      ClinStorageList <- StorePermutations(ClinCoefGenes,colnames(as.data.frame(ClinData[i])))
      #Initializing these lists again because it was part of a change I made somewhere and afterward the code worked so not touching this
      ClinORList <- list()
      ClinFitAUCList <- list()
      ClinicalCutoffList <- list()
      ClinCoefficientList <- list()
      ClinChiResultList <- list()
      StorageOfClinStorage <- list()
        for (j in 1:length(ClinStorageList)) {
          print("j is")
          print(j)
          #To do: change groupnames to a variable so it's easy to alter
          bloop <<- ClinGroup
          floop <<- as.data.frame(clinicaltables[i])
          if (ncol(as.data.frame(ClinGroup))==2) {
            FinalLasso <- LassoModeling(genedata=as.data.frame(clinicaltables[i]),genenames=c("ModelSums",datanames),groupnames="status",id=as.numeric(as.character(unlist(ClinStorageList[j]))))
          } else {
            FinalLasso <- LassoModeling(genedata=as.data.frame(clinicaltables[i]),genenames=c("ModelSums",datanames),groupnames=colnames(as.data.frame(clinicaltables[i]))[length(datanames)+2],id=as.numeric(as.character(unlist(ClinStorageList[j]))))
          }
          ClinCutoffList <- FitModel(as.data.frame(FinalLasso[1]),id=as.numeric(as.character(unlist(ClinStorageList[j]))),senspec=maxsenspec,senspecvalue=minsenspec)
          ClinCutoffLasso <- cbind(as.data.frame(FinalLasso[1]),id=as.numeric(as.character(unlist(ClinCutoffList[4]))))
          ClinNameStorage <- paste(colnames(ClinData)[unlist(ClinStorageList[j])],collapse=",")
          ClinORList <- c(ClinORList,list(ClinCutoffList[[1]]))
          ClinFitAUCList <- c(ClinFitAUCList,list(cbind(as.data.frame(ClinCutoffList[2]),ClinNameStorage)))
          ClinicalCutoffList <- c(ClinicalCutoffList,list(ClinCutoffList[[3]]))
          ClinCoefficientList <- c(ClinCoefficientList,list(rbind(as.data.frame(unlist(FinalLasso[2])),ClinNameStorage)))
          
          ##################
          if (length(unique(ClinCutoffList[[4]]))==1) {
            ClinChiList <- list(c(rep("NA",8)))
            print("Error: This entry is skipped due to all values being either above or below cutoff.")
          } else {
            ClinChiList <- chisq.processing(ClinCutoffLasso)
          }
          ####################
          #ClinChiList <- chisq.processing(ClinCutoffLasso)
          ClinChiResultList <- c(ClinChiResultList,ClinChiList)
        }
      if (maxsenspec=="none") {
        ClinTrainPredictionTable <- GeneratePredTable(ClinFitAUCList,ClinChiResultList,ClinicalCutoffList)
      } else {
        ClinTrainPredictionTable <- GenerateMaxPredTable(ClinFitAUCList,ClinChiResultList,ClinicalCutoffList)
      }
      LoopList <- list(ClinTrainPredictionTable,ClinCoefficientList,ClinORList,ClinicalCutoffList,ClinStorageList)
      GiantClinTrainList <- c(GiantClinTrainList,LoopList)
    }
  } else {
    #This is the validation cohort
    for (i in 1:length(clinicaltables)) {
      cat("This model is",i)
      mylength <- length(datanames) + 1
      MakeTable <- as.data.frame(clinicaltables[i])
      ClinData <- as.matrix(as.data.frame(lapply(MakeTable[,c(1:mylength)],as.numeric)))
      ClinGroup <- as.numeric(as.data.frame(clinicaltables[[i]])$'ClinicalList[[2]]')
      #
      ClinORList <- list()
      ClinFitAUCList <- list()
      ClinChiResultList <- list()
      #
      for (j in 1:length(cutofflist[[5*i]])) {
        FinalLasso <- LassoModeling(genedata=MakeTable,genenames=c("ModelSums",datanames),groupnames="ClinicalList..2..",id=as.numeric(as.character(unlist(cutofflist[[5*i]][j]))),validation=TRUE,coefvalues=as.numeric(unlist(cutofflist[[2+(5*(i-1))]][[j]])))
        #Adding check for whether a max sens/max spec was used
        if (is.null(names(unlist(cutofflist[[4+(5*(i-1))]][[j]])))) {
          #Maximizing both sensitivity and specificity
          ClinCutoffList <- FitModel(as.data.frame(FinalLasso[1]),validation=TRUE,cutoffvalues=as.numeric(unlist(cutofflist[[4+(5*(i-1))]][[j]]))[3],senspec=maxsenspec)
        } else {
          #Maximizing either sens or spec
          ClinCutoffList <- FitModel(as.data.frame(FinalLasso[1]),validation=TRUE,cutoffvalues=as.numeric(unlist(cutofflist[[4+(5*(i-1))]][[j]]))[1],senspec=maxsenspec)
        }
        ClinCutoffLasso <- cbind(as.data.frame(FinalLasso[1]),id=as.numeric(as.character(unlist(ClinCutoffList[4]))))
        ClinNameStorage <- unlist(cutofflist[[2+(5*(i-1))]][[j]])[[length(unlist(cutofflist[[2+(5*(i-1))]][[j]]))]]
        ClinORList <- c(ClinORList,list(ClinCutoffList[[1]]))
        ClinFitAUCList <- c(ClinFitAUCList,list(cbind(as.data.frame(ClinCutoffList[2]),ClinNameStorage)))
        if (length(unique(ClinCutoffList[[4]]))>1) {
          ClinChiList <- chisq.processing(ClinCutoffLasso)
        } else {
          ClinChiList <- "All values were either above or below cutoff."
          print("Error: This entry is skipped due to all values being either above or below cutoff.")
        }
        ClinChiResultList <- c(ClinChiResultList,ClinChiList)
      }
      ClinValPredictionTable <- GeneratePredTable(ClinFitAUCList,ClinChiResultList)
      LoopList <- list(ClinValPredictionTable,ClinORList)
      GiantClinTrainList <- c(GiantClinTrainList,LoopList)
    }
  }
  #BLOOP
  return(GiantClinTrainList)
}
 
#DeLong test functionality
#Where the first argument is your categorical variable and the second argument is your continuous variable. Third argument should be the name of your table
#How to use (example)
#traindata <- LoadData("mydatafile.txt")
#RunDeLong("Cl_vs_CS_GS","PSA_at_initial_urine_collection","Cl_vs_CS_GS","ProCUrE",traindata)

RunDeLong <- function(catvariable1,contvariable1,catvariable2,contvariable2,dataframe) {
  myroc1 <- roc(dataframe[,catvariable1],dataframe[,contvariable1])
  myroc2 <- roc(dataframe[,catvariable2],dataframe[,contvariable2])
  roc.test(myroc1,myroc2,method="delong")
}