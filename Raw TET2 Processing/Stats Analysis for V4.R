#Additional analysis for Clinical Epigenetics submission Feb 2019

#Import all functions from LassoFuncsv5.R
source("LassoFuncsv5.R")

#Loading and pre-processing of data
NineGenesExp <- read.table("NineGenes.txt",sep="\t",header=TRUE)
DataNames <- c("ASB2.51676","ETNK2.55224","KCNJ15.3772","MEIS2.4212","NRG1.3084","NTN1.9423","NUDT10.170685","PDE4A.5141","SRPX.8406")
GroupNames <- "Sample.Type"
traindata <- NineGeneExp
HoldingList <- StratifyData(traindata,DataNames,GroupNames)
FinalDat <- HoldingList[[1]]
FinalGroup <- HoldingList[[2]]
corr <- cor(FinalDat, use="complete.obs")

#Lasso (just to see what happens)
if (length(GroupNames)==2) {
LassoCoefs <- LassoRidge(id="lasso",family="cox")
RidgeCoefs <- LassoRidge(id="ridge",family="cox")
} else {
LassoCoefs <- LassoRidge(id="lasso",family="binomial")
RidgeCoefs <- LassoRidge(id="lasso",family="binomial")
}


#Univariate logistic regressions
lapply(DataNames,
function(var) {
formula    <- as.formula(paste("Sample.Type ~", var))
res.logist <- glm(formula, data = traindata, family = binomial)
summary(res.logist)
})


#Multivariate logistic regressions

MVID <- noquote(paste(DataNames,collapse="+"))
res.logist <- glm(as.formula(paste("Sample.Type ~",MVID)),data=traindata,family=binomial)
summary(res.logist)


#Calculating AUCs, sensitivities, and specificities 

#ASB2
ASB2_fit <- glm(formula=Dich.Sample~ASB2.51676,data=NineGeneExp,family="binomial")
prob_ASB2=predict(ASB2_fit,type=c("response"))
pred_ASB2 <- prediction(prob_ASB2,NineGeneExp$Dich.Sample)
perf_ASB2 <- performance(pred_ASB2,measure="tpr",x.measure="fpr")
auc_ASB2 <- performance(pred_ASB2,"auc")
Model.AUC_ASB2 <- as.numeric(auc@y.values)
CutoffValues_ASB2 <- print(opt.cut(perf_ASB2,pred_ASB2))

#ETNK2
ETNK2_fit <- glm(formula=Dich.Sample~ETNK2.55224,data=NineGeneExp,family="binomial")
prob_ETNK2=predict(ETNK2_fit,type=c("response"))
pred_ETNK2 <- prediction(prob_ETNK2,NineGeneExp$Dich.Sample)
perf_ETNK2 <- performance(pred_ETNK2,measure="tpr",x.measure="fpr")
auc_ETNK2 <- performance(pred_ETNK2,"auc")
Model.AUC_ETNK2 <- as.numeric(auc@y.values)
CutoffValues_ETNK2 <- print(opt.cut(perf_ETNK2,pred_ETNK2))

#KCNJ15
KCNJ15_fit <- glm(formula=Dich.Sample~KCNJ15.3772,data=NineGeneExp,family="binomial")
prob_KCNJ15=predict(KCNJ15_fit,type=c("response"))
pred_KCNJ15 <- prediction(prob_KCNJ15,NineGeneExp$Dich.Sample)
perf_KCNJ15 <- performance(pred_KCNJ15,measure="tpr",x.measure="fpr")
auc_KCNJ15 <- performance(pred_KCNJ15,"auc")
Model.AUC_KCNJ15 <- as.numeric(auc@y.values)
CutoffValues_KCNJ15 <- print(opt.cut(perf_KCNJ15,pred_KCNJ15))

#MEIS2
MEIS2_fit <- glm(formula=Dich.Sample~MEIS2.4212,data=NineGeneExp,family="binomial")
prob_MEIS2=predict(MEIS2_fit,type=c("response"))
pred_MEIS2 <- prediction(prob_MEIS2,NineGeneExp$Dich.Sample)
perf_MEIS2 <- performance(pred_MEIS2,measure="tpr",x.measure="fpr")
auc_MEIS2 <- performance(pred_MEIS2,"auc")
Model.AUC_MEIS2 <- as.numeric(auc@y.values)
CutoffValues_MEIS2 <- print(opt.cut(perf_MEIS2,pred_MEIS2))

#NRG1
NRG1_fit <- glm(formula=Dich.Sample~NRG1.3084,data=NineGeneExp,family="binomial")
prob_NRG1=predict(NRG1_fit,type=c("response"))
pred_NRG1 <- prediction(prob_NRG1,NineGeneExp$Dich.Sample)
perf_NRG1 <- performance(pred_NRG1,measure="tpr",x.measure="fpr")
auc_NRG1 <- performance(pred_NRG1,"auc")
Model.AUC_NRG1 <- as.numeric(auc@y.values)
CutoffValues_NRG1 <- print(opt.cut(perf_NRG1,pred_NRG1))

#NTN1
NTN1_fit <- glm(formula=Dich.Sample~NTN1.9423,data=NineGeneExp,family="binomial")
prob_NTN1=predict(NTN1_fit,type=c("response"))
pred_NTN1 <- prediction(prob_NTN1,NineGeneExp$Dich.Sample)
perf_NTN1 <- performance(pred_NTN1,measure="tpr",x.measure="fpr")
auc_NTN1 <- performance(pred_NTN1,"auc")
Model.AUC_NTN1 <- as.numeric(auc@y.values)
CutoffValues_NTN1 <- print(opt.cut(perf_NTN1,pred_NTN1))

#NUDT10
NUDT10_fit <- glm(formula=Dich.Sample~NUDT10.170685,data=NineGeneExp,family="binomial")
prob_NUDT10=predict(NUDT10_fit,type=c("response"))
pred_NUDT10 <- prediction(prob_NUDT10,NineGeneExp$Dich.Sample)
perf_NUDT10 <- performance(pred_NUDT10,measure="tpr",x.measure="fpr")
auc_NUDT10 <- performance(pred_NUDT10,"auc")
Model.AUC_NUDT10 <- as.numeric(auc@y.values)
CutoffValues_NUDT10 <- print(opt.cut(perf_NUDT10,pred_NUDT10))

#PDE4A
PDE4A_fit <- glm(formula=Dich.Sample~PDE4A.5141,data=NineGeneExp,family="binomial")
prob_PDE4A=predict(PDE4A_fit,type=c("response"))
pred_PDE4A <- prediction(prob_PDE4A,NineGeneExp$Dich.Sample)
perf_PDE4A <- performance(pred_PDE4A,measure="tpr",x.measure="fpr")
auc_PDE4A <- performance(pred_PDE4A,"auc")
Model.AUC_PDE4A <- as.numeric(auc@y.values)
CutoffValues_PDE4A <- print(opt.cut(perf_PDE4A,pred_PDE4A))

#SRPX
SRPX_fit <- glm(formula=Dich.Sample~SRPX.8406,data=NineGeneExp,family="binomial")
prob_SRPX=predict(SRPX_fit,type=c("response"))
pred_SRPX <- prediction(prob_SRPX,NineGeneExp$Dich.Sample)
perf_SRPX <- performance(pred_SRPX,measure="tpr",x.measure="fpr")
auc_SRPX <- performance(pred_SRPX,"auc")
Model.AUC_SRPX <- as.numeric(auc@y.values)
CutoffValues_SRPX <- print(opt.cut(perf_SRPX,pred_SRPX))

#Using color values from viridis magma, plotting "#000004FF" "#1D1147FF" "#51127CFF" "#822681FF" "#B63679FF" "#E65164FF" "#FB8861FF" "#FEC287FF" "#FCFDBFFF"

library(viridis)

plot(perf_ASB2,col="red")
plot(perf_ETNK2,col="blue",add=TRUE)
plot(perf_KCNJ15,col="orange",add=TRUE)
plot(perf_MEIS2,col="yellow",add=TRUE)
plot(perf_NRG1,col="green",add=TRUE)
plot(perf_NTN1,col="brown",add=TRUE)
plot(perf_NUDT10,col="darkblue",add=TRUE)
plot(perf_PDE4A,col="purple",add=TRUE)
plot(perf_SRPX,col="pink",add=TRUE)
abline(0,1)

##############################################################################################################################

#Decided to play around with these cutoffs (maxsenspec="none") for GS, recurrence, and stage

NineGeneExp$Dich.Sample <- NineGeneExp$Sample.Type
NineGeneExp$Dich.Sample <- as.character(NineGeneExp$Dich.Sample)
NineGeneExp$Dich.Sample[which(NineGeneExp$Dich.Sample=="Primary solid Tumor")] <- "1"
NineGeneExp$Dich.Sample[which(NineGeneExp$Dich.Sample=="Solid Tissue Normal")] <- "0"
NineGeneExp$Dich.Sample <- as.numeric(as.character(NineGeneExp$Dich.Sample))

NineGeneExp$Biochemical.recurrence <- as.character(NineGeneExp$Biochemical.recurrence)
NineGeneExp$Biochemical.recurrence[which(NineGeneExp$Biochemical.recurrence=="Unknown")] <- ""
NineGeneExp$Biochemical.recurrence <- as.factor(NineGeneExp$Biochemical.recurrence)

#Ran Lasso.

##############################################################################################################################

#Now trying a 2/3, 1/3 approach so we have a validation set.

mylist <- GetSplitDataFrame(NineGeneExp,NineGeneExp$Dich.Sample,0.66)
traindata <- as.data.frame(mylist[[1]])
testdata <- as.data.frame(mylist[[2]])

#Ran LassoFuncs 

#Dichotomize data based on cutoff values

traindata$ASB2_Dich <- ifelse(traindata$ASB2.51676<130.1889,1,0)
traindata$ETNK2_Dich <- ifelse(traindata$ETNK2.55224<161.7002,1,0)
traindata$KCNJ15_Dich <- ifelse(traindata$KCNJ15.3772<14.1392,1,0)
traindata$MEIS2_Dich <- ifelse(traindata$MEIS2.4212<373.0769,1,0)
traindata$NRG1_Dich <- ifelse(traindata$NRG1.3084<4.825,1,0)
traindata$NTN1_Dich <- ifelse(traindata$NTN1.9423<151.1045,1,0)
traindata$NUDT10_Dich <- ifelse(traindata$NUDT10.170685<175.4788,1,0)
traindata$PDE4A_Dich <- ifelse(traindata$PDE4A.5141<192.0698,1,0)
traindata$SRPX_Dich <- ifelse(traindata$SRPX.8406<42.936,1,0)

#Do Kaplan-Meier survival curves for high (0) vs low (1) expression for each gene in training and test
#First we need to update the training and test data so it includes days to recurrence

traindata <- LoadData("traindata.txt")
testdata <- LoadData("testdata.txt")

#Create SurvOBj

traindata$SurvObj <- with(traindata, Surv(Time_to_BCR,Biochemical.recurrence==1))
testdata$SurvObj <- with(testdata, Surv(Time_to_BCR,Biochemical.recurrence==1))

mysurvfit <- survfit(SurvObj~ASB2_Dich,data=testdata,conf.type="log-log")
ggsurvplot(mysurvfit)
summary(coxph(SurvObj~ASB2_Dich,data=traindata))


#Alas, all data using the cutpoints from TvN sucks. Going to try X-tile analysis.

write.table(traindata,"train.txt",sep="\t",col.names=TRUE)
write.table(testdata,"test.txt",sep="\t",col.names=TRUE)

#Testdata has only ten events and this may be skewing the analysis. Trying to redo with a new test and train set.
#Also must remember to remove normal cases before porting into X-tile.

Tumor9 <- NineGeneExp2[which(NineGeneExp2$Sample.Type=="Primary solid Tumor"),]
Normal9 <- NineGeneExp2[which(NineGeneExp2$Sample.Type=="Solid Tissue Normal"),]
Tumor9$Biochemical.recurrence <- as.character(Tumor9$Biochemical.recurrence)
Tumor9$Biochemical.recurrence[which(Tumor9$Biochemical.recurrence=="Unknown")] <- ""
Tumor9$Biochemical.recurrence[which(Tumor9$Biochemical.recurrence=="N")] <- "0"
Tumor9$Biochemical.recurrence[which(Tumor9$Biochemical.recurrence=="Y")] <- "1"
Tumor9$Biochemical.recurrence <- as.numeric(as.character(Tumor9$Biochemical.recurrence))

Normal9$Biochemical.recurrence <- as.character(Normal9$Biochemical.recurrence)
Normal9$Biochemical.recurrence <- rep("",35)
Normal9$Biochemical.recurrence <- as.numeric(Normal9$Biochemical.recurrence)

mylist <- GetSplitDataFrame(Tumor9,Tumor9$Biochemical.recurrence,0.66)
traindata <- as.data.frame(mylist[[1]])
testdata <- as.data.frame(mylist[[2]])

#Randomly grab some normal samples and put 2/3 in train and 1/3 in test
set.seed(1)
trainnorms <- sample.int(35,23)

trainnormdata <- Normal9[trainnorms,]
testnormdata <- Normal9[-trainnorms,]

traindata <- as.data.frame(rbind(traindata,trainnormdata))
testdata <- as.data.frame(rbind(testdata,testnormdata))

testdata$Dich.Sample <- testdata$Sample.Type
testdata$Dich.Sample <- as.character(testdata$Dich.Sample)
testdata$Dich.Sample[which(testdata$Dich.Sample=="Primary solid Tumor")] <- "1"
testdata$Dich.Sample[which(testdata$Dich.Sample=="Solid Tissue Normal")] <- "0"
testdata$Dich.Sample <- as.numeric(as.character(testdata$Dich.Sample))

#Reran the AUCs and everything.
#Now going to try adjusting for X-tile export

traindata$Category <- rep("training",258)
testdata$Category <- rep("test",200)

Xtile <- as.data.frame(rbind(traindata,testdata))

write.table(Xtile,"TCGA_Data.txt",sep="\t",col.names=TRUE)
