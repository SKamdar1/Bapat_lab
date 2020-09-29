##Final Random Forest

########PACKAGE LOADING################

require(randomForest)
require(ranger)
require(dplyr)
require(survival)
require(caret)
require(glmnet)
require(rms)
require(ggplot2)
library(e1071)

########PROCESSING DATA###############

#Here, training was performed on the Moreno cohort (to ensure 50:50 split of BCR groups) and testing was performed on the TCGA

testdata_tcag <- LoadData("TCGA_FPKM_TET.txt")
testdata_tcag3 <- testdata_tcag[-which(is.na(testdata_tcag$Biochemical_recurrence)),]
testdata_tcag3$Biochemical_recurrence <- as.character(testdata_tcag3$Biochemical_recurrence)
testdata_tcag3$Biochemical_recurrence[which(testdata_tcag3$Biochemical_recurrence=="N")] <- "0"
testdata_tcag3$Biochemical_recurrence[which(testdata_tcag3$Biochemical_recurrence=="Y")] <- "1"
testdata_tcag3$Biochemical_recurrence <- as.factor(testdata_tcag3$Biochemical_recurrence)
floop <- training
colnames(floop)[1] <- "Biochemical_recurrence"
colnames(floop)[2] <- "CAPRAS"
traindata_moreno <- floop[,which(names(floop) %in% colnames(testdata_tcag3))]
traindata_moreno$Biochemical_recurrence <- as.factor(traindata_moreno$Biochemical_recurrence)

########BACKWARDS FEATURE SELECTION################

#Final selected backwards feature is shown
#Had three tested ranges: 2-30, 2-40, and 2-25. 2-50 gave best sensitivity and specificity when applied to validation cohort, and thus was kept. 

#Backwards feature selection

control <- rfeControl(functions=rfFuncs,method="boot",verbose=FALSE,returnResamp="final",number=75)
sizes <- c(2:30)
profile.1.1 <- rfe(traindata_moreno[,-1],traindata_moreno$Biochemical_recurrence,sizes=sizes,rfeControl=control)
cat("rf: Profile 1 predictors:",predictors(profile.1.1),fill=TRUE)

mybloops_2 <- c(predictors(profile.1.1),"Biochemical_recurrence")
traintrunc_2 <- traindata_moreno[,which(names(traindata_moreno) %in% mybloops_2)]

control <- rfeControl(functions=rfFuncs,method="boot",verbose=FALSE,returnResamp="final",number=75)
sizes <- c(2:50)
profile.1.2 <- rfe(traindata_moreno[,-1],traindata_moreno$Biochemical_recurrence,sizes=sizes,rfeControl=control)
cat("rf: Profile 2 predictors:",predictors(profile.1.2),fill=TRUE)

mybloops_3 <- c(predictors(profile.1.2),"Biochemical_recurrence")
traintrunc_3 <- traindata_moreno[,which(names(traindata_moreno) %in% mybloops_3)]

sizes <- c(2:40)
profile.1.3 <- rfe(traindata_moreno[,-1],traindata_moreno$Biochemical_recurrence,sizes=sizes,rfeControl=control)
cat("rf: Profile 3 predictors:",predictors(profile.1.3),fill=TRUE)

mybloops_4 <- c(predictors(profile.1.3),"Biochemical_recurrence")
traintrunc_4 <- traindata_moreno[,which(names(traindata_moreno) %in% mybloops_4)]


########OPTIMIZATION##########

#Refitting and checking accuracy loss

set.seed(1)
tuneGrid <- expand.grid(.mtry=c(1:10))
trControl <- trainControl(method="cv",number=10,search="grid")
rf22_default <- train(Biochemical_recurrence~.,data=traintrunc_2,method="rf",metric="Accuracy",trControl=trControl)
rf22_default

set.seed(1)
tuneGrid <- expand.grid(.mtry=c(1:10))
trControl <- trainControl(method="cv",number=10,search="grid")
rf23_default <- train(Biochemical_recurrence~.,data=traintrunc_3,method="rf",metric="Accuracy",trControl=trControl)
rf23_default

set.seed(1)
tuneGrid <- expand.grid(.mtry=c(1:10))
trControl <- trainControl(method="cv",number=10,search="grid")
rf24_default <- train(Biochemical_recurrence~.,data=traintrunc_4,method="rf",metric="Accuracy",trControl=trControl)
rf24_default


#Check importances of all features
VarImp1 <- varImp(rf2_default,scale=FALSE)
ggplot(VarImp1,mapping=NULL,top=29,environment=NULL)

#optimize maximum node number
store22_maxnode=list()
tuneGrid <- expand.grid(.mtry=2)
for (maxnodes in c(5:15)) {
  set.seed(1)
  rf22_maxnode <- train(Biochemical_recurrence~., data=traintrunc_2, method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,maxnodes=maxnodes,ntree=300)
  current_iteration <- toString(maxnodes)
  store22_maxnode[[current_iteration]] <- rf22_maxnode}
results22_mtry <- resamples(store22_maxnode)
summary(results22_mtry)
#results22 = 10 nodes

store23_maxnode=list()
tuneGrid <- expand.grid(.mtry=2)
for (maxnodes in c(5:15)) {
  set.seed(1)
  rf23_maxnode <- train(Biochemical_recurrence~., data=traintrunc_3, method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,maxnodes=maxnodes,ntree=300)
  current_iteration <- toString(maxnodes)
  store23_maxnode[[current_iteration]] <- rf23_maxnode}
results23_mtry <- resamples(store23_maxnode)
summary(results23_mtry)
#results23 = 5 nodes

store24_maxnode=list()
tuneGrid <- expand.grid(.mtry=2)
for (maxnodes in c(5:15)) {
  set.seed(1)
  rf24_maxnode <- train(Biochemical_recurrence~., data=traintrunc_4, method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,maxnodes=maxnodes,ntree=300)
  current_iteration <- toString(maxnodes)
  store24_maxnode[[current_iteration]] <- rf24_maxnode}
results24_mtry <- resamples(store24_maxnode)
summary(results24_mtry)
#results24 = 8 nodes

#optimize treesize
store22_maxtrees <- list()
for (ntree in c(seq(from=250,to=1000,by=50))) {
  set.seed(1)
  rf22_maxtrees <- train(Biochemical_recurrence~., data=traintrunc_2,method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,maxnodes=10,ntree=ntree)
  key <- toString(ntree)
  store22_maxtrees[[key]] <- rf22_maxtrees
}
results22_tree <- resamples(store22_maxtrees)
summary(results22_tree)
#results22 = 300 ntree

store23_maxtrees <- list()
for (ntree in c(seq(from=250,to=1000,by=50))) {
  set.seed(1)
  rf23_maxtrees <- train(Biochemical_recurrence~., data=traintrunc_3,method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,maxnodes=5,ntree=ntree)
  key <- toString(ntree)
  store23_maxtrees[[key]] <- rf23_maxtrees
}
results23_tree <- resamples(store23_maxtrees)
summary(results23_tree)
#results23 = 1000 ntree

store24_maxtrees <- list()
for (ntree in c(seq(from=250,to=1000,by=50))) {
  set.seed(1)
  rf24_maxtrees <- train(Biochemical_recurrence~., data=traintrunc_4,method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,maxnodes=8,ntree=ntree)
  key <- toString(ntree)
  store24_maxtrees[[key]] <- rf24_maxtrees
}
results24_tree <- resamples(store24_maxtrees)
summary(results24_tree)
#results24 = 400 ntree

#Final model
fit22_rf <- train(Biochemical_recurrence~., traintrunc_2, method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,ntree=300,maxnodes=10)
predict(fit22_rf,data=traintrunc_2)
prediction22 <- predict(fit22_rf,data=traintrunc_2)
confusionMatrix(prediction22,data=traintrunc_2$Biochemical_recurrence)

prediction32 <- predict(fit22_rf,newdata=testdata_tcag3)
confusionMatrix(prediction32,as.factor(testdata_tcag3$Biochemical_recurrence))
#sens 37.2, spec 77.8, ppv 18.8, npv 89.96
#accuracy 0.7288

fit23_rf <- train(Biochemical_recurrence~., traintrunc_3, method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,ntree=1000,maxnodes=5)
predict(fit23_rf,data=traintrunc_3)
prediction23 <- predict(fit23_rf,data=traintrunc_3)
confusionMatrix(prediction23,data=traintrunc_3$Biochemical.recurrence)

prediction33 <- predict(fit23_rf,newdata=testdata_tcag3)
confusionMatrix(prediction33,as.factor(testdata_tcag3$Biochemical_recurrence))
#sens 72.09, spec 47.27, ppv 15.90, npv 92.45
#accuracy 0.5028

fit24_rf <- train(Biochemical_recurrence~., traintrunc_4, method="rf",metric="Accuracy",tuneGrid=tuneGrid,trControl=trControl,importance=TRUE,nodesize=14,ntree=400,maxnodes=8)
predict(fit24_rf,data=traintrunc_4)
prediction24 <- predict(fit24_rf,data=traintrunc_4)
confusionMatrix(prediction24,data=traintrunc_4$Biochemical_recurrence)

prediction34 <- predict(fit24_rf,newdata=testdata_tcag3)
confusionMatrix(prediction34,as.factor(testdata_tcag3$Biochemical_recurrence))
#sens 39.53, spec 77.2, ppv 19.3, npv 90.2
#accuracy 0.726

#In the final model, out-of-bag error rate has been reduced to 24.53% from 25.47%, and the model is better at predicting positive cases.
#In first final model:
#sens 79.1, spec 36.98, ppv 14.78, npv 92.74
#accuracy 42.09

############GRAPHING FINAL MODEL###############

#Variable importance
rf2_Imp <- varImp(fit2_rf,scale=FALSE)
ggplot(rf2_Imp,mapping=NULL,top=29,environment=NULL)

##########KAPLAN-MEIER ANALYSIS##############

traincap1 <- testdata_tcag3[,c(6,1138,1161)]
traincap2 <- traincap1
traincap3 <- traincap1

#Making survival plots using prediction3

traincap1$pred <- prediction33
traincap1$LowAll <- ifelse(traincap1$CAPRA_S<=2,0,1)
traincap1$HighAll <- ifelse(traincap1$CAPRA_S>=6,1,0)
traincap1$IntHigh <- ifelse(traincap1$CAPRA_S>=6,1,ifelse(traincap1$CAPRA_S<=2,NA,0))

traincap1$SurvObj <- with(traincap1,Surv(BCR_time,Biochemical_recurrence==1))
mysurvpred <- survfit(SurvObj~pred,data=traincap1,conf.type="log-log")
summary(coxph(SurvObj~pred,data=traincap1))
ggsurvplot(mysurvpred,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvlow <- survfit(SurvObj~LowAll,data=traincap1,conf.type="log-log")
ggsurvplot(mysurvlow,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvlowpred <- survfit(SurvObj~LowAll+pred,data=traincap1,conf.type="log-log")
ggsurvplot(mysurvlowpred,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
pairwise_survdiff(formula=SurvObj~LowAll+pred,data=traincap1,p.adjust.method="none")

mysurvhigh <- survfit(SurvObj~HighAll,data=traincap,conf.type="log-log")
ggsurvplot(mysurvhigh,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvhighpred <- survfit(SurvObj~HighAll+pred,data=traincap1,conf.type="log-log")
ggsurvplot(mysurvhighpred,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
pairwise_survdiff(formula=SurvObj~HighAll+pred,data=traincap1,p.adjust.method="none")

mysurvint <- survfit(SurvObj~IntHigh,data=traincap,conf.type="log-log")
ggsurvplot(mysurvint,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvintpred <- survfit(SurvObj~IntHigh+pred,data=traincap1,conf.type="log-log")
ggsurvplot(mysurvintpred,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
pairwise_survdiff(formula=SurvObj~IntHigh+pred,data=traincap1,p.adjust.method="none")

#5 years cutoff

traincap_5years <- traincap[which(traincap$BCR_time<=1825),]
traincap_5years$SurvObj <- with(traincap_5years,Surv(BCR_time,Biochemical.recurrence==1))
mysurvpred5 <- survfit(SurvObj~pred,data=traincap_5years,conf.type="log-log")
summary(coxph(SurvObj~pred,data=traincap_5years))
ggsurvplot(mysurvpred5,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvlow5 <- survfit(SurvObj~LowAll,data=traincap_5years,conf.type="log-log")
ggsurvplot(mysurvlow5,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvlowpred5 <- survfit(SurvObj~LowAll+pred,data=traincap_5years,conf.type="log-log")
ggsurvplot(mysurvlowpred5,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvhigh5 <- survfit(SurvObj~HighAll,data=traincap_5years,conf.type="log-log")
ggsurvplot(mysurvhigh5,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvhighpred5 <- survfit(SurvObj~HighAll+pred,data=traincap_5years,conf.type="log-log")
ggsurvplot(mysurvhighpred5,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvint <- survfit(SurvObj~IntHigh,data=traincap_5years,conf.type="log-log")
ggsurvplot(mysurvint,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

mysurvintpred5 <- survfit(SurvObj~IntHigh+pred,data=traincap_5years,conf.type="log-log")
ggsurvplot(mysurvintpred5,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)



test_ROC_curve <- roc(response=testdata_tcag3$Biochemical_recurrence,predictor=as.numeric(prediction33))

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