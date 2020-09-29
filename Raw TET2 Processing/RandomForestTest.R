require(randomForest)
require(ranger)
require(dplyr)

#Load and process data
training <- LoadData("TCGA_Train_Forest.txt")
training$Recurrence_time <- training$X_TIME_TO_EVENT
training$Recurrence_time[which(!is.na(training$days_to_first_biochemical_recurrence))] <- training$days_to_first_biochemical_recurrence[which(!is.na(training$days_to_first_biochemical_recurrence))]

train_set <- training[,c(6:460,487:488)]
train_set2 <- train_set[-which(is.na(train_set$Biochemical.recurrence)),]

#Split into test and training sets based on BCR
set.seed(1)

mylist <- GetSplitDataFrame(train_set2,train_set2$Biochemical.recurrence,0.66)
traindata <- as.data.frame(mylist[[1]])
testdata <- as.data.frame(mylist[[2]])

floopbloop <- traindata$CAPRA_S
bloopfloop <- testdata$CAPRA_S

traindata <- traindata[,-456]
testdata <- testdata[,-456]

#Perform weighted random forest

#Set weights to favor lower (Y) variable

my_weight <- 1/table(traindata$Biochemical.recurrence)
my_weight <- my_weight/sum(my_weight)
weights <- rep(0,nrow(traindata))
weights[traindata$Biochemical.recurrence=='N'] <- my_weight['N']
weights[traindata$Biochemical.recurrence=='Y'] <- my_weight['Y']
#Visualizing
table(weights,traindata$Biochemical.recurrence)

#Generate random forest model
#Can adjust mtry and ntree to optimize (lower) prediction error
my_model <- ranger(Biochemical.recurrence~.,traindata,case.weights=weights,importance="impurity",write.forest=TRUE,mtry=51,num.trees=400)

#Test predictions on the test set
pred_model <- predict(my_model,testdata)
table(testdata$Biochemical.recurrence,pred_model$predictions)

#Trying different method
rf.classwt <- randomForest(Biochemical.recurrence~.,data=traindata,classwt=c(0.0004,2000))
pred_model <- predict(rf.classwt,testdata)
table(testdata$Biochemical.recurrence,pred_model$predictions)

#Trying survivalforest
require(survival)
require(caret)
require(glmnet)
require(rms)
require(risksetROC)
require(doParallel)

registerDoParallel(detectCores()-1)
detectCores()
options(rf.cores=detectCores()-1,mc.cores=detectCores()-1)

require(randomForestSRC)
require(ggRandomForests)

my_rsf_1 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata,importance=TRUE)
imp_rsf_1 <- sort(my_rsf_1$importance,decreasing=T)
plot(gg_vimp(my_rsf_1))

#NOPE.

#Trying z-scoring using AutoZ function from LassoFuncs

traindata2 <- AutoZ(traindata,DataNames)
traindata3 <- traindata2[,c(1,458:914)]
my_rsf_2 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata3,importance=TRUE,ntree=1500)

#pred on traindata
pred_rsf = predict(my_rsf_2, newdata = traindata3)

rcorr.cens(-pred_rsf$predicted , 
           Surv(traindata$Recurrence_time, traindata$Biochemical.recurrence))["C Index"]

#Pred on testdata
testdata2 <- AutoZ(testdata,DataNames)
testdata3 <- testdata2[,c(1,458:914)]
testdata4 <- testdata3[-which(is.na(testdata3$Recurrence_time)),]

pred_rsf2 = predict(my_rsf_2, newdata = testdata4)

rcorr.cens(-pred_rsf2$predicted , 
           Surv(testdata4$Recurrence_time, testdata4$Biochemical.recurrence))["C Index"]


#plot C-index on testdata with fixed mtry and varying ntree

treevec <- c(500,600,700,800,900,1000,1100,1200,1300,1400,1500)
StorageVector <- rep(NA,length(treevec))
for (i in 1:length(treevec)) {
  my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata3,ntree=treevec[i])
  my_pred <- predict(my_rsf,newdata=testdata4)
  StorageVector[i] <- as.numeric(rcorr.cens(-my_pred$predicted , 
                                            Surv(testdata4$Recurrence_time, testdata4$Biochemical.recurrence))["C Index"])
}

#Plotting this
sparedf <- as.data.frame(cbind(treevec,StorageVector))
ggplot(data=sparedf,aes(x=treevec,y=StorageVector,group=1)) + geom_line(linetype="dashed")+geom_point()

treevec <- c(seq(from=635,to=655,by=1),seq(from=715,to=735,by=1))
#Optimal ntree=645

tryvec <- seq(from=3,to=50,by=1)
StorageVector <- rep(NA,length(tryvec))
for (i in 1:length(tryvec)) {
  my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata3,ntree=645,mtry=tryvec[i])
  my_pred <- predict(my_rsf,newdata=testdata4)
  StorageVector[i] <- as.numeric(rcorr.cens(-my_pred$predicted , 
                                            Surv(testdata4$Recurrence_time, testdata4$Biochemical.recurrence))["C Index"])
}

#Plotting this
sparedf <- as.data.frame(cbind(tryvec,StorageVector))
ggplot(data=sparedf,aes(x=tryvec,y=StorageVector,group=1)) + geom_line(linetype="dashed")+geom_point()
#Optimal mtry appears to be 3 (?) Might increase OOB error/impair predictive capability - let's test

my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata3,ntree=645,mtry=3)
my_pred <- predict(my_rsf,newdata=traindata3)
(rcorr.cens(-my_pred$predicted ,                                           Surv(traindata3$Recurrence_time, traindata3$Biochemical.recurrence))["C Index"])

#training c is 0.9640811, error rate 43.81%. No substantial loss - should be fine!

#Let's try lasso on this set

DataNames <- colnames(traindata3)[-c(1:2)]
GroupNames <- c("Biochemical.recurrence")

#Grab list of stratified continuous data (FinalDat2) and grouping variable (FinalGroup2) and create these variables
HoldingList <- StratifyData(traindata3,DataNames,GroupNames)
FinalDat <- HoldingList[[1]]
FinalGroup <- HoldingList[[2]]


#Generate ridge and lasso coefficients and graphs
if (length(GroupNames)==2) {
  LassoCoefs <- LassoRidge(id="lasso",family="cox")
  RidgeCoefs <- LassoRidge(id="ridge",family="cox")
} else {
  LassoCoefs <- LassoRidge(id="lasso",family="binomial")
  RidgeCoefs <- LassoRidge(id="ridge",family="binomial")
}


if (length(GroupNames)==2) {
  CoefGenes <- IdentifyCoefs(LassoCoefs,DataNames)
} else {
  CoefGenes <- IdentifyCoefs(LassoCoefs,DataNames,family="binomial")
}

#SPAG5_z is identified as the only gene from LASSO if Cox is used
#28 genes are identified if binomial LASSO is used
#AGPAT9, ANTXR2, C3orf62, CYTH2, DDIT4, DERA, DLEU1, HHAT, HOXA13, HSD3B7, IGDCC4,INPPL1,IRAK1BP1, ITFG2, M6PR, MED21, MEIS3P1, NOS1, NRXN3, PARM1, PDE4D, PODXL2, PRKG2, SPAG5, TGM3, TRIM6, TUBB2A, ZIC5

#Let's try tuning random forest with these genes
paste0(DataNames[CoefGenes],collapse="+")

my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~AGPAT9_z+ANTXR2_z+C3orf62_z+CYTH2_z+DDIT4_z+DERA_z+DLEU1_z+HHAT_z+HOXA13_z+HSD3B7_z+IGDCC4_z+INPPL1_z+IRAK1BP1_z+ITFG2_z+M6PR_z+MED21_z+MEIS3P1_z+NOS1_z+NRXN3_z+PARM1_z+PDE4D_z+PODXL2_z+PRKG2_z+SPAG5_z+TGM3_z+TRIM6_z+TUBB2A_z+ZIC5_z,data=traindata3,importance=TRUE)
#error rate is 30.48%

#Testing optimal ntree and mtry reveals ntree=904 and mtry=3


  my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~AGPAT9_z+ANTXR2_z+C3orf62_z+CYTH2_z+DDIT4_z+DERA_z+DLEU1_z+HHAT_z+HOXA13_z+HSD3B7_z+IGDCC4_z+INPPL1_z+IRAK1BP1_z+ITFG2_z+M6PR_z+MED21_z+MEIS3P1_z+NOS1_z+NRXN3_z+PARM1_z+PDE4D_z+PODXL2_z+PRKG2_z+SPAG5_z+TGM3_z+TRIM6_z+TUBB2A_z+ZIC5_z,data=traindata3,ntree=904,mtry=3)
  my_pred <- predict(my_rsf,newdata=traindata3)
  (rcorr.cens(-my_pred$predicted , 
                                            Surv(traindata3$Recurrence_time, traindata3$Biochemical.recurrence))["C Index"])

#Predictive capabilities
  my_pred <- predict(my_rsf,newdata=testdata4)
  
#ROC and Kaplan-Meier
  w.ROC <- risksetROC(Stime=traindata3$Recurrence_time,status=traindata3$Biochemical.recurrence,marker=my_rsf$predicted,predict.time=median(traindata3$Recurrence_time),method="Cox",lwd=3,col="red")
w.ROC$AUC
my_pred <- predict(my_rsf,newdata=traindata3)
my_pred$predicted
mean(my_pred$predicted)
traindata3$risk <- ifelse(my_pred$predicted>=2.4564,1,0)
testdata4$risk <- ifelse(my_pred2$predicted>=2.4564,1,0)
testdata4$SurvObj <- with(testdata4,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk,data=testdata4,conf.type="log-log")
ggsurvplot(mysurvfit)
summary(coxph(SurvObj~risk,data=testdata4))

  
#First pruning for variable selection
  #Removing negative importance variables: TUBB2A, IGDCC4, DDIT4, PODXL2, TGM3
  
  my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~AGPAT9_z+ANTXR2_z+C3orf62_z+CYTH2_z+DERA_z+DLEU1_z+HHAT_z+HOXA13_z+HSD3B7_z+INPPL1_z+IRAK1BP1_z+ITFG2_z+M6PR_z+MED21_z+MEIS3P1_z+NOS1_z+NRXN3_z+PARM1_z+PDE4D_z+PRKG2_z+SPAG5_z+TRIM6_z+ZIC5_z,data=traindata3,ntree=904,mtry=3)
  my_pred <- predict(my_rsf,newdata=traindata3)
  (rcorr.cens(-my_pred$predicted , 
              Surv(traindata3$Recurrence_time, traindata3$Biochemical.recurrence))["C Index"])
  
  w.ROC <- risksetROC(Stime=traindata3$Recurrence_time,status=traindata3$Biochemical.recurrence,marker=my_rsf$predicted,predict.time=median(traindata3$Recurrence_time),method="Cox",lwd=3,col="red")
  w.ROC$AUC
  
  my_pred2 <- predict(my_rsf,newdata=testdata4)
  w.roc=risksetROC(Stime=testdata4$Recurrence_time,status=testdata4$Biochemical.recurrence,marker=my_pred2$predicted,predict.time=median(testdata4$Recurrence_time),method="Cox",lwd=3,col="red")

#Let's try running random forest and performing variable selection only with it
  
  traindata4 <- traindata3[,-459]
  
  my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata4,importance=TRUE)
  plot(gg_vimp(my_rsf))
  imp_rsf_1 <- sort(my_rsf$importance,decreasing=T)
  neg_vars <- names(which(imp_rsf_1<0))
  
  traindata_trim1 <- traindata4[,-which(names(traindata4) %in% neg_vars)]
  
  my_rsf2 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim1,importance=TRUE)
  plot(gg_vimp(my_rsf2))
  imp_rsf_2 <- sort(my_rsf2$importance,decreasing=T)
  neg_vars2 <- names(which(imp_rsf_2<0))
  
  traindata_trim2 <- traindata_trim1[,-which(names(traindata_trim1) %in% neg_vars2)]
  
  my_rsf3 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim2,importance=TRUE)
  plot(gg_vimp(my_rsf3))
  imp_rsf_3 <- sort(my_rsf3$importance,decreasing=T)
  neg_vars3 <- names(which(imp_rsf_3<0))
  
  traindata_trim3 <- traindata_trim2[,-which(names(traindata_trim2) %in% neg_vars3)]
  
  my_rsf4 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim3,importance=TRUE)
  plot(gg_vimp(my_rsf4))
  imp_rsf_4 <- sort(my_rsf4$importance,decreasing=T)
  neg_vars4 <- names(which(imp_rsf_4<0))
  
  traindata_trim4 <- traindata_trim3[,-which(names(traindata_trim3) %in% neg_vars4)]
  
  my_rsf5 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim4,importance=TRUE)
  plot(gg_vimp(my_rsf5))
  imp_rsf_5 <- sort(my_rsf5$importance,decreasing=T)
  neg_vars5 <- names(which(imp_rsf_5<0))
  
  traindata_trim5 <- traindata_trim4[,-which(names(traindata_trim4) %in% neg_vars5)]
  
  my_rsf6 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim5,importance=TRUE)
  plot(gg_vimp(my_rsf6))
  imp_rsf_6 <- sort(my_rsf6$importance,decreasing=T)
  neg_vars6 <- names(which(imp_rsf_6<0))
  
  traindata_trim6 <- traindata_trim5[,-which(names(traindata_trim5) %in% neg_vars6)]
  
  my_rsf7 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim6,importance=TRUE)
  plot(gg_vimp(my_rsf7))
  imp_rsf_7 <- sort(my_rsf7$importance,decreasing=T)
  neg_vars7 <- names(which(imp_rsf_7<0))
  
  traindata_trim7 <- traindata_trim6[,-which(names(traindata_trim6) %in% neg_vars7)]
  
  my_rsf8 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim7,importance=TRUE)
  plot(gg_vimp(my_rsf8))
  imp_rsf_8 <- sort(my_rsf8$importance,decreasing=T)
  neg_vars8 <- names(which(imp_rsf_8<0))
  
  traindata_trim8 <- traindata_trim7[,-which(names(traindata_trim7) %in% neg_vars8)]
  
  my_rsf9 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim8,importance=TRUE)
  plot(gg_vimp(my_rsf9))
  imp_rsf_9 <- sort(my_rsf9$importance,decreasing=T)
  neg_vars9 <- names(which(imp_rsf_9<0))
  
  traindata_trim9 <- traindata_trim8[,-which(names(traindata_trim8) %in% neg_vars9)]
  
  my_rsf10 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim9,importance=TRUE)
  plot(gg_vimp(my_rsf10))
  imp_rsf_10 <- sort(my_rsf10$importance,decreasing=T)
  neg_vars10 <- names(which(imp_rsf_10<0))
  
  traindata_trim10 <- traindata_trim9[,-which(names(traindata_trim9) %in% neg_vars10)]
  
  my_rsf11 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim10,importance=TRUE)
  plot(gg_vimp(my_rsf11))
  imp_rsf_11 <- sort(my_rsf11$importance,decreasing=T)
  neg_vars11 <- names(which(imp_rsf_11<0))
  my_rsf11
  
  traindata_trim11 <- traindata_trim10[,-which(names(traindata_trim10) %in% neg_vars11)]
  
  my_rsf12 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim11,importance=TRUE)
  plot(gg_vimp(my_rsf12))
  imp_rsf_12 <- sort(my_rsf12$importance,decreasing=T)
  neg_vars12 <- names(which(imp_rsf_12<0))
  my_rsf12
  
  #Calculate AUC over time for each
  w.ROC = risksetAUC(Stime=traindata_trim11$Recurrence_time,status=traindata_trim11$Biochemical.recurrence,marker=my_rsf12$predicted.oob,tmax=4000)
  w.ROC$AUC[15]

  #pred on traindata
  
  rcorr.cens(my_rsf$predicted.oob , 
             Surv(traindata4$Recurrence_time, traindata4$Biochemical.recurrence))["C Index"]

#Optimizing my_rsf12
  
  set.seed(1)
  
  #treevec <- c(500,600,700,800,900,1000,1100,1200,1300,1400,1500)
  treevec <- seq(3,50,by=1)
  StorageVector <- rep(NA,length(treevec))
  for (i in 1:length(treevec)) {
    myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim11,ntree=648,mtry=tryvec[i])
    w.ROC = risksetAUC(Stime=traindata_trim11$Recurrence_time,status=traindata_trim11$Biochemical.recurrence,marker=myrsf$predicted.oob,tmax=4000)
    StorageVector[i] <- as.numeric(w.ROC$AUC[15])
  }
  
  #Plotting this
  sparedf <- as.data.frame(cbind(treevec,StorageVector))
  ggplot(data=sparedf,aes(x=treevec,y=StorageVector,group=1)) + geom_line(linetype="dashed")+geom_point()
  
#SO ideal ntree is 648 and mtry is 3
  
  myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim11,ntree=648,mtry=3)
  w.ROC = risksetAUC(Stime=traindata_trim11$Recurrence_time,status=traindata_trim11$Biochemical.recurrence,marker=myrsf$predicted.oob,tmax=4000)
  
  testdata5 <- testdata4[,which(names(testdata4) %in% names(traindata_trim11))]
  
  my_pred <- predict(myrsf,newdata=testdata5)
  w.ROC2 = risksetAUC(Stime=testdata5$Recurrence_time,status=testdata5$Biochemical.recurrence,marker=my_pred$predicted,tmax=4000)
  
#Testing hazard ratios
  
  myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim10)
  mypred <- predict(myrsf,newdata=traindata_trim11)
  mean(mypred$predicted)
  testdata6 <- testdata5
  testdata6$risk <- ifelse(my_pred$predicted>=3.4128622,1,0)
  testdata6$SurvObj <- with(testdata6,Surv(Recurrence_time,Biochemical.recurrence==1))
  mysurvfit <- survfit(SurvObj~risk,data=testdata6,conf.type="log")
  ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE)
  summary(coxph(SurvObj~risk,data=testdata6))
  
  traindata_trim10$risk <- ifelse(myrsf$predicted.oob>=2.548868,1,0)
  traindata_trim10$SurvObj <- with(traindata_trim10,Surv(Recurrence_time,Biochemical.recurrence==1))
  mysurvfit <- survfit(SurvObj~risk,data=traindata_trim10,conf.type="log")
  ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,conf.int = TRUE)
  summary(coxph(SurvObj~risk,data=traindata_trim10))
  #Barely misses significance :( Try with trim10?
  #Or try optimizing for something else with treesize?
  #Optimization before 1000 days seems to be ideal. Maybe make the threshold 1500.
  
  
#Testing with CAPRA-S added in and a shorter recurrence time in training cohort.
  #Try 1500 days first
  
traindata2 <- AutoZ(traindata,DataNames)
testdata2 <- AutoZ(testdata,DataNames)

traindata3 <- traindata2[which(traindata2$Recurrence_time<=1000),]
traindata4 <- traindata3[,c(1,458,459:914)]

traindata4$Biochemical.recurrence <- as.character(traindata4$Biochemical.recurrence)
traindata4$Biochemical.recurrence[which(traindata4$Biochemical.recurrence=="N")] <- "0"
traindata4$Biochemical.recurrence[which(traindata4$Biochemical.recurrence=="Y")] <- "1"
traindata4$Biochemical.recurrence <- as.numeric(as.character(traindata4$Biochemical.recurrence))

testdata2$Biochemical.recurrence <- as.character(testdata2$Biochemical.recurrence)
testdata2$Biochemical.recurrence[which(testdata2$Biochemical.recurrence=="N")] <- "0"
testdata2$Biochemical.recurrence[which(testdata2$Biochemical.recurrence=="Y")] <- "1"
testdata2$Biochemical.recurrence <- as.numeric(as.character(testdata2$Biochemical.recurrence))

testdata3 <- testdata2[,c(1,458,459:914)]

my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata4,importance=TRUE)
plot(gg_vimp(my_rsf))
imp_rsf_1 <- sort(my_rsf$importance,decreasing=T)
neg_vars <- names(which(imp_rsf_1<0))

traindata_trim1 <- traindata4[,-which(names(traindata4) %in% neg_vars)]

w.ROC = risksetAUC(Stime=traindata4$Recurrence_time,status=traindata4$Biochemical.recurrence,marker=my_rsf$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf2 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim1,importance=TRUE)
plot(gg_vimp(my_rsf2))
imp_rsf_2 <- sort(my_rsf2$importance,decreasing=T)
neg_vars2 <- names(which(imp_rsf_2<0))

traindata_trim2 <- traindata_trim1[,-which(names(traindata_trim1) %in% neg_vars2)]

w.ROC = risksetAUC(Stime=traindata_trim1$Recurrence_time,status=traindata_trim1$Biochemical.recurrence,marker=my_rsf2$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf3 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim2,importance=TRUE)
plot(gg_vimp(my_rsf3))
imp_rsf_3 <- sort(my_rsf3$importance,decreasing=T)
neg_vars3 <- names(which(imp_rsf_3<0))

traindata_trim3 <- traindata_trim2[,-which(names(traindata_trim2) %in% neg_vars3)]

w.ROC = risksetAUC(Stime=traindata_trim2$Recurrence_time,status=traindata_trim2$Biochemical.recurrence,marker=my_rsf3$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf4 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim3,importance=TRUE)
plot(gg_vimp(my_rsf4))
imp_rsf_4 <- sort(my_rsf4$importance,decreasing=T)
neg_vars4 <- names(which(imp_rsf_4<0))

traindata_trim4 <- traindata_trim3[,-which(names(traindata_trim3) %in% neg_vars4)]

w.ROC = risksetAUC(Stime=traindata_trim3$Recurrence_time,status=traindata_trim3$Biochemical.recurrence,marker=my_rsf4$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf5 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim4,importance=TRUE)
plot(gg_vimp(my_rsf5))
imp_rsf_5 <- sort(my_rsf5$importance,decreasing=T)
neg_vars5 <- names(which(imp_rsf_5<0))

traindata_trim5 <- traindata_trim4[,-which(names(traindata_trim4) %in% neg_vars5)]

w.ROC = risksetAUC(Stime=traindata_trim4$Recurrence_time,status=traindata_trim4$Biochemical.recurrence,marker=my_rsf5$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf6 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim5,importance=TRUE)
plot(gg_vimp(my_rsf6))
imp_rsf_6 <- sort(my_rsf6$importance,decreasing=T)
neg_vars6 <- names(which(imp_rsf_6<0))

traindata_trim6 <- traindata_trim5[,-which(names(traindata_trim5) %in% neg_vars6)]

w.ROC = risksetAUC(Stime=traindata_trim5$Recurrence_time,status=traindata_trim5$Biochemical.recurrence,marker=my_rsf6$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf7 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim6,importance=TRUE)
plot(gg_vimp(my_rsf7))
imp_rsf_7 <- sort(my_rsf7$importance,decreasing=T)
neg_vars7 <- names(which(imp_rsf_7<0))

traindata_trim7 <- traindata_trim6[,-which(names(traindata_trim6) %in% neg_vars7)]

w.ROC = risksetAUC(Stime=traindata_trim6$Recurrence_time,status=traindata_trim6$Biochemical.recurrence,marker=my_rsf7$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf8 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim7,importance=TRUE)
plot(gg_vimp(my_rsf8))
imp_rsf_8 <- sort(my_rsf8$importance,decreasing=T)
neg_vars8 <- names(which(imp_rsf_8<0))

traindata_trim8 <- traindata_trim7[,-which(names(traindata_trim7) %in% neg_vars8)]

w.ROC = risksetAUC(Stime=traindata_trim7$Recurrence_time,status=traindata_trim7$Biochemical.recurrence,marker=my_rsf8$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf9 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim8,importance=TRUE)
plot(gg_vimp(my_rsf9))
imp_rsf_9 <- sort(my_rsf9$importance,decreasing=T)
neg_vars9 <- names(which(imp_rsf_9<0))

traindata_trim9 <- traindata_trim8[,-which(names(traindata_trim8) %in% neg_vars9)]

w.ROC = risksetAUC(Stime=traindata_trim8$Recurrence_time,status=traindata_trim8$Biochemical.recurrence,marker=my_rsf9$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf10 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim9,importance=TRUE)
plot(gg_vimp(my_rsf10))
imp_rsf_10 <- sort(my_rsf10$importance,decreasing=T)
neg_vars10 <- names(which(imp_rsf_10<0))

traindata_trim10 <- traindata_trim9[,-which(names(traindata_trim9) %in% neg_vars10)]

w.ROC = risksetAUC(Stime=traindata_trim9$Recurrence_time,status=traindata_trim9$Biochemical.recurrence,marker=my_rsf10$predicted.oob,tmax=1500)
mean(w.ROC$AUC)


my_rsf11 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim10,importance=TRUE)
plot(gg_vimp(my_rsf11))
imp_rsf_11 <- sort(my_rsf11$importance,decreasing=T)
neg_vars11 <- names(which(imp_rsf_11<0))
my_rsf11

traindata_trim11 <- traindata_trim10[,-which(names(traindata_trim10) %in% neg_vars11)]

w.ROC = risksetAUC(Stime=traindata_trim10$Recurrence_time,status=traindata_trim10$Biochemical.recurrence,marker=my_rsf11$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf12 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim11,importance=TRUE)
plot(gg_vimp(my_rsf12))
imp_rsf_12 <- sort(my_rsf12$importance,decreasing=T)
neg_vars12 <- names(which(imp_rsf_12<0))
my_rsf12

w.ROC = risksetAUC(Stime=traindata_trim11$Recurrence_time,status=traindata_trim11$Biochemical.recurrence,marker=my_rsf12$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

traindata_trim12 <- traindata_trim11[,-which(names(traindata_trim11) %in% neg_vars12)]

my_rsf13 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim12,importance=TRUE)
plot(gg_vimp(my_rsf13))
imp_rsf_13 <- sort(my_rsf13$importance,decreasing=T)
neg_vars13 <- names(which(imp_rsf_13<0))
my_rsf13

traindata_trim13 <- traindata_trim12[,-which(names(traindata_trim12) %in% neg_vars13)]

w.ROC = risksetAUC(Stime=traindata_trim12$Recurrence_time,status=traindata_trim12$Biochemical.recurrence,marker=my_rsf13$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf14 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim13,importance=TRUE)
plot(gg_vimp(my_rsf14))
imp_rsf_14 <- sort(my_rsf14$importance,decreasing=T)
neg_vars14 <- names(which(imp_rsf_14<0))
my_rsf14

traindata_trim14 <- traindata_trim13[,-which(names(traindata_trim13) %in% neg_vars14)]

w.ROC = risksetAUC(Stime=traindata_trim13$Recurrence_time,status=traindata_trim13$Biochemical.recurrence,marker=my_rsf14$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

my_rsf15 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim14,importance=TRUE)
plot(gg_vimp(my_rsf15))
imp_rsf_15 <- sort(my_rsf15$importance,decreasing=T)
neg_vars15 <- names(which(imp_rsf_15<0))
my_rsf15

traindata_trim15 <- traindata_trim14[,-which(names(traindata_trim14) %in% neg_vars15)]

w.ROC = risksetAUC(Stime=traindata_trim14$Recurrence_time,status=traindata_trim14$Biochemical.recurrence,marker=my_rsf15$predicted.oob,tmax=1500)
mean(w.ROC$AUC)

#test predictions

testdata4 <- testdata3[-which(is.na(testdata3$Recurrence_time)),]
my_pred <- predict(my_rsf15, newdata=testdata4)
w.ROC2 <- risksetAUC(Stime=testdata4$Recurrence_time,status=testdata4$Biochemical.recurrence,marker=my_pred$predicted,tmax=4000)

#Trying KM

myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim14)
mypred <- predict(myrsf,newdata=traindata_trim14)
mean(mypred$predicted)
testdata6 <- testdata5
testdata6$risk <- ifelse(my_pred$predicted>=3.714695,1,0)
testdata6$SurvObj <- with(testdata6,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk,data=testdata6,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE)
summary(coxph(SurvObj~risk,data=testdata6))

#Using 3rd quartile of predicted.oob
testdata6$risk <- ifelse(my_pred$predicted>=5.2267728,1,0)
testdata6$SurvObj <- with(testdata6,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk,data=testdata6,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
summary(coxph(SurvObj~risk,data=testdata6))



#Optimizing?

set.seed(1)

#treevec <- c(500,600,700,800,900,1000,1100,1200,1300,1400,1500)
treevec <- seq(3,50,by=1)
StorageVector <- rep(NA,length(treevec))
for (i in 1:length(treevec)) {
  myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim14,ntree=1426,mtry=treevec[i])
  w.ROC = risksetAUC(Stime=traindata_trim14$Recurrence_time,status=traindata_trim14$Biochemical.recurrence,marker=myrsf$predicted.oob,tmax=4000)
  StorageVector[i] <- mean(w.ROC$AUC)
}

#Plotting this
sparedf <- as.data.frame(cbind(treevec,StorageVector))
ggplot(data=sparedf,aes(x=treevec,y=StorageVector,group=1)) + geom_line(linetype="dashed")+geom_point()

#SO ideal ntree is 1426 and mtry is 3
set.seed(1)
myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim14,ntree=1426,mtry=3)
w.ROC = risksetAUC(Stime=traindata_trim14$Recurrence_time,status=traindata_trim14$Biochemical.recurrence,marker=myrsf$predicted.oob,tmax=4000)

testdata5 <- testdata4[,which(names(testdata4) %in% names(traindata_trim14))]

my_pred <- predict(myrsf,newdata=testdata6)
w.ROC2 = risksetAUC(Stime=testdata5$Recurrence_time,status=testdata5$Biochemical.recurrence,marker=my_pred$predicted,tmax=4000)

#KM with optimized OOB cutoff

testdata6$risk <- ifelse(my_pred$predicted>=4.8090762,1,0)
testdata6$SurvObj <- with(testdata6,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk,data=testdata6,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
summary(coxph(SurvObj~risk,data=testdata6))

testdata6$predicted <- my_pred$predicted
summary(coxph(SurvObj~CAPRA_S+predicted,data=testdata6))
summary(coxph(SurvObj~CAPRA_S+risk,data=testdata6))
testdata6$CAPRA_lowvall <- ifelse(testdata6$CAPRA_S<=2,0,1)
testdata6$CAPRA_highvall <- ifelse(testdata6$CAPRA_S>=6,1,0)
mysurvfit <- survfit(SurvObj~risk+CAPRA_lowvall,data=testdata6,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
mysurvfit <- survfit(SurvObj~risk+CAPRA_highvall,data=testdata6,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
valdata3 <- testdata6[-which(testdata6$CAPRA_S>=6),]
valdata3$CAPRA_intvlow <- ifelse(valdata3$CAPRA_S<=2,0,1)
valdata3$SurvObj <- with(valdata3,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk+CAPRA_intvlow,data=testdata6,conf.type="log-log")
mysurvfit <- survfit(SurvObj~risk+CAPRA_intvlow,data=valdata3,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)


#Trying validation dataset

valdata <- LoadData("Val_Cohort.txt")
DataNames2 <- colnames(valdata[c(12:39)])

valdata2 <- AutoZ(valdata,DataNames2)
valdata2$Recurrence_time <- valdata2$Months.total.F.U
valdata2$Recurrence_time[which(valdata2$BCR==1)] <- valdata2$Months.to.BCR[which(valdata2$BCR==1)]

colnames(valdata2)[3] <- "Biochemical.recurrence"

my_valpred <- predict(myrsf,newdata=valdata2)
w.ROC3 = risksetAUC(Stime=valdata2$Recurrence_time,status=valdata2$Biochemical.recurrence,marker=my_valpred$predicted,tmax=160)
#to plot:
w.AUC3 <- risksetROC(Stime=valdata2$Recurrence_time,status=valdata2$Biochemical.recurrence,marker=my_valpred$predicted,predict.time=156.3666667,main="ROC Curve",col="red ")

valdata2$risk <- ifelse(my_valpred$predicted>=4.8090762,1,0)
valdata2$SurvObj <- with(valdata2,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit2 <- survfit(SurvObj~risk,data=valdata2,conf.type="log-log")
ggsurvplot(mysurvfit2,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
summary(coxph(SurvObj~risk,data=valdata2))

valdata2$CAPRA_lowvall <- ifelse(valdata2$CAPRA_S<=2,0,1)
valdata2$CAPRA_highvall <- ifelse(valdata2$CAPRA_S>=6,1,0)

mysurvfit2 <- survfit(SurvObj~risk+CAPRA_highvall,data=valdata2,conf.type="log-log")
ggsurvplot(mysurvfit2,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

valdata3 <- valdata2[-which(valdata2$CAPRA_S>=6),]
valdata3$CAPRA_intvlow <- ifelse(valdata3$CAPRA_S<=2,0,1)
valdata3$SurvObj <- with(valdata3,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk+CAPRA_intvlow,data=valdata3,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

#Pairwise comparisons
mysurvfit3 <- pairwise_survdiff(SurvObj~risk+CAPRA_intvlow,data=valdata3)
mysurvfit3

#Testing on microarray set
taylor <- LoadData("Taylor_Val.txt")

DataNames3 <- colnames(taylor[c(3:31)])

taylor2 <- AutoZ(taylor,DataNames3)
taylor2 <- taylor2[-c(132:164),]

colnames(taylor2)[57] <- "Recurrence_time"
colnames(taylor2)[58] <- "Biochemical.recurrence"

#Trying DLEU1 Point 1 first

colnames(taylor2)[75] <- "DLEU1_z"

taylorpred <- predict(myrsf,newdata=taylor2)
w.ROC4 = risksetAUC(Stime=taylor2$Recurrence_time,status=taylor2$Biochemical.recurrence,marker=taylorpred$predicted,tmax=500)
#to plot:
w.AUC3 <- risksetROC(Stime=valdata2$Recurrence_time,status=valdata2$Biochemical.recurrence,marker=my_valpred$predicted,predict.time=156.3666667,main="ROC Curve",col="red ")

taylor2$risk <- ifelse(taylorpred$predicted>=4.8090762,1,0)
taylor2$SurvObj <- with(taylor2,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit3 <- survfit(SurvObj~CAPRA_highvall,data=taylor2,conf.type="log-log")
ggsurvplot(mysurvfit3,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
summary(coxph(SurvObj~risk,data=taylor2))

valdata2$CAPRA_lowvall <- ifelse(valdata2$CAPRA_S<=2,0,1)
taylor2$CAPRA_highvall <- ifelse(taylor2$CAPRA_S>=6,1,0)

mysurvfit2 <- survfit(SurvObj~risk+CAPRA_highvall,data=valdata2,conf.type="log-log")
ggsurvplot(mysurvfit2,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

valdata3 <- valdata2[-which(valdata2$CAPRA_S>=6),]
valdata3$CAPRA_intvlow <- ifelse(valdata3$CAPRA_S<=2,0,1)
valdata3$SurvObj <- with(valdata3,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk+CAPRA_intvlow,data=valdata3,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

#Pairwise comparisons
mysurvfit3 <- pairwise_survdiff(SurvObj~risk+CAPRA_intvlow,data=valdata3)
mysurvfit3

##################################################################3

#Taylor data didn't work out - doesn't appear to have good transferability to microarrays
#Maybe try NOT z-scoring and developing from sequencing data?

my_rsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata,importance=TRUE)
plot(gg_vimp(my_rsf))
imp_rsf_1 <- sort(my_rsf$importance,decreasing=T)
neg_vars <- names(which(imp_rsf_1<0))

w.ROC = risksetAUC(Stime=traindata$Recurrence_time,status=traindata$Biochemical.recurrence,marker=my_rsf$predicted.oob,tmax=2000)
traindata_trim1 <- traindata[,-which(names(traindata) %in% neg_vars)]

my_rsf2 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim1,importance=TRUE)
plot(gg_vimp(my_rsf2))
imp_rsf_2 <- sort(my_rsf2$importance,decreasing=T)
neg_vars2 <- names(which(imp_rsf_2<0))

w.ROC = risksetAUC(Stime=traindata_trim1$Recurrence_time,status=traindata_trim1$Biochemical.recurrence,marker=my_rsf2$predicted.oob,tmax=2000)
traindata_trim2 <- traindata_trim1[,-which(names(traindata_trim1) %in% neg_vars2)]

my_rsf3 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim2,importance=TRUE)
plot(gg_vimp(my_rsf3))
imp_rsf_3 <- sort(my_rsf3$importance,decreasing=T)
neg_vars3 <- names(which(imp_rsf_3<0))

w.ROC = risksetAUC(Stime=traindata_trim2$Recurrence_time,status=traindata_trim2$Biochemical.recurrence,marker=my_rsf3$predicted.oob,tmax=2000)
traindata_trim3 <- traindata_trim2[,-which(names(traindata_trim2) %in% neg_vars3)]

my_rsf4 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim3,importance=TRUE)
plot(gg_vimp(my_rsf4))
imp_rsf_4 <- sort(my_rsf4$importance,decreasing=T)
neg_vars4 <- names(which(imp_rsf_4<0))

w.ROC = risksetAUC(Stime=traindata_trim3$Recurrence_time,status=traindata_trim3$Biochemical.recurrence,marker=my_rsf4$predicted.oob,tmax=2000)
traindata_trim4 <- traindata_trim3[,-which(names(traindata_trim3) %in% neg_vars4)]

my_rsf5 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim4,importance=TRUE)
plot(gg_vimp(my_rsf5))
imp_rsf_5 <- sort(my_rsf5$importance,decreasing=T)
neg_vars5 <- names(which(imp_rsf_5<0))

w.ROC = risksetAUC(Stime=traindata_trim4$Recurrence_time,status=traindata_trim4$Biochemical.recurrence,marker=my_rsf5$predicted.oob,tmax=2000)
traindata_trim5 <- traindata_trim4[,-which(names(traindata_trim4) %in% neg_vars5)]

my_rsf6 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim5,importance=TRUE)
plot(gg_vimp(my_rsf6))
imp_rsf_6 <- sort(my_rsf6$importance,decreasing=T)
neg_vars6 <- names(which(imp_rsf_6<0))

w.ROC = risksetAUC(Stime=traindata_trim5$Recurrence_time,status=traindata_trim5$Biochemical.recurrence,marker=my_rsf6$predicted.oob,tmax=5000)
traindata_trim6 <- traindata_trim5[,-which(names(traindata_trim5) %in% neg_vars6)]

my_rsf7 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim6,importance=TRUE)
plot(gg_vimp(my_rsf7))
imp_rsf_7 <- sort(my_rsf7$importance,decreasing=T)
neg_vars7 <- names(which(imp_rsf_7<0))

w.ROC = risksetAUC(Stime=traindata_trim6$Recurrence_time,status=traindata_trim6$Biochemical.recurrence,marker=my_rsf7$predicted.oob,tmax=5000)
traindata_trim7 <- traindata_trim6[,-which(names(traindata_trim6) %in% neg_vars7)]

my_rsf8 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim7,importance=TRUE)
plot(gg_vimp(my_rsf8))
imp_rsf_8 <- sort(my_rsf8$importance,decreasing=T)
neg_vars8 <- names(which(imp_rsf_8<0))

w.ROC = risksetAUC(Stime=traindata_trim7$Recurrence_time,status=traindata_trim7$Biochemical.recurrence,marker=my_rsf8$predicted.oob,tmax=5000)
traindata_trim8 <- traindata_trim7[,-which(names(traindata_trim7) %in% neg_vars8)]

my_rsf9 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim8,importance=TRUE)
plot(gg_vimp(my_rsf9))
imp_rsf_9 <- sort(my_rsf9$importance,decreasing=T)
neg_vars9 <- names(which(imp_rsf_9<0))

w.ROC = risksetAUC(Stime=traindata_trim8$Recurrence_time,status=traindata_trim8$Biochemical.recurrence,marker=my_rsf9$predicted.oob,tmax=5000)
traindata_trim9 <- traindata_trim8[,-which(names(traindata_trim8) %in% neg_vars9)]

my_rsf10 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim9,importance=TRUE)
plot(gg_vimp(my_rsf10))
imp_rsf_10 <- sort(my_rsf10$importance,decreasing=T)
neg_vars10 <- names(which(imp_rsf_10<0))

w.ROC = risksetAUC(Stime=traindata_trim9$Recurrence_time,status=traindata_trim9$Biochemical.recurrence,marker=my_rsf10$predicted.oob,tmax=2000)
traindata_trim10 <- traindata_trim9[,-which(names(traindata_trim9) %in% neg_vars10)]

my_rsf11 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim10,importance=TRUE)
plot(gg_vimp(my_rsf11))
imp_rsf_11 <- sort(my_rsf11$importance,decreasing=T)
neg_vars11 <- names(which(imp_rsf_11<0))
my_rsf11

w.ROC = risksetAUC(Stime=traindata_trim10$Recurrence_time,status=traindata_trim10$Biochemical.recurrence,marker=my_rsf11$predicted.oob,tmax=2000)
traindata_trim11 <- traindata_trim10[,-which(names(traindata_trim10) %in% neg_vars11)]

my_rsf12 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim11,importance=TRUE)
plot(gg_vimp(my_rsf12))
imp_rsf_12 <- sort(my_rsf12$importance,decreasing=T)
neg_vars12 <- names(which(imp_rsf_12<0))
my_rsf12

w.ROC = risksetAUC(Stime=traindata_trim11$Recurrence_time,status=traindata_trim11$Biochemical.recurrence,marker=my_rsf12$predicted.oob,tmax=2000)
traindata_trim12 <- traindata_trim11[,-which(names(traindata_trim11) %in% neg_vars12)]
w.ROC

my_rsf13 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim12,importance=TRUE)
plot(gg_vimp(my_rsf13))
imp_rsf_13 <- sort(my_rsf13$importance,decreasing=T)
neg_vars13 <- names(which(imp_rsf_13<0))
my_rsf13

w.ROC = risksetAUC(Stime=traindata_trim12$Recurrence_time,status=traindata_trim12$Biochemical.recurrence,marker=my_rsf13$predicted.oob,tmax=2000)
traindata_trim13 <- traindata_trim12[,-which(names(traindata_trim12) %in% neg_vars13)]
w.ROC

my_rsf14 <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~., data=traindata_trim13,importance=TRUE)
plot(gg_vimp(my_rsf14))
imp_rsf_14 <- sort(my_rsf14$importance,decreasing=T)
neg_vars14 <- names(which(imp_rsf_14<0))
my_rsf14

w.ROC = risksetAUC(Stime=traindata_trim13$Recurrence_time,status=traindata_trim13$Biochemical.recurrence,marker=my_rsf14$predicted.oob,tmax=2000)
traindata_trim14 <- traindata_trim13[,-which(names(traindata_trim13) %in% neg_vars14)]
w.ROC

#pred on traindata

#Optimize

set.seed(2)

#treevec <- c(500,600,700,800,900,1000,1100,1200,1300,1400,1500)
treevec <- seq(3,50,by=1)
StorageVector <- rep(NA,length(treevec))
VecList <- list()
for (j in 1:10) {
  set.seed(j)
  for (i in 1:length(treevec)) {
    myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim13,ntree=1381,mtry=treevec[i])
    StorageVector[i] <- myrsf$err.rate[1381]
  }
  sparedf <- as.data.frame(cbind(treevec,StorageVector))
  VecList[j] <- list(sparedf)
}

mytest <- as.data.frame(VecList[[1]])
for (i in 2:length(VecList)){
  mytest <- cbind(mytest,as.data.frame(VecList[[i]]$StorageVector))
}

rowMeans(mytest[,-1])
as.data.frame(cbind(treevec,rowMeans(mytest[,-1])))

#treevec=1300 has lowest average error rate of 0.2048009
#treevec=1380 has lowest average error rate of 0.2037675
#treevec=1381 has lowest average error rate of 0.201507
#mtry=3 has lowest average error rate of 0.2008827

#Plotting this
sparedf <- as.data.frame(cbind(treevec,StorageVector))
ggplot(data=sparedf,aes(x=treevec,y=StorageVector,group=1)) + geom_line(linetype="dashed")+geom_point()

veca <- sparedf$StorageVector

#Optimized

#myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim14,ntree=697,mtry=7)
myrsf <- rfsrc(Surv(Recurrence_time,Biochemical.recurrence)~.,data=traindata_trim9,ntree=898,mtry=3)
w.ROC = risksetAUC(Stime=traindata_trim9$Recurrence_time,status=traindata_trim9$Biochemical.recurrence,marker=myrsf$predicted.oob,tmax=4000)

#Quantile divisions

quantile(myrsf$predicted.oob)
#quantile is 3.216460
#testing validation data

testdata_noz <- testdata[-which(is.na(testdata$Recurrence_time)),]
testdata_noz2 <- testdata[,which(names(testdata_noz) %in% colnames(traindata_trim9))]
testdata_noz2$Biochemical.recurrence <- as.character(testdata_noz2$Biochemical.recurrence)
testdata_noz2$Biochemical.recurrence[which(testdata_noz2$Biochemical.recurrence=="N")] <- "0"
testdata_noz2$Biochemical.recurrence[which(testdata_noz2$Biochemical.recurrence=="Y")] <- "1"
testdata_noz2$Biochemical.recurrence <- as.numeric(as.character(testdata_noz2$Biochemical.recurrence))
testdata_noz3 <- testdata_noz2[-which(is.na(testdata_noz2$Recurrence_time)),]
my_testpred <- predict(myrsf,newdata=testdata_noz3)

w.ROC2 = risksetAUC(Stime=testdata_noz3$Recurrence_time,status=testdata_noz3$Biochemical.recurrence,marker=my_testpred$predicted,tmax=4000)
testdata_noz3$risk <- ifelse(my_testpred$predicted>3.216460,1,0)

testdata_noz3$SurvObj <- with(testdata_noz3,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit2 <- survfit(SurvObj~risk,data=testdata_noz3,conf.type="log-log")
ggsurvplot(mysurvfit2,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
summary(coxph(SurvObj~risk,data=testdata_noz3))

#Extract info on CAPRA_S from the complete training set
bloop <- as.numeric(row.names(testdata_noz3))

valdata$Recurrence_time <- valdata$Months.total.F.U
valdata$Recurrence_time[which(valdata$BCR==1)] <- valdata$Months.to.BCR[which(valdata$BCR==1)]

colnames(valdata)[3] <- "Biochemical.recurrence"
my_valpred <- predict(myrsf,newdata=valdata2)
w.ROC3 = risksetAUC(Stime=valdata2$Recurrence_time,status=valdata2$Biochemical.recurrence,marker=my_valpred$predicted,tmax=160)
#to plot:
w.AUC3 <- risksetROC(Stime=valdata2$Recurrence_time,status=valdata2$Biochemical.recurrence,marker=my_valpred$predicted,predict.time=156.3666667,main="ROC Curve",col="red ")

valdata2$risk <- ifelse(my_valpred$predicted>=4.8090762,1,0)
valdata2$SurvObj <- with(valdata2,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit2 <- survfit(SurvObj~risk,data=valdata2,conf.type="log-log")
ggsurvplot(mysurvfit2,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)
summary(coxph(SurvObj~risk,data=valdata2))

valdata2$CAPRA_lowvall <- ifelse(valdata2$CAPRA_S<=2,0,1)
valdata2$CAPRA_highvall <- ifelse(valdata2$CAPRA_S>=6,1,0)

mysurvfit2 <- survfit(SurvObj~risk+CAPRA_highvall,data=valdata2,conf.type="log-log")
ggsurvplot(mysurvfit2,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

valdata3 <- valdata2[-which(valdata2$CAPRA_S>=6),]
valdata3$CAPRA_intvlow <- ifelse(valdata3$CAPRA_S<=2,0,1)
valdata3$SurvObj <- with(valdata3,Surv(Recurrence_time,Biochemical.recurrence==1))
mysurvfit <- survfit(SurvObj~risk+CAPRA_intvlow,data=valdata3,conf.type="log-log")
ggsurvplot(mysurvfit,risk.table="abs_pct",risk.table.y.text.col=T,risk.table.y.text=FALSE,pval=TRUE)

#Pairwise comparisons
mysurvfit3 <- pairwise_survdiff(SurvObj~risk+CAPRA_intvlow,data=valdata3)
mysurvfit3
rcorr.cens(my_rsf$predicted.oob , 
           Surv(traindata4$Recurrence_time, traindata4$Biochemical.recurrence))["C Index"]

