library(survivalROC)

traindata$ModelSums <- -2.13349+(0.301827*traindata$RGS11_z)+(0.435413*traindata$SPAG5_z)
traindata$BCR <- as.numeric(traindata$Biochemical.recurrence)
traindata$BCR[which(traindata$BCR==1)] <- 0
traindata$BCR[which(traindata$BCR==2)] <- 1
traindata <- traindata[-which(is.na(traindata$BCR)),]
traindata$modcap <- -1.45983+(0.733329*traindata$ModelSums)+(0.183707*traindata$CAPRA_S)
traindata$modcapage <- -1.35768 + (0.859355*traindata$ModelSums) + (0.168511*traindata$CAPRA_S) + (0.003623*traindata$Age)

testdata$modcap <- -1.45983+(0.733329*testdata$ModelSums)+(0.183707*testdata$CAPRA_S)

#Need to adjust for different prediction times (right now set to date of last BCR event)
#set to time of last event in train and use diff times in test?
traincap <- survivalROC(traindata$BCR_time,traindata$BCR,traindata$CAPRA_S,predict.time=2000,method="NNE",lambda=0.001)
trainmodel <- survivalROC(traindata$BCR_time,traindata$BCR,traindata$ModelSums,predict.time=2000,method="NNE",lambda=0.001)
trainage <- survivalROC(traindata$BCR_time,traindata$BCR,traindata$Age,predict.time=2000,method="NNE",lambda=0.001)
trainmodcap <- survivalROC(traindata$BCR_time,traindata$BCR,traindata$modcap,predict.time=2000,method="NNE",lambda=0.001)
trainmodcapage <- survivalROC(traindata$BCR_time,traindata$BCR,traindata$modcapage,predict.time=2000,method="NNE",lambda=0.001)


library(viridis)

plot(traincap$FP,traincap$TP,type="n")
lines(traincap$FP,traincap$TP,type="s",col="pink")
lines(trainmodel$FP,trainmodel$TP,type="s",col="red")
lines(trainage$FP,trainage$TP,type="s",col="orange")
lines(trainmodcap$FP,trainmodcap$TP,type="s",col="green")
lines(trainmodcapage$FP,trainmodcapage$TP,type="s",col="purple")
abline(0,1)

#Establish new cutpoints

#Find the optimal cutoff point maximizing sensitivity and specificity, based on a performance and prediction object
opt.cut_surv = function(survdat){
  max_sum <- survdat$TP + (1-survdat$FP)
  cut.ind <- which(max_sum==max(max_sum))
  return(cut.ind)
}

#Alternate function to find the optimal cutoff point with sensitivity or specificity set to some minimum value (default=90) while maximizing the other.
opt.cut.senspec_surv = function(survdat,sens=TRUE,minval=0.90){
  cutoffs=as.data.frame(cbind(survdat$cut.values,survdat$TP,(1-survdat$FP)))
  colnames(cutoffs) <- c("cutpoint","sens","spec")
  if (sens==TRUE) {
    cutoffs <- cutoffs[order(cutoffs$spec,decreasing=TRUE),]
    cutoffchoice <- subset(cutoffs,sens>minval)
  } else {
    cutoffs <- cutoffs[order(cutoffs$sens,decreasing=TRUE),]
    cutoffchoice <- subset(cutoffs,spec>minval)
  }
  return(cutoffchoice[1,])
}

#For each model, we want to:

#a) get the cutpoint (maxing both, maxing sens, or maxing spec - so three cutpoints are possible)
#b) apply this cutpoint to test data in different subdivisions
#c) plot out the ROCs
#d) Return ROCs, sens, and spec values

my_times <- c(12,24,36,60,84,183)
my_iterative_list <- list(trainmodel,traincap,trainage,trainmodcap,trainmodcapage)

find_cutoffs <- function(data,maxsenspec) {
  if (maxsenspec=="none"){
    mycut <- opt.cut_surv(data)
    cutoffs <- as.data.frame(cbind(data$cut.values[mycut],data$TP[mycut],(1-data$FP[mycut])))
    colnames(cutoffs) <- c("cutpoint","sens","spec")
  } else if (maxsenspec=="sens") {
    cutoffs <- opt.cut.senspec_surv(data)
  } else {
    cutoffs <- opt.cut.senspec_surv(data,sens=FALSE)
  }
  return(cutoffs)
}

time_AUC_analysis <- function(data,predict_time,time_to_event,censor,modelname,cutdata=cutoffs) {
  test_result_auc <- survivalROC(data[,time_to_event],data[,censor],data[,modelname],predict.time=predict_time,method="NNE",lambda=0.001)
  test_result <- survivalROC(data[,time_to_event],data[,censor],data[,modelname],predict.time=predict_time,method="NNE",lambda=0.001,cut.values = as.numeric(cutdata[1]))
  time_AUC_tests <- list(test_result_auc,test_result)
  return(time_AUC_tests)
}

my_rocobject_list <- list(trainmodel,traincap,trainage,trainmodcap,trainmodcapage)

FindSenspecTable <- function(data,survROCobject,cutparameter,time_name,event_name,model_name,times=my_times) {
  #Initialize storage table and list
  SenspecTransferTable <- data.frame(matrix(data=rep("NA"),nrow=length(my_times),ncol=3),stringsAsFactors = FALSE)
  overall_AUC_list <- list()
  #For the my_times list, get sens/spec/AUC values for either:
  #a)All possible cutoffs generated (overall_AUC_list)
  #b)Using cutoffs derived from the training set (SenspecTransferTable)
  cutoffs <- find_cutoffs(survROCobject,cutparameter)
  for (i in 1:length(times)) {
    time_result <- time_AUC_analysis(data,times[i],time_name,event_name,model_name,cutdata=cutoffs)
    overall_AUC <- as.data.frame(cbind(time_result[1][[1]]$TP,time_result[1][[1]]$FP))
    just_AUC <- time_result[1][[1]]$AUC
    colnames(overall_AUC) <- c("TP","FP")
    senspec_time <- c(time_result[2][[1]]$TP[2],(1-time_result[2][[1]]$FP[2]),just_AUC)
    SenspecTransferTable[i,] <- as.numeric(senspec_time)
    overall_AUC_list <- c(overall_AUC_list,overall_AUC)
  }
  colnames(SenspecTransferTable) <- c("sens","spec","AUC")
  row.names(SenspecTransferTable) <- times
  list_of_lists <- list(SenspecTransferTable,overall_AUC_list)
  return(list_of_lists)
}


SenspecTransferTable <- data.frame(matrix(data=rep("NA"),nrow=length(my_times),ncol=3),stringsAsFactors = FALSE)
overall_AUC_list <- list()

cutoffs <- find_cutoffs(traincap,"sens")
for (i in 1:length(my_times)){
  time_result <- time_AUC_analysis(testdata,my_times[i],"time_to_BCR","Biochemical.recurrence","CAPRA_S",cutdata=2)
  overall_AUC <- as.data.frame(cbind(time_result[1][[1]]$TP,time_result[1][[1]]$FP))
  just_AUC <- time_result[1][[1]]$AUC
  colnames(overall_AUC) <- c("TP","FP")
  senspec_time <- c(time_result[2][[1]]$TP[2],(1-time_result[2][[1]]$FP[2]),just_AUC)
  SenspecTransferTable[i,] <- as.numeric(senspec_time)
  overall_AUC_list <- c(overall_AUC_list,overall_AUC)
}

colnames(SenspecTransferTable) <- c("sens","spec","AUC")
row.names(SenspecTransferTable) <- my_times

none_cap_table <- SenspecTransferTable
none_cap_AUC <- overall_AUC_list


sens_model_table
sens_model_AUC
sens_cap_table
sens_cap_AUC
sens_age_table
sens_age_AUC
sens_modcap_table
sens_modcap_AUC
sens_modcapage_table
sens_modcapage_AUC

sink('Time To Event AUCs Validation.csv')

cat("2G Model","\n")
write.csv(sens_model_table)
cat("\n")

cat("CAPRA-S","\n")
write.csv(sens_cap_table)
cat("\n")

cat("Age","\n")
write.csv(sens_age_table)
cat("\n")

cat("2G Model + CAPRA-S","\n")
write.csv(sens_modcap_table)
cat("\n")

cat("2G Model + CAPRA-S + Age","\n")
write.csv(sens_modcapage_table)
cat("\n")

sink()

my_times <- c(12,24,36,60,84,150)
my_iterative_list <- list(trainmodel,traincap,trainage,trainmodcap,trainmodcapage)

SenspecTransferTable <- data.frame(matrix(data=rep("NA"),nrow=length(my_times),ncol=3),stringsAsFactors = FALSE)
overall_AUC_list <- list()

cutoffs <- find_cutoffs(trainmodcapage,"sens")
for (i in 1:length(my_times)){
  time_result <- time_AUC_analysis(taylordata,my_times[i],"SurvTime","Metastasis","modcapage",cutdata=cutoffs)
  overall_AUC <- as.data.frame(cbind(time_result[1][[1]]$TP,time_result[1][[1]]$FP))
  just_AUC <- time_result[1][[1]]$AUC
  colnames(overall_AUC) <- c("TP","FP")
  senspec_time <- c(time_result[2][[1]]$TP[2],(1-time_result[2][[1]]$FP[2]),just_AUC)
  SenspecTransferTable[i,] <- as.numeric(senspec_time)
  overall_AUC_list <- c(overall_AUC_list,overall_AUC)
}

colnames(SenspecTransferTable) <- c("sens","spec","AUC")
row.names(SenspecTransferTable) <- my_times

sens_modcapage_mets_taylor_table <- SenspecTransferTable
sens_modcapage_mets_taylor_AUC <- overall_AUC_list

sink('Time To Event AUCs Test Taylor Cohort.csv')

cat("BCR","\n")

cat("2G Model","\n")
write.csv(sens_model_taylor_table)
cat("\n")

cat("CAPRA-S","\n")
write.csv(sens_cap_taylor_table)
cat("\n")

cat("Age","\n")
write.csv(sens_age_taylor_table)
cat("\n")

cat("2G Model + CAPRA-S","\n")
write.csv(sens_modcap_taylor_table)
cat("\n")

cat("2G Model + CAPRA-S + Age","\n")
write.csv(sens_modcapage_taylor_table)
cat("\n")

cat("\n")
cat("\n")
cat("____________________________________________________","\n")

cat("Metastasis","\n")

cat("2G Model","\n")
write.csv(sens_model_mets_taylor_table)
cat("\n")

cat("CAPRA-S","\n")
write.csv(sens_cap_mets_taylor_table)
cat("\n")

cat("Age","\n")
write.csv(sens_age_mets_taylor_table)
cat("\n")

cat("2G Model + CAPRA-S","\n")
write.csv(sens_modcap_mets_taylor_table)
cat("\n")

cat("2G Model + CAPRA-S + Age","\n")
write.csv(sens_modcapage_mets_taylor_table)
cat("\n")

sink()

#2 years
plot((unlist(sens_model_mets_taylor_AUC[12])),unlist(sens_model_mets_taylor_AUC[11]),type="n")
lines((unlist(sens_model_mets_taylor_AUC[12])),unlist(sens_model_mets_taylor_AUC[11]),type="s",col="red")
lines(unlist(sens_cap_mets_taylor_AUC[12]),unlist(sens_cap_mets_taylor_AUC[11]),type="s",col="pink")
lines(unlist(sens_age_mets_taylor_AUC[12]),unlist(sens_age_mets_taylor_AUC[11]),type="s",col="orange")
lines(unlist(sens_modcap_mets_taylor_AUC[12]),unlist(sens_modcap_mets_taylor_AUC[11]),type="s",col="green")
lines(unlist(sens_modcapage_mets_taylor_AUC[12]),unlist(sens_modcapage_mets_taylor_AUC[11]),type="s",col="purple")
abline(0,1)
