#Plotting scatters using ggplot2

if(!require("ggplot2")) {
  install.packages("ggplot2")
}

if(!require("gridExtra")) {
  install.packages("gridExtra")
}

if(!require("dplyr")) {
  install.packages("dplyr")
}


library(ggplot2)
library(gridExtra)
library(dplyr)
library(plotrix)

UNC5D <- read.table("UNC5DMatched.txt",sep="\t",header=TRUE)
SORCS3 <- read.table("SORCS3Matched.txt",sep="\t",header=TRUE)
NELL2 <- read.table("NELL2Matched.txt",sep="\t",header=TRUE)
GRAMD1B <- read.table("GRAMD1BMatched.txt",sep="\t",header=TRUE)
GNAO1 <- read.table("GNAO1Matched.txt",sep="\t",header=TRUE)
ELFN1 <- read.table("ELFN1Matched.txt",sep="\t",header=TRUE)
ASB2 <- read.table("ASB2Matched.txt",sep="\t",header=TRUE)

ASB2Summary <- ASB2 %>% group_by(Sample.Type) %>% summarise_all(.funs=c(mean="mean",std.error))

ggplot(ASB2,aes(x=factor(Sample.Type),y=ASB2_cg01956154,color=factor(Sample.Type))) +
  geom_point(aes(y=ASB2_cg01956154,color=factor(Sample.Type)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=ASB2_cg01956154_mean),color="black",size=2,data=ASB2Summary) +
  geom_errorbar(aes(y=ASB2_cg01956154_mean,ymin=ASB2_cg01956154_mean-cg1_se, ymax=ASB2_cg01956154_mean+cg1_se), color="black",width=0.2,data=ASB2Summary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("ASB2_cg01956154") +
  labs(x = "Group",y="Methylation beta value")

UNC5DSummary <- UNC5D %>% group_by(Class) %>% summarize(Gene_mean = mean(UNC5D),Gene_se=sqrt(var(UNC5D)/length(UNC5D)))

UNC5DPlot <- ggplot(UNC5D,aes(x=factor(Class),y=UNC5D,color=factor(Class))) +
  geom_point(aes(y=UNC5D,color=factor(Class)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=Gene_mean),color="black",size=2,data=UNC5DSummary) +
  geom_errorbar(aes(y=Gene_mean,ymin=Gene_mean-Gene_se, ymax=Gene_mean+Gene_se), color="black",width=0.2,data=UNC5DSummary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("UNC5D") +
  labs(x = "Group",y="Methylation beta value")

SORCS3Summary <- SORCS3 %>% group_by(Class) %>% summarize(Gene_mean = mean(SORCS3),Gene_se=sqrt(var(SORCS3)/length(SORCS3)))

SORCS3Plot <- ggplot(SORCS3,aes(x=factor(Class),y=SORCS3,color=factor(Class))) +
  geom_point(aes(y=SORCS3,color=factor(Class)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=Gene_mean),color="black",size=2,data=SORCS3Summary) +
  geom_errorbar(aes(y=Gene_mean,ymin=Gene_mean-Gene_se, ymax=Gene_mean+Gene_se), color="black",width=0.2,data=SORCS3Summary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("SORCS3") +
  labs(x = "Group",y="Methylation beta value")

NELL2Summary <- NELL2 %>% group_by(Class) %>% summarize(Gene_mean = mean(NELL2),Gene_se=sqrt(var(NELL2)/length(NELL2)))

NELL2Plot <- ggplot(NELL2,aes(x=factor(Class),y=NELL2,color=factor(Class))) +
  geom_point(aes(y=NELL2,color=factor(Class)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=Gene_mean),color="black",size=2,data=NELL2Summary) +
  geom_errorbar(aes(y=Gene_mean,ymin=Gene_mean-Gene_se, ymax=Gene_mean+Gene_se), color="black",width=0.2,data=NELL2Summary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("NELL2") +
  labs(x = "Group",y="Methylation beta value")

GRAMD1BSummary <- GRAMD1B %>% group_by(Class) %>% summarize(Gene_mean = mean(GRAMD1B),Gene_se=sqrt(var(GRAMD1B)/length(GRAMD1B)))

GRAMD1BPlot <- ggplot(GRAMD1B,aes(x=factor(Class),y=GRAMD1B,color=factor(Class))) +
  geom_point(aes(y=GRAMD1B,color=factor(Class)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=Gene_mean),color="black",size=2,data=GRAMD1BSummary) +
  geom_errorbar(aes(y=Gene_mean,ymin=Gene_mean-Gene_se, ymax=Gene_mean+Gene_se), color="black",width=0.2,data=GRAMD1BSummary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("GRAMD1B") +
  labs(x = "Group",y="Methylation beta value")

GNAO1Summary <- GNAO1 %>% group_by(Class) %>% summarize(Gene_mean = mean(GNAO1),Gene_se=sqrt(var(GNAO1)/length(GNAO1)))

GNAO1Plot <- ggplot(GNAO1,aes(x=factor(Class),y=GNAO1,color=factor(Class))) +
  geom_point(aes(y=GNAO1,color=factor(Class)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=Gene_mean),color="black",size=2,data=GNAO1Summary) +
  geom_errorbar(aes(y=Gene_mean,ymin=Gene_mean-Gene_se, ymax=Gene_mean+Gene_se), color="black",width=0.2,data=GNAO1Summary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("GNAO1") +
  labs(x = "Group",y="Methylation beta value")

ASB2Summary <- ASB2 %>% group_by(Class) %>% summarize(ASB2_mean = mean(ASB2),ASB2_se=sqrt(var(ASB2)/length(ASB2)))

ASB2Plot <- ggplot(ASB2,aes(x=factor(Class),y=ASB2,color=factor(Class))) +
  geom_point(aes(y=ASB2,color=factor(Class)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=ASB2_mean),color="black",size=2,data=ASB2Summary) +
  geom_errorbar(aes(y=ASB2_mean,ymin=ASB2_mean-ASB2_se, ymax=ASB2_mean+ASB2_se), color="black",width=0.2,data=ASB2Summary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("ASB2") +
  labs(x = "Group",y="Methylation beta value")

ELFN1Summary <- ELFN1 %>% group_by(Class) %>% summarize(Gene_mean = mean(ELFN1),Gene_se=sqrt(var(ELFN1)/length(ELFN1)))

ELFN1Plot <- ggplot(ELFN1,aes(x=factor(Class),y=ELFN1,color=factor(Class))) +
  geom_point(aes(y=ELFN1,color=factor(Class)), position=position_jitter(width=0.05, height=0.0), alpha=0.6) +
  geom_point(aes(y=Gene_mean),color="black",size=2,data=ELFN1Summary) +
  geom_errorbar(aes(y=Gene_mean,ymin=Gene_mean-Gene_se, ymax=Gene_mean+Gene_se), color="black",width=0.2,data=ELFN1Summary) +
  ylim(0,1) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("ELFN1") +
  labs(x = "Group",y="Methylation beta value")

#Boxplots for PMR

ggplot(DataTable, aes(x = Type, y = V2, color = Type)) +
geom_boxplot(notch=TRUE) +
theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("PMR for Matched GS7 Tumor vs Normal Samples") +
labs(x="Tissue Type",y="Percent of Methylated Reference (PMR)")

ggplot(SRPX,aes(x=Sample.Type,y=SRPX_cg03509565,fill=Sample.Type)) +
   geom_boxplot(notch=TRUE) +
   theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
   ggtitle("SRPX_cg03509565") +
   labs(x="Group",y="Methylation beta value")

#Paired scatterplot for PMR

ggplot(DataTable2, aes(x=Type, y=V2)) +
geom_point(size=4,aes(colour=factor(Type))) +
geom_line(aes(group = IDs)) +
  theme(panel.background=element_blank(),axis.line = element_line(colour = "black"),legend.position="none") +
  ggtitle("Paired PMRs for Matched GS7 Tumor vs Normal Samples") +
  labs(x="Tissue Type",y="Percent of Methylated Reference (PMR)")
