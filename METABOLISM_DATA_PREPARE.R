rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/")
#================================================
#分别挑出两个队列的代谢数据
#================================================
library(data.table)
library(readxl)
set.seed(12345)
# rawdata <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/input_data/filtered_data.xlsx")
# rawdata <- subset(rawdata,rawdata$MS2_name!="NA"|rawdata$MS1_name!="NA")
# save(rawdata,file="metabolism_rawdata.RData")
load("input_data/metabolism_rawdata.RData")
rawdata <- rawdata %>% column_to_rownames(var="MS2_name")

MGS_meta <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt")

MGS <- read.table("input_data/MGS_metabolism_link.txt",header = T)

#内标代谢物
IS <- c("IS1","IS2","IS3","IS4","IS5","IS6")

#====bigmeta====
bigmeta_data <- rawdata[,colnames(rawdata)%in%MGS$metab_sample]
colnames(bigmeta_data) <- MGS$mgs_sample[match(colnames(bigmeta_data),MGS$metab_sample)]
colnames(bigmeta_data) <- gsub("-",".",colnames(bigmeta_data))
bigmeta_data <- bigmeta_data[,colnames(bigmeta_data)%in%MGS_meta$Row.names]# 415
bigmeta_data <- bigmeta_data[!rownames(bigmeta_data)%in%IS,]

#====renxx====
renxx_data <- rawdata[,!colnames(rawdata)%in%MGS$metab_sample]
renxx_data <- renxx_data[,-(1:7)]
renxx_data <- renxx_data[,!grepl("QC",colnames(renxx_data))]
colnames(renxx_data) <- gsub("X","",colnames(renxx_data))
colnames(renxx_data) <- gsub("_J","",colnames(renxx_data))
renxx_data <- renxx_data[!rownames(renxx_data)%in%IS,]

renxx_sample <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
renxx_data <- renxx_data[,colnames(renxx_data)%in%renxx_sample$sample] # 95


#====数据预处理====
#====bigmeta====
bigmeta <- bigmeta_data[rowSums(bigmeta_data!=0)>ncol(bigmeta_data)*0.95,]
#==================================
#插补
#==================================
#install.packages("DMwR2")
library(DMwR2)
raw <- bigmeta
raw[raw==0]=NA
#看了一下函数给的参考data,是行为样本列为变量
raw_knn <- knnImputation(as.data.frame(t(raw)),k=10,scale = T, meth ='weighAvg',distData = NULL)
raw_knn <- as.data.frame(t(raw_knn))

#====================================
#mean-centering
#====================================
##均值归一化
for(i in 1:ncol(raw_knn)){
  a <- raw_knn[,i]
  a <- mean(a)
  raw_knn[,i]<- sapply(raw_knn[,i], function(x)y <- x/a)
}

#===========================================================================
##log2转换
#===========================================================================
raw_knn <- log2(raw_knn)
boxplot(raw_knn,las=2)
sc <- as.data.frame(t(raw_knn))
save(sc,file="./sc_BIGMeta2500618.RData")

#==========renxx====================
#师兄队列的样本
#==============================
renxx <- renxx_data
renxx <- renxx[rowSums(renxx!=0)/ncol(renxx)>0.95,]

#==================================
#插补
#==================================
#install.packages("DMwR2")
library(DMwR2)
raw <- renxx
raw[raw==0]=NA
raw_knn <- knnImputation(as.data.frame(t(raw)),k=10,scale = T, meth ='mean',distData = NULL)
raw_knn <- as.data.frame(t(raw_knn))

#====================================
#mean-centering
#====================================
##均值归一化
for(i in 1:ncol(raw_knn)){
  a <- raw_knn[,i]
  a <- mean(a)
  raw_knn[,i]<- sapply(raw_knn[,i], function(x)y <- x/a)
}

#raw_knn <- scale(as.data.frame(t(raw_knn)))
#===========================================================================
##log2转换
#===========================================================================
raw_knn <- log2(raw_knn)
boxplot(raw_knn,las=2)
sc_rxx <- as.data.frame(t(raw_knn))
# sc_rxx <- as.data.frame(scale(sc_rxx))

save(sc_rxx,file="sc_renxx250618.RData")

#====数据预处理2====
#====renxx====
renxx <- renxx_data
renxx <- renxx[rowSums(renxx!=0)/ncol(renxx)>0.95,]
renxx[renxx==0] <- min(renxx[renxx>0])/2
sc_rxx <- as.data.frame(t(log2(renxx)))
save(sc_rxx2,file="input_data/sc_renxx2.RData")
