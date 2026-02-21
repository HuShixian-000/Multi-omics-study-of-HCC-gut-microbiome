rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/big_adonis/")
library(vegan)
library(ggplot2)
library(ggsci)
library(data.table)
library(foreach)
library(dplyr)
#source("/home/liuyuyao/bin/Microbiome.function.R")
#library(microbiome)
library(crayon)
library(nlme)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggprism)
set.seed(12345)
#setwd( "/groups/ProHuShiX/home/liuyuyao/big_adonis/V2_250106")

#==========load function===============
rmOutlier=function(x){
  tmp.sd=sd(x,na.rm=T)
  tmp.rm=c(x[x<median(x)-1*tmp.sd],x[x>median(x)+1*tmp.sd])
  x[x %in% tmp.rm]=NA
  return(x)
}

count_to_composition=function(x,samples_row=T){
  if (!samples_row){
    x=as.data.frame(t(x))
  } 
  x=as.data.frame(apply(x,1,function(x){
    x=x/sum(x)
    return(x)
  }))
  return(as.data.frame(t(x)))
}

transform_and_filter_taxa=function(x, samples_row=T, method="asin", missing_filter=0){
  x[x=="NA"]=0
  x[x=="NaN"]=0
  #if samples are in columns transpose
  if (!samples_row){
    
    x=as.data.frame(t(x))
    
  } 
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>100){
    
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 100") 
    
  }
  
  x_filt=x[,((colSums(x !=0) / nrow(x)) *100 )>missing_filter]
  my_num_removed=ncol(x)-ncol(x_filt)
  print (paste(my_num_removed, "species removed due to many missing values"))
  
  if (method=="asin"){
    print ("ASIN")
    x_filt=x_filt/100
    x_filt=asin(sqrt(x_filt))
  } else if (method=="log"){
    print ("LOG10")
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0])/2
    x_filt=x_filt+my_min
    x_filt=log10(x_filt)
  }else if (method=="clr"){
    print ("CLR")
    #Adapted from Alexander Kurilshikov 
    x_filt=x_filt/100
    #replace 0 by the half of the smallest value observed
    my_min=min(x[x>0])/2
    x=x+my_min
    #Calculate geometric mean
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x), na.rm=na.rm) / length(x))
    }
    Gmean_core = apply(x, 1, gm_mean)
    data_prepared = cbind(Gmean_core,x)
    d <- t(apply(data_prepared, 1, function(b) {
      log(b/b[1])[-1]
    }))
    x=d
    x_filt=x[,colnames(x) %in%colnames(x_filt)]
  }
  return(as.data.frame(x_filt))
}
library(glmnet)

lasso_fit=function(feature,outcome){
  if(length(colnames(feature))==0){
    result=data.frame(Protein=colnames(outcome),Feature="No",Lasso.beta=0)
  }else if(length(colnames(feature))==1){
    model=lm(outcome[,1]~ feature[,1])
    beta = summary(model)$coefficients[2]
    result=data.frame(Protein=colnames(outcome),Feature=colnames(feature),Lasso.beta=beta)
  }else{
    cv=cv.glmnet(as.matrix(feature),as.matrix(outcome), alpha = 1, nfolds = 5, type.measure="mse",standardize=T)
    beta <- (coef(cv, s = "lambda.min"))
    beta=as.data.frame(beta[-1,1])
    beta$variable=rownames(beta)
    colnames(beta)[1]="beta"
    result=data.frame(Protein=colnames(outcome),Feature=beta$variable,Lasso.beta=beta$beta)
  }
}

#=======MGS==================
data <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/species_raw_filter.txt")
uncorrect=as.data.frame(t(data))

##去掉低检出率
uncorrect=uncorrect[,colSums(uncorrect!=0) > (nrow(uncorrect) * 0.1) ]
##把菌踢掉之后再算一遍相对丰度
correct=apply(uncorrect,1,function(x){
  x=x/sum(x)
  return(x)
})
uncorrect <- as.data.frame(t(correct))
data <- uncorrect

data_clr <- transform_and_filter_taxa(data,samples_row = T,method = "clr",missing_filter = 0)

load("V2_250106/MGS_meta.RData")
MGS_meta[MGS_meta$sample=="BM424-0020283","BMI"] <- 21.887
MGS_meta[MGS_meta$sample=="BM424-0020167","BMI"] <- 18.75
MGS_meta[MGS_meta$sample=="BM424-0020264","BMI"] <- 27.99
MGS_meta[MGS_meta$sample=="BM424-0020274","BMI"] <- 25
MGS_meta[MGS_meta$sample=="BM424-0020454","BMI"] <- 19
MGS_meta$sample <- gsub("-",".",MGS_meta$sample)
MGS_meta <- MGS_meta[MGS_meta$sample%in%row.names(data),]

##剩下的人都没有酒精肝了，删掉
MGS_meta$Alcoholic_liver <- NULL
clinical <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/MGSinfo_0514.xlsx")
MGS_meta$AFP <- NULL
MGS_meta$Antibiotics <- NULL
MGS_meta$Cirrhosis <- NULL
# 保留Health
MGS_meta <- merge(MGS_meta,clinical[,c("sample","AFP","Cirrhosis","Antibiotics")],by="sample",all.x = T)
MGS_meta[MGS_meta$group=="Health",]$Antibiotics <- "No"

save(MGS_meta,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/442MGS_clinical.RData")

#开始循环
row.names(MGS_meta) <- MGS_meta[,1]
mgs_ad = foreach(i=1:ncol(MGS_meta),.combine = rbind) %do% {
  tmp.cov=MGS_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.taxa=data[rownames(data) %in% rownames(tmp.cov),]
  tmp.taxa=tmp.taxa[order(rownames(tmp.taxa)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  #tmp.cov[,1] <- as.factor(tmp.cov[,1])
  tmp.ad=adonis2(tmp.taxa ~ tmp.cov[,1], permutations = 999, method = "bray")
  cat(green(colnames(MGS_meta)[i],"\n"))
  return.string=data.frame(Cov=colnames(MGS_meta)[i],DF=tmp.ad[1,1],R2=tmp.ad[1,3],Pvalue=tmp.ad[1,5])
}
mgs_ad$Dataset="MGS"
write.table(mgs_ad,"mgs_ad1.txt",sep="\t",row.names = F,col.names = T,quote = F)
#save(MGS_meta,file="MGS_meta.RData")
test <- read.table("./mgs_ad.txt",header=T)

#========Pathway==============
# load("/groups/ProHuShiX/home/liuyuyao/humann/trim_mgs_path.RData")
mgs_path <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/humann_MGS_trim.txt")
mgs_path <- as.data.frame(t(mgs_path))
mgs_path_clr <-  transform_and_filter_taxa(mgs_path,samples_row = T,method = "clr",missing_filter = 0)

metacyc_ad = foreach(i=1:ncol(MGS_meta),.combine = rbind) %do% {
  tmp.cov=MGS_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.taxa=mgs_path [rownames(mgs_path ) %in% rownames(tmp.cov),]
  tmp.taxa=tmp.taxa[order(rownames(tmp.taxa)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  #tmp.cov[,1] <- as.factor(tmp.cov[,1])
  tmp.ad=adonis2(tmp.taxa ~ tmp.cov[,1] ,  permutations = 999, method = "bray")
  cat(green(colnames(MGS_meta)[i],"\n"))
  return.string=data.frame(Cov=colnames(MGS_meta)[i],DF=tmp.ad[1,1],R2=tmp.ad[1,3],Pvalue=tmp.ad[1,5])
}
metacyc_ad$Dataset="MetaCyc"
write.table(metacyc_ad,"metacyc_ad1.txt",sep="\t",row.names = F,col.names = T,quote = F)


#=======VF==============
VF <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/big_adonis/V2_250106/VF_merged_new.txt",header=T,check.names = F)
VF <- VF[,colnames(VF)%in%row.names(data)]

VF=as.data.frame(t(VF))

##去掉低检出率
VF=VF[,colSums(VF!=0) > (nrow(VF) * 0.1) ]
##把菌踢掉之后再算一遍相对丰度
correct=apply(VF,1,function(x){
  x=x/sum(x)
  return(x)
})
VF <- as.data.frame(t(correct))


VF_clr <- transform_and_filter_taxa(VF,samples_row = T,method = "clr",missing_filter = 0)


VF_ad = foreach(i=1:ncol(MGS_meta),.combine = rbind) %do% {
  tmp.cov=MGS_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.taxa=VF[rownames(VF) %in% rownames(tmp.cov),]
  tmp.taxa=tmp.taxa[order(rownames(tmp.taxa)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  #tmp.cov[,1] <- as.factor(tmp.cov[,1])
  tmp.ad=adonis2(tmp.taxa ~ tmp.cov[,1] ,  permutations = 999, method = "bray")
  cat(green(colnames(MGS_meta)[i],"\n"))
  return.string=data.frame(Cov=colnames(MGS_meta)[i],DF=tmp.ad[1,1],R2=tmp.ad[1,3],Pvalue=tmp.ad[1,5])
}
VF_ad$Dataset="VF"
write.table(VF_ad,"VF_ad1.txt",sep="\t",row.names = F,col.names = T,quote = F)

#========ARG==============
#load("/groups/ProHuShiX/home/liuyuyao/ARG/trim_ARG.RData")
ARG <- read.table("/groups/ProHuShiX/home/share/BIGMetaG_result/ARGs/ARGs_merged.txt",header=T,check.names = F,row.names = 1)
colnames(ARG) <- gsub("-",".",colnames(ARG))
ARG <- ARG[,colnames(ARG)%in%rownames(data)]
ARG <- ARG[rowSums(ARG>0)>0.1*ncol(ARG),]
ARG <- apply(ARG,2,function(x){
  y<-x/sum(x)
  return(y)
})
ARG <- as.data.frame(t(ARG))

#先修改一下名字，不然不能clr
ARG_name <- data.frame(name=colnames(ARG))
ARG_name$num <- paste("ARG",1:nrow(ARG_name),sep = "")

colnames(ARG) <- paste("ARG",1:ncol(ARG),sep = "")
ARG_name$Des <- gsub("^.*\\|","",ARG_name$name)
ARG_name$Des <- gsub("_.*$","",ARG_name$Des)
ARG_name$Species <- NA

ARG_clr=transform_and_filter_taxa(ARG,samples_row = T,method = "clr",missing_filter = 0)

ARG_ad = foreach(i=1:ncol(MGS_meta),.combine = rbind) %do% {
  tmp.cov=MGS_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.taxa=ARG[rownames(ARG) %in% rownames(tmp.cov),]
  tmp.taxa=tmp.taxa[order(rownames(tmp.taxa)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  #tmp.cov[,1] <- as.factor(tmp.cov[,1])
  tmp.ad=adonis2(tmp.taxa ~ tmp.cov[,1] ,  permutations = 999, method = "bray")
  cat(green(colnames(MGS_meta)[i],"\n"))
  return.string=data.frame(Cov=colnames(MGS_meta)[i],DF=tmp.ad[1,1],R2=tmp.ad[1,3],Pvalue=tmp.ad[1,5])
}
ARG_ad$Dataset="ARG"
write.table(ARG_ad,"ARG_ad1.txt",sep="\t",row.names = F,col.names = T,quote = F)

#=======Metabolism=================
library(readxl)
set.seed(12345)
rawdata <- read_excel("/groups/ProHuShiX/home/liuyuyao/metabolism_1.30/filtered_data.xlsx")
load("/groups/ProHuShiX/home/liuyuyao/metabolism_1.30/WGCNA/Health+HCC_V2/HCC+Health_only_sc_old.RData")
bigmeta <- as.data.frame(rawdata[,colnames(rawdata)%in%row.names(sc)])


#==================================
#插补
#==================================
#install.packages("DMwR2")
library(DMwR2)
raw <- bigmeta
raw[raw==0]=NA
#看了一下函数给的参考data,是行为样本列为变量
#raw_knn <- knnImputation(as.data.frame(t(raw)),k=5,scale = T, meth ='weighAvg',distData = NULL)
#看了下在所有代谢物都能检出的样本小于5
#raw_knn <- as.data.frame(t(raw_knn))
raw[raw==0]=1e-10

#====================================
#mean-centering
#====================================
##均值归一化
for(i in 1:ncol(raw)){
  a <- raw[,i]
  a <- mean(a)
  raw[,i]<- sapply(raw[,i], function(x)y <- x/a)
}

#===========================================================================
##log2转换
#===========================================================================
raw<- log2(raw)
sc <- as.data.frame(t(raw))
save(sc,file="sc_BIGMeta_all.RData")

load("/groups/ProHuShiX/home/liuyuyao/big_adonis/metab_meta.RData")

metab_meta[metab_meta$sample=="X1199_J","BMI"] <- 21.887
metab_meta[metab_meta$sample=="X1626_J","BMI"] <- 18.75
metab_meta[metab_meta$sample=="X3408_J","BMI"] <- 27.99
metab_meta[metab_meta$sample=="X3542_J","BMI"] <- 25
metab_meta <- metab_meta[metab_meta$sample%in%row.names(sc),]


metab_ad = foreach(i=1:ncol(metab_meta),.combine = rbind) %do% {
  tmp.cov=metab_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.sc=sc[rownames(sc) %in% rownames(tmp.cov),]
  tmp.sc=tmp.sc[order(rownames(tmp.sc)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  #tmp.cov[,1] <- as.factor(tmp.cov[,1])
  tmp.ad=adonis(tmp.sc ~ tmp.cov[,1] ,  permutations = 999, method = "euclidian")
  cat(green(colnames(metab_meta)[i],"\n"))
  return.string=data.frame(Cov=colnames(metab_meta)[i],DF=tmp.ad$aov.tab[1,1],R2=tmp.ad$aov.tab[1,5],Pvalue=tmp.ad$aov.tab[1,6])
}
metab_ad $Dataset="Metabolism"
write.table(metab_ad,"metab_all_ad.txt",sep="\t",row.names = F,col.names = T,quote = F)
###代谢不用画

#=======用林福师兄的green genes画16s==========
data_16S_genus=read.table("/groups/ProHuShiX/home/liuyuyao/BIGMeta_16S/translocate/PCA/Greengenes2.16S.genus.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t",check.names = F)
colnames(data_16S_genus) <- gsub("963_Trelative_abundance","963_T",colnames(data_16S_genus))
colnames(data_16S_genus) <- gsub("2142","2120",colnames(data_16S_genus))
colnames(data_16S_genus) <- gsub("3775","3775/RH29",colnames(data_16S_genus))


load("/groups/ProHuShiX/home/liuyuyao/big_adonis/16s_meta.RData")
row.names(s16_meta) <- gsub("3755","3755/RH29",row.names(s16_meta))
row.names(s16_meta) <- gsub("2142","2120",row.names(s16_meta))

#准备16S的信息表
group_16s <- as.data.frame(colnames(data_16S_genus))
colnames(group_16s) <- "sample"
group_16s$sampleID <- gsub("_.*$","",group_16s$sample)
group_16s$type <- ifelse(group_16s$sample%like%"T","Tumor",ifelse(group_16s$sample%like%"P","Peri","Liver"))


#踢掉非HCC和两个只有癌旁的
library(readxl)
no_HCC <- read_excel("/groups/ProHuShiX/home/liuyuyao/BIGMeta_16S/tissue_0826/大样本队列非HCC患者.xlsx",sheet = "16s中非HCC患者")
group_16s <- group_16s[!group_16s$sampleID%in%no_HCC$术中标本编号,]
group_16s <- group_16s[!group_16s$sampleID%in%c("1944","1466"),]

data_16S_genus <- data_16S_genus[,colnames(data_16S_genus)%in%group_16s$sample]
data_16S_genus=as.data.frame(t(data_16S_genus))
##去污染
contamination=read.table("/groups/ProHuShiX/home/liuyuyao/BIGMeta_16S/translocate/PCA/Known.contaminates.list.txt",header = F,stringsAsFactors = F)
return.string <- data.frame()
contaminated=foreach(i=1:(ncol(data_16S_genus)),.combine = rbind) %do%  {
  tmp.feature=colnames(data_16S_genus)[i]
  for(n in 1:nrow(contamination)){
    tmp.conta=contamination$V1[n]
    if(tmp.feature %like% tmp.conta){
      return.string=data.frame(taxa=tmp.feature,Contaminated=tmp.conta)
    }
  }
  return.string=return.string
}
contaminated=unique(contaminated$taxa)
data_16S_genus=data_16S_genus[,!colnames(data_16S_genus) %in% contaminated]
data_16S_genus=data_16S_genus[,!colnames(data_16S_genus) %like% "uncultured"]
data_16S_genus=data_16S_genus[,!colnames(data_16S_genus) %like% "Mitochondria"]
#去掉没有注释到属的
data_16S_genus <- as.data.frame(t(data_16S_genus))
data_16S_genus<- data_16S_genus %>%dplyr::mutate(genus=str_extract(row.names(data_16S_genus),"g__.*"))
data_16S_genus=data_16S_genus[!data_16S_genus$genus=="g__",]
row.names(data_16S_genus) <- data_16S_genus$genus
data_16S_genus$genus <- NULL
row.names(data_16S_genus) <- gsub("g__","",row.names(data_16S_genus))

# remove low count(剔样本？)
data_16S_genus <- as.data.frame(t(data_16S_genus))
count=data.frame(count=rowSums(data_16S_genus))
data_16S_genus=data_16S_genus[rownames(data_16S_genus) %in% rownames(count)[count$count>1000],]

# remove low presence and abundance
data_16S_genus=data_16S_genus[,colSums(data_16S_genus>0)>(nrow(data_16S_genus)*0.2)]

#最后剩下的属太少了，这些属在有些样本中都检测不到，先踢掉
count=data.frame(count=rowSums(data_16S_genus))
data_16S_genus=data_16S_genus[rownames(data_16S_genus) %in% rownames(count)[count$count>0],]

# re-calucalte proportion
data_16S_genus <- count_to_composition(data_16S_genus)
range(colMeans(data_16S_genus))


#========重新整理16s信息表=====================
HCCinfo02APR2024 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical/HCCinfo02APR2024.xlsx")
HCCinfo02APR2024 <- HCCinfo02APR2024[HCCinfo02APR2024$术中标本编号%in%row.names(s16_tumor),]


s16_meta <- HCCinfo02APR2024[,c("术中标本编号" ,  "性别"  , "年龄"  ,"身高" ,                     
                                "体重" , "吸烟史" ,                  
                                "饮酒史" ,  "酒精肝" ,                 
                                 "脂肪肝"  ,"总胆红素"  , "谷丙转氨酶 ALT" ,           
                                "谷草转氨酶 AST"     , "HBsAg"   ,   "HBV-DNA" ,                  
                                "AFP"   ,  "影像简表+病理BCLC" ,  "肿瘤个数" ,                 
                                "肿瘤大小" ,                  "最早抗生素时间距离手术天数")]



ATB <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical/HCCinfo02APR2024.xlsx", 
             sheet = "抗生素")

ATB$time <- gsub("JAN","-01-",ATB$医嘱时间)
ATB$time <- gsub("FEB","-02-",ATB$time)
ATB$time <- gsub("MAR","-03-",ATB$time)
ATB$time <- gsub("APR","-04-",ATB$time)
ATB$time <- gsub("MAY","-05-",ATB$time)
ATB$time <- gsub("JUN","-06-",ATB$time)
ATB$time <- gsub("JUL","-07-",ATB$time)
ATB$time <- gsub("AUG","-08-",ATB$time)
ATB$time <- gsub("SEP","-09-",ATB$time)
ATB$time <- gsub("OCT","-10-",ATB$time)
ATB$time <- gsub("NOV","-11-",ATB$time)
ATB$time <- gsub("DEC","-12-",ATB$time)
ATB$time <- gsub(":.*$","",ATB$time)

ATB$date <- as.Date(ATB$time,"%d-%m-%Y")
#把ATB分为有标本号的和只有住院号的
ATB1 <- subset(ATB,is.na(ATB$tid)==F)
ATB2 <- subset(ATB,is.na(ATB$tid)==T)


library(readxl)
biaoben <- read_excel("/groups/ProHuShiX/home/liuyuyao/metaG_clinical/标本非trial17MAY2023-宏基因组.xlsx")

library(tidyverse)
library(stringi)
biaoben <- biaoben%>%select("术中标本编号","ltname","lthosid","粪便采集日期")
colnames(biaoben) <- c("sample","name","patientid","select_time")
biaoben$select_time <- as.Date(biaoben$select_time,"%Y-%m-%d")
ATB1 <- merge(ATB1,biaoben,by.y = "sample",by.x="tid",all.x=T)
library(lubridate)
ATB1$time_result <- ATB1$date-ATB1$select_time
ATB1$result <- ifelse(ATB1$time_result>0,"No",ifelse(is.na(ATB1$time_result)==T,NA,"Yes"))

ATB2 <- merge(ATB2,biaoben,by.y = "patientid",by.x="住院号",all.x=T)
library(lubridate)
ATB2$time_result <- ATB2$date-ATB2$select_time
ATB2$result <- ifelse(ATB2$time_result>0,"No",ifelse(is.na(ATB2$time_result)==T,NA,"Yes"))
##ATB2没有能merge上的

ATB2$tid <- "-"
ATB1 <- ATB1%>%arrange(tid,time_result)
ATB1 <- ATB1[!duplicated(ATB1$tid),]


s16_meta <- merge(s16_meta,ATB1[,c("tid","result")],by.x="术中标本编号",by.y="tid",all.x=T)
s16_meta$BMI <- s16_meta$体重/((s16_meta$身高/100)^2)
s16_meta$身高 <- NULL
s16_meta$体重 <- NULL
s16_meta$酒精肝 <- NULL
s16_meta$`HBV-DNA` <- NULL
s16_meta[s16_meta=="-3"] <- NA
s16_meta[s16_meta=="."] <- NA
s16_meta[s16_meta=="-1"] <- NA
s16_meta[s16_meta=="不详"] <- NA
s16_meta$最早抗生素时间距离手术天数 <- NULL
s16_meta$总胆红素 <- NULL
s16_meta$肿瘤大小 <- NULL
s16_meta$肿瘤个数 <- NULL
colnames(s16_meta) <- c( "sample" ,"Gender" ,"Age","Smoke",           
                         "Drink" ,  "Fatty_liver" ,   "ALT"  , 
                         "AST" , "HBsAg", "AFP"   , "BCLC",
                          "Antibiotics","BMI"  )

s16_meta[s16_meta$sample=="1199","BMI"] <- 21.887
s16_meta[s16_meta$sample=="1626","BMI"] <- 18.75
s16_meta[s16_meta$sample=="3408","BMI"] <- 27.99
s16_meta[s16_meta$sample=="3542","BMI"] <- 25
s16_meta[s16_meta$sample=="1788","AFP"] <- 4639.2
s16_meta[s16_meta$sample=="1789","AFP"] <- 2765.36
s16_meta[s16_meta$sample=="1801","AFP"] <- 5.21
s16_meta[s16_meta$sample=="1830","AFP"] <- 205.06
s16_meta[s16_meta$sample=="1871","AFP"] <- 851.27
s16_meta[s16_meta$sample=="2043","AFP"] <- 21318.08
s16_meta[s16_meta$sample=="3560","AFP"] <- 371.81
s16_meta[s16_meta$sample=="3694","AFP"] <- 1.9
s16_meta[s16_meta$sample=="3744","AFP"] <- 1.9
s16_meta[s16_meta$sample=="3784","AFP"] <- 3.3
s16_meta[s16_meta$sample=="3794","AFP"] <- 1.9

nrow(s16_meta[is.na(s16_meta$Antibiotics)==T,])
#382
library(readxl)
 Cirrhosis <- read_excel("Cirrhosis_res.xlsx")
 Cirrhosis$sample <- gsub("2142","2120",Cirrhosis$sample)
s16_meta <- merge(s16_meta,Cirrhosis[,c("sample","Cirrhosis")],by="sample",all.x=T)




#========16S肿瘤===============
s16_tumor <- data_16S_genus[row.names(data_16S_genus)%in%group_16s[group_16s$type=="Tumor",]$sample,]
s16_tumor$sample <- gsub("_.*$","",row.names(s16_tumor))
s16_tumor <- aggregate(s16_tumor[,1:(ncol(s16_tumor)-1)],by=list(s16_tumor$sample),mean)

row.names(s16_tumor) <- s16_tumor$Group.1
s16_tumor$Group.1 <- NULL

s16_meta <- s16_meta[!duplicated(s16_meta$sample),]
row.names(s16_meta) <- s16_meta[,1]
s16_tumor_ad = foreach(i=1:ncol(s16_meta),.combine = rbind) %do% {
  tmp.cov=s16_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.taxa=s16_tumor[rownames(s16_tumor) %in% rownames(tmp.cov),]
  tmp.taxa=tmp.taxa[order(rownames(tmp.taxa)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  #tmp.cov[,1] <- as.factor(tmp.cov[,1])
  tmp.ad=adonis(tmp.taxa ~ tmp.cov[,1] ,  permutations = 999, method = "bray")
  cat(green(colnames(s16_meta)[i],"\n"))
  return.string=data.frame(Cov=colnames(s16_meta)[i],DF=tmp.ad$aov.tab[1,1],R2=tmp.ad$aov.tab[1,5],Pvalue=tmp.ad$aov.tab[1,6])
}
s16_tumor_ad[15,1] <- "group"
s16_tumor_ad$Dataset="16s_tumor"
write.table(s16_tumor_ad,"16s_tumor_ad.txt",sep="\t",row.names = F,col.names = T,quote = F)
save(s16_meta,file="s16_meta.RData")

#========16S非肿瘤===============
s16_nontumor <- data_16S_genus[row.names(data_16S_genus)%in%group_16s[!group_16s$type=="Tumor",]$sample,]
s16_nontumor$sample <- gsub("_.*$","",row.names(s16_nontumor))
s16_nontumor <- aggregate(s16_nontumor[,1:(ncol(s16_nontumor)-1)],by=list(s16_nontumor$sample),mean)

row.names(s16_nontumor) <- s16_nontumor$Group.1
s16_nontumor$Group.1 <- NULL


s16_nontumor_ad = foreach(i=1:ncol(s16_meta),.combine = rbind) %do% {
  tmp.cov=s16_meta[,i,drop=F]
  tmp.cov=na.omit(tmp.cov)
  tmp.taxa=s16_nontumor[rownames(s16_nontumor) %in% rownames(tmp.cov),]
  tmp.cov <- tmp.cov[row.names(tmp.cov)%in%row.names(tmp.taxa),,drop=F]
  tmp.taxa=tmp.taxa[order(rownames(tmp.taxa)),]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  #tmp.cov[,1] <- as.factor(tmp.cov[,1])
  tmp.ad=adonis(tmp.taxa ~ tmp.cov[,1] ,  permutations = 999, method = "bray")
  cat(green(colnames(s16_meta)[i],"\n"))
  return.string=data.frame(Cov=colnames(s16_meta)[i],DF=tmp.ad$aov.tab[1,1],R2=tmp.ad$aov.tab[1,5],Pvalue=tmp.ad$aov.tab[1,6])
}
s16_nontumor_ad[15,1] <- "group"
s16_nontumor_ad$Dataset="16s_peri"
write.table(s16_nontumor_ad,"16s_nontumor_ad.txt",sep="\t",row.names = F,col.names = T,quote = F)

