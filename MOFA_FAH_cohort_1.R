library(data.table)
library(vegan)
library(ggplot2)
library(edgeR)
library(limma)
library(circlize)
library(biomaRt)
library(cowplot)
library(dplyr)
library(foreach)
library(tidyverse)
library(MOFA2)

#============load function================
lasso_fit=function(feature,outcome){
  if(length(colnames(feature))==0){
    result=data.frame(Protein=colnames(outcome),Feature="No",Lasso.beta=0)
  }else if(length(colnames(feature))==1){
    model=lm(outcome[,1]~ feature[,1])
    beta = summary(model)$coefficients[2]
    result=data.frame(Protein=colnames(outcome),Feature=colnames(feature),Lasso.beta=beta)
  }else{
    cv=cv.glmnet(as.matrix(feature),as.matrix(outcome), alpha = 1, nfolds = 10, type.measure="mse",standardize=T)
    beta <- (coef(cv, s = "lambda.min"))
    beta=as.data.frame(beta[-1,1])
    beta$variable=rownames(beta)
    colnames(beta)[1]="beta"
    result=data.frame(Protein=colnames(outcome),Feature=beta$variable,Lasso.beta=beta$beta)
  }
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


#=======整理具有匹配组学的样本=======================
setwd("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/原始文件MOFAbig")
data_phenotype=readxl::read_xlsx('./MGSinfo_250310.xlsx')
data_phenotype$sampleID[data_phenotype$sampleID=="2120"]=2142


data_div=read.table("./Balance.score1.txt",sep = "\t",stringsAsFactors = F,header = T)
data_div$Row_names=str_replace_all(data_div$Row_names,"\\.","-")

data_div=data_div %>% dplyr::select(Row_names,balance_value)

data_phenotype=left_join(data_div,data_phenotype,by=c('Row_names'='sample'))


data_phenotype=data_phenotype %>% dplyr::filter(group=='HCC')




#data_coupling=read.table("./MOFA_Hu0106/MOFA_big0112/所有匹配宏基因组的RNA样本.csv",header = T,stringsAsFactors = F,sep = ",")
#和我整理的版本有些样本差别，用我的样本版本
load('./匹配宏基因组的转录组临床信息与表达值1203.RData')

data_coupling=rbind(clinical_p,clinical_t %>% dplyr::select(-group3))


data_rnaNorm=read.table("./Data.RNA.all.rpkm.txt",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
raw_data=data_rnaNorm
data_rnaNorm=raw_data
data_coupling=data_coupling[data_coupling$patient %in% data_phenotype$sampleID,] #删除了一个患者3743的3个样本 从502到499

length(unique(data_coupling$patient))
data_rnaNorm=data_rnaNorm[rownames(data_rnaNorm) %in% data_coupling$sampleID,]
data_coupling=data_coupling[data_coupling$sampleID %in% rownames(data_rnaNorm),]
data_rnaNorm=data_rnaNorm[order(rownames(data_rnaNorm)),]
data_coupling=data_coupling[order(data_coupling$sampleID),]
table(data_coupling$sampleID==rownames(data_rnaNorm))
# data_rnaNorm=data_rnaNorm[,colSums(data_rnaNorm>0)>(nrow(data_rnaNorm)*0.8)]
data_rnaNorm=log(data_rnaNorm+0.1)
data_rnaNorm=data_rnaNorm[,!colnames(data_rnaNorm) %like% "MT-"]

data_coupling$group=as.factor(data_coupling$group)
data_coupling <- within(data_coupling, group <- relevel(group, ref = "g5840014"))
stopifnot(rownames(data_rnaNorm)==data_coupling$sampleID)
data_rnaNorm=apply(data_rnaNorm,2,function(x){
  resids=resid(lm(x~data_coupling$group))
  return(resids)
})
data_rnaNorm=as.data.frame(data_rnaNorm)
data_coupling=data_coupling[data_coupling$sampleID %in% rownames(data_rnaNorm),]

# rnaseq, non-tumor
data_coupling_nontumor=data_coupling %>% dplyr::filter(sampleID %in% clinical_p$sampleID)

# rnaseq, tumor
data_coupling_tumor=data_coupling %>% dplyr::filter(sampleID %in% clinical_t$sampleID)

# metabolites
meta_data=read.table("./Mean.txt",header = T,check.names = F,stringsAsFactors = F,fill = T,sep = "\t")
meta_annot=read.table("./Metabolite Mapping.txt",header = T,sep = "\t",fill=T)
meta_coupling=read.table("./sample_clientID.txt",header = T,stringsAsFactors = F)
meta_coupling$sample=gsub("_J","",meta_coupling$sample)
colnames(meta_data)=gsub("_J","",colnames(meta_data))
meta_coupling$sampleID[meta_coupling$group!="HCC"]=meta_coupling$sample[meta_coupling$group!="HCC"]
meta_data=meta_data[meta_data$id %in% meta_annot$id,]
rownames(meta_data)=meta_data$id
meta_data=meta_data[,colnames(meta_data) %in% meta_coupling$sample]
meta_data=as.data.frame(t(meta_data))
meta_data[meta_data==0]=min(meta_data[meta_data>0])/2
meta_data=log2(meta_data)

# metagenomic path
mgs_path=read.table("./pathabundance_all.tsv",header = T,stringsAsFactors = F,row.names = 1,sep="\t",fill = T,comment.char = "", check.names = FALSE)
colnames(mgs_path)=gsub("_Abundance","",colnames(mgs_path))
colnames(mgs_path)=gsub("\\.","-",colnames(mgs_path))
mgs_path=mgs_path[!rownames(mgs_path) %like% "\\|",]
mgs_path=mgs_path[rowSums(mgs_path>0)>(ncol(mgs_path)*0.1),]
mgs_path=mgs_path[!rownames(mgs_path) %in% c("UNMAPPED","UNINTEGRATED"),]
mgs_path=as.data.frame(apply(mgs_path,2,function (x){
  x=x/sum(x)
  return(x)
}))
mgs_path_clr=transform_and_filter_taxa(mgs_path,samples_row = F,method = "clr",missing_filter = 0)



tmp.path=mgs_path_clr[rownames(mgs_path_clr) %in% data_phenotype$Row_names,]
rownames(tmp.path)==data_phenotype$Row.names
tmp.path=as.data.frame(tmp.path)
rownames(tmp.path)=data_phenotype$sampleID

tmp.meta=meta_data[rownames(meta_data) %in% meta_coupling$sample,]
tmp.meta=tmp.meta[order(rownames(tmp.meta)),]
table(meta_coupling=meta_coupling[order(meta_coupling$sample),])
rownames(tmp.meta)==meta_coupling$sample
rownames(tmp.meta)=meta_coupling$sampleID

tmp.score=data_phenotype[,c("sampleID","balance_value")] %>% as.data.frame()
tmp.meta=tmp.meta[rownames(tmp.meta) %in% tmp.score$sampleID,] ## 442 个患者 和代谢取了交集 只有129个患者了
tmp.score=tmp.score[tmp.score$sampleID %in% rownames(tmp.meta),]
tmp.meta=tmp.meta[order(rownames(tmp.meta)),]
tmp.score=tmp.score[order(tmp.score$sampleID),]
rownames(tmp.score)=tmp.score$sampleID
tmp.score$sampleID=NULL
stopifnot(rownames(tmp.score)==rownames(tmp.meta))


tmp.nontumor=data_rnaNorm[rownames(data_rnaNorm) %in% data_coupling_nontumor$sampleID,]
tmp.nontumor=tmp.nontumor[order(rownames(tmp.nontumor)),]
data_coupling_nontumor=data_coupling_nontumor[order(data_coupling_nontumor$sampleID),]
table(rownames(tmp.nontumor)==data_coupling_nontumor$sampleID)
data_coupling_nontumor$patient.multi=make.unique(as.character(data_coupling_nontumor$patient))
rownames(tmp.nontumor)=data_coupling_nontumor$patient.multi

tmp.tumor=data_rnaNorm[rownames(data_rnaNorm) %in% data_coupling_tumor$sampleID,]
tmp.tumor=tmp.tumor[order(rownames(tmp.tumor)),]
data_coupling_tumor=data_coupling_tumor[order(data_coupling_tumor$sampleID),]
table(rownames(tmp.tumor)==data_coupling_tumor$sampleID)
data_coupling_tumor$patient.multi=make.unique(as.character(data_coupling_tumor$patient))
rownames(tmp.tumor)=data_coupling_tumor$patient.multi

tmp.patient=intersect(rownames(tmp.score),intersect(data_coupling_tumor$patient,data_coupling_nontumor$patient)) #129 个患者都有RNA样本




data_coupling_tumor=data_coupling_tumor[data_coupling_tumor$patient %in% tmp.patient,] #现在只有180个样本
data_coupling_nontumor=data_coupling_nontumor[data_coupling_nontumor$patient %in% tmp.patient,] #现在只有270个样本
tmp.score=tmp.score[rownames(tmp.score) %in% tmp.patient,,drop=F]
tmp.tumor=tmp.tumor[rownames(tmp.tumor) %in% data_coupling_tumor$patient.multi,]
tmp.nontumor=tmp.nontumor[rownames(tmp.nontumor) %in% data_coupling_nontumor$patient.multi,]

tmp.score=tmp.score[order(rownames(tmp.score)),,drop=F]
data_coupling_tumor=data_coupling_tumor[order((data_coupling_tumor$patient)),,drop=F]
data_coupling_nontumor=data_coupling_nontumor[order((data_coupling_nontumor$patient)),,drop=F]
tmp.tumor=tmp.tumor[order(rownames(tmp.tumor)),]
tmp.nontumor=tmp.nontumor[order(rownames(tmp.nontumor)),]

stopifnot(((rownames(tmp.score) == unique(data_coupling_tumor$patient))))
stopifnot(((rownames(tmp.score) == unique(data_coupling_nontumor$patient))))

stopifnot(((rownames(tmp.nontumor) == data_coupling_nontumor$patient.multi)))
stopifnot(((rownames(tmp.tumor) == data_coupling_tumor$patient.multi)))









var_meta=foreach(i=1:ncol(tmp.meta),.combine = rbind) %do%  {
  tmp.feature=colnames(tmp.meta)[i]
  tmp.var=var(tmp.meta[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_path=foreach(i=1:ncol(tmp.path),.combine = rbind) %do%  {
  tmp.feature=colnames(tmp.path)[i]
  tmp.var=var(tmp.path[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_tumor=foreach(i=1:ncol(tmp.tumor),.combine = rbind) %do%  {
  tmp.feature=colnames(tmp.tumor)[i]
  tmp.var=var(tmp.tumor[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_nontumor=foreach(i=1:ncol(tmp.nontumor),.combine = rbind) %do%  {
  tmp.feature=colnames(tmp.nontumor)[i]
  tmp.var=var(tmp.nontumor[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}





quantile(var_meta$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_path$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_tumor$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_nontumor$var, probs = c(0,0.25,0.5,0.85,1)) 

##根据变异系数进行筛选纳入的通路和基因

tmp.meta=tmp.meta[,colnames(tmp.meta) %in% var_meta$feature[var_meta$var>quantile(var_meta$var,0.85)]]
tmp.path=tmp.path[,colnames(tmp.path) %in% var_path$feature[var_path$var>quantile(var_path$var,0.85)]]
tmp.tumor=tmp.tumor[,colnames(tmp.tumor) %in% var_tumor$feature[var_tumor$var>quantile(var_tumor$var,0.85)]]
tmp.nontumor=tmp.nontumor[,colnames(tmp.nontumor) %in% var_nontumor$feature[var_nontumor$var>quantile(var_nontumor$var,0.85)]]



tmp.meta1=merge(data_coupling_tumor,tmp.meta,by.x="patient",by.y="row.names",all=F)
tmp.path1=merge(data_coupling_tumor,tmp.path,by.x="patient",by.y="row.names",all=F)
data_coupling_tumor=data_coupling_tumor[order(data_coupling_tumor$patient.multi),]
tmp.tumor=tmp.tumor[order(rownames(tmp.tumor)),]
stopifnot(((data_coupling_tumor$patient.multi == rownames(tmp.tumor))))
tmp.tumor=tmp.tumor[rownames(tmp.tumor) %in% tmp.meta1$patient.multi,]
tmp.path1=tmp.path1[tmp.path1$patient.multi %in% tmp.meta1$patient.multi,]
tmp.meta1=tmp.meta1[order(tmp.meta1$patient.multi),]
tmp.path1=tmp.path1[order(tmp.path1$patient.multi),]

stopifnot(tmp.meta1$patient.multi == rownames(tmp.tumor))
stopifnot(tmp.path1$patient.multi == rownames(tmp.tumor))


rownames(tmp.meta1)=tmp.meta1$patient.multi;tmp.meta1=tmp.meta1[,-1:-6];tmp.meta1$patient.multi=NULL
rownames(tmp.path1)=tmp.path1$patient.multi;tmp.path1=tmp.path1[,-1:-6];tmp.path1$patient.multi=NULL
MOFA_data1=list(as.matrix(t(tmp.path1)),as.matrix(t(tmp.meta1)),as.matrix(t(tmp.tumor)))

write.table(data_coupling_tumor,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/big_mofa_rnasample.txt',sep = "\t",row.names = T)

#=========reload data===================
rm(list=ls())
gc()
setwd("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/原始文件MOFAbig")
data_phenotype=readxl::read_xlsx('./MGSinfo_250310.xlsx')
data_phenotype$sampleID[data_phenotype$sampleID=="2120"]=2142


data_div=read.table("./Balance.score1.txt",sep = "\t",stringsAsFactors = F,header = T)
data_div$Row_names=str_replace_all(data_div$Row_names,"\\.","-")

data_div=data_div %>% dplyr::select(Row_names,balance_value)

data_phenotype=left_join(data_div,data_phenotype,by=c('Row_names'='sample'))


data_phenotype=data_phenotype %>% dplyr::filter(group=='HCC')




#data_coupling=read.table("./MOFA_Hu0106/MOFA_big0112/所有匹配宏基因组的RNA样本.csv",header = T,stringsAsFactors = F,sep = ",")
#和我整理的版本有些样本差别，用我的样本版本
load('./匹配宏基因组的转录组临床信息与表达值1203.RData')

data_coupling=rbind(clinical_p,clinical_t %>% dplyr::select(-group3))


data_rnaNorm=read.table("./Data.RNA.all.rpkm.txt",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
raw_data=data_rnaNorm
data_rnaNorm=raw_data
sample_list <- read.table("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/big_mofa_rnasample.txt",header=T)
data_coupling=sample_list #直接用已经挑好的有匹配的样本

length(unique(data_coupling$patient))
data_rnaNorm=data_rnaNorm[rownames(data_rnaNorm) %in% data_coupling$sampleID,]
data_coupling=data_coupling[data_coupling$sampleID %in% rownames(data_rnaNorm),]
# data_coupling <- data_coupling[data_coupling$Sample%in%sample_list$Sample,]
data_rnaNorm=data_rnaNorm[rownames(data_rnaNorm) %in% data_coupling$sampleID,]
data_rnaNorm=data_rnaNorm[order(rownames(data_rnaNorm)),]
data_coupling=data_coupling[order(data_coupling$sampleID),]
table(data_coupling$sampleID==rownames(data_rnaNorm))
data_rnaNorm=data_rnaNorm[,colSums(data_rnaNorm>0)>(nrow(data_rnaNorm)*0.8)]
data_rnaNorm=log(data_rnaNorm+0.1)
data_rnaNorm=data_rnaNorm[,!colnames(data_rnaNorm) %like% "MT-"]


data_coupling$group=as.factor(data_coupling$group)
data_coupling <- within(data_coupling, group <- relevel(group, ref = "g5840014"))
stopifnot(rownames(data_rnaNorm)==data_coupling$sampleID)
data_rnaNorm=apply(data_rnaNorm,2,function(x){
  resids=resid(lm(x~data_coupling$group))
  return(resids)
})
data_rnaNorm=as.data.frame(data_rnaNorm)
data_coupling=data_coupling[data_coupling$sampleID %in% rownames(data_rnaNorm),]

# rnaseq, non-tumor
# data_coupling_nontumor=data_coupling %>% dplyr::filter(sampleID %in% clinical_p$sampleID)


# rnaseq, tumor
data_coupling_tumor=data_coupling %>% dplyr::filter(sampleID %in% clinical_t$sampleID)

# metabolites
# meta_data=read.table("./Mean.txt",header = T,check.names = F,stringsAsFactors = F,fill = T,sep = "\t")
# meta_annot=read.table("./Metabolite Mapping.txt",header = T,sep = "\t",fill=T)
# meta_coupling=read.table("./sample_clientID.txt",header = T,stringsAsFactors = F)
# meta_coupling$sample=gsub("_J","",meta_coupling$sample)
# colnames(meta_data)=gsub("_J","",colnames(meta_data))

load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/淑芬结果/sc_BIGMeta250618.RData")
meta_data <- sc
row.names(meta_data) <- gsub("\\.","-",row.names(meta_data))
##修改行名为sampleID
for (i in 1:nrow(meta_data)){
  x <- row.names(meta_data)[i]
  row.names(meta_data)[i] <- ifelse(nrow(data_phenotype[data_phenotype$Row_names==x,])==0,x,data_phenotype[data_phenotype$Row_names==x,]$sampleID)
}

# meta_coupling <- meta_coupling[meta_coupling$sampleID%in%sample_list$patient,]
# meta_coupling$sampleID[meta_coupling$group!="HCC"]=meta_coupling$sample[meta_coupling$group!="HCC"]
meta_data=meta_data[row.names(meta_data)%in%sample_list$patient,]
# rownames(meta_data)=meta_data$id
# meta_data=meta_data[,colnames(meta_data) %in% meta_coupling$sample]
# meta_data=as.data.frame(t(meta_data))
# meta_data[meta_data==0]=min(meta_data[meta_data>0])/2
# meta_data=log2(meta_data)

# metagenomic path
mgs_path=read.table("./pathabundance_all.tsv",header = T,stringsAsFactors = F,row.names = 1,sep="\t",fill = T,comment.char = "", check.names = FALSE)
colnames(mgs_path)=gsub("_Abundance","",colnames(mgs_path))
colnames(mgs_path)=gsub("\\.","-",colnames(mgs_path))
score <- read.table("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/淑芬结果/Balance.score1.txt",header=T)
score$Row.names <- gsub("\\.","-",score$Row.names)
## 取大样本队列最终纳入的样本
mgs_path <- mgs_path[,colnames(mgs_path)%in%score$Row.names]
mgs_path=mgs_path[!rownames(mgs_path) %like% "\\|",]
mgs_path=mgs_path[rowSums(mgs_path>0)>(ncol(mgs_path)*0.1),]
mgs_path=mgs_path[!rownames(mgs_path) %in% c("UNMAPPED","UNINTEGRATED"),]
mgs_path=as.data.frame(apply(mgs_path,2,function (x){
  x=x/sum(x)
  return(x)
}))
mgs_path_clr=transform_and_filter_taxa(mgs_path,samples_row = F,method = "clr",missing_filter = 0)

##修改列名为sampleID
for (i in 1:nrow(mgs_path_clr)){
  x <- row.names(mgs_path_clr)[i]
  row.names(mgs_path_clr)[i] <- ifelse(nrow(data_phenotype[data_phenotype$Row_names==x,])==0,x,data_phenotype[data_phenotype$Row_names==x,]$sampleID)
}


tmp.path=mgs_path_clr[rownames(mgs_path_clr) %in% sample_list$patient,]
tmp.path <- tmp.path[order(rownames(tmp.path)),]
tmp.meta=meta_data
tmp.meta=tmp.meta[order(rownames(tmp.meta)),]
# table(meta_coupling=meta_coupling[order(meta_coupling$sample),])
# rownames(tmp.meta)==meta_coupling$sample
# rownames(tmp.meta)=meta_coupling$sampleID

tmp.score=data_phenotype[,c("sampleID","balance_value")] %>% as.data.frame()
tmp.meta=tmp.meta[rownames(tmp.meta) %in% tmp.score$sampleID,] ## 442 个患者 和代谢取了交集 只有129个患者了
tmp.score=tmp.score[tmp.score$sampleID %in% rownames(tmp.meta),]
tmp.meta=tmp.meta[order(rownames(tmp.meta)),]
tmp.score=tmp.score[order(tmp.score$sampleID),]
rownames(tmp.score)=tmp.score$sampleID
tmp.score$sampleID=NULL
stopifnot(rownames(tmp.score)==rownames(tmp.meta))


# tmp.nontumor=data_rnaNorm[rownames(data_rnaNorm) %in% data_coupling_nontumor$sampleID,]
# tmp.nontumor=tmp.nontumor[order(rownames(tmp.nontumor)),]
# data_coupling_nontumor=data_coupling_nontumor[order(data_coupling_nontumor$sampleID),]
# table(rownames(tmp.nontumor)==data_coupling_nontumor$sampleID)
# data_coupling_nontumor$patient.multi=make.unique(as.character(data_coupling_nontumor$patient))
# rownames(tmp.nontumor)=data_coupling_nontumor$patient.multi

tmp.tumor=data_rnaNorm
tmp.tumor=tmp.tumor[order(rownames(tmp.tumor)),]
data_coupling_tumor=data_coupling_tumor[order(data_coupling_tumor$sampleID),]
table(rownames(tmp.tumor)==data_coupling_tumor$sampleID)
data_coupling_tumor$patient.multi=make.unique(as.character(data_coupling_tumor$patient))
rownames(tmp.tumor)=data_coupling_tumor$patient.multi

# tmp.patient=intersect(rownames(tmp.score),intersect(data_coupling_tumor$patient,data_coupling_nontumor$patient)) #129 个患者都有RNA样本




# data_coupling_tumor=data_coupling_tumor[data_coupling_tumor$patient %in% tmp.patient,] #现在只有180个样本
# data_coupling_nontumor=data_coupling_nontumor[data_coupling_nontumor$patient %in% tmp.patient,] #现在只有270个样本
# tmp.score=tmp.score[rownames(tmp.score) %in% tmp.patient,,drop=F]
# tmp.tumor=tmp.tumor[rownames(tmp.tumor) %in% data_coupling_tumor$patient.multi,]
# # tmp.nontumor=tmp.nontumor[rownames(tmp.nontumor) %in% data_coupling_nontumor$patient.multi,]
# 
# tmp.score=tmp.score[order(rownames(tmp.score)),,drop=F]
# data_coupling_tumor=data_coupling_tumor[order((data_coupling_tumor$patient)),,drop=F]
# # data_coupling_nontumor=data_coupling_nontumor[order((data_coupling_nontumor$patient)),,drop=F]
# tmp.tumor=tmp.tumor[order(rownames(tmp.tumor)),]
# # tmp.nontumor=tmp.nontumor[order(rownames(tmp.nontumor)),]
# 
# stopifnot(((rownames(tmp.score) == unique(data_coupling_tumor$patient))))
# # stopifnot(((rownames(tmp.score) == unique(data_coupling_nontumor$patient))))
# 
# # stopifnot(((rownames(tmp.nontumor) == data_coupling_nontumor$patient.multi)))
# stopifnot(((rownames(tmp.tumor) == data_coupling_tumor$patient.multi)))









var_meta=foreach(i=1:ncol(tmp.meta),.combine = rbind) %do%  {
  tmp.feature=colnames(tmp.meta)[i]
  tmp.var=var(tmp.meta[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_path=foreach(i=1:ncol(tmp.path),.combine = rbind) %do%  {
  tmp.feature=colnames(tmp.path)[i]
  tmp.var=var(tmp.path[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_tumor=foreach(i=1:ncol(tmp.tumor),.combine = rbind) %do%  {
  tmp.feature=colnames(tmp.tumor)[i]
  tmp.var=var(tmp.tumor[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
# var_nontumor=foreach(i=1:ncol(tmp.nontumor),.combine = rbind) %do%  {
#   tmp.feature=colnames(tmp.nontumor)[i]
#   tmp.var=var(tmp.nontumor[,tmp.feature])
  
#   return.string=data.frame(feature=tmp.feature,var=tmp.var)
# }





quantile(var_meta$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_path$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_tumor$var, probs = c(0,0.25,0.5,0.85,1)) 
# quantile(var_nontumor$var, probs = c(0,0.25,0.5,0.85,1)) 

##根据变异系数进行筛选纳入的通路和基因

tmp.meta=tmp.meta[,colnames(tmp.meta) %in% var_meta$feature[var_meta$var>quantile(var_meta$var,0.85)]]
tmp.path=tmp.path[,colnames(tmp.path) %in% var_path$feature[var_path$var>quantile(var_path$var,0.85)]]
tmp.tumor=tmp.tumor[,colnames(tmp.tumor) %in% var_tumor$feature[var_tumor$var>quantile(var_tumor$var,0.85)]]
# tmp.nontumor=tmp.nontumor[,colnames(tmp.nontumor) %in% var_nontumor$feature[var_nontumor$var>quantile(var_nontumor$var,0.85)]]



tmp.meta1=merge(data_coupling_tumor,tmp.meta,by.x="patient",by.y="row.names",all=F)
tmp.path1=merge(data_coupling_tumor,tmp.path,by.x="patient",by.y="row.names",all=F)
data_coupling_tumor=data_coupling_tumor[order(data_coupling_tumor$patient.multi),]
tmp.tumor=tmp.tumor[order(rownames(tmp.tumor)),]
stopifnot(((data_coupling_tumor$patient.multi == rownames(tmp.tumor))))
tmp.tumor=tmp.tumor[rownames(tmp.tumor) %in% tmp.meta1$patient.multi,]
tmp.path1=tmp.path1[tmp.path1$patient.multi %in% tmp.meta1$patient.multi,]
tmp.meta1=tmp.meta1[order(tmp.meta1$patient.multi),]
tmp.path1=tmp.path1[order(tmp.path1$patient.multi),]

stopifnot(tmp.meta1$patient.multi == rownames(tmp.tumor))
stopifnot(tmp.path1$patient.multi == rownames(tmp.tumor))


rownames(tmp.meta1)=tmp.meta1$patient.multi;tmp.meta1=tmp.meta1[,-1:-6];tmp.meta1$patient.multi=NULL
rownames(tmp.path1)=tmp.path1$patient.multi;tmp.path1=tmp.path1[,-1:-6];tmp.path1$patient.multi=NULL
MOFA_data1=list(as.matrix(t(tmp.path1)),as.matrix(t(tmp.meta1)),as.matrix(t(tmp.tumor)))

save(MOFA_data1,file="D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/MOFA_data1.RData")

#==========开始跑MOFA================
MOFAobject1 <- create_mofa(MOFA_data1)

data_opts <- get_default_data_options(MOFAobject1)
data_opts$scale_views=F
model_opts <- get_default_model_options(MOFAobject1)
model_opts$num_factors=30
train_opts <- get_default_training_options(MOFAobject1)
train_opts$convergence_mode="fast"

MOFAobject1 <- prepare_mofa(
  object = MOFAobject1,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/","model.hdf5")
library(reticulate)
use_python("D:/miniconda3/envs/mofa_env/python.exe", required = TRUE)
# use_condaenv("mofa_env", required = TRUE)
# py_config()
MOFAobject1.trained <- run_mofa(MOFAobject1, outfile, use_basilisk = F)


model1 <- load_model(outfile)
plot_data_overview(model1)
Nsamples = sum(model1@dimensions$N)
sample_metadata <- data.frame(
  sample = samples_names(model1)[[1]]
)

sample_metadata=merge(sample_metadata,data_coupling_tumor,by.x="sample",by.y="patient.multi",all=F)
sample_metadata=merge(sample_metadata,data_phenotype,by.x="patient",by.y="sampleID",all=F)
samples_names(model1)[[1]]==sample_metadata$sample
sample_metadata$group=NULL
samples_metadata(model1) <- sample_metadata



## 计算MOFA的Factor和肠型评分的关系----------
MOFA_variance=(model1@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(model1, x="view", y="factor",max_r2 =15)
# write.table(MOFA_variance,file = './新版肠型评分结果整理/big_mofa_variance.txt',sep = "\t")
test <- read.table("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/整图/big_mofa_variance.txt",header=T)

factors <- get_factors(model1, as.data.frame = T)
factors$factor=as.character(factors$factor)
factors$factor=gsub("Factor","",factors$factor)
factors$factor=as.numeric(factors$factor)

sample_metadata$Gender=as.factor(sample_metadata$Gender)
sample_metadata$Smoke=as.factor(sample_metadata$Smoke)
sample_metadata$Drink=as.factor(sample_metadata$Drink)




MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lmerTest::lmer(value ~ Age+Gender+BMI+Smoke+Drink+balance_value + (1|patient),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="balance_value"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="balance_value"])
}
#MOFA_correlation=MOFA_correlation[MOFA_correlation$Factor %in% c(1:10),] # in case that no FDR passed
MOFA_correlation$FDR=p.adjust(MOFA_correlation$Pvalue,method = "BH")
# write.table(MOFA_correlation,file = './新版肠型评分结果整理/big_mofa_correlation.txt',sep = "\t")


### version2 
MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Drink','AFP','TumorNumber','TumorSize')
  tmp.data=tmp.data[c('value',feature_name,'patient')]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    mm=lmerTest::lmer(as.formula(paste0("value~",j,"+(1|patient)")) ,data = tmp.data)
    mm=as.data.frame(summary(mm)$coef)
    return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                             Pvalue=mm$`Pr(>|t|)`[2])
  }
  
  
  return.string=cor1
}



i=1
tmp.data=factors[factors$factor==i,]
tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
tmp.data$group=NULL
feature_name=c('Age','Gender','BMI','Smoke','Drink','AFP','TumorNumber','TumorSize','balance_value')
tmp.data=tmp.data[c('value',feature_name,'patient')]


p1 <- ggplot(data = tmp.data,aes(y=value,x=balance_value))+geom_point(color="#DA5F48",fill="#DA5F48",alpha=0.8,size=3)+
  geom_smooth(method = 'lm',color="#C45641",fill="#DA7E6D")+
  ggpubr::stat_cor()+theme_bw()+
  theme(axis.title.x =element_text(size=60/.pt), 
        axis.title.y=element_text(size=60/.pt),
        axis.text.x = element_text(size=50/.pt,colour = 'black'),
        axis.text.y = element_text(size=50/.pt,colour = 'black'),
        legend.text = element_text(size=50/.pt,colour = 'black'),
        legend.title = element_text(size=50/.pt,colour = 'black'),
        plot.title=element_text(hjust=0.5,size=50/.pt))

pdf("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/plot/FigSC_bigmeta_score_corelation.pdf",width = 5,height = 5)
p1
dev.off()

## 计算前5Factor转录组的富集结果-正相关--------
library(reactome.db)
library(ReactomePA)
library(clusterProfiler)
library(MOFAdata)
genes = bitr(features_names(model1)[["view_3"]], fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db",drop=F)
genes=genes[!duplicated(genes$SYMBOL),]


genes$ENSEMBL[is.na(genes$ENSEMBL)]=genes$SYMBOL[is.na(genes$ENSEMBL)]
data(MSigDB_v6.0_C2_human)
features_names(model1)[["view_3"]] <- toupper(genes$ENSEMBL)
enrich_path=enrichment.parametric <- run_enrichment(model1,
                                                    view = "view_3", factors = 1:5,
                                                    feature.sets = MSigDB_v6.0_C2_human,
                                                    sign = "positive",
                                                    statistical.test = "parametric",set.statistic = "rank.sum"
)$pval.adj %>%
  as.data.frame() %>% 
  rownames_to_column('path') %>% pivot_longer(.,colnames(.)[-1],names_to = 'factor',
                                              values_to = 'padj')


enrich_path$type=ifelse(enrich_path$padj<0.05,"sig",'nosig')
enrich_path <- enrich_path %>% dplyr::filter(str_detect(path,'REACT')) ##只保留REACTOME信号通路

#enrich_path <- enrich_path %>% dplyr::filter(factor=='Factor1') %>% dplyr::arrange(padj)

##提取通路富集的基因，确认哪些基因需要画图
gene_weight=plot_enrichment_detailed(run_enrichment(model1,
                                                    view = "view_3", factors = 1:5,
                                                    feature.sets = MSigDB_v6.0_C2_human,
                                                    sign = "positive",
                                                    statistical.test = "parametric",set.statistic = "rank.sum"
), 
factor = 1
)$data

gene_weight <- gene_weight %>% dplyr::filter(str_detect(pathway_long_name,"REACT"))
gene_weight$symbol=genes$SYMBOL[match(gene_weight$feature,genes$ENSEMBL)]
gene_weight <- gene_weight %>% dplyr::arrange(pathway,feature)

gene_weight <- gene_weight %>% dplyr::arrange(pathway,desc(feature.statistic))

library(clusterProfiler)

gene_name=unique(gene_weight$symbol[gene_weight$feature.statistic>0])

#write.table(data.frame(gene=unique(gene_weight$symbol[gene_weight$feature.statistic>0])),file = './gene.txt',row.names = F,quote = F) 





factor_gene=factors %>% dplyr::filter(factor==1) %>% dplyr::select(-group)

factor_gene=left_join(factor_gene,tmp.tumor %>% rownames_to_column('sample'),by=c('sample'='sample'))

library(foreach)

gene_name=str_replace_all(gene_name,"-","_")
colnames(factor_gene)=str_replace_all(colnames(factor_gene),"-","_")

gene_correlation = foreach(i=gene_name,.combine = rbind) %do%  {
  tmp.data=factor_gene %>% dplyr::select(c(sample,value,i))
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  # mm=lmerTest::lmer(as.formula(paste0('value ~ Age+Gender+BMI+Smoke+Drink+',i,'+(1|patient)')),data = tmp.data)
  mm=lmerTest::lmer(as.formula(paste0(i,' ~ value+Age+Gender+BMI+Smoke+Drink+(1|patient)')),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=='value'],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="value"])
}
gene_correlation$padjust=p.adjust(gene_correlation$Pvalue,method = 'BH') 
library(writexl)
write_xlsx(gene_correlation,"./整图/bigmeta_output/bigmeta_gene_cor.xlsx")

write.table(gene_correlation,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/big_mofa_gene_correlation',row.names = F,quote = F,sep = "\t") 


## 用来整附表
gene_correlation = foreach(i=gene_name,.combine = rbind) %do%  {
  tmp.data=factor_gene %>% dplyr::select(c(sample,value,i))
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  # mm=lmerTest::lmer(as.formula(paste0('value ~',i, '+Age+Gender+BMI+Smoke+Drink+(1|patient)')),data = tmp.data)
  mm=lmerTest::lmer(as.formula(paste0(i,' ~ value+Age+Gender+BMI+Smoke+Drink+(1|patient)')),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  mm=mm[2,]
  mm$gene=i
  return.string=mm
}
gene_correlation$FDR <- p.adjust(gene_correlation$`Pr(>|t|)`,"BH")
gene_correlation$gene <- gsub("_","-",gene_correlation$gene)
write_xlsx(gene_correlation,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/bigmeta_gene_cor.xlsx")

#=========写个for循环批量画图并合并===============
target <- c("HLA_DOB","HLA_DQA1","HLA_DQA2","CD3D","CD3E","ICOS")
plist <- list()
for (i in 1:length(target)){
  tmp.data <- factor_gene[,c("value",target[i])]
  colnames(tmp.data)[2] <- "gene"
  plist[[i]] <- ggplot(data = tmp.data,aes(x=value,y=gene))+geom_point(color="#DA5F48",fill="#DA5F48",alpha=0.8,size=7)+
    geom_smooth(method = 'lm',color="#C45641",fill="#DA7E6D")+
    # geom_smooth(method = 'lm',color='#b93f33',fill='#b93f33',size=1.5)+ggpubr::stat_cor()+theme_classic()+
    theme(panel.background =element_blank(),legend.position ="none",
          axis.line=element_line(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text=element_text(colour="black"),
          axis.title = element_text(colour="black"))+
    labs(title = gsub("_","-",target[i]),y="Gene expression",x="value")
  
}
library(gridExtra)
library(patchwork)
pdf("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/plot/FigSD_bigmeta_gene_cor_ar.pdf",width=35,height=6)
plist[[1]]+plist[[2]]+plist[[3]]+plist[[4]]+plist[[5]]+plist[[6]]+plot_layout(ncol=6,nrow=1)
dev.off()

save(data_coupling_tumor,data_phenotype,meta_data,tmp.meta,tmp.meta1,
     tmp.path,tmp.path1,file = "D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/mofa_big_0730.RData")

#=========整理主图RData============
rm(list=ls())
load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/mofa_big_0730.RData")
# load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/mofa_big_0714.rdata")

outfile = file.path("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/","model.hdf5")

model1 <- load_model(outfile)
plot_data_overview(model1)
Nsamples = sum(model1@dimensions$N)
sample_metadata <- data.frame(
  sample = samples_names(model1)[[1]]
)
sample_metadata=merge(sample_metadata,data_coupling_tumor,by.x="sample",by.y="patient.multi",all=F)
sample_metadata=merge(sample_metadata,data_phenotype,by.x="patient",by.y="sampleID",all=F)
samples_names(model1)[[1]]==sample_metadata$sample
sample_metadata$group=NULL
samples_metadata(model1) <- sample_metadata

## 计算MOFA的Factor和肠型评分的关系----------
MOFA_variance=(model1@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(model1, x="view", y="factor",max_r2 =15)
# out <- cbind(row.names(MOFA_variance),MOFA_variance)
# out <- as.data.frame(out)
# library(writexl)
# write_xlsx(out,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/Bigmeta_factor解释度.xlsx")
# write.table(MOFA_variance,file = './新版肠型评分结果整理/big_mofa_variance.txt',sep = "\t")
# test <- read.table("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/整图/big_mofa_variance.txt",header=T)

factors <- get_factors(model1, as.data.frame = T)
factors$factor=as.character(factors$factor)
factors$factor=gsub("Factor","",factors$factor)
factors$factor=as.numeric(factors$factor)

sample_metadata$Gender=as.factor(sample_metadata$Gender)
sample_metadata$Smoke=as.factor(sample_metadata$Smoke)
sample_metadata$Drink=as.factor(sample_metadata$Drink)

MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lmerTest::lmer(value ~ Age+Gender+BMI+Smoke+Drink+balance_value + (1|patient),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="balance_value"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="balance_value"])
}
#MOFA_correlation=MOFA_correlation[MOFA_correlation$Factor %in% c(1:10),] # in case that no FDR passed
MOFA_correlation$FDR=p.adjust(MOFA_correlation$Pvalue,method = "BH")

## 用来整附表
MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lmerTest::lmer(value ~ balance_value + Age+Gender+BMI+Smoke+Drink+(1|patient),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  mm <- mm[2,]
  mm$Factor <- i
  return.string=mm
}
#MOFA_correlation=MOFA_correlation[MOFA_correlation$Factor %in% c(1:10),] # in case that no FDR passed
MOFA_correlation$FDR=p.adjust(MOFA_correlation$`Pr(>|t|)`,method = "BH")
write_xlsx(MOFA_correlation,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/Bigmeta_factor和评分相关性.xlsx")



#临床表型相关性热图----------
sample_metadata$BCLC[sample_metadata$BCLC %in% c('0','A')]='0/A'
sample_metadata$BCLC[sample_metadata$BCLC %in% c('B','C')]='B/C'

sample_metadata$TumorNumber[sample_metadata$TumorNumber==">3"]='4'

sample_metadata$TumorNumber=as.character(sample_metadata$TumorNumber)



sample_metadata$MVI=as.factor(sample_metadata$MVI)


#========第二版本热图===================
MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Drink','AFP','TumorNumber','TumorSize')
  basic=c('Age','Gender','BMI','Smoke','Drink')
  tmp.data=tmp.data[c('value',feature_name,'patient')]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    if (j %in% basic){
      basic1 <- basic[!basic==j]
      mm=lmerTest::lmer(as.formula(paste0("value~",j,"+(1|patient)+",basic1[1],"+",basic1[2],"+",basic1[3],"+",basic1[4])) ,data = tmp.data)
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
    }
    else {
      mm=lmerTest::lmer(as.formula(paste0("value~",j,"+(1|patient)+Age+Gender+BMI+Smoke+Drink")) ,data = tmp.data)
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
    }
  }
  
  
  return.string=cor1
}



i=1
tmp.data=factors[factors$factor==i,]
tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
tmp.data$group=NULL
feature_name=c('Age','Gender','BMI','Smoke','Drink','AFP','TumorNumber','TumorSize','balance_value')
tmp.data=tmp.data[c('value',feature_name,'patient')]


p1 <- ggplot(data = tmp.data,aes(y=value,x=balance_value))+geom_point(color='#d73027',alpha=0.7,size=5)+geom_smooth(method = 'lm',color='#d73027')+
  ggpubr::stat_cor()+theme_bw()+
  theme(axis.title.x =element_text(size=60/.pt), 
        axis.title.y=element_text(size=60/.pt),
        axis.text.x = element_text(size=50/.pt,colour = 'black',angle = 45,hjust = 1),
        axis.text.y = element_text(size=50/.pt,colour = 'black'),
        legend.text = element_text(size=50/.pt,colour = 'black'),
        legend.title = element_text(size=50/.pt,colour = 'black'),
        plot.title=element_text(hjust=0.5,size=50/.pt))


p1


Plot_heatmap=MOFA_correlation2 %>% dplyr::select(-Pvalue) %>% pivot_wider(.,names_from = feature,values_from = Beta) %>% column_to_rownames('Factor')
Plot_heatmap_pvalue <- MOFA_correlation2 %>% dplyr::select(-Beta) %>% pivot_wider(.,names_from = feature,values_from = Pvalue) %>% column_to_rownames('Factor')

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.8, 0,0.8), c("#4475b4","#eeeeee","#d73027")) # 设置连续颜色

Part1=Plot_heatmap 
Part1_p <- Plot_heatmap_pvalue


## version2 
MOFA_correlation=MOFA_correlation %>% dplyr::filter(Factor %in% seq(1,15))

Part2=MOFA_correlation %>% dplyr::select(-c(Pvalue,FDR))  %>%  column_to_rownames('Factor')
Part2_p=MOFA_correlation %>% dplyr::select(-c(Beta,FDR))  %>%  column_to_rownames('Factor')


















H1 <- Heatmap(t(Part1),cluster_rows = F,cluster_columns = F,col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if ((t(Part1)[i,j]>=0.5 | t(Part1)[i,j]<(-0.5))  & t(Part1_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part1)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='white'))
                }else if((t(Part1)[i,j]>=0.5 | t(Part1)[i,j]<(-0.5))  & t(Part1_p)[i,j]>0.05){
                  grid.text(sprintf("%.2f", t(Part1)[i, j]), x, y, gp = gpar(fontsize = 10,col='white'))
                  
                }else if( t(Part1_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part1)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='black'))
                }else if( t(Part1_p)[i,j]>0.05) {
                  grid.text(sprintf("%.2f", t(Part1)[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                  
                }
              }
)
H1



col_fun2 = colorRamp2(c(-0.1, 0,0.15), c("#9a382e","#f4f5bc","#385531")) # 设置连续颜色
H2 <- Heatmap(t(Part2),cluster_rows = F,cluster_columns = F,col = col_fun2,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (t(Part2)[i,j]>=0.1 & t(Part2_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part2)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='white'))
                  
                  
                }else if(t(Part2)[i,j]>=0.1 & t(Part2_p)[i,j]>=0.05){
                  grid.text(sprintf("%.2f", t(Part2)[i, j]), x, y, gp = gpar(fontsize = 10,col='white'))
                  
                }else if(t(Part2)[i,j]<(-0.1) & t(Part2_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part2)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='black'))
                }else{
                  grid.text(sprintf("%.2f", t(Part2)[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                  
                }},
)



H2




Plot_heatmap2=t(MOFA_variance) %>% as.data.frame()
Plot_heatmap2=Plot_heatmap2[1:15]
rownames(Plot_heatmap2)=c('Microbiome','Metabolome','Transcriptome')


col_fun2 = colorRamp2(c(0,3), c("white","#000a85"))

p2 <- Heatmap(Plot_heatmap2,cluster_rows = F,cluster_columns = F,col = col_fun2,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (Plot_heatmap2[i, j]>3) {
                  
                  grid.text(sprintf("%.2f", Plot_heatmap2[i, j]), x, y, gp = gpar(fontsize = 10,col='white'))  
                  
                }else{
                  grid.text(sprintf("%.2f", Plot_heatmap2[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                }
                
              },
              row_gap = unit(4,"mm"),
              rect_gp = gpar(col = "white", lwd = 2)
)

ht_list = H1 %v% p2%v% H2


#pdf('./新版肠型评分结果整理//bigdata_heatmap_v1.pdf',width = 8)
draw(ht_list)
#dev.off()



## 保存version2
saveRDS(ht_list,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/bigdata_heatmap_v2.rds')

## 计算前5Factor转录组的富集结果-负相关-----------
library(reactome.db)
library(ReactomePA)
library(clusterProfiler)
library(MOFAdata)
data(MSigDB_v6.0_C2_human) # MSigDB_v6.0_C2_human;reactomeGS; MSigDB_v6.0_C5_human
library(org.Hs.eg.db)

#====整附表==============
weights_df <- get_weights(model1,
                          views = "view_3",     # 指定视图
                          factors = "Factor1",  # 指定 factor
                          as.data.frame = TRUE) # 返回 data.frame 格式

weights_df1 <- weights_df[weights_df$feature%like%"ENSG",]
weights_df2 <- weights_df[!weights_df$feature%like%"ENSG",]

genes = bitr(weights_df1$feature, toType="SYMBOL", fromType="ENSEMBL", OrgDb="org.Hs.eg.db",drop=F)
weights_df1 <- merge(weights_df1,genes,by.x="feature",by.y="ENSEMBL",all.x=T)
weights_df1 <- weights_df1[!duplicated(weights_df1$feature,fromLast=T),]
# table(duplicated(weights_df1$feature))
# test <- weights_df1[duplicated(weights_df1$feature),]
weights_df1$feature <- as.character(weights_df1$feature)
for (i in 1:nrow(weights_df1)){
  weights_df1$feature[i] <- ifelse(is.na(weights_df1$SYMBOL[i])==T,weights_df1$feature[i],weights_df1$SYMBOL[i]) 
}
weights_df1$SYMBOL <- NULL
weights_df <- rbind(weights_df1,weights_df2)
write_xlsx(weights_df,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/bigmeta_gene_weight.xlsx")

#=======通路富集===================
genes = bitr(features_names(model1)[["view_3"]], fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db",drop=F)
genes=genes[!duplicated(genes$SYMBOL),]


genes$ENSEMBL[is.na(genes$ENSEMBL)]=genes$SYMBOL[is.na(genes$ENSEMBL)]
data(MSigDB_v6.0_C2_human)
features_names(model1)[["view_3"]] <- toupper(genes$ENSEMBL)

enrich_path=enrichment.parametric <- run_enrichment(model1,
                                                    view = "view_3", factors = 1:5,
                                                    feature.sets = MSigDB_v6.0_C2_human,
                                                    sign = "positive",
                                                    statistical.test = "parametric",set.statistic = "rank.sum"
)$pval.adj %>%
  as.data.frame() %>% 
  rownames_to_column('path') %>% pivot_longer(.,colnames(.)[-1],names_to = 'factor',
                                              values_to = 'padj')

enrich_path$type=ifelse(enrich_path$padj<0.05,"sig",'nosig')
enrich_path <- enrich_path %>% dplyr::filter(str_detect(path,'REACT')) ##只保留REACTOME信号通路

enrich_path2=run_enrichment(model1,
                            view = "view_3", factors = 1:5,
                            feature.sets = MSigDB_v6.0_C2_human,
                            sign = "negative",
                            statistical.test = "parametric",set.statistic = "rank.sum"
)$pval.adj %>%
  as.data.frame() %>% 
  rownames_to_column('path') %>% pivot_longer(.,colnames(.)[-1],names_to = 'factor',
                                              values_to = 'padj')
enrich_path2$type=ifelse(enrich_path2$padj<0.05,"sig",'nosig')
enrich_path2 <- enrich_path2 %>% dplyr::filter(str_detect(path,'REACT'))




## 合并起来绘制火山图--------

enrich_path_all=rbind(enrich_path %>% dplyr::mutate(group='UP'),enrich_path2 %>% dplyr::mutate(group='Down'))
write_xlsx(enrich_path_all,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/bigmeta_GSVA_pathway.xlsx")

enrich_path_all$factor=factor(enrich_path_all$factor,levels = paste0("Factor",seq(1,5)))
enrich_path_all$size=-log(enrich_path_all$padj)
enrich_path_all$padj=ifelse(enrich_path_all$group=='UP',-log10(enrich_path_all$padj),log10(enrich_path_all$padj))

new_pallete=c('#6e7ca5','#7c9680','#fe9014','#fcd364',
              '#ffcad3','#2482BC','#7ca6be','#92B096','#8B516E','#D65741','#8298A3','#52AEB5')


enrich_path_all_label=enrich_path_all %>% group_by(factor) %>% dplyr::arrange(desc(padj),.by_group = T) #

enrich_path_all_label_up=enrich_path_all_label%>%
  group_by(factor) %>% dplyr::filter(padj>0) %>% 
  top_n(.,3,padj)
enrich_path_all_label_down=enrich_path_all_label%>%
  group_by(factor)  %>% dplyr::filter(padj<0) %>% 
  top_n(.,-3,padj)

enrich_path_all_label_all=rbind(enrich_path_all_label_down,enrich_path_all_label_up)
enrich_path_all_label_all$match=paste0(enrich_path_all_label_all$path,enrich_path_all_label_all$factor,
                                       enrich_path_all_label_all$padj)

enrich_path_all$match=paste0(enrich_path_all$path,enrich_path_all$factor,enrich_path_all$padj)

enrich_path_all$label=ifelse(enrich_path_all$match %in% enrich_path_all_label_all$match,enrich_path_all$path,"")

enrich_path_all$label=str_remove_all(enrich_path_all$label,"KEGG_") %>% 
  str_remove_all(.,"REACTOME_") %>% str_remove_all(.,"NABA_")

enrich_path_all$label=str_to_title(enrich_path_all$label)
enrich_path_all_sig=enrich_path_all %>% dplyr::filter(type=='sig')
enrich_path_all_nosig=enrich_path_all %>% dplyr::filter(type=='nosig')

p1_ren <- ggplot()+
  geom_jitter(aes(x=factor ,y=padj,size=size),alpha=0.7,color='grey',data = enrich_path_all_nosig,width = 0.3)+
  
  geom_jitter(aes(x=factor,y=padj,size=size,color=factor),alpha=0.7,data = enrich_path_all_sig,width = 0.3)+
  geom_tile(aes(x=factor,y=0,fill=factor),height=1,color='black',data = enrich_path_all%>% distinct(factor,.keep_all = T))+
  
  scale_fill_manual(values = c('Factor1'='#6e7ca5','Factor2'='#7c9680','Factor3'='#fe9014','Factor4'='#8B516E','Factor5'='#2482BC'))+
  scale_color_manual(values = c('Factor1'='#6e7ca5','Factor2'='#7c9680','Factor3'='#fe9014','Factor4'='#8B516E','Factor5'='#2482BC'))+
  geom_text(data=enrich_path_all%>% distinct(factor,.keep_all = T), # 绘制中心分组标记图文本注释
            aes(x=factor, 
                y=0, 
                label=factor),
            size = 5,
            color ="white")+
  ggrepel::geom_text_repel(aes(x=factor,y=padj,label = label),data = enrich_path_all_sig,size=2,max.overlaps = 50)+
  # ylim(c(-50,50))+
  theme_minimal() + 
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())
p1_ren
# ggsave('./新版肠型评分结果整理/rna_big.pdf',width = 10,height = 5)

##改一下颜色
p1_ren <- ggplot()+
  geom_jitter(aes(x=factor ,y=padj,size=size),alpha=0.7,color='grey',data = enrich_path_all_nosig,width = 0.3)+
  
  geom_jitter(aes(x=factor,y=padj,size=size,color=factor),alpha=0.7,data = enrich_path_all_sig,width = 0.3)+
  geom_tile(aes(x=factor,y=0,fill=factor),height=2,color='black',data = enrich_path_all%>% distinct(factor,.keep_all = T))+
  
  scale_fill_manual(values = c('Factor1'='#E31A1C','Factor2'='#1F78B4','Factor3'='#FF7F00','Factor4'='#6A3D9A','Factor5'='#33A02C'))+
  scale_color_manual(values = c('Factor1'='#E31A1C','Factor2'='#1F78B4','Factor3'='#FF7F00','Factor4'='#6A3D9A','Factor5'='#33A02C'))+
  geom_text(data=enrich_path_all%>% distinct(factor,.keep_all = T), # 绘制中心分组标记图文本注释
            aes(x=factor, 
                y=0, 
                label=factor),
            size = 5,
            color ="white")+
  ggrepel::geom_text_repel(aes(x=factor,y=padj,label = label),data = enrich_path_all_sig,size=2,max.overlaps = 50)+ylim(c(-50,50))+
  theme_minimal() + ylim(c(-25,25))+
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())
p1_ren


## 修改jitter
p1_ren <- ggplot()+
  geom_jitter(aes(x=factor ,y=padj,size=size),alpha=0.7,color='grey',data = enrich_path_all_nosig,width = 0.3)+
  
  geom_jitter(aes(x=factor,y=padj,size=size,color=factor),alpha=0.7,data = enrich_path_all_sig,width = 0.3)+
  geom_tile(aes(x=factor,y=0,fill=factor),height=1,color='black',data = enrich_path_all%>% distinct(factor,.keep_all = T))+
  
  scale_fill_manual(values = c('Factor1'='#E31A1C','Factor2'='#1F78B4','Factor3'='#FF7F00','Factor4'='#6A3D9A','Factor5'='#33A02C'))+
  scale_color_manual(values = c('Factor1'='#E31A1C','Factor2'='#1F78B4','Factor3'='#FF7F00','Factor4'='#6A3D9A','Factor5'='#33A02C'))+
  geom_text(data=enrich_path_all%>% distinct(factor,.keep_all = T), # 绘制中心分组标记图文本注释
            aes(x=factor, 
                y=0, 
                label=factor),
            size = 5,
            color ="white")+
  ggrepel::geom_text_repel(aes(x=factor,y=padj,label = label),data = enrich_path_all_sig,size=2,max.overlaps = 50)+
  #ylim(c(-50,50))+
  theme_minimal() + 
  #ylim(c(-25,25))+
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())
p1_ren

saveRDS(p1_ren,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/Big_rna.rds')

## 重新画火山图
saveRDS(p1_ren,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/Big_rna_260206.rds')
pdf("Figure拼图/FigC_left_260208.pdf",width = 7.5,height = 5)
p1_ren
dev.off()

#==========画代谢图=============

library(ggridges)

meta_data2=meta_data

meta_merge <- merge(tmp.meta1,factors %>% dplyr::filter(factor==1),by.x="row.names",by.y="sample",all=F)


meta_merge <- meta_merge %>% dplyr::select(-c(group,factor,Row.names))

meta_merge <- meta_merge %>%  dplyr::select(value,everything()) %>% 
  dplyr::mutate(group=ifelse(value>mean(value),"high","low"))
meta_merge$value=NULL
meta_merge <- meta_merge %>% dplyr::select(group,everything()) %>% 
  pivot_longer(.,colnames(.)[-1],names_to = 'meta',values_to = 'value')

# meta_merge$meta=meta_annot$MS2_name[match(meta_merge$meta,meta_annot$id)]

meta_merge_try=meta_merge %>% dplyr::filter(meta %in% c('MethylIndole-3-acetate','indolin-2-one'))
library(ggridges)
meta_merge_try$group=factor(meta_merge_try$group,levels = c('low','high'))
ptry_big <- ggplot(data=meta_merge_try,aes(x=value,y=meta,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+ggtitle('Tryptophan metabolism')+
  scale_fill_manual(values = c("#97d0ab","#4579b2"))+scale_color_manual(values = c("#97d0ab","#4579b2"))+
  theme_classic()
ptry_big

meta_merge_acid=meta_merge %>% dplyr::filter(meta %in% c('Taurocholic acid','Taurodeoxycholic acid','Glycohyocholic acid'))
library(ggridges)
meta_merge_acid$group=factor(meta_merge_acid$group,levels = c('low','high'))
meta_merge_acid$meta <- factor(meta_merge_acid$meta,levels = rev(c('Taurocholic acid','Taurodeoxycholic acid','Glycohyocholic acid')))
pacid_big <- ggplot(data=meta_merge_acid,aes(x=value,y=meta,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+ggtitle('Bile acid metabolism')+
  scale_fill_manual(values = c("#8c4944","#756dac"))+scale_color_manual(values = c("#8c4944","#756dac"))+
  theme_classic()
pacid_big
# ggsave(filename = './新版肠型评分结果整理/ptry_ren.pdf',width = 10)
# ggsave(filename = './新版肠型评分结果整理/pacid_ren.pdf',width = 10)

###改一下配色
meta_merge_try=meta_merge %>% dplyr::filter(meta %in% c('MethylIndole-3-acetate','indolin-2-one'))
library(ggridges)
meta_merge_try$group=factor(meta_merge_try$group,levels = c('low','high'))
ptry_big <- ggplot(data=meta_merge_try,aes(x=value,y=meta,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+ggtitle('Tryptophan metabolism')+
  scale_fill_manual(values = c( "#A1D99B","#238B45"))+scale_color_manual(values = c("#A1D99B","#238B45"))+
  theme_classic()
ptry_big

meta_merge_acid=meta_merge %>% dplyr::filter(meta %in% c('Taurocholic acid','Taurodeoxycholic acid','Glycohyocholic acid'))
library(ggridges)
meta_merge_acid$group=factor(meta_merge_acid$group,levels = c('low','high'))
pacid_big <- ggplot(data=meta_merge_acid,aes(x=value,y=meta,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+ggtitle('Bile acid metabolism')+
  scale_fill_manual(values = c("#A7A3C4","#504394"))+scale_color_manual(values = c("#A7A3C4","#504394"))+
  theme_classic()
pacid_big



save(ptry_big,pacid_big,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/Pmetabolism_big.RData')

## metapath

MOFA_weights_all <- get_weights(model1, as.data.frame = T)

meta_path=MOFA_weights_all[MOFA_weights_all$view=="view_1",] %>% dplyr::filter(factor %in% c(paste0('Factor',seq(1,15))))


meta_path_label_up=meta_path %>%
  group_by(factor) %>% dplyr::filter(value>0) %>% 
  top_n(.,3,value)
meta_path_label_down=meta_path %>%
  group_by(factor) %>% dplyr::filter(value<0) %>% 
  top_n(.,-3,value)

meta_path_label=rbind(meta_path_label_up,meta_path_label_down)
meta_path_label$match=paste0(meta_path_label$feature,meta_path_label$factor,meta_path_label$value)

meta_path$match=paste0(meta_path$feature,meta_path$factor,meta_path$value)

meta_path$group=ifelse(meta_path$match %in% meta_path_label$match,"yes","no")
meta_path_label$feature=str_remove_all(meta_path_label$feature,":.*")


library(ggrepel)
library(ggbeeswarm)
# meta_path <- meta_path %>% dplyr::filter(factor %in% paste0('Factor',seq(1,5)))
# meta_path_label=meta_path_label%>% dplyr::filter(factor %in% paste0('Factor',seq(1,5)))
# Pmeta <- ggplot()+geom_beeswarm(data = meta_path,aes(x=factor,y=value,color=value),cex = 0.2,priority = 'descending',alpha=0.7,size=2)+
#   scale_color_viridis_c(option = 'B')+
#   geom_text_repel(aes(x=factor,y=value,label = feature),data = meta_path_label,show.legend = F,size=3,
#                   max.overlaps = Inf)+
#   geom_hline(yintercept = 0,lty=4,col="#666666",lwd=0.3)+
#   coord_flip()+ggthemes::theme_base()
# Pmeta
# ggsave('./metapath.pdf')
# save(Pmeta,file = './Pmeta.RData')

meta_path_f <- tmp.path1 %>% dplyr::select(c('PROPFERM-PWY: superpathway of L-alanine fermentation (Stickland reaction)',
                                             'PWY-8188: L-alanine degradation VI (reductive Stickland reaction)',
                                             'PWY-8189: L-alanine degradation V (oxidative Stickland reaction)')) %>% 
  rownames_to_column('sample') %>% pivot_longer(.,colnames(.)[-1],names_to = 'path',values_to = 'value')
meta_path_f <-left_join(meta_path_f,factors %>% dplyr::filter(factor=='1'),by=c('sample'='sample'))
meta_path_f$group=ifelse(meta_path_f$value.y>mean(meta_path_f$value.y),"high",'low')


Pmeta_big_p1 <- ggplot(data=meta_path_f,aes(x=value.x,y=path,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+
  scale_color_manual(values = c( "#F79927","#3F5763"))+
  scale_fill_manual(values = c( "#F79927","#3F5763"))+theme_classic()




meta_path_f2 <- tmp.path1 %>% dplyr::select(c('LIPA-CORESYN-PWY: lipid A-core biosynthesis (E. coli K-12)',
                                              'LPSSYN-PWY: superpathway of lipopolysaccharide biosynthesis')) %>% 
  rownames_to_column('sample') %>% pivot_longer(.,colnames(.)[-1],names_to = 'path',values_to = 'value')
meta_path_f2 <-left_join(meta_path_f2,factors %>% dplyr::filter(factor=='1'),by=c('sample'='sample'))
meta_path_f2$group=ifelse(meta_path_f2$value.y>mean(meta_path_f2$value.y),"high",'low')


Pmeta_big_p2 <- ggplot(data=meta_path_f2,aes(x=value.x,y=path,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+
  scale_color_manual(values = c( "#e58579","#8AB1D2"))+
  scale_fill_manual(values = c( "#e58579","#8AB1D2"))+theme_classic()

Pmeta_big_p1
Pmeta_big_p2
# ggsave(filename = './新版肠型评分结果整理/Pmeta_ren_p1.pdf',width = 10)
# ggsave(filename = './新版肠型评分结果整理/Pmeta_ren_p2.pdf',width = 10)



####改一下颜色
meta_path_f$group <- factor(meta_path_f$group,levels = c("low","high"))
Pmeta_big_p1 <- ggplot(data=meta_path_f,aes(x=value.x,y=path,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+
  scale_color_manual(values = c("#A6CEE3", "#1F78B4"))+
  scale_fill_manual(values = c(  "#A6CEE3", "#1F78B4"))+theme_classic()

Pmeta_big_p1

meta_path_f2$group <- factor(meta_path_f2$group,levels = c("low","high"))
Pmeta_big_p2 <- ggplot(data=meta_path_f2,aes(x=value.x,y=path,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+
  scale_color_manual(values = c( "#f9cd9c","#e18182"))+
  scale_fill_manual(values = c( "#f9cd9c","#e18182"))+theme_classic()
Pmeta_big_p2



# save(Pmeta_big_p1,Pmeta_big_p2,file = './新版肠型评分结果整理/Pmeta_density_big.RData')

save(Pmeta_big_p1,Pmeta_big_p2,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/Pmeta_density_big.RData')
library(patchwork)
ptry_big+pacid_big+Pmeta_big_p1+Pmeta_big_p2+plot_layout(ncol=2,nrow=2,widths = c(1,1,1,1),heights = c(1,1,1,1))


##=======整理代谢物附表==========
weights_df <- get_weights(model1,
                          views = "view_2",     # 指定视图
                          factors = "Factor1",  # 指定 factor
                          as.data.frame = TRUE) # 返回 data.fr
# weights_df <- merge(weights_df,meta_annot[,c("id","MS2_name")],by.x="feature",by.y="id",all.x=T)
# weights_df <- weights_df[,c("MS2_name","factor","value" ,"view","feature")]
write_xlsx(weights_df,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/Bigmeta_meta_weight.xlsx")

##=======整理菌代谢通路附表==========
weights_df <- get_weights(model1,
                          views = "view_1",     # 指定视图
                          factors = "Factor1",  # 指定 factor
                          as.data.frame = TRUE) # 返回 data.fr
write_xlsx(weights_df,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/Bigmeta_path_weight.xlsx")


#========把RData导出pdf，拼图===========
setwd("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure")
Big_rna <- readRDS("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/Big_rna.rds")
pdf("FigC_left.pdf",width = 7.5,height = 5)
Big_rna
dev.off()

bigdata_heatmap_v2 <- readRDS("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/bigdata_heatmap_v2.rds")
pdf("FigB_left.pdf",width = 10,height = 8)
bigdata_heatmap_v2
dev.off()

load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/Pmeta_density_big.RData")
load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/Figure/Pmetabolism_big.RData")

library(patchwork)
pdf("FigG_left.pdf",width=18,height = 10)
ptry_big+pacid_big+Pmeta_big_p1+Pmeta_big_p2+plot_layout(ncol=2,nrow=2,widths = c(1,1,1,1),heights = c(1,1,1,1))
dev.off()


#===========新版热图2601=======================
rm(list=ls())
load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/mofa_big_0730.RData")
# load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/mofa_big_0714.rdata")

outfile = file.path("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/bigmeta_output_250730/","model.hdf5")

model1 <- load_model(outfile)
plot_data_overview(model1)
Nsamples = sum(model1@dimensions$N)
sample_metadata <- data.frame(
  sample = samples_names(model1)[[1]]
)
sample_metadata=merge(sample_metadata,data_coupling_tumor,by.x="sample",by.y="patient.multi",all=F)
sample_metadata=merge(sample_metadata,data_phenotype,by.x="patient",by.y="sampleID",all=F)
samples_names(model1)[[1]]==sample_metadata$sample
sample_metadata$group=NULL
samples_metadata(model1) <- sample_metadata

## 计算MOFA的Factor和肠型评分的关系----------
MOFA_variance=(model1@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(model1, x="view", y="factor",max_r2 =15)
# out <- cbind(row.names(MOFA_variance),MOFA_variance)
# out <- as.data.frame(out)
# library(writexl)
# write_xlsx(out,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/Bigmeta_factor解释度.xlsx")
# write.table(MOFA_variance,file = './新版肠型评分结果整理/big_mofa_variance.txt',sep = "\t")
# test <- read.table("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/整图/big_mofa_variance.txt",header=T)

factors <- get_factors(model1, as.data.frame = T)
factors$factor=as.character(factors$factor)
factors$factor=gsub("Factor","",factors$factor)
factors$factor=as.numeric(factors$factor)

sample_metadata$Gender=as.factor(sample_metadata$Gender)
sample_metadata$Smoke=as.factor(sample_metadata$Smoke)
sample_metadata$Drink=as.factor(sample_metadata$Drink)

MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lmerTest::lmer(value ~ Age+Gender+BMI+Smoke+Drink+balance_value + (1|patient),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="balance_value"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="balance_value"])
}
#MOFA_correlation=MOFA_correlation[MOFA_correlation$Factor %in% c(1:10),] # in case that no FDR passed
MOFA_correlation$FDR=p.adjust(MOFA_correlation$Pvalue,method = "BH")





#临床表型相关性热图----------
sample_metadata$BCLC[sample_metadata$BCLC %in% c('0','A')]='0/A'
sample_metadata$BCLC[sample_metadata$BCLC %in% c('B','C')]='B/C'

sample_metadata$TumorNumber[sample_metadata$TumorNumber==">3"]='4'

sample_metadata$TumorNumber=as.numeric(sample_metadata$TumorNumber)



sample_metadata$MVI=as.factor(sample_metadata$MVI)

#========第二版本热图===================
MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  print(i)
  # i=6
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  # tmp.data$AFP <- log2(tmp.data$AFP)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Drink','AFP','TumorNumber','TumorSize')
  basic=c('Age','Gender','BMI','Smoke','Drink')
  tmp.data=tmp.data[c('value',feature_name,'patient')]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    # j=feature_name[8]
    if (j %in% basic){
      print(j)
      basic1 <- basic[!basic==j]
      mm=lmerTest::lmer(as.formula(paste0("value~",j,"+(1|patient)+",basic1[1],"+",basic1[2],"+",basic1[3],"+",basic1[4])) ,data = tmp.data)
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
      # print(return.string)
    }
    else {
      print(j)
      # mm=lmerTest::lmer(as.formula(paste0(j,"~value+(1|patient)+Age+Gender+BMI+Smoke+Drink")) ,data = tmp.data) ##这里可能会报错
      mm=lm(as.formula(paste0(j,"~value+Age+Gender+BMI+Smoke+Drink")) ,data = tmp.data) 
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
      print(return.string)
    }
    # else {
    #   print(j)
    # 
    #   # 尝试混合模型，如果出错则使用线性模型
    #   result <- tryCatch({
    #     # 尝试混合线性模型
    #     tmp.data1 <- tmp.data
    #     tmp.data1$value=scale(tmp.data1$value)
    #     mm <- lmerTest::lmer(as.formula(paste0(j,"~value+(1|patient)+Age+Gender+BMI+Smoke+Drink")),
    #                          data = tmp.data1)
    #     mm <- as.data.frame(summary(mm)$coef)
    #     mm
    #   }, error = function(e) {
    #     # 如果混合模型出错，使用线性模型
    #     print(paste0("lm:Factor ",i,"-",j))
    #     mm <- lm(as.formula(paste0(j, "~value+Age+Gender+BMI+Smoke+Drink")),
    #              data = tmp.data)
    #     mm <- as.data.frame(summary(mm)$coef)
    #     mm
    #   })
    # 
    #   return.string=data.frame(Factor=i,feature=j,Beta=result$Estimate[2],
    #                            Pvalue=result$`Pr(>|t|)`[2])
    #   # print(return.string)
    # }
  }
  
  # print(cor1)
  return.string=cor1
}

# lmerTest::lmer(as.formula(paste0(j,"~value+(1|patient)")), data = tmp.data1)


## AFP的beta scale一下
MOFA_correlation2_part <- MOFA_correlation2[MOFA_correlation2$feature=="AFP",]
MOFA_correlation2_part$Beta <- scale(MOFA_correlation2_part$Beta , center = FALSE, scale = TRUE)
MOFA_correlation2 <- MOFA_correlation2[!MOFA_correlation2$feature=="AFP",]
MOFA_correlation2 <- rbind(MOFA_correlation2,MOFA_correlation2_part)

i=1
tmp.data=factors[factors$factor==i,]
tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
tmp.data$group=NULL
feature_name=c('Age','Gender','BMI','Smoke','Drink','AFP','TumorNumber','TumorSize','balance_value')
tmp.data=tmp.data[c('value',feature_name,'patient')]


p1 <- ggplot(data = tmp.data,aes(y=value,x=balance_value))+geom_point(color='#d73027',alpha=0.7,size=5)+geom_smooth(method = 'lm',color='#d73027')+
  ggpubr::stat_cor()+theme_bw()+
  theme(axis.title.x =element_text(size=60/.pt), 
        axis.title.y=element_text(size=60/.pt),
        axis.text.x = element_text(size=50/.pt,colour = 'black',angle = 45,hjust = 1),
        axis.text.y = element_text(size=50/.pt,colour = 'black'),
        legend.text = element_text(size=50/.pt,colour = 'black'),
        legend.title = element_text(size=50/.pt,colour = 'black'),
        plot.title=element_text(hjust=0.5,size=50/.pt))


p1


Plot_heatmap=MOFA_correlation2 %>% dplyr::select(-Pvalue) %>% pivot_wider(.,names_from = feature,values_from = Beta) %>% column_to_rownames('Factor')
Plot_heatmap_pvalue <- MOFA_correlation2 %>% dplyr::select(-Beta) %>% pivot_wider(.,names_from = feature,values_from = Pvalue) %>% column_to_rownames('Factor')

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-1, 0,1), c("#4475b4","#eeeeee","#d73027")) # 设置连续颜色

Part1=Plot_heatmap 
Part1_p <- Plot_heatmap_pvalue


## version2 
MOFA_correlation=MOFA_correlation %>% dplyr::filter(Factor %in% seq(1,15))

Part2=MOFA_correlation %>% dplyr::select(-c(Pvalue,FDR))  %>%  column_to_rownames('Factor')
Part2_p=MOFA_correlation %>% dplyr::select(-c(Beta,FDR))  %>%  column_to_rownames('Factor')


















H1 <- Heatmap(t(Part1),cluster_rows = F,cluster_columns = F,col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if ((t(Part1)[i,j]>=0.5 | t(Part1)[i,j]<(-0.5))  & t(Part1_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part1)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='white'))
                }else if((t(Part1)[i,j]>=0.5 | t(Part1)[i,j]<(-0.5))  & t(Part1_p)[i,j]>0.05){
                  grid.text(sprintf("%.2f", t(Part1)[i, j]), x, y, gp = gpar(fontsize = 10,col='white'))
                  
                }else if( t(Part1_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part1)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='black'))
                }else if( t(Part1_p)[i,j]>0.05) {
                  grid.text(sprintf("%.2f", t(Part1)[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                  
                }
              }
)
H1



col_fun2 = colorRamp2(c(-0.1, 0,0.15), c("#9a382e","#f4f5bc","#385531")) # 设置连续颜色
H2 <- Heatmap(t(Part2),cluster_rows = F,cluster_columns = F,col = col_fun2,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (t(Part2)[i,j]>=0.1 & t(Part2_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part2)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='white'))
                  
                  
                }else if(t(Part2)[i,j]>=0.1 & t(Part2_p)[i,j]>=0.05){
                  grid.text(sprintf("%.2f", t(Part2)[i, j]), x, y, gp = gpar(fontsize = 10,col='white'))
                  
                }else if(t(Part2)[i,j]<(-0.1) & t(Part2_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part2)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='black'))
                }else{
                  grid.text(sprintf("%.2f", t(Part2)[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                  
                }},
)



H2




Plot_heatmap2=t(MOFA_variance) %>% as.data.frame()
Plot_heatmap2=Plot_heatmap2[1:15]
rownames(Plot_heatmap2)=c('Microbiome','Metabolome','Transcriptome')


col_fun2 = colorRamp2(c(0,3), c("white","#000a85"))

p2 <- Heatmap(Plot_heatmap2,cluster_rows = F,cluster_columns = F,col = col_fun2,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (Plot_heatmap2[i, j]>3) {
                  
                  grid.text(sprintf("%.2f", Plot_heatmap2[i, j]), x, y, gp = gpar(fontsize = 10,col='white'))  
                  
                }else{
                  grid.text(sprintf("%.2f", Plot_heatmap2[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                }
                
              },
              row_gap = unit(4,"mm"),
              rect_gp = gpar(col = "white", lwd = 2)
)

ht_list = H1 %v% p2%v% H2


pdf('D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/bigdata_heatmap_v1.pdf',width = 8)
draw(ht_list)
dev.off()



## 保存version2
saveRDS(ht_list,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/bigdata_heatmap_v2.rds')


pdf('D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/bigdata_heatmap_AFP.pdf',width = 8)
draw(ht_list)
dev.off()



## 保存version2
saveRDS(ht_list,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/bigdata_heatmap_AFP.rds')
