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

#=======load function=====================
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

#=========整理最后需要跑的样本名单===================
setwd("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/原始文件MOFAren")
library(readxl)
data_phenotype=readxl::read_xlsx("./RenxxInfo_250310.xlsx")


data_div <- read.table('D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/淑芬结果/renxx_BS_250916.txt',header=T)

data_phenotype=merge(data_phenotype,data_div[,c("sample","balance_score")],by.x="sample",by.y="sample",all=F)

##not_keep 1769
# data_phenotype=data_phenotype[data_phenotype$sample!="1769",]

data_score=data_phenotype


# final patients included
patients=data_phenotype$sample

# import rnaseq data
data_rna=read.csv("./Merge_118patient_gene_rpkm.csv",header = T,row.names = 1)



data_rna=data_rna[!duplicated(data_rna$Description),]
rownames(data_rna)=data_rna$Description
data_rna$Description=NULL
raw_data=data_rna
data_rnaNorm=raw_data
data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %like% "T"]
colnames(data_rnaNorm)=gsub("HCC","",colnames(data_rnaNorm))
colnames(data_rnaNorm)=gsub(".P.T","",colnames(data_rnaNorm))
colnames(data_rnaNorm)=gsub("T.*","",colnames(data_rnaNorm))
data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %in% patients]
data_rnaNorm=as.data.frame(t(data_rnaNorm))
# data_rnaNorm=data_rnaNorm[,colSums(data_rnaNorm>0)>(nrow(data_rnaNorm)*0.8)]
# data_rnaNorm=log(data_rnaNorm+0.1)
# data_rnaNorm=data_rnaNorm[,!colnames(data_rnaNorm) %like% "MT-"]
data_rnaNorm=as.data.frame(t(data_rnaNorm))

# import metacyc
data_path=read.table("./Metacyc_Renxx.txt",header = T,stringsAsFactors = F,row.names = 1,sep="\t",fill = T,comment.char = "", check.names = FALSE)
colnames(data_path)=gsub("_Abundance","",colnames(data_path))
colnames(data_path)=gsub("\\.","-",colnames(data_path))
data_path=data_path[!rownames(data_path) %like% "\\|",]
data_path=data_path[rowSums(data_path>0)>(ncol(data_path)*0.1),]
data_path=data_path[!rownames(data_path) %in% c("UNMAPPED","UNINTEGRATED"),]
data_path=as.data.frame(apply(data_path,2,function (x){
  x=x/sum(x)
  return(x)
}))
data_path_clr=transform_and_filter_taxa(data_path,samples_row = F,method = "clr",missing_filter = 0)
data_path_clr=as.data.frame(t(data_path_clr))
data_path_clr=data_path_clr[,colnames(data_path_clr) %in% patients]

load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/淑芬结果/sc_renxx250618.RData")
data_meta <- sc_rxx
# import meta
# data_meta=data_metabolism_rxx
data_meta=as.data.frame(t(data_meta))
data_meta=data_meta[,colnames(data_meta) %in% patients]

# overlap samples
samples=intersect(colnames(data_rnaNorm),intersect(colnames(data_meta),colnames(data_path_clr)))
data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %in% samples]
data_meta=data_meta[,colnames(data_meta) %in% samples]
data_path_clr=data_path_clr[,colnames(data_path_clr) %in% samples]

data_rnaNorm=data_rnaNorm[,order(colnames(data_rnaNorm))]
data_meta=data_meta[,order(colnames(data_meta))]
data_path_clr=data_path_clr[,order(colnames(data_path_clr))]

stopifnot(colnames(data_path_clr)==colnames(data_meta))
stopifnot(colnames(data_path_clr)==colnames(data_rnaNorm))

# top 10%-5% variance selection
data_path_clr=as.data.frame(t(data_path_clr))
data_meta=as.data.frame(t(data_meta))
data_rnaNorm=as.data.frame(t(data_rnaNorm))

var_meta=foreach(i=1:ncol(data_meta),.combine = rbind) %do%  {
  tmp.feature=colnames(data_meta)[i]
  tmp.var=var(data_meta[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_path=foreach(i=1:ncol(data_path_clr),.combine = rbind) %do%  {
  tmp.feature=colnames(data_path_clr)[i]
  tmp.var=var(data_path_clr[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_tumor=foreach(i=1:ncol(data_rnaNorm),.combine = rbind) %do%  {
  tmp.feature=colnames(data_rnaNorm)[i]
  tmp.var=var(data_rnaNorm[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
quantile(var_meta$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_path$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_tumor$var, probs = c(0,0.25,0.5,0.85,1)) 

data_meta=data_meta[,colnames(data_meta) %in% var_meta$feature[var_meta$var>quantile(var_meta$var,0.85)]]
data_path_clr=data_path_clr[,colnames(data_path_clr) %in% var_path$feature[var_path$var>quantile(var_path$var,0.85)]]
data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %in% var_tumor$feature[var_tumor$var>quantile(var_tumor$var,0.85)]]

sample <- as.data.frame(samples)
write.table(sample,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250730/samples.txt",sep="\t",quote = F,col.names = T,row.names = F)

#=========只纳入最终纳入的样本==================
rm(list=ls())
gc()

setwd("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/原始文件MOFAren")
library(readxl)
data_phenotype=readxl::read_xlsx("./RenxxInfo_250310.xlsx")


data_div <- read.table('D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/淑芬结果/renxx_BS_250916.txt',header=T)

data_phenotype=merge(data_phenotype,data_div[,c("sample","balance_score")],by.x="sample",by.y="sample",all=F)

##not_keep 1769
# data_phenotype=data_phenotype[data_phenotype$sample!="1769",]

data_score=data_phenotype


# final patients included
patients=data_phenotype$sample

# import rnaseq data
data_rna=read.csv("./Merge_118patient_gene_rpkm.csv",header = T,row.names = 1)



data_rna=data_rna[!duplicated(data_rna$Description),]
rownames(data_rna)=data_rna$Description
data_rna$Description=NULL
raw_data=data_rna
data_rnaNorm=raw_data
data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %like% "T"]
colnames(data_rnaNorm)=gsub("HCC","",colnames(data_rnaNorm))
colnames(data_rnaNorm)=gsub(".P.T","",colnames(data_rnaNorm))
colnames(data_rnaNorm)=gsub("T.*","",colnames(data_rnaNorm))

samples <- read.table("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250730/samples.txt",header=T)
data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %in% samples$samples]
data_rnaNorm=as.data.frame(t(data_rnaNorm))
data_rnaNorm=data_rnaNorm[,colSums(data_rnaNorm>0)>(nrow(data_rnaNorm)*0.8)]
data_rnaNorm=log(data_rnaNorm+0.1)
data_rnaNorm=data_rnaNorm[,!colnames(data_rnaNorm) %like% "MT-"]
data_rnaNorm=as.data.frame(t(data_rnaNorm))

# import metacyc
# data_path=read.table("./Metacyc_Renxx.txt",header = T,stringsAsFactors = F,row.names = 1,sep="\t",fill = T,comment.char = "", check.names = FALSE)
data_path <- read_tsv("./Renxx_HCC_Health_pathabundance_all.tsv")
data_path <- as.data.frame(data_path)
row.names(data_path) <- data_path[,1]
data_path$`# Pathway` <- NULL
colnames(data_path)=gsub("_Abundance.*$","",colnames(data_path))
colnames(data_path)=gsub("_pair.*$","",colnames(data_path))
# colnames(data_path)=gsub("\\.","-",colnames(data_path))
data_path <- data_path[,colnames(data_path)%in%data_score$sample]

data_path=data_path[!rownames(data_path) %like% "\\|",]
data_path=data_path[rowSums(data_path>0)>(ncol(data_path)*0.1),]
data_path=data_path[!rownames(data_path) %in% c("UNMAPPED","UNINTEGRATED"),]
data_path=as.data.frame(apply(data_path,2,function (x){
  x=x/sum(x)
  return(x)
}))
data_path_clr=transform_and_filter_taxa(data_path,samples_row = F,method = "clr",missing_filter = 0)
data_path_clr=as.data.frame(t(data_path_clr))
data_path_clr=data_path_clr[,colnames(data_path_clr) %in% samples$samples]

load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/淑芬结果/sc_renxx250618.RData")
data_meta <- sc_rxx
# import meta
# data_meta=data_metabolism_rxx
data_meta=as.data.frame(t(data_meta))
data_meta=data_meta[,colnames(data_meta) %in% samples$samples]

# overlap samples
# samples=intersect(colnames(data_rnaNorm),intersect(colnames(data_meta),colnames(data_path_clr)))
# data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %in% samples]
# data_meta=data_meta[,colnames(data_meta) %in% samples]
# data_path_clr=data_path_clr[,colnames(data_path_clr) %in% samples]

data_rnaNorm=data_rnaNorm[,order(colnames(data_rnaNorm))]
data_meta=data_meta[,order(colnames(data_meta))]
data_path_clr=data_path_clr[,order(colnames(data_path_clr))]

stopifnot(colnames(data_path_clr)==colnames(data_meta))
stopifnot(colnames(data_path_clr)==colnames(data_rnaNorm))

# top 10%-5% variance selection
data_path_clr=as.data.frame(t(data_path_clr))
data_meta=as.data.frame(t(data_meta))
data_rnaNorm=as.data.frame(t(data_rnaNorm))

var_meta=foreach(i=1:ncol(data_meta),.combine = rbind) %do%  {
  tmp.feature=colnames(data_meta)[i]
  tmp.var=var(data_meta[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_path=foreach(i=1:ncol(data_path_clr),.combine = rbind) %do%  {
  tmp.feature=colnames(data_path_clr)[i]
  tmp.var=var(data_path_clr[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
var_tumor=foreach(i=1:ncol(data_rnaNorm),.combine = rbind) %do%  {
  tmp.feature=colnames(data_rnaNorm)[i]
  tmp.var=var(data_rnaNorm[,tmp.feature])
  
  return.string=data.frame(feature=tmp.feature,var=tmp.var)
}
quantile(var_meta$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_path$var, probs = c(0,0.25,0.5,0.85,1)) 
quantile(var_tumor$var, probs = c(0,0.25,0.5,0.85,1)) 

data_meta=data_meta[,colnames(data_meta) %in% var_meta$feature[var_meta$var>quantile(var_meta$var,0.85)]]
data_path_clr=data_path_clr[,colnames(data_path_clr) %in% var_path$feature[var_path$var>quantile(var_path$var,0.85)]]
data_rnaNorm=data_rnaNorm[,colnames(data_rnaNorm) %in% var_tumor$feature[var_tumor$var>quantile(var_tumor$var,0.85)]]

MOFA_data1=list(as.matrix(t(data_path_clr)),as.matrix(t(data_meta)),as.matrix(t(data_rnaNorm)))
MOFAobject1 <- create_mofa(MOFA_data1)

data_opts <- get_default_data_options(MOFAobject1)
data_opts$scale_views=F
model_opts <- get_default_model_options(MOFAobject1)
model_opts$nu=m_factors=20
train_opts <- get_default_training_options(MOFAobject1)
train_opts$convergence_mode="fast"

MOFAobject1 <- prepare_mofa(
  object = MOFAobject1,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/","model.hdf5")

library(reticulate)
use_python("D:/miniconda3/envs/mofa_env/python.exe", required = TRUE)
MOFAobject1.trained <- run_mofa(MOFAobject1, outfile, use_basilisk = F)

model1 <- load_model(outfile)



plot_data_overview(model1)
Nsamples = sum(model1@dimensions$N)
sample_metadata <- data.frame(
  sample = samples_names(model1)[[1]]
)

sample_metadata <- left_join(sample_metadata,data_phenotype,by=c('sample'='sample'))

data_score=data_score[data_score$sample %in% samples$samples,]
samples_names(model1)[[1]]==data_score$sample
data_score$group=NULL
samples_metadata(model1) <- data_score

MOFA_variance=(model1@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(model1, x="view", y="factor",max_r2 =15)

#save.image(file = './mofa_ren_0630.rdata')

write.table(MOFA_variance,file = '../renxx_output_250924/output/ren_mofa_variance.txt',sep = "\t")


factors <- get_factors(model1, as.data.frame = T)
factors$factor=as.character(factors$factor)
factors$factor=gsub("Factor","",factors$factor)
factors$factor=as.numeric(factors$factor)

sample_metadata$Gender=as.factor(sample_metadata$Gender)
sample_metadata$Smoke=as.factor(sample_metadata$Smoke)
sample_metadata$Alcohol=as.factor(sample_metadata$Alcohol)



## 矫正了混杂因素之后肠型评分和factor的关联
MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lm(value ~ Age+Gender+BMI+Smoke+Alcohol+balance_score,data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="balance_score"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="balance_score"])
}
write.table(MOFA_correlation,file = '../renxx_output_250924/output/ren_mofa_correlation.txt',sep = "\t")


sample_metadata$BCLC[sample_metadata$BCLC %in% c('0','A')]='0/A'
sample_metadata$BCLC[sample_metadata$BCLC %in% c('B','C')]='B/C'

sample_metadata$`Tumor number`[sample_metadata$`Tumor number`==">3"]='4'

sample_metadata$`Tumor number`= as.character(sample_metadata$`Tumor number`)

sample_metadata=sample_metadata %>% dplyr::rename(Tumor_number=`Tumor number`)
sample_metadata$MVI=as.factor(sample_metadata$MVI)

i=2
## version 1 na.omit
MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
  tmp.data=tmp.data[c('value',feature_name)]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    tmp.data=na.omit(tmp.data)
    mm=lm(as.formula(paste0("value~",j)) ,data = tmp.data)
    mm=as.data.frame(summary(mm)$coef)
    return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                             Pvalue=mm$`Pr(>|t|)`[2])
  }
  
  
  return.string=cor1
}


## version 1
i=2
tmp.data=factors[factors$factor==i,]
colnames(sample_metadata) <- gsub("balance_score","balance_value",colnames(sample_metadata))
tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
tmp.data$group=NULL
feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
tmp.data=tmp.data[c('value',feature_name)]
tmp.data=na.omit(tmp.data)



p1 <- ggplot(data = tmp.data,aes(y=value,x=balance_value))+geom_point(color='#195881',alpha=0.7,size=5)+geom_smooth(method = 'lm',color='#195881')+
  ggpubr::stat_cor()+theme_bw()+
  theme(axis.title.x =element_text(size=60/.pt), 
        axis.title.y=element_text(size=60/.pt),
        axis.text.x = element_text(size=50/.pt,colour = 'black',angle = 45,hjust = 1),
        axis.text.y = element_text(size=50/.pt,colour = 'black'),
        legend.text = element_text(size=50/.pt,colour = 'black'),
        legend.title = element_text(size=50/.pt,colour = 'black'),
        plot.title=element_text(hjust=0.5,size=50/.pt))


## version1
saveRDS(p1,file = './新版肠型评分结果整理/ren_cor_heatmap_v1.rds')

save(MOFA_data1,file="D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/MOFA_data1.RData")
colnames(data_phenotype) <- gsub("balance_score","balance_value",colnames(data_phenotype))
colnames(data_score) <- gsub("balance_score","balance_value",colnames(data_score))
save(data_meta,data_path_clr,data_rnaNorm,data_phenotype,data_score,
     file="D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/mofa_rxx_250924.RData")


#=========整理附图====================
rm(list=ls())
gc()

outfile = file.path("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/","model.hdf5")
load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/mofa_rxx_250924.RData")
model1 <- load_model(outfile)

plot_data_overview(model1)
Nsamples = sum(model1@dimensions$N)
sample_metadata <- data.frame(
  sample = samples_names(model1)[[1]]
)

sample_metadata <- left_join(sample_metadata,data_phenotype,by=c('sample'='sample'))

# data_score=data_score[data_score$sample %in% samples$samples,]
samples_names(model1)[[1]]==data_score$sample
data_score$group=NULL
samples_metadata(model1) <- data_score

MOFA_variance=(model1@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(model1, x="view", y="factor",max_r2 =15)

#save.image(file = './mofa_ren_0630.rdata')
# MOFA_variance <- as.data.frame(MOFA_variance)
# MOFA_variance$factor <- row.names(MOFA_variance)
# # write.table(MOFA_variance,file = './新版肠型评分结果整理/ren_mofa_variance.txt',sep = "\t")
# library(writexl)
# write_xlsx(MOFA_variance,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/ren_mofa_variance.xlsx")

factors <- get_factors(model1, as.data.frame = T)
factors$factor=as.character(factors$factor)
factors$factor=gsub("Factor","",factors$factor)
factors$factor=as.numeric(factors$factor)

sample_metadata$Gender=as.factor(sample_metadata$Gender)
sample_metadata$Smoke=as.factor(sample_metadata$Smoke)
sample_metadata$Alcohol=as.factor(sample_metadata$Alcohol)



## 矫正了混杂因素之后肠型评分和factor的关联
MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lm(value ~ Age+Gender+BMI+Smoke+Alcohol+balance_value,data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="balance_value"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="balance_value"])
}
# write.table(MOFA_correlation,file = './新版肠型评分结果整理/ren_mofa_correlation.txt',sep = "\t")


sample_metadata$BCLC[sample_metadata$BCLC %in% c('0','A')]='0/A'
sample_metadata$BCLC[sample_metadata$BCLC %in% c('B','C')]='B/C'

sample_metadata$`Tumor number`[sample_metadata$`Tumor number`==">3"]='4'

sample_metadata$`Tumor number`= as.character(sample_metadata$`Tumor number`)

sample_metadata=sample_metadata %>% dplyr::rename(Tumor_number=`Tumor number`)
sample_metadata$MVI=as.factor(sample_metadata$MVI)

i=2
## version 1 na.omit
MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
  tmp.data=tmp.data[c('value',feature_name)]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    tmp.data=na.omit(tmp.data)
    mm=lm(as.formula(paste0("value~",j)) ,data = tmp.data)
    mm=as.data.frame(summary(mm)$coef)
    return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                             Pvalue=mm$`Pr(>|t|)`[2])
  }
  
  
  return.string=cor1
}


## version 1
i=2
tmp.data=factors[factors$factor==i,]
tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
tmp.data$group=NULL
feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
tmp.data=tmp.data[c('value',feature_name)]
tmp.data=na.omit(tmp.data)



p1 <- ggplot(data = tmp.data,aes(y=value,x=balance_value))+geom_point(color="#4B3DA2",fill="#4B3DA2",alpha=0.8,size=3)+
  geom_smooth(method = 'lm',color="#332583",fill="#9d94bd")+
  ggpubr::stat_cor()+theme_bw()+
  theme(axis.title.x =element_text(size=60/.pt), 
        axis.title.y=element_text(size=60/.pt),
        axis.text.x = element_text(size=50/.pt,colour = 'black'),
        axis.text.y = element_text(size=50/.pt,colour = 'black'),
        legend.text = element_text(size=50/.pt,colour = 'black'),
        legend.title = element_text(size=50/.pt,colour = 'black'),
        plot.title=element_text(hjust=0.5,size=50/.pt))

pdf("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/plot/FigSC_renxx_score_corelation.pdf",width = 5,height = 5)
p1
dev.off()



## 计算前5Factor转录组的富集结果-正相关--------
library(reactome.db)
library(ReactomePA)
library(clusterProfiler)
library(MOFAdata)
data(MSigDB_v6.0_C2_human) # MSigDB_v6.0_C2_human;reactomeGS; MSigDB_v6.0_C5_human
library(org.Hs.eg.db)
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

enrich_path <- enrich_path %>% dplyr::filter(factor=='Factor4') %>% dplyr::arrange(padj)



##提取通路富集的基因，确认哪些基因需要画图
gene_weight=plot_enrichment_detailed(run_enrichment(model1,
                                                    view = "view_3", factors = 1:5,
                                                    feature.sets = MSigDB_v6.0_C2_human,
                                                    sign = "positive",
                                                    statistical.test = "parametric",set.statistic = "rank.sum"
), 
factor = 2
)$data

gene_weight <- gene_weight %>% dplyr::filter(str_detect(pathway_long_name,"REACT"))
gene_weight$symbol=genes$SYMBOL[match(gene_weight$feature,genes$ENSEMBL)]
gene_weight <- gene_weight %>% dplyr::arrange(pathway,feature)

gene_weight <- gene_weight %>% dplyr::arrange(pathway,desc(feature.statistic))

library(clusterProfiler)

gene_name=unique(gene_weight$symbol[gene_weight$feature.statistic>0])

# write.table(data.frame(gene=unique(gene_weight$symbol[gene_weight$feature.statistic>0])),file = './gene.txt',row.names = F,quote = F) 





factor_gene=factors %>% dplyr::filter(factor==2) %>% dplyr::select(-group)

factor_gene=left_join(factor_gene,data_rnaNorm %>% rownames_to_column('sample'),by=c('sample'='sample'))

library(foreach)

gene_name=str_replace_all(gene_name,"-","_")
colnames(factor_gene)=str_replace_all(colnames(factor_gene),"-","_")

gene_correlation = foreach(i=gene_name,.combine = rbind) %do%  {
  tmp.data=factor_gene %>% dplyr::select(c(sample,value,i))
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  # mm=lm(as.formula(paste0('value ~ Age+Gender+BMI+Smoke+Alcohol+',i)),data = tmp.data)
  mm=lm(as.formula(paste0(i,' ~ value+Age+Gender+BMI+Smoke+Alcohol')),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="value"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="value"])
}
gene_correlation$padj <- p.adjust(gene_correlation$Pvalue,"BH")


# 用来整附表
gene_correlation = foreach(i=gene_name,.combine = rbind) %do%  {
  tmp.data=factor_gene %>% dplyr::select(c(sample,value,i))
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  # mm=lm(as.formula(paste0('value ~',i,'+Age+Gender+BMI+Smoke+Alcohol')),data = tmp.data)
  mm=lm(as.formula(paste0(i,' ~ value+Age+Gender+BMI+Smoke+Alcohol')),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  mm=mm[2,]
  mm$gene=i
  return.string=mm
}
gene_correlation$FDR <- p.adjust(gene_correlation$`Pr(>|t|)`,"BH")
gene_correlation$gene <- gsub("_","-",gene_correlation$gene)
write_xlsx(gene_correlation,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/renxx_gene_cor.xlsx")



#=========写个for循环批量画图并合并===============
target <- c("HLA_DOB","HLA_DQA1","HLA_DQA2","CD3D","CD3E","ICOS")
plist <- list()
for (i in 1:length(target)){
  tmp.data <- factor_gene[,c("value",target[i])]
  colnames(tmp.data)[2] <- "gene"
  plist[[i]] <- ggplot(data = tmp.data,aes(x=value,y=gene))+geom_point(color="#4B3DA2",fill="#4B3DA2",alpha=0.8,size=7)+
    geom_smooth(method = 'lm',color="#332583",fill="#9d94bd")+
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
pdf("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/plot/FigSD_renxx_gene_cor_ar.pdf",width=35,height=6)
plist[[1]]+plist[[2]]+plist[[3]]+plist[[4]]+plist[[5]]+plist[[6]]+plot_layout(ncol=6,nrow=1)
dev.off()


#==========开始画主图================
rm(list=ls())
outfile = file.path("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/","model.hdf5")
load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/mofa_rxx_250924.RData")
model1 <- load_model(outfile)


plot_data_overview(model1)
Nsamples = sum(model1@dimensions$N)
sample_metadata <- data.frame(
  sample = samples_names(model1)[[1]]
)

sample_metadata <- left_join(sample_metadata,data_phenotype,by=c('sample'='sample'))

# data_score=data_score[data_score$sample %in% samples$samples,]
samples_names(model1)[[1]]==data_score$sample
data_score$group=NULL
samples_metadata(model1) <- data_score

MOFA_variance=(model1@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(model1, x="view", y="factor",max_r2 =15)

#save.image(file = './mofa_ren_0630.rdata')
# MOFA_variance <- as.data.frame(MOFA_variance)
# MOFA_variance$factor <- row.names(MOFA_variance)
# # write.table(MOFA_variance,file = './新版肠型评分结果整理/ren_mofa_variance.txt',sep = "\t")
# library(writexl)
# write_xlsx(MOFA_variance,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/ren_mofa_variance.xlsx")

factors <- get_factors(model1, as.data.frame = T)
factors$factor=as.character(factors$factor)
factors$factor=gsub("Factor","",factors$factor)
factors$factor=as.numeric(factors$factor)

sample_metadata$Gender=as.factor(sample_metadata$Gender)
sample_metadata$Smoke=as.factor(sample_metadata$Smoke)
sample_metadata$Alcohol=as.factor(sample_metadata$Alcohol)



## 矫正了混杂因素之后肠型评分和factor的关联
MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lm(value ~ Age+Gender+BMI+Smoke+Alcohol+balance_value,data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="balance_value"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="balance_value"])
}
MOFA_correlation$FDR <- p.adjust(MOFA_correlation$Pvalue,"BH")
# write.table(MOFA_correlation,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/ren_mofa_correlation.txt',sep = "\t")


sample_metadata$BCLC[sample_metadata$BCLC %in% c('0','A')]='0/A'
sample_metadata$BCLC[sample_metadata$BCLC %in% c('B','C')]='B/C'

sample_metadata$`Tumor number`[sample_metadata$`Tumor number`==">3"]='4'

sample_metadata$`Tumor number`= as.character(sample_metadata$`Tumor number`)

sample_metadata=sample_metadata %>% dplyr::rename(Tumor_number=`Tumor number`)
sample_metadata$MVI=as.factor(sample_metadata$MVI)

i=2
## version 1 na.omit
MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
  tmp.data=tmp.data[c('value',feature_name)]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    tmp.data=na.omit(tmp.data)
    mm=lm(as.formula(paste0("value~",j)) ,data = tmp.data)
    mm=as.data.frame(summary(mm)$coef)
    return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                             Pvalue=mm$`Pr(>|t|)`[2])
  }
  
  
  return.string=cor1
}


## version 1
i=2
tmp.data=factors[factors$factor==i,]
tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
tmp.data$group=NULL
feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
tmp.data=tmp.data[c('value',feature_name)]
tmp.data=na.omit(tmp.data)



p1 <- ggplot(data = tmp.data,aes(y=value,x=balance_value))+geom_point(color="#4B3DA2",fill="#4B3DA2",alpha=0.8,size=3)+
  geom_smooth(method = 'lm',color="#332583",fill="#9d94bd")+
  ggpubr::stat_cor()+theme_bw()+
  theme(axis.title.x =element_text(size=60/.pt), 
        axis.title.y=element_text(size=60/.pt),
        axis.text.x = element_text(size=50/.pt,colour = 'black'),
        axis.text.y = element_text(size=50/.pt,colour = 'black'),
        legend.text = element_text(size=50/.pt,colour = 'black'),
        legend.title = element_text(size=50/.pt,colour = 'black'),
        plot.title=element_text(hjust=0.5,size=50/.pt))
p1

# 用来整附表
# MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
#   tmp.data=factors[factors$factor==i,]
#   tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
#   mm=lm(value ~ balance_value + Age+Gender+BMI+Smoke+Alcohol,data = tmp.data)
#   mm=as.data.frame(summary(mm)$coef)
#   mm <- mm[2,]
#   mm$Factor <- i
#   return.string=mm
# }
# #MOFA_correlation=MOFA_correlation[MOFA_correlation$Factor %in% c(1:10),] # in case that no FDR passed
# MOFA_correlation$FDR=p.adjust(MOFA_correlation$`Pr(>|t|)`,method = "BH")
# write_xlsx(MOFA_correlation,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/renxx_factor和评分相关性.xlsx")

#========开始画热图===========

MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','Cirrhosis','BCLC','balance_value')
  tmp.data=tmp.data[c('value',feature_name)]
  #tmp.data=na.omit(tmp.data)
  ##用混合模型计算
  basic=c('Age','Gender','BMI','Smoke','Alcohol')
  # tmp.data=tmp.data[c('value',feature_name,'patient')]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    if (j %in% basic){
      basic1 <- basic[!basic==j]
      mm=lm(as.formula(paste0("value~",j,"+",basic1[1],"+",basic1[2],"+",basic1[3],"+",basic1[4])) ,data = tmp.data)
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
    }
    else {
      mm=lm(as.formula(paste0("value~",j,"+Age+Gender+BMI+Smoke+Alcohol")) ,data = tmp.data)
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
    }
  }
  
  return.string=cor1
}

Plot_heatmap=MOFA_correlation2 %>% dplyr::select(-Pvalue) %>% pivot_wider(.,names_from = feature,values_from = Beta) %>% column_to_rownames('Factor')
Plot_heatmap <- Plot_heatmap[,c(1:8)]  ##把多的几列数据踢走
Plot_heatmap_pvalue <- MOFA_correlation2 %>% dplyr::select(-Beta) %>% pivot_wider(.,names_from = feature,values_from = Pvalue) %>% column_to_rownames('Factor')
Plot_heatmap_pvalue <- Plot_heatmap_pvalue[,c(1:8)]   ##把多的几列数据踢走

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
                  print(c(i, j, t(Part2_p)[i,j]))
                  }else if( t(Part2_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part2)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='black'))
                    print(c(i, j, t(Part2_p)[i,j]))
                    }else{
                  grid.text(sprintf("%.2f", t(Part2)[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                  
                }},
)



H2
draw(H2)



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

save(H1,H2,p2,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/renxx_heatmap.RData')
saveRDS(ht_list,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/renxx_heatmap_v2.rds')

#======整附表基因weight===============
weights_df <- get_weights(model1,
                          views = "view_3",     # 指定视图
                          factors = "Factor2",  # 指定 factor
                          as.data.frame = TRUE) # 返回 data.frame 格式

write_xlsx(weights_df,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/renxx_gene_weight.xlsx")

## 计算前5Factor转录组的富集结果-正相关--------
MOFA_weights <- get_weights(model1, as.data.frame = T)
MOFA_weights=MOFA_weights[MOFA_weights$factor=="Factor2",]
table(MOFA_weights$view)
med_path=MOFA_weights[MOFA_weights$view=="view_1",]
med_meta=MOFA_weights[MOFA_weights$view=="view_2",]
med_gene=MOFA_weights[MOFA_weights$view=="view_3",]

library(reactome.db)
library(ReactomePA)
library(clusterProfiler)
library(MOFAdata)
data(MSigDB_v6.0_C2_human) # MSigDB_v6.0_C2_human;reactomeGS; MSigDB_v6.0_C5_human
library(org.Hs.eg.db)
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

enrich_path <- enrich_path %>% dplyr::filter(factor=='Factor4') %>% dplyr::arrange(padj)



##提取通路富集的基因，确认哪些基因需要画图
gene_weight=plot_enrichment_detailed(run_enrichment(model1,
                                                    view = "view_3", factors = 1:5,
                                                    feature.sets = MSigDB_v6.0_C2_human,
                                                    sign = "positive",
                                                    statistical.test = "parametric",set.statistic = "rank.sum"
), 
factor = 2
)$data

gene_weight <- gene_weight %>% dplyr::filter(str_detect(pathway_long_name,"REACT"))
gene_weight$symbol=genes$SYMBOL[match(gene_weight$feature,genes$ENSEMBL)]
gene_weight <- gene_weight %>% dplyr::arrange(pathway,feature)

gene_weight <- gene_weight %>% dplyr::arrange(pathway,desc(feature.statistic))





library(clusterProfiler)

gene_name=unique(gene_weight$symbol[gene_weight$feature.statistic>0])

#write.table(data.frame(gene=unique(gene_weight$symbol[gene_weight$feature.statistic>0])),file = './gene.txt',row.names = F,quote = F) 





factor_gene=factors %>% dplyr::filter(factor==2) %>% dplyr::select(-group)

factor_gene=left_join(factor_gene,data_rnaNorm %>% rownames_to_column('sample'),by=c('sample'='sample'))

library(foreach)

gene_name=str_replace_all(gene_name,"-","_")
colnames(factor_gene)=str_replace_all(colnames(factor_gene),"-","_")

gene_correlation = foreach(i=gene_name,.combine = rbind) %do%  {
  tmp.data=factor_gene %>% dplyr::select(c(sample,value,i))
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lm(as.formula(paste0('value ~ Age+Gender+BMI+Smoke+Alcohol+',i)),data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)==i],Pvalue=mm$`Pr(>|t|)`[rownames(mm)==i])
}
gene_correlation$padjust=p.adjust(gene_correlation$Pvalue,method = 'BH') 
# write.table(gene_correlation,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/ren_mofa_gene_correlation',row.names = F,quote = F,sep = "\t")





## 火山图

library(RColorBrewer)
new_pallete=brewer.pal(5, "Paired") 

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
enrich_path <- enrich_path %>% dplyr::filter(str_detect(path,'REACT'))
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


enrich_path_all=rbind(enrich_path %>% dplyr::mutate(group='UP'),enrich_path2 %>% dplyr::mutate(group='Down'))
write_xlsx(enrich_path_all,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/附表/V2_250730/renxx_GSVA_pathway.xlsx")




enrich_path_all$factor=factor(enrich_path_all$factor,levels = paste0("Factor",seq(1,5)))




enrich_path_all$size=-log(enrich_path_all$padj)
enrich_path_all$padj=ifelse(enrich_path_all$group=='UP',-log10(enrich_path_all$padj),log10(enrich_path_all$padj))




#enrich_path_all=enrich_path_all %>% dplyr::filter(padj!=0)


# 
#  enrich_path_all=enrich_path_all %>% dplyr::rename(p_val_adj=size,cluster=factor,gene=path,avg_log2FC=padj)
# 
#  enrich_path_all$p_val=enrich_path_all$p_val_adj

#jjVolcano(enrich_path_all,topGeneN = 3,log2FC.cutoff = log(0.05),adjustP.cutoff = 0.05,)



#new_pallete=c('6e7ca5','7c9680','fe9014','fcd364',
#           'ffcad3','2482BC','7ca6be','92B096','8B516E','D65741','8298A3','52AEB5')


enrich_path_all_label=enrich_path_all %>% group_by(factor) %>% dplyr::arrange(desc(padj),.by_group = T) 

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

##调整一下，展示factor2中的适应性免疫通路
# enrich_path_all$label_new <- ifelse(enrich_path_all$path=="REACTOME_ADAPTIVE_IMMUNE_SYSTEM"&enrich_path_all$factor=="Factor2",
                                    # str_to_title("ADAPTIVE_IMMUNE_SYSTEM"),enrich_path_all$label)
# enrich_path_all$label_new <- gsub("Interferon_signaling","",enrich_path_all$label_new)

enrich_path_all_sig=enrich_path_all %>% dplyr::filter(type=='sig')
enrich_path_all_nosig=enrich_path_all %>% dplyr::filter(type=='nosig')

new_pallete <- c("#66C2A4","#FB6A4A","#1D91C0","#F768A1","#810F7C")
p1 <- ggplot()+
  geom_jitter(aes(x=factor,y=padj,size=size),alpha=0.7,color='grey',data = enrich_path_all_nosig,width = 0.3)+
  
  geom_jitter(aes(x=factor,y=padj,size=size,color=factor),alpha=0.7,data = enrich_path_all_sig,width = 0.3)+
  geom_tile(aes(x=factor,y=0,fill=factor),height=1,color='black',data = enrich_path_all%>% distinct(factor,.keep_all = T))+
  
  scale_fill_manual(values = new_pallete)+
  scale_color_manual(values =new_pallete)+
  geom_text(data=enrich_path_all %>% distinct(factor,.keep_all = T) , # 绘制中心分组标记图文本注释
            aes(x=factor, 
                y=0, 
                label=factor),
            color ="white")+
  ggrepel::geom_text_repel(aes(x=factor,y=padj,label = label),data = enrich_path_all_sig,size=2,max.overlaps = Inf)+
  theme_minimal() +
  # ylim(c(-32,32))+
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank())
p1
saveRDS(p1,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/rna_plot_ren_260206.rds')

Big_rna+p1+plot_layout(ncol=2,nrow=1,heights =c(1,1),widths = c(1,1))
pdf("Figure拼图/FigC_right_260208.pdf",width = 7.5,height = 5)
p1
dev.off()


#=======特征代谢图============
## metabolism
library(ggridges)



meta_merge <- merge(data_meta,factors %>% dplyr::filter(factor==2),by.x="row.names",by.y="sample",all=F)


meta_merge <- meta_merge %>% dplyr::select(-c(group,factor,Row.names))

meta_merge <- meta_merge %>%  dplyr::select(value,everything()) %>% 
  dplyr::mutate(group=ifelse(value>mean(value),"high","low"))
meta_merge$value=NULL
meta_merge <- meta_merge %>% dplyr::select(group,everything()) %>% 
  pivot_longer(.,colnames(.)[-1],names_to = 'meta',values_to = 'value')

meta_merge_try=meta_merge %>% dplyr::filter(meta %in% c('indolin-2-one','MethylIndole-3-acetate','Indoxyl sulfate'))
library(ggridges)
meta_merge_try$group=factor(meta_merge_try$group,levels = c('low','high'))
ptry_ren <- ggplot(data=meta_merge_try,aes(x=value,y=meta,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+ggtitle('Tryptophan metabolism')+
  scale_fill_manual(values = c("#A1D99B","#238B45"))+scale_color_manual(values = c("#A1D99B","#238B45"))+
  theme_classic()
ptry_ren

# meta_merge_acid=meta_merge %>% dplyr::filter(meta %in% c('Nor-Desoxycholic acid','Nutriacholic acid','Deoxycholic acid','Glycocholic acid',
#                                                          'Glycohyocholic acid','Cholic acid','Taurodeoxycholic acid','Taurocholic acid'))

meta_merge_acid=meta_merge %>% dplyr::filter(meta %in% c('Taurodeoxycholic acid','Taurocholic acid','Taurochenodeoxycholate'))
meta_merge_acid$meta <- gsub('Taurochenodeoxycholate',"Taurochenodesoxycholic acid",meta_merge_acid$meta)

library(ggridges)
meta_merge_acid$group=factor(meta_merge_acid$group,levels = c('low','high'))
meta_merge_acid$meta <- factor(meta_merge_acid$meta,levels = rev(c('Taurocholic acid','Taurodeoxycholic acid',"Taurochenodesoxycholic acid")))
pacid_ren <- ggplot(data=meta_merge_acid,aes(x=value,y=meta,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+ggtitle('Bile acid metabolism')+
  scale_fill_manual(values = c("#A7A3C4","#504394"))+scale_color_manual(values = c("#A7A3C4","#504394"))+
  theme_classic()
pacid_ren
save(ptry_ren,pacid_ren,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/Pmetabolism_ren.RData')

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


# library(ggrepel)
# library(ggbeeswarm)
# meta_path <- meta_path %>% dplyr::filter(factor %in% paste0('Factor',seq(1,5)))
# meta_path_label=meta_path_label%>% dplyr::filter(factor %in% paste0('Factor',seq(1,5)))
# Pmeta <- ggplot()+geom_beeswarm(data = meta_path,aes(x=factor,y=value,color=value),cex = 0.2,priority = 'descending',alpha=0.7,size=2)+
#   scale_color_viridis_c(option = 'B')+
#   geom_text_repel(aes(x=factor,y=value,label = feature),data = meta_path_label,show.legend = F,size=3,
#                   max.overlaps = Inf)+
#   geom_hline(yintercept = 0,lty=4,col="#666666",lwd=0.3)+
#   coord_flip()+ggthemes::theme_base()
# ggsave('./metapath.pdf')
# save(Pmeta,file = './Pmeta.RData')

meta_path_f <- data_path_clr %>% dplyr::select(c('PROPFERM-PWY: superpathway of L-alanine fermentation (Stickland reaction)',
                                                 'PWY-8188: L-alanine degradation VI (reductive Stickland reaction)',
                                                 'PWY-8189: L-alanine degradation V (oxidative Stickland reaction)')) %>% 
  rownames_to_column('sample') %>% pivot_longer(.,colnames(.)[-1],names_to = 'path',values_to = 'value')
meta_path_f <-left_join(meta_path_f,factors %>% dplyr::filter(factor=='2'),by=c('sample'='sample'))
meta_path_f$group=ifelse(meta_path_f$value.y>mean(meta_path_f$value.y),"high",'low')

meta_path_f$group <- factor(meta_path_f$group,levels = c("low","high"))
Pmeta_ren_p1 <- ggplot(data=meta_path_f,aes(x=value.x,y=path,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+
  scale_color_manual(values = c( "#A6CEE3", "#1F78B4"))+
  scale_fill_manual(values = c( "#A6CEE3", "#1F78B4"))+theme_classic()
Pmeta_ren_p1



meta_path_f2 <- data_path_clr %>% dplyr::select(c('LIPA-CORESYN-PWY: lipid A-core biosynthesis (E. coli K-12)',
                                                  'LPSSYN-PWY: superpathway of lipopolysaccharide biosynthesis')) %>% 
  rownames_to_column('sample') %>% pivot_longer(.,colnames(.)[-1],names_to = 'path',values_to = 'value')
meta_path_f2 <-left_join(meta_path_f2,factors %>% dplyr::filter(factor=='2'),by=c('sample'='sample'))
meta_path_f2$group=ifelse(meta_path_f2$value.y>mean(meta_path_f2$value.y),"high",'low')

meta_path_f2$group <- factor(meta_path_f2$group,levels = c("low","high"))
Pmeta_ren_p2 <- ggplot(data=meta_path_f2,aes(x=value.x,y=path,fill=group,color=group))+geom_density_ridges(alpha=0.4,scale=T)+
  scale_color_manual(values = c("#f9cd9c","#e18182"))+
  scale_fill_manual(values = c( "#f9cd9c","#e18182"))+theme_classic()

Pmeta_ren_p1
Pmeta_ren_p2
save(Pmeta_ren_p1,Pmeta_ren_p2,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/Pmeta_density_ren.RData')


library(patchwork)
ptry_ren+pacid_ren+Pmeta_ren_p1+Pmeta_ren_p2+plot_layout(ncol=2,nrow=2,widths = c(1,1,1,1),heights = c(1,1,1,1))


##=======整理代谢物附表==========
weights_df <- get_weights(model1,
                          views = "view_2",     # 指定视图
                          factors = "Factor2",  # 指定 factor
                          as.data.frame = TRUE) # 返回 data.fr
write_xlsx(weights_df,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/renxx_meta_weight.xlsx")

##=======整理菌代谢通路附表==========
weights_df <- get_weights(model1,
                          views = "view_1",     # 指定视图
                          factors = "Factor2",  # 指定 factor
                          as.data.frame = TRUE) # 返回 data.fr
write_xlsx(weights_df,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/renxx_path_weight.xlsx")


#========把RData导出pdf，拼图===========
setwd("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure")
rna_plot_ren <- readRDS("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/rna_plot_ren.rds")
pdf("FigC_right.pdf",width = 7.5,height = 5)
rna_plot_ren
dev.off()

renxx_heatmap_v2 <- readRDS("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/renxx_heatmap_v2.rds")
ht_list = H1 %v% p2%v% H2
pdf("FigB_right.pdf",width = 10,height = 8)
renxx_heatmap_v2
dev.off()

load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/Pmeta_density_ren.RData")
load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/Figure/Pmetabolism_ren.RData")
library(patchwork)

pdf("FigG_right.pdf",width=18,height = 10)
ptry_ren+pacid_ren+Pmeta_ren_p1+Pmeta_ren_p2+plot_layout(ncol=2,nrow=2,widths = c(1,1,1,1),heights = c(1,1,1,1))
dev.off()



#==========新版热图2601================
rm(list=ls())
outfile = file.path("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/","model.hdf5")
load("D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/mofa_rxx_250924.RData")
model1 <- load_model(outfile)


plot_data_overview(model1)
Nsamples = sum(model1@dimensions$N)
sample_metadata <- data.frame(
  sample = samples_names(model1)[[1]]
)

sample_metadata <- left_join(sample_metadata,data_phenotype,by=c('sample'='sample'))

# data_score=data_score[data_score$sample %in% samples$samples,]
samples_names(model1)[[1]]==data_score$sample
data_score$group=NULL
samples_metadata(model1) <- data_score

MOFA_variance=(model1@cache$variance_explained$r2_per_factor[[1]])
plot_variance_explained(model1, x="view", y="factor",max_r2 =15)

#save.image(file = './mofa_ren_0630.rdata')
# MOFA_variance <- as.data.frame(MOFA_variance)
# MOFA_variance$factor <- row.names(MOFA_variance)
# # write.table(MOFA_variance,file = './新版肠型评分结果整理/ren_mofa_variance.txt',sep = "\t")
# library(writexl)
# write_xlsx(MOFA_variance,"D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/ren_mofa_variance.xlsx")

factors <- get_factors(model1, as.data.frame = T)
factors$factor=as.character(factors$factor)
factors$factor=gsub("Factor","",factors$factor)
factors$factor=as.numeric(factors$factor)

sample_metadata$Gender=as.factor(sample_metadata$Gender)
sample_metadata$Smoke=as.factor(sample_metadata$Smoke)
sample_metadata$Alcohol=as.factor(sample_metadata$Alcohol)



## 矫正了混杂因素之后肠型评分和factor的关联
MOFA_correlation = foreach(i=1:nrow(MOFA_variance),.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  mm=lm(value ~ Age+Gender+BMI+Smoke+Alcohol+balance_value,data = tmp.data)
  mm=as.data.frame(summary(mm)$coef)
  
  return.string=data.frame(Factor=i,Beta=mm$Estimate[rownames(mm)=="balance_value"],Pvalue=mm$`Pr(>|t|)`[rownames(mm)=="balance_value"])
}
MOFA_correlation$FDR <- p.adjust(MOFA_correlation$Pvalue,"BH")
# write.table(MOFA_correlation,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/renxx_output_250924/output/ren_mofa_correlation.txt',sep = "\t")


sample_metadata$BCLC[sample_metadata$BCLC %in% c('0','A')]='0/A'
sample_metadata$BCLC[sample_metadata$BCLC %in% c('B','C')]='B/C'

sample_metadata$`Tumor number`[sample_metadata$`Tumor number`==">3"]='4'

sample_metadata$`Tumor number`= as.character(sample_metadata$`Tumor number`)

sample_metadata=sample_metadata %>% dplyr::rename(Tumor_number=`Tumor number`)
sample_metadata$MVI=as.factor(sample_metadata$MVI)

i=2
## version 1 na.omit
MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
  tmp.data=tmp.data[c('value',feature_name)]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    tmp.data=na.omit(tmp.data)
    mm=lm(as.formula(paste0("value~",j)) ,data = tmp.data)
    mm=as.data.frame(summary(mm)$coef)
    return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                             Pvalue=mm$`Pr(>|t|)`[2])
  }
  
  
  return.string=cor1
}


## AFP的beta scale一下
MOFA_correlation2_part <- MOFA_correlation2[MOFA_correlation2$feature=="AFP",]
MOFA_correlation2_part$Beta <- scale(MOFA_correlation2_part$Beta)
MOFA_correlation2 <- MOFA_correlation2[!MOFA_correlation2$feature=="AFP",]
MOFA_correlation2 <- rbind(MOFA_correlation2,MOFA_correlation2_part)


## version 1
i=2
tmp.data=factors[factors$factor==i,]
tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
tmp.data$group=NULL
feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize','balance_value')
tmp.data=tmp.data[c('value',feature_name)]
tmp.data=na.omit(tmp.data)



p1 <- ggplot(data = tmp.data,aes(y=value,x=balance_value))+geom_point(color="#4B3DA2",fill="#4B3DA2",alpha=0.8,size=3)+
  geom_smooth(method = 'lm',color="#332583",fill="#9d94bd")+
  ggpubr::stat_cor()+theme_bw()+
  theme(axis.title.x =element_text(size=60/.pt), 
        axis.title.y=element_text(size=60/.pt),
        axis.text.x = element_text(size=50/.pt,colour = 'black'),
        axis.text.y = element_text(size=50/.pt,colour = 'black'),
        legend.text = element_text(size=50/.pt,colour = 'black'),
        legend.title = element_text(size=50/.pt,colour = 'black'),
        plot.title=element_text(hjust=0.5,size=50/.pt))
p1


#========开始画热图===========

MOFA_correlation2 = foreach(i=1:15,.combine = rbind) %do%  {
  tmp.data=factors[factors$factor==i,]
  tmp.data=merge(tmp.data,sample_metadata,by.x="sample",by.y="sample",all=F)
  tmp.data$group=NULL
  feature_name=c('Age','Gender','BMI','Smoke','Alcohol','AFP','Tumor_number','TumorSize')
  tmp.data=tmp.data[c('value',feature_name)]
  #tmp.data=na.omit(tmp.data)
  ##用混合模型计算
  basic=c('Age','Gender','BMI','Smoke','Alcohol')
  # tmp.data=tmp.data[c('value',feature_name,'patient')]
  ##用混合模型计算
  cor1=foreach(j=feature_name,.combine = rbind) %do%{
    if (j %in% basic){
      basic1 <- basic[!basic==j]
      mm=lm(as.formula(paste0("value~",j,"+",basic1[1],"+",basic1[2],"+",basic1[3],"+",basic1[4])) ,data = tmp.data)
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
    }
    else {
      mm=lm(as.formula(paste0(j,"~value+Age+Gender+BMI+Smoke+Alcohol")) ,data = tmp.data)
      mm=as.data.frame(summary(mm)$coef)
      return.string=data.frame(Factor=i,feature=j,Beta=mm$Estimate[2],
                               Pvalue=mm$`Pr(>|t|)`[2])
    }
  }
  
  return.string=cor1
}


## AFP的beta scale一下
MOFA_correlation2_part <- MOFA_correlation2[MOFA_correlation2$feature=="AFP",]
MOFA_correlation2_part$Beta <- scale(MOFA_correlation2_part$Beta , center = FALSE, scale = TRUE)
MOFA_correlation2 <- MOFA_correlation2[!MOFA_correlation2$feature=="AFP",]
MOFA_correlation2 <- rbind(MOFA_correlation2,MOFA_correlation2_part)
scale(MOFA_correlation2_part$Beta , center = FALSE, scale = TRUE)

Plot_heatmap=MOFA_correlation2 %>% dplyr::select(-Pvalue) %>% pivot_wider(.,names_from = feature,values_from = Beta) %>% column_to_rownames('Factor')
Plot_heatmap <- Plot_heatmap[,c(1:8)]  ##把多的几列数据踢走
Plot_heatmap_pvalue <- MOFA_correlation2 %>% dplyr::select(-Beta) %>% pivot_wider(.,names_from = feature,values_from = Pvalue) %>% column_to_rownames('Factor')
Plot_heatmap_pvalue <- Plot_heatmap_pvalue[,c(1:8)]   ##把多的几列数据踢走

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
                  print(c(i, j, t(Part2_p)[i,j]))
                }else if( t(Part2_p)[i,j]<0.05){
                  grid.text(paste0(sprintf("%.2f", t(Part2)[i, j]),"\n","*") , x, y, gp = gpar(fontsize = 10,col='black'))
                  print(c(i, j, t(Part2_p)[i,j]))
                }else{
                  grid.text(sprintf("%.2f", t(Part2)[i, j]), x, y, gp = gpar(fontsize = 10,col='black'))
                  
                }},
)



H2
draw(H2)



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

save(H1,H2,p2,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/renxx_heatmap.RData')
saveRDS(ht_list,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/renxx_heatmap.rds')


pdf('D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/renxx_heatmap.pdf',width = 8)
draw(ht_list)
dev.off()


save(H1,H2,p2,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/renxx_heatmap_AFP.RData')
saveRDS(ht_list,file = 'D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/renxx_heatmap_AFP.rds')


pdf('D:/研一/任务/生信学习/宏基因组—12.7/最后整理2504/MOFA整图/新版热图2601/renxx_heatmap_AFP.pdf',width = 8)
draw(ht_list)
dev.off()
