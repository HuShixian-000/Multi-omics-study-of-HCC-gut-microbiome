rm(list=ls())
set.seed(12345)
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/")

#==========load function===============
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/functions.R")

#============MGS(筛平均丰度)======================
data_mgs_genus=read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/旧的16S代码/胡教授给的/Greengenes2.mgs.genus.ready.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
#单独挑出最后用来分析的HCC和Health
data <- read.table("../Final_input/species_raw_filter.txt")
sample <- gsub("\\.","-",colnames(data))
data_mgs_genus <- data_mgs_genus[,colnames(data_mgs_genus)%in%sample]
data_mgs_genus=data_mgs_genus[rowSums(data_mgs_genus!=0)>(ncol(data_mgs_genus)*0.1),]
#踢掉没有注释出来的属
data_mgs_genus<- data_mgs_genus %>%dplyr::mutate(genus=str_extract(row.names(data_mgs_genus),"g__.*"))
test=data_mgs_genus[data_mgs_genus$genus=="g__",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus=="g__",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus %like% "uncultured",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus%like% "Mitochondria",]
row.names(data_mgs_genus) <- data_mgs_genus$genus
data_mgs_genus$genus <- NULL
row.names(data_mgs_genus) <- gsub("g__","",row.names(data_mgs_genus))
data_mgs_genus=as.data.frame(apply(data_mgs_genus,2,function(x){
  x=x/sum(x)
  return(x)
}))

data_mgs_genus=as.data.frame(t(data_mgs_genus))
range(colMeans(data_mgs_genus))
data_mgs_genus=data_mgs_genus[,colMeans(data_mgs_genus)>0.001]
data_mgs_genus=data_mgs_genus[!duplicated(rownames(data_mgs_genus)),]


library(readxl)
MGSinfo_240829 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/MGSinfo_240829.xlsx")
data_mgs_genus <- merge(data_mgs_genus,MGSinfo_240829[,c("sampleID","sample")],by.x="row.names",by.y="sample")  ##改名字这步顺便把健康人踢掉了
row.names(data_mgs_genus) <- data_mgs_genus$sampleID
data_mgs_genus$sampleID <- NULL
data_mgs_genus$Row.names <- NULL
write.table(data_mgs_genus,"/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/final_16s_output/trim.mgs.greengenes.txt",
            row.names = T,col.names = T,sep="\t",quote=F)

#==========16S data===============
data_16S_genus=read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/旧的16S代码/胡教授给的/Greengenes2.16S.genus.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t",check.names = F)
colnames(data_16S_genus) <- gsub("963_Trelative_abundance","963_T",colnames(data_16S_genus))
colnames(data_16S_genus) <- gsub("2142","2120",colnames(data_16S_genus))
colnames(data_16S_genus) <- gsub("3775","3775/RH29",colnames(data_16S_genus))

#准备16S的信息表
group_16s <- as.data.frame(colnames(data_16S_genus))
colnames(group_16s) <- "sample"
group_16s$sampleID <- gsub("_.*$","",group_16s$sample)
group_16s$type <- ifelse(group_16s$sample%like%"T","Tumor",ifelse(group_16s$sample%like%"P","Peri","Liver"))

#踢掉非HCC和两个只有癌旁的
library(readxl)
no_HCC <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/input_data/大样本队列非HCC患者.xlsx",sheet = "16s中非HCC患者")
group_16s <- group_16s[!group_16s$sampleID%in%no_HCC$术中标本编号,]
group_16s <- group_16s[!group_16s$sampleID%in%c("1944","1466"),]

#重复送测样本剔除第二次送测的： 1072_P_2/1217_T1_2 1217_T2_2 1223_T1_2 897_T1_2 
group_16s <- group_16s[!group_16s$sample%in%c("1072_P_2","1217_T1_2","1217_T2_2","1223_T1_2","897_T1_2"),]

#踢掉1个年龄小于18岁的
group_16s <- group_16s[!group_16s$sampleID%in%c("2088"),]

#踢掉rxx队列重复的
rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon1.txt")
group_16s <- group_16s[!group_16s$sampleID%in%rxx$sample,]

#踢掉4个后面提了MGS样本的
#group_16s <- group_16s[!group_16s$sampleID%in%c("1943","1811","2088","3743"),]

data_16S_genus <- data_16S_genus[,colnames(data_16S_genus)%in%group_16s$sample]
data_16S_genus=as.data.frame(t(data_16S_genus))


# remove low count(会剔样本)
count=data.frame(count=rowSums(data_16S_genus))#行为样本
data_16S_genus=data_16S_genus[rownames(data_16S_genus) %in% row.names(count[count$count>10000,,drop=F]),]

##去污染
contamination=read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/input_data/Known.contaminates.list.txt",header = F,stringsAsFactors = F)
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
data_16S_genus=data_16S_genus[,!colnames(data_16S_genus) %like% "Chloroplast"]
#去掉没有注释到属的
data_16S_genus <- as.data.frame(t(data_16S_genus))
data_16S_genus<- data_16S_genus %>%dplyr::mutate(genus=str_extract(row.names(data_16S_genus),"g__.*"))
data_16S_genus=data_16S_genus[!data_16S_genus$genus=="g__",]
row.names(data_16S_genus) <- data_16S_genus$genus
data_16S_genus$genus <- NULL
row.names(data_16S_genus) <- gsub("g__","",row.names(data_16S_genus))

# 取mgs也能检测到的属
data_16S_genus <- as.data.frame(t(data_16S_genus))
# data_16S <- data_16S_genus
# data_16S_genus=data_16S_genus[,colnames(data_16S_genus) %in% colnames(data_mgs_genus)]


# remove low presence and abundance
data_16S_genus=data_16S_genus[,colSums(data_16S_genus>0)>(nrow(data_16S_genus)*0.2)]
# data_16S=data_16S[,colSums(data_16S>0)>(nrow(data_16S)*0.2)]

#最后剩下的属太少了，这些属在有些样本中都检测不到，先踢掉
count=data.frame(count=rowSums(data_16S_genus))

data_16S_genus=data_16S_genus[rownames(data_16S_genus) %in% rownames(count)[count$count>0],]

# re-calucalte proportion
data_16S_genus <- count_to_composition(data_16S_genus,samples_row = T)
range(colMeans(data_16S_genus))
# 0.0006057024 0.0514067182
# data_16S <- count_to_composition(data_16S,samples_row = T)
# range(colMeans(data_16S))
#踢掉低检出丰度
data_16S_genus <- data_16S_genus[,colMeans(data_16S_genus)>0.001]

##clr转化
data_clr <- transform_and_filter_taxa(data_16S_genus,samples_row = T,method = "clr",missing_filter = 0)
# data_16S_clr <- transform_and_filter_taxa(data_16S,samples_row = T,method = "clr",missing_filter = 0)
write.table(data_clr,"./16S_250509/final_16s_output/data_clr1.txt",sep="\t",quote=F,col.names = T,row.names = T)
# write.table(data_16S_clr,"data/data_16S_clr.txt",sep="\t",quote=F,col.names = T,row.names = T)

group_16s <- group_16s[group_16s$sample%in%row.names(data_16S_genus),]
write.table(group_16s,"./16S_250509/final_16s_output/group_16s1.txt",sep="\t",quote=F,col.names = T,row.names = T)

#====重跑开始处====
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/")
#======colorset=============
col1 <- c("#2482BC","#f99170","#e73133")
col2 <- c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
col3 <- c("#fcd364","#b66d32")
col4 <- c("#b39658","#b66d32")

#=========每个患者随机抽取1个肿瘤样本并整理样本信息表=============
data_16S <- read.table("./final_16s_output/data_16S_genus1.txt",header=T)
data_clr <- read.table("./final_16s_output/data_clr1.txt",header = T)

group_16s <- read.table("./final_16s_output/group_16s1.txt",header=T)
# group_16s1 <- read.table("./final_16s_output/group_16s1.txt",header=T)


MGSinfo_250310 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/MGSinfo_250310.xlsx")
MGSinfo_250310 <- MGSinfo_250310[!MGSinfo_250310$sampleID%in%c("1943","1811","2088","3743"),]
MGSinfo_250310$BMI[is.na(MGSinfo_250310$BMI)==T] <- median(MGSinfo_250310$BMI,na.rm=T)

score <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt",header=T)
score$sampleID <- gsub("2142","2120",score$sampleID)
score$sampleID <- gsub("1378574","2138",score$sampleID)

#加载批次分组
batch <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/input_data/pici_tissue.txt",sep="\t",header=T)
colnames(batch)[2] <- "batch"
batch$sample <- gsub("2142","2120",batch$sample)
batch$sample <- gsub("3775","3775/RH29",batch$sample)
row.names(batch) <- batch$sample

group_16s <- group_16s[group_16s$sampleID%in%MGSinfo_250310$sampleID,]

group_tumor <- group_16s[group_16s$type=="Tumor",]

set.seed(1711)
group_tumor_all <- group_tumor
group_tumor <- group_tumor %>%
  group_by(sampleID) %>%      # 按患者分组
  sample_n(1) %>%              # 随机抽取1个样本
  ungroup()                    # 取消分组

group_tumor <- merge(group_tumor,score[,c("sampleID","Cluster","balance_value")],by="sampleID",all.x=T)
group_tumor <- merge(group_tumor,MGSinfo_250310[,c("sampleID","Age","Gender","BMI","Smoke","Drink")],by="sampleID",all.x=T)
group_tumor$Cluster_new <- ifelse(group_tumor$Cluster=="4","E.coli","non_E.coli")
group_tumor <- merge(group_tumor,batch,by="sample")
group_tumor$batch <- as.factor(group_tumor$batch)
group_tumor[group_tumor=="Yes"] <- 1
group_tumor[group_tumor=="No"] <- 0
group_tumor[is.na(group_tumor)==T] <- 0.5
group_tumor1 <- group_tumor

group_nontumor <- group_16s[group_16s$type!="Tumor",]
# nontumor 直接取的去重部分
group_nontumor <- group_nontumor[!duplicated(group_nontumor$sampleID),]
group_nontumor <- merge(group_nontumor,score[,c("sampleID","Cluster","balance_value")],by="sampleID",all.x=T)
group_nontumor <- merge(group_nontumor,MGSinfo_250310[,c("sampleID","Age","Gender","BMI","Smoke","Drink")],by="sampleID",all.x=T)
group_nontumor$Cluster_new <- ifelse(group_nontumor$Cluster=="4","E.coli","non_E.coli")
group_nontumor <- merge(group_nontumor,batch,by="sample")
group_nontumor$batch <- as.factor(group_nontumor$batch)
group_nontumor[group_nontumor=="Yes"] <- 1
group_nontumor[group_nontumor=="No"] <- 0
group_nontumor[is.na(group_nontumor)==T] <- 0.5
group_nontumor1 <- group_nontumor

# write.table(group_tumor,"group_tumor.test1.txt",col.names = T,row.names = F,sep="\t",quote = F)
# write.table(group_nontumor,"group_nontumor.test1.txt",col.names = T,row.names = F,sep="\t",quote = F)

#=========Fig3A alpha diversity========================
data_tumor_16S <- data_16S[row.names(data_16S)%in%group_tumor$sample,]

##################肿瘤
diversity<-vegan::diversity

shannon_tumor<- data.frame(Richness=specnumber(data_tumor_16S),
                           shannon=diversity(data_tumor_16S, index="shannon"),
                           simpson=diversity(data_tumor_16S, index="simpson"),
                           inverse_simpson=diversity(data_tumor_16S, index="invsimpson"),
                           Evenness=diversity(data_tumor_16S, index="shannon")/log(specnumber(data_tumor_16S)))
group_tumor <- merge(group_tumor,shannon_tumor,by.x="sample",by.y="row.names")

group_tumor$batch <- as.factor(group_tumor$batch)
group_tumor <- na.omit(group_tumor)
group_tumor$Enterotype <- ifelse(group_tumor$Cluster=="4","ETE","ETnonE")

summary(lm(Richness~Enterotype+batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))
group_tumor$ad_richness <- resid(lm(Richness~batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))

summary(lm(shannon~Enterotype+batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))
group_tumor$ad_shannon <- resid(lm(shannon~batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))

group_tumor$ad_simpson <- resid(lm(simpson~batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))

group_tumor$ad_inverse_simpson <- resid(lm(inverse_simpson~batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))

group_tumor$ad_Evenness <- resid(lm(Evenness~batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))


##########癌旁
data_nontumor_16S <- data_16S[row.names(data_16S)%in%group_nontumor$sample,]

shannon_nontumor<- data.frame(Richness=specnumber(data_nontumor_16S),
                              shannon=diversity(data_nontumor_16S, index="shannon"),
                              simpson=diversity(data_nontumor_16S, index="simpson"),
                              inverse_simpson=diversity(data_nontumor_16S, index="invsimpson"),
                              Evenness=diversity(data_nontumor_16S, index="shannon")/log(specnumber(data_nontumor_16S)))
group_nontumor <- merge(group_nontumor,shannon_nontumor,by.x="sample",by.y="row.names")

group_nontumor$batch <- as.factor(group_nontumor$batch)
group_nontumor <- na.omit(group_nontumor)
group_nontumor$Enterotype <- ifelse(group_nontumor$Cluster=="4","ETE","ETnonE")


group_nontumor$ad_richness <- resid(lm(Richness~batch+Age+Gender+BMI+Smoke+Drink,data=group_nontumor))

summary(lm(shannon~Enterotype+batch+Age+Gender+BMI+Smoke+Drink,data=group_nontumor))
group_nontumor$ad_shannon <- resid(lm(shannon~batch+Age+Gender+BMI+Smoke+Drink,data=group_nontumor))

group_nontumor$ad_simpson <- resid(lm(simpson~batch+Age+Gender+BMI+Smoke+Drink,data=group_nontumor))

group_nontumor$ad_inverse_simpson <- resid(lm(inverse_simpson~batch+Age+Gender+BMI+Smoke+Drink,data=group_nontumor))

group_nontumor$ad_Evenness <- resid(lm(Evenness~batch+Age+Gender+BMI+Smoke+Drink,data=group_nontumor))


#############画图
group_tumor1 <- group_tumor[,c("sample","sampleID", "type","Cluster","Cluster_new"  , "ad_richness","ad_shannon","ad_simpson" ,"ad_inverse_simpson",
                               "ad_Evenness"  )]
group_tumor1$Enterotype <- ifelse(group_tumor1$Cluster=="4","ETE","ETnonE")
group_tumor1$Type <- "Tumor"

group_nontumor1 <- group_nontumor[,c("sample","sampleID", "type","Cluster","Cluster_new" ,"ad_richness","ad_shannon","ad_simpson","ad_inverse_simpson",
                                     "ad_Evenness"  )]
group_nontumor1$Enterotype <- ifelse(group_nontumor1$Cluster=="4","ETE","ETnonE")
group_nontumor1$Type <- "Peri"

plot_adiv <- rbind(group_tumor1,group_nontumor1)
plot_adiv$Type <- factor(plot_adiv$Type,levels = c("Peri","Tumor"))
plot_adiv$Enterotype <- factor(plot_adiv$Enterotype,levels = c("ETnonE","ETE"))


#====Fig3A_alpha_diversity====
library(ggh4x)
pdf("../../Fig3/Fig3A_adiversity.pdf",width = 6,height = 5)
for (i in colnames(plot_adiv)[6:10]){
  tmp.data <- plot_adiv[,c("Type","Enterotype",i)]
  colnames(tmp.data)[3] <- "diversity"
  p <- ggplot(tmp.data, aes(x=Enterotype, y=diversity, color=Enterotype)) + 
    geom_violin(trim=FALSE,aes(fill=Enterotype,alpha=0.8))+
    geom_boxplot(width=0.15, fill="white")+
    facet_grid2(. ~ Type, scales = "free", space = "free",strip = strip_themed(
      background_x = list(
        "Peri" = element_rect(fill = col1[1], color = "white"),
        "Tumor" = element_rect(fill = col1[3], color = "white")),
      text_x = list(
        "Peri" = element_text(color = "white", face = "bold"),
        "Tumor" = element_text(color = "white", face = "bold")))) +
    theme_classic()+
    stat_compare_means(comparisons =list(c("ETE","ETnonE")))+
    theme(strip.text.y = element_text(angle = 0),
          legend.position = "none")+
    scale_color_manual(values = col4)+
    scale_fill_manual(values=col3)+labs(y=i)
  print(p)
}
dev.off()

#====α diversity & balance value====
cor.test(group_tumor$ad_shannon,group_tumor$balance_value,method="spearman")
ggplot(group_tumor,aes(balance_value,ad_shannon)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +  stat_cor(method = 'spearman', label.y=-0.505)+theme_bw()
group_tumor$Enterotype <- factor(group_tumor$Enterotype,levels = c("ETnonE","ETE"))
pdf("../../Fig3/Fig3S_tumor_adiver&BS.pdf",width = 6,height = 5)
ggplot(data = group_tumor,aes(x=balance_value,y=ad_shannon))+
  geom_point(aes(color=Enterotype))+geom_smooth(method = 'lm',color=col4[2],fill=col3[1])+
  scale_color_manual(values = col3)+
  stat_cor(method = 'spearman', label.y=0.6)+
  theme(panel.background =element_blank(),legend.position.inside = c(0.15,0.15),
        axis.line=element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()


cor.test(group_nontumor$ad_shannon,group_nontumor$balance_value,method="spearman")
ggplot(group_nontumor,aes(balance_value,ad_shannon)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +  stat_cor(method = 'spearman', label.y=-0.505)+theme_bw()
group_nontumor$Enterotype <- factor(group_nontumor$Enterotype,levels = c("ETnonE","ETE"))
pdf("../../Fig3/Fig3S_nontumor_adiver&BS.pdf",width = 8,height = 5)
ggplot(data = group_nontumor,aes(x=balance_value,y=ad_shannon))+
  geom_point(aes(color=Enterotype))+geom_smooth(method = 'lm',color=col4[2],fill=col3[1])+
  scale_color_manual(values = col3)+
  stat_cor(method = 'spearman', label.y=0.6)+
  theme(panel.background =element_blank(),legend.position.inside = c(0.15,0.15),
        axis.line=element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

#==========Fig3B_beta_diversity==================
library(ape)
data_tumor <- data_16S[row.names(data_16S)%in%group_tumor$sample,]

dist_BC <- vegdist(data_tumor,method = "bray")

pcoa <- pcoa(dist_BC,correction = "none", rn = NULL) 
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]
#groups<-group_tumor1[match(rownames(pcoa$vectors),group_tumor1$sample),]
groups<-group_tumor[match(rownames(pcoa$vectors),group_tumor$sample),]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups[,c("Cluster" , "Cluster_new","balance_value","Age","Gender","BMI","Smoke","Drink","sampleID")])
colnames(pcoadata) <- c("sample","PC1","PC2","Cluster" , "Cluster_new","balance_value","Age","Gender","BMI","Smoke","Drink","sampleID")

#pcoadata$Cluster <- as.factor(pcoadata$Cluster)
pcoadata$Enterotype <- ifelse(pcoadata$Cluster=="4","ETE","ETnonE")
#Adonis test
set.seed(12345)
otu.table <- adonis2(dist_BC~Enterotype,data = pcoadata,permutations = 999)
otu.table.bs <- adonis2(dist_BC~balance_value,data = pcoadata,permutations = 999)

batch <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/input_data/pici_tissue.txt",sep="\t",header=T)
colnames(batch)[2] <- "batch"
batch$sample <- gsub("2142","2120",batch$sample)
batch$sample <- gsub("3775","3775/RH29",batch$sample)
row.names(batch) <- batch$sample
pcoadata <- merge(pcoadata,batch,by="sample")
pcoadata$batch <- as.factor(pcoadata$batch)
pcoadata$sampleID <- gsub("_.*$","",pcoadata$sample)
pcoadata$ad_PC1 <- resid(lm(PC1~batch+Age+Gender+BMI+Smoke+Drink,data=pcoadata))
pcoadata$ad_PC2 <- resid(lm(PC2~batch+Age+Gender+BMI+Smoke+Drink,data=pcoadata))
pcoadata$Enterotype <- factor(pcoadata$Enterotype,levels = c("ETnonE","ETE"))

otu.table1 <- otu.table
pcoadata1 <- pcoadata

#=======画肿瘤第一种曲线密度图=============
library(ggplot2)
library(ape)
library(tidyverse)
library(multcomp)
library(extrafont)
library(patchwork)
library(vegan)
library(ggtern)
library(reshape2)
library(ggsci)
library(scales)
library(patchwork)

bdiversity1 <- ggplot(pcoadata1,aes(ad_PC1, ad_PC2, color = Enterotype))+
  geom_point(size = 2)+#ggtitle("France")+
  stat_ellipse(aes(ad_PC1,ad_PC2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )"),
       color="Enterotype") +
  scale_color_manual(values =c(col3))+
  #guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position =c(.15, .15),
        legend.title = element_text(size =12,family ="sans",color="black"),
        legend.text = element_text(size =9,family ="sans",color="black"),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))+
  annotate("text",x=-0.4,y=0.2,label = paste("PERMANOVA:\np=",otu.table1[1,5],
                                             "\nR2 = ", round(otu.table1[1,3],4),sep = ""
  ))+ggtitle("Tumor")
bdiversity1

bdiversity1 <- ggExtra::ggMarginal(bdiversity1, type = "density", groupColour = F, groupFill = T,
                                   xparams = list(alpha = 0.5,color=NA),
                                   yparams = list(alpha = 0.5,color=NA))

#ggsave("../../Fig3/test.pdf",device = cairo_pdf,width =5.5, height =5,bdiversity1)

#==== Tumor β diversity & balance value=====
bdiversity1.bs <- ggplot(pcoadata1,aes(ad_PC1,ad_PC2, color = balance_value))+
  geom_point(size = 2)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )"),
       color="balance_value") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.2,9.1))+
  theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position ="top",
        # legend.title = element_text(size =12,family ="bold",color="black"),
        #  legend.text = element_text(size =9,family ="sans",color="black"),
        legend.key = element_blank(), 
        # axis.text.x = element_text(size = 10,family ="sans",color="black"),
        #axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))+
  annotate("text",x=-0.3,y=0.1,label=paste("PERMANOVA:\npvalue=",signif(otu.table.bs[1,5],2),
                                           "\nR2=",signif(otu.table.bs[1,3],2)))+ggtitle("Tumor")
bdiversity1.bs
ggsave("../../Fig3/tumor_betadiver_BS.pdf",device = cairo_pdf,width =5.5, height =5,bdiversity1.bs)

#========画非肿瘤的PCA图密度图===========
data_nontumor <- data_16S[row.names(data_16S)%in%group_nontumor$sample,]

dist_BC <- vegdist(data_nontumor,method = "bray")

pcoa <- pcoa(dist_BC,correction = "none", rn = NULL) 
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]
groups<-group_nontumor[match(rownames(pcoa$vectors),group_nontumor$sample),]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups[,c("Cluster" , "Cluster_new","balance_value","Age","Gender","BMI","Smoke","Drink","sampleID")])
colnames(pcoadata) <- c("sample","PC1","PC2","Cluster" , "Cluster_new","balance_value","Age","Gender","BMI","Smoke","Drink","sampleID")

pcoadata$Cluster <- as.factor(pcoadata$Cluster)
pcoadata$Enterotype <- ifelse(pcoadata$Cluster=="3","ETE","ETnonE")
#Adonis test
otu.table <- adonis2(dist_BC~Enterotype,data = pcoadata,permutations = 999)
otu.table.bs2 <- adonis2(dist_BC~balance_value,data = pcoadata,permutations = 999)

pcoadata <- merge(pcoadata,batch,by="sample")
pcoadata$batch <- as.factor(pcoadata$batch)
pcoadata$sampleID <- gsub("_.*$","",pcoadata$sample)
pcoadata$ad_PC1 <- resid(lm(PC1~batch+Age+Gender+BMI+Smoke+Drink,data=pcoadata))
pcoadata$ad_PC2 <- resid(lm(PC2~batch+Age+Gender+BMI+Smoke+Drink,data=pcoadata))

tmp_comparisons <- list(c("ETE","ETnonE"))
pcoadata$Enterotype <- factor(pcoadata$Enterotype,levels = c("ETnonE","ETE"))

pcoadata2 <- pcoadata
otu.table2 <- otu.table

#=======画非肿瘤第一种曲线密度图=============
library(ggplot2)
library(ape)
library(tidyverse)
library(multcomp)
library(extrafont)
library(patchwork)
library(vegan)
library(ggtern)
library(reshape2)
library(ggsci)
library(scales)
library(patchwork)

bdiversity2 <- ggplot(pcoadata2,aes(ad_PC1, ad_PC2, color = Enterotype))+
  geom_point(size = 2)+#ggtitle("France")+
  stat_ellipse(aes(ad_PC1,ad_PC2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )"),
       color="Enterotype") +
  scale_color_manual(values =c(col3))+
  #guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position =c(.15, .15),
        legend.title = element_text(size =12,family ="bold",color="black"),
        legend.text = element_text(size =9,family ="sans",color="black"),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))+
  annotate("text",x=-0.4,y=0.25,label = paste("PERMANOVA:\np=",otu.table2[1,5],
                                              "\nR2 = ", round(otu.table2[1,3],4),sep = ""
  ))+ggtitle("Peri")
bdiversity2

bdiversity2 <- ggExtra::ggMarginal(bdiversity2, type = "density", groupColour = F, groupFill = T,
                                   xparams = list(alpha = 0.5,color=NA),
                                   yparams = list(alpha = 0.5,color=NA))

bdiversity_all <- ggarrange(bdiversity2,bdiversity1, ncol = 2, nrow = 1)
ggsave("../../Fig3/Fig3B_betadiversity.pdf",device = cairo_pdf,width =11, height =5,bdiversity_all)

#====Peri β diversity & balance value=====
bdiversity2.bs <- ggplot(pcoadata1,aes(ad_PC1,ad_PC2, color = balance_value))+
  geom_point(size = 2)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )"),
       color="balance_value") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.6,9.45))+
  theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position ="top",
        # legend.title = element_text(size =12,family ="bold",color="black"),
        #  legend.text = element_text(size =9,family ="sans",color="black"),
        legend.key = element_blank(), 
        # axis.text.x = element_text(size = 10,family ="sans",color="black"),
        #axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))+
  annotate("text",x=-0.3,y=0.1,label=paste("PERMANOVA:\npvalue=",signif(otu.table.bs2[1,5],2),
                                           "\nR2=",signif(otu.table.bs2[1,3],2)))+ggtitle("Peri")
bdiversity2.bs
ggsave("../../Fig3/nontumor_betadiver_BS.pdf",device = cairo_pdf,width =5.5, height =5,bdiversity2.bs)


#========Fig3F_prevalance================
# group_tumor <- read.table("group_tumor.ready.txt",header=T)
# group_nontumor <- read.table("group_nontumor.ready.txt",header=T)
s16_tumor_E.coli <- data_tumor_16S[row.names(data_tumor_16S)%in%group_tumor[group_tumor$Cluster_new=="E.coli","sample"],]
s16_tumor_nonE.coli<- data_tumor_16S[row.names(data_tumor_16S)%in%group_tumor[group_tumor$Cluster_new=="non_E.coli","sample"],]

res_tumor <- list()
i=1
for ( i in 1:ncol(s16_tumor_E.coli)){
  bac <- colnames(s16_tumor_E.coli)[i]
  a <- s16_tumor_E.coli[,i,drop=F]
  colnames(a) <- "bac"
  a$bac <- ifelse(a$bac==0,"0","1") # 只要a菌出现，则计数1（同卡检出率的逻辑）
  b <- s16_tumor_nonE.coli[,i,drop=F]
  colnames(b) <- "bac"
  b$bac <- ifelse(b$bac==0,"0","1")
  a.freq <- as.data.frame(table(a$bac))
  row.names(a.freq) <- a.freq$Var1
  a.freq$Var1 <- NULL
  colnames(a.freq) <- "E.coli"
  b.freq <- as.data.frame(table(b$bac))
  row.names(b.freq) <- b.freq$Var1
  b.freq$Var1 <- NULL
  colnames(b.freq) <- "non_E.coli"
  if(nrow(a.freq)>1){
    freq <- cbind(a.freq,b.freq)
    freq <- as.data.frame(t(freq))
    pval <- chisq.test(freq)[["p.value"]]
    res_tumor[[i]] <- data.frame(genus=colnames(s16_tumor_E.coli)[i],pvalue=pval,
                                 ratio=(a.freq[2,1]/nrow(a))/(b.freq[2,1]/nrow(b)),
                                 prevalance.E.coli=a.freq[2,1]/nrow(a),# a菌在ETE中的流行度
                                 prevalance.nonE.coli=b.freq[2,1]/nrow(b))} # a菌在ETnonE中的流行度
  else{
    pre <- ifelse(row.names(a.freq)==1,1,0) # 如果a.frq只有一行说明该菌在所有ETE样本中的prevalence为1（不存在0的情况因为筛了检出率）
    a.freq[2,] <- 0
    row.names(a.freq)[2] <- 0
    freq <- cbind(a.freq,b.freq)
    freq <- as.data.frame(t(freq)) # 不排除b.freq中出现a菌在nonETE中prevalence为1的情况，但实际并没有，所以就没有报错
    pval <- chisq.test(freq)[["p.value"]]
    res_tumor[[i]] <- data.frame(genus=colnames(s16_tumor_E.coli)[i],pvalue=pval,
                                 ratio=pre/(b.freq[2,1]/nrow(b)),
                                 prevalance.E.coli=pre,
                                 prevalance.nonE.coli=b.freq[2,1]/nrow(b))
  }
}
res_tumor <- do.call("rbind",res_tumor)
res_tumor$p.adj <- p.adjust(res_tumor$pvalue,"BH")

res_tumor$significance <- ifelse(res_tumor$pvalue<0.05,"pvalue<0.05","pvalue>=0.05")
#res_tumor$label <- ifelse(res_tumor$genus%in%intersect(colnames(data_tumor),res_tumor[res_tumor$p.adj<0.05,"genus"]),
#                          res_tumor$genus,NA)

#target <- read.table("/groups/ProHuShiX/home/liuyuyao/BIGMeta_16S/translocate/PCA/prevalance_abundance_overlap.txt")
#target <- c(  "Lawsonibacter"  ,          "Parabacteroides_B_862066",
#             "Blautia_A_141780"    ,   "Parasutterella",         "Haemophilus_D_735815"   )
res_tumor$label <- ifelse(res_tumor$pvalue<0.05,
                          res_tumor$genus,NA)
res_tumor$group <- ifelse(res_tumor$pvalue<0.05&res_tumor$ratio>1,"Higher colonization rate\nin ETE",
                          ifelse(res_tumor$pvalue<0.05&res_tumor$ratio<1,"Higher colonization rate\nin ETnonE",
                                 "No significant colonization rate\ndifference"))

res_tumor$group <- factor(res_tumor$group,levels = c("Higher colonization rate\nin ETE","Higher colonization rate\nin ETnonE",
                                                     "No significant colonization rate\ndifference"))
library(ggrepel)
pdf("../../Fig3/Fig3E_prevalance_tumor.pdf",width = 5,height=5)
ggplot(data=res_tumor, aes(x=prevalance.nonE.coli, y=prevalance.E.coli,label=label)) + 
  geom_point(aes(color = group,size=ratio), alpha=0.75) + 
  xlab("Intratumoral bacteria in ETnonE") + 
  ylab("Intratumoral bacteria in ETE") + 
  scale_color_manual(values = c(col3[2],col3[1],"grey")) +
  #scale_y_continuous(limits=c(-30,60.6), expand = c(0,0))+
  #scale_x_continuous(limits=c(-30,60.6), expand = c(0,0))+
  annotate("segment", x = 0, xend = 1, y = 0, yend = 1, colour = col3[2], linewidth=0.5)+
  geom_text_repel(min.segment.length = 0, seed = 42,
                  box.padding = 0.2,fontface="bold")+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position =c(.85, .6),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))
dev.off()

write_xlsx(res_tumor,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/Fig3F_data250711.xlsx")

#========FigS3F_prevalance&BS================
i=7
res_glm_ad <- list()
for (i in seq_len(ncol(data_tumor_16S))) {
  genus <- colnames(data_tumor_16S)[i]
  # 二分类处理 若该菌出现则计数为1
  # glm with BS
  tmp.data <- data.frame(sample=rownames(data_tumor_16S),prevalance = as.numeric(ifelse(data_tumor_16S[, i] == 0, "0", "1")))
  tmp.data <- merge(tmp.data,group_tumor,by="sample")
  mod2<-glm(prevalance~balance_value+Age+Gender+BMI+Drink+Smoke,data = tmp.data,family = "binomial")
  a <- summary(mod2)$coefficients %>% as.data.frame()
  a<-  as.data.frame(a[2,])
  a$genus <- colnames(data_tumor_16S)[i]
  res_glm_ad[[i]] <-a
}
res_glm_ad <- do.call("rbind",res_glm_ad)
res_glm_ad$adj.p <- p.adjust(res_glm_ad$`Pr(>|z|)`,"BH")

#========FigS3G_abundance&BS================
i=7
res_lm_ad <- list()
for (i in seq_len(ncol(data_tumor_16S))) {
  genus <- colnames(data_tumor_16S)[i]
  # 二分类处理 若该菌出现则计数为1
  # glm with BS
  tmp.data <- data.frame(sample=rownames(data_tumor_16S), abundance = data_tumor_16S[, i])
  tmp.data <- merge(tmp.data,group_tumor,by="sample")
  mod2<-lm(abundance~balance_value+Age+Gender+BMI+Drink+Smoke,data = tmp.data)
  a <- summary(mod2)$coefficients %>% as.data.frame()
  a<-  as.data.frame(a[2,])
  a$genus <- colnames(data_tumor_16S)[i]
  res_lm_ad[[i]] <-a
}
res_lm_ad <- do.call("rbind",res_lm_ad)
res_lm_ad$adj.p <- p.adjust(res_lm_ad$`Pr(>|t|)`,"BH")

write_xlsx(list(prevalance_BS_glm=res_glm_ad,abundance_BS_lm=res_lm_ad),path="/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/final_16s_output/Fig3FG_supply.xlsx")

#===========Masslin肿瘤==================
library(Maaslin2)
#group_tumor <- na.omit(group_tumor)
a <- colnames(data_tumor)
data_tumor_16S_maaslin <- cbind(row.names(data_tumor),data_tumor)
colnames(data_tumor_16S_maaslin) <- c("ID",a)
write_tsv(data_tumor_16S_maaslin,"./Maaslin_input/data_tumor_16S.tsv")
group_tumor_maaslin <- group_tumor
colnames(group_tumor_maaslin) <- c("ID",colnames(group_tumor)[-1])
group_tumor_maaslin$batch <- gsub("1","one",group_tumor_maaslin$batch)
group_tumor_maaslin$batch <- gsub("2","two",group_tumor_maaslin$batch)
group_tumor_maaslin$batch <- gsub("3","three",group_tumor_maaslin$batch)
group_tumor_maaslin$batch <- gsub("4","four",group_tumor_maaslin$batch)
group_tumor_maaslin$batch <- gsub("5","five",group_tumor_maaslin$batch)
write_tsv(group_tumor_maaslin,"./Maaslin_input/group_tumor.tsv")
test <- read_tsv("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/Maaslin_input/data_tumor_16S.tsv")

fit_data <- Maaslin2(
  '/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/Maaslin_input/data_tumor_16S.tsv', 
  '/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/Maaslin_input/group_tumor.tsv',
  'Maaslin_batch_patient',
  fixed_effects = c('Cluster_new','batch',"Age",'Gender','BMI','Drink','Smoke'),
  #random_effects = c('sampleID'),
  min_abundance = 0,min_prevalence = 0,
  reference=c("Cluster_new,non_E.coli",
              "batch,three"),
  standardize = FALSE, cores=1)

Maaslin_batch_patient <- read_tsv("Maaslin_batch_patient/significant_results.tsv") %>% 
  mutate(statistic = coef / stderr) 



#======Fig3G_tumor_enriched=============
Maaslin_batch_patient <- read_tsv("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/Maaslin_batch_patient/all_results.tsv") %>% 
  mutate(statistic = coef / stderr) 
library(data.table)
test <- Maaslin_batch_patient[Maaslin_batch_patient$metadata=="Cluster_new",]

diff_plot <- test

#diff_plot <- diff_plot[abs(diff_plot$statistic)>1.5,]
diff_plot$change <- ifelse(diff_plot$pval<0.05&diff_plot$statistic>0,"Up",ifelse(diff_plot$pval<0.05&diff_plot$statistic<0,"Down","Stable"))
#diff_plot <- diff_plot%>%arrange(-statistic)
#diff_plot$feature <- factor(diff_plot$feature,levels=diff_plot$feature)

#target <- read.table("/groups/ProHuShiX/home/liuyuyao/BIGMeta_16S/translocate/PCA/prevalance_abundance_overlap.txt")
diff_plot$label <- ifelse(diff_plot$pval<0.05,
                          diff_plot$feature,NA)

diff_plot$change <- factor(diff_plot$change,levels = c("Up","Down","Stable"))
pdf("../../Fig3/Fig3F_tumor_enriched.pdf",width = 5,height = 5)
ggplot(diff_plot,aes(x=statistic,y=-log10(pval),label=label))+
  geom_point(aes(color=change))+
  scale_color_manual(values=c(col3[2],col3[1],"grey"),label=c("ETE enriched","ETnonE enriched","Non significance"))+
  labs(y=expression(-Log[10](italic(q)~value)))+
  theme_bw()+
  theme(panel.background =element_blank(),legend.position = c(0.15,0.15),
        axis.line=element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        #axis.text=element_text(face="bold",colour="black",size=20),
        #axis.title = element_text(family = "Cuti",face="bold",colour="black",size=20),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5))+
  geom_vline(xintercept=c(-2,2),lty=3,col="black",lwd=0.8)+
  geom_hline(yintercept=-log10(0.05),lty=3,col="black",lwd=0.8)+
  #geom_text(data=mettle,
  # aes(x= logFC,y= -log10(P.Value),label= name), check_overlap = TRUE,
  #nudge_y = 0.15,fontface="bold")+
  geom_text_repel(min.segment.length = 0, seed = 42,
                  box.padding = 0.2,fontface="bold")
#xlim(-2.5,5)
dev.off()
#save(diff_plot,file="Fig2/Fig2.RData")

write_xlsx(diff_plot,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/Fig3G_data250711.xlsx")

#========Fig3D source tracker============
#============数据预处理==========
#============MGS======================
#============MGS(不筛平均丰度)======================
data_mgs_genus=read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/旧的16S代码/胡教授给的/Greengenes2.mgs.genus.ready.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
#单独挑出最后用来分析的HCC和Health
data <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/species_raw_filter.txt")
sample <- gsub("\\.","-",colnames(data))
data_mgs_genus <- data_mgs_genus[,colnames(data_mgs_genus)%in%sample]
data_mgs_genus=data_mgs_genus[rowSums(data_mgs_genus!=0)>(ncol(data_mgs_genus)*0.1),]
#踢掉没有注释出来的属
data_mgs_genus<- data_mgs_genus %>%dplyr::mutate(genus=str_extract(row.names(data_mgs_genus),"g__.*"))
test=data_mgs_genus[data_mgs_genus$genus=="g__",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus=="g__",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus %like% "uncultured",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus%like% "Mitochondria",]
row.names(data_mgs_genus) <- data_mgs_genus$genus
data_mgs_genus$genus <- NULL
row.names(data_mgs_genus) <- gsub("g__","",row.names(data_mgs_genus))
data_mgs_genus=as.data.frame(apply(data_mgs_genus,2,function(x){
  x=x/sum(x)
  return(x)
}))

data_mgs_genus=as.data.frame(t(data_mgs_genus))
range(colMeans(data_mgs_genus))
# data_mgs_genus=data_mgs_genus[,colMeans(data_mgs_genus)>0.001]
data_mgs_genus=data_mgs_genus[!duplicated(rownames(data_mgs_genus)),]

library(readxl)
MGSinfo_240829 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/MGSinfo_240829.xlsx")
data_mgs_genus <- merge(data_mgs_genus,MGSinfo_240829[,c("sampleID","sample")],by.x="row.names",by.y="sample")  ##改名字这步顺便把健康人踢掉了
row.names(data_mgs_genus) <- data_mgs_genus$sampleID
data_mgs_genus$sampleID <- NULL
data_mgs_genus$Row.names <- NULL

#====16S data====
data_16S_genus <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/final_16s_output/data_16S_genus1.txt",header=T)

#==========source tracker============
library(gridExtra)
library(Rcpp)
library(RcppArmadillo)
library(vegan)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
require("FEAST")
#devtools::install_github("cozygene/FEAST")
library(FEAST)
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/for_sourcetracker/FEAST.R")
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/for_sourcetracker/Infer_LatentVariables.R")
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/for_sourcetracker/RcppExports.R")
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/for_sourcetracker/utils_plot.R")
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/for_sourcetracker/utils.R")

setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/")
# group_tumor <- read.table("group_tumor.test.txt",header=T)
# group_tumor[group_tumor=="Yes"] <- 1
# group_tumor[group_tumor=="No"] <- 0
# group_tumor[is.na(group_tumor)==T] <- 0.5
# group_tumor <- na.omit(group_tumor)
# data_mgs_genus <- read.table("/groups/ProHuShiX/home/liuyuyao/bigmeta_final/16S_gg250509/data/mgs.greengenes.txt",header=T)

group_tumor <- read.table("final_16s_output/group_tumor.test.txt",header = T)

data_mgs_genus <- data_mgs_genus[,colnames(data_mgs_genus)%in%colnames(data_16S_genus)]
data_mgs_genus <- count_to_composition(data_mgs_genus,samples_row=T)

tmp2=data_16S_genus[row.names(data_16S_genus)%in%group_tumor$sample,]
group_tumor1 <- group_tumor[group_tumor$sample%in%row.names(tmp2),]
data_mgs_feast <- data_mgs_genus[row.names(data_mgs_genus)%in%group_tumor1$sampleID,]
data_mgs_feast <- data_mgs_feast[match(group_tumor1$sampleID,row.names(data_mgs_feast)),]

data_mgs_feast <- as.data.frame(t(data_mgs_feast))
tmp2 <- as.data.frame(t(tmp2))

data_feast_otu <- merge(tmp2,data_mgs_feast,by="row.names",all=T) # mgs和16S合并
data_feast_otu[is.na(data_feast_otu)==T] <- 0
row.names(data_feast_otu) <- data_feast_otu$Row.names
data_feast_otu$Row.names <- NULL
data_feast_otu <- as.data.frame(t(data_feast_otu))
#data_feast_otu=rbind(tmp2,data_mgs_feast)

data_feast_otu=data_feast_otu*200000000  # 10e6-10e8
data_feast_otu <- mutate_all(data_feast_otu, function(x) (as.integer(x)))

tmp2 <- as.data.frame(t(tmp2))
data_mgs_feast <- as.data.frame(t(data_mgs_feast))
data_feast_file_tumor <- data.frame(SampleID=rownames(tmp2),Env=c("Tumor"),SourceSink=c("Sink"))
data_feast_file_mgs <- data_frame(SampleID=rownames(data_mgs_feast),Env=c("Fecal"),SourceSink=c("Source"))
data_feast_file_merge <- cbind(data_feast_file_tumor,data_feast_file_mgs)
colnames(data_feast_file_merge) <- c("SampleID","Env", "SourceSink", "SampleID.1" ,"Env.1" ,"SourceSink.1")
data_feast_file_tumor <- data_feast_file_merge[,c("SampleID","Env","SourceSink","SampleID.1")]
colnames(data_feast_file_tumor)[4] <- "id"
data_feast_file_mgs <- data_feast_file_merge[,c("SampleID.1","Env.1","SourceSink.1","SampleID.1")]
colnames(data_feast_file_mgs) <- c("SampleID","Env","SourceSink","id")

data_feast_file <- rbind(data_feast_file_tumor,data_feast_file_mgs)

data_feast_otu=as.data.frame(t(data_feast_otu))

#=====Feast data====
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/sourcetracker_data")
# write.table(data_feast_otu,file = "feast/data_feast_otu_Tumor1.txt",row.names = T,quote = F,sep = "\t")
# write.table(data_feast_file,file = "feast/data_feast_file_Tumor1.txt",row.names = F,quote = F,sep = "\t")
write.table(data_feast_otu,file = "feast/data_feast_otu_Tumor2e8.txt",row.names = T,quote = F,sep = "\t")

# feast in tumor
data_feast_file <- Load_metadata(metadata_path = "feast/data_feast_file_Tumor1.txt")
# data_feast_otu <- Load_CountMatrix(CountMatrix_path = "feast/data_feast_otu_Tumor1.txt")
data_feast_otu <- Load_CountMatrix(CountMatrix_path = "feast/data_feast_otu_Tumor2e8.txt")

#====source tracker====
set.seed(12345)
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/sourcetracker_data/")
FEAST_output <- FEAST(C =data_feast_otu, metadata = data_feast_file, different_sources_flag = 1, dir_path = "./Feast_out1",outfile="Tumor")
#====重跑sourceTracker可从此====
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/sourcetracker_data/")
FEAST_output = read.table("./Feast_out1/Tumor_source_contributions_matrix.txt",header = T,row.names = 1,sep = "\t",check.names = F)
PlotSourceContribution(SinkNames = rownames(FEAST_output)[1:4],
                       SourceNames = colnames(FEAST_output), dir_path = "./Feast_out1",
                       mixing_proportions = FEAST_output, Plot_title = "Tumor",Same_sources_flag = 0, N = 4)

setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/sourcetracker_data/")
feast_tumor=foreach(i=1:nrow(FEAST_output),.combine = rbind) %do%  {
  tmp.sample=rownames(FEAST_output)[i]
  tmp.portion=(FEAST_output[i,][!is.na(FEAST_output[i,])])
  return.string=data.frame(SampleID=tmp.sample,SourceFecal=tmp.portion[1],SourceUnknown=tmp.portion[2])
}
hist(feast_tumor$SourceFecal)

feast_tumor$SourceFecal=log2(feast_tumor$SourceFecal)
#feast_tumor$ID=gsub("tumor_","",feast_tumor$SampleID)
feast_tumor$ID=gsub("_Tumor","",feast_tumor$SampleID)
feast_tumor=merge(feast_tumor,group_tumor,by.x="ID",by.y="sample",all.x=T)
feast_tumor=feast_tumor[feast_tumor$SourceFecal!=-Inf,]
feast_tumor$balance_value <- as.numeric(feast_tumor$balance_value)
feast_tumor$Cluster <- as.factor(feast_tumor$Cluster)
feast_tumor$Cluster_new <- as.factor(feast_tumor$Cluster_new)

summary(lm(SourceFecal~Cluster_new+batch+Age+Gender+BMI+Drink+Smoke,data=feast_tumor))
feast_tumor$ad_SourceFecal <- resid(lm(SourceFecal~batch+Age+Gender+BMI+Drink+Smoke,data=feast_tumor))

ggboxplot(feast_tumor,x="Cluster_new",y="ad_SourceFecal",color =  "Cluster_new",)+
  stat_compare_means(comparisons = list(c("non_E.coli","E.coli")))+
  scale_color_jco()

cor.test(feast_tumor$ad_SourceFecal,feast_tumor$balance_value,method = "spearman")
ggplot(feast_tumor,aes(balance_value,ad_SourceFecal)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +  stat_cor(method = 'spearman', label.y=-0.505)+theme_bw()

pdf("../../../Fig3/Fig3S_ad_sorcetracker&BS.pdf",width = 5,height = 5)
ggplot(data = feast_tumor,aes(x=balance_value,y=ad_SourceFecal))+
  geom_point(aes(color=Enterotype))+geom_smooth(method = 'lm',color=col4[2],fill=col3[1])+
  scale_color_manual(values = col3)+
  stat_cor(method = 'spearman', label.y=0.6)+
  theme(panel.background =element_blank(),legend.position = c(0.15,0.15),
        axis.line=element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

pdf("../../../Fig3/Fig3S_sorcetracker&BS.pdf",width = 5,height = 5)
ggplot(data = feast_tumor,aes(x=balance_value,y=SourceFecal))+
  geom_point(aes(color=Enterotype))+geom_smooth(method = 'lm',color=col4[2],fill=col3[1])+
  scale_color_manual(values = col3)+
  stat_cor(method = 'spearman', label.y=0.6)+
  theme(panel.background =element_blank(),legend.position = c(0.15,0.15),
        axis.line=element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()


cor.test(feast_tumor$SourceFecal,feast_tumor$balance_value,method = "spearman")
ggplot(feast_tumor,aes(balance_value,SourceFecal)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +  stat_cor(method = 'spearman', label.y=-0.505)+theme_bw()

feast_tumor$Enterotype <- ifelse(feast_tumor$Cluster=="4","ETE","ETnonE")
feast_tumor$Enterotype <- factor(feast_tumor$Enterotype,levels=c("ETnonE","ETE"))
pdf("../../../Fig3/Fig3C_ad_sorce_tracker250630.pdf",width = 4,height = 5)
ggplot(feast_tumor, aes(x=Enterotype, y=ad_SourceFecal, color=Enterotype)) + 
  geom_violin(trim=FALSE,aes(fill=Enterotype,alpha=0.8))+
  geom_boxplot(width=0.15, fill="white")+theme_classic()+
  stat_compare_means(comparisons = list(c("ETE","ETnonE")))+
  theme(
    legend.position = "none")+
  scale_color_manual(values = col4)+
  scale_fill_manual(values=col3)
dev.off()

pdf("../../../Fig3/Fig3C_sorce_tracker250630.pdf",width = 4,height = 5)
ggplot(feast_tumor, aes(x=Enterotype, y=SourceFecal, color=Enterotype)) + 
  geom_violin(trim=FALSE,aes(fill=Enterotype,alpha=0.8))+
  geom_boxplot(width=0.15, fill="white")+theme_classic()+
  stat_compare_means(comparisons = list(c("ETE","ETnonE")))+
  theme(
    legend.position = "none")+
  scale_color_manual(values = col4)+
  scale_fill_manual(values=col3)
dev.off()


#=======Fig3E_Heterogeneity肿瘤Gut-Liver距离===========
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/")
#============MGS数据预处理======================
data_mgs_genus=read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/旧的16S代码/胡教授给的/Greengenes2.mgs.genus.ready.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
#单独挑出最后用来分析的HCC和Health
data <- read.table("../../Final_input/species_raw_filter.txt")
sample <- gsub("\\.","-",colnames(data))
data_mgs_genus <- data_mgs_genus[,colnames(data_mgs_genus)%in%sample]
data_mgs_genus=data_mgs_genus[rowSums(data_mgs_genus!=0)>(ncol(data_mgs_genus)*0.1),]
#踢掉没有注释出来的属
data_mgs_genus<- data_mgs_genus %>%dplyr::mutate(genus=str_extract(row.names(data_mgs_genus),"g__.*"))
test=data_mgs_genus[data_mgs_genus$genus=="g__",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus=="g__",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus %like% "uncultured",]
data_mgs_genus=data_mgs_genus[!data_mgs_genus$genus%like% "Mitochondria",]
row.names(data_mgs_genus) <- data_mgs_genus$genus
data_mgs_genus$genus <- NULL
row.names(data_mgs_genus) <- gsub("g__","",row.names(data_mgs_genus))
data_mgs_genus=as.data.frame(apply(data_mgs_genus,2,function(x){
  x=x/sum(x)
  return(x)
}))

data_mgs_genus=as.data.frame(t(data_mgs_genus))
#data_mgs_genus=data_mgs_genus[,colMeans(data_mgs_genus)>0.001]
rownames(data_mgs_genus)=gsub("\\.","-",rownames(data_mgs_genus))
rownames(data_mgs_genus)=gsub("_count","",rownames(data_mgs_genus))
data_mgs_genus=data_mgs_genus[!duplicated(rownames(data_mgs_genus)),]

library(readxl)
MGSinfo_240829 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/MGSinfo_240829.xlsx")
data_mgs_genus <- merge(data_mgs_genus,MGSinfo_240829[,c("sampleID","sample")],by.x="row.names",by.y="sample")  ##改名字这步顺便把健康人踢掉了
row.names(data_mgs_genus) <- data_mgs_genus$sampleID
data_mgs_genus$sampleID <- NULL
data_mgs_genus$Row.names <- NULL

#=========计算Gut-Liver dissimilarity==========
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/")
# group_tumor <- read.table("group_tumor.ready.txt",header=T)
# group_tumor <- na.omit(group_tumor)
# data_mgs_genus <- read.table("/groups/ProHuShiX/home/liuyuyao/bigmeta_final/16S_gg250509/data/trim.mgs.greengenes.txt",header=T)
# 目前的版本是mgs只剔除<10%检出率，16S不取MGS也能检测到的属

data_16S_genus <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/final_16s_output/data_16S_genus1.txt",header=T)

# data_mgs_genus <- data_mgs_genus[,colnames(data_mgs_genus)%in%colnames(data_16S_genus)]
# data_mgs_genus <- count_to_composition(data_mgs_genus,samples_row=T)
data_tumor <- data_16S_genus[row.names(data_16S_genus)%in%group_tumor$sample,]

data_mgs_genus <- data_mgs_genus[row.names(data_mgs_genus)%in%group_tumor$sampleID,]
data_tumor <- as.data.frame(t(data_tumor))
data_mgs_genus <- as.data.frame(t(data_mgs_genus))

data_RDA <- merge(data_tumor,data_mgs_genus,by="row.names",all=T)
data_RDA[is.na(data_RDA)==T] <- 0
row.names(data_RDA) <- data_RDA$Row.names
data_RDA$Row.names <- NULL
data_RDA <- as.data.frame(t(data_RDA))
# data_RDA <- rbind(data_tumor,data_mgs_genus[row.names(data_mgs_genus)%in%group_tumor$sampleID,])

data_mgs_genus <- as.data.frame(t(data_mgs_genus))
RDA_group_mgs <- as.data.frame(row.names(data_mgs_genus[row.names(data_mgs_genus)%in%group_tumor$sampleID,]))
colnames(RDA_group_mgs) <- "sampleID"
RDA_group_mgs$group <- "gut"
RDA_group_mgs <- merge(RDA_group_mgs,group_tumor,by="sampleID")
RDA_group_mgs <- RDA_group_mgs[!duplicated(RDA_group_mgs$sampleID),]
RDA_group_mgs$sample <- RDA_group_mgs$sampleID

data_tumor <- as.data.frame(t(data_tumor))
RDA_group_16S <- as.data.frame(row.names(data_tumor))
colnames(RDA_group_16S) <- "sample"
RDA_group_16S$group <- "tumor"
RDA_group_16S <- merge(RDA_group_16S,group_tumor,by="sample")
RDA_group <- rbind(RDA_group_16S,RDA_group_mgs)

data_RDA_clr <- transform_and_filter_taxa(data_RDA,samples_row = T,method = "clr",missing_filter = 0)
dist_EU <- vegdist(data_RDA_clr,method = "euclidean")
distmat <- as.matrix(dist_EU)

# 提取Gut和Liver(mgs和tumor两两之间计算的矩阵部分)
i=1
dist_inter=foreach(i=1:length(RDA_group_16S$sample),.combine = rbind) %do% {
  cat(i,"\n")
  tmp_sample=RDA_group_16S$sample[i]
  # 提取一个sample[i]的2×2矩阵
  tmp_dist=distmat[rownames(distmat) %in%c(tmp_sample,gsub("_.*$","",tmp_sample)),colnames(distmat)%in%c(tmp_sample,gsub("_.*$","",tmp_sample))]
  
  return.string=(data.frame(as.table(tmp_dist))[lower.tri(tmp_dist, diag = F), ])
}

dist_inter <- merge(dist_inter,RDA_group_16S,by.x="Var2",by.y="sample")

dist_inter$batch <- as.factor(dist_inter$batch)
dist_inter$sampleID <- as.factor(dist_inter$sampleID)
dist_inter$Cluster_new <- factor(dist_inter$Cluster_new,levels=c("non_E.coli","E.coli"))

summary(lm(Freq~Cluster_new+batch+Age+Gender+BMI+Drink+Smoke,data=dist_inter))
dist_inter$ad_Freq <- resid(lm(Freq~batch+Age+Gender+BMI+Drink+Smoke,data=dist_inter))

# pdf("/groups/ProHuShiX/home/liuyuyao/Figure_250105/16s_diversity/ET_heterogenity.pdf",width=5,height=5)
# ggboxplot(dist_inter, x = "Cluster_new", y = "ad_Freq",
#           color = "Cluster_new", palette = "jco")+
#   stat_compare_means(comparisons = list(c("E.coli","non_E.coli")))
# dev.off()

col3 <- c("#fcd364","#b35c31")  ##ETnonE,ETE
col4 <- c("#BD9E4B","#b35c31")  ##画小提琴图加粗的线的颜色 
pdf("../../Fig3/ET_heterogenity.pdf",width=4,height=5)
ggplot(dist_inter, aes(x=Cluster_new, y=ad_Freq, color=Cluster_new)) + 
  geom_violin(trim=FALSE,aes(fill=Cluster_new,alpha=0.8))+
  geom_boxplot(width=0.15, fill="white")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("E.coli","non_E.coli")))+
  theme(legend.position = "none")+
  scale_color_manual(values = col4)+
  scale_fill_manual(values=col3)
dev.off()

cor.test(dist_inter$ad_Freq,dist_inter$balance_value,method="spearman")
ggplot(dist_inter,aes(balance_value, ad_Freq)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +  stat_cor(method = 'spearman', label.y=-0.505)+theme_bw()

pdf("../../Fig3/Fig3S_ad_GutLiver_dissimilarity&BS.pdf",width = 5,height = 5)
ggplot(data = dist_inter,aes(x=balance_value,y=ad_Freq))+
  geom_point(aes(color=Cluster_new))+geom_smooth(method = 'lm',color=col4[2],fill=col3[1])+
  scale_color_manual(values = col3)+
  stat_cor(method = 'spearman', label.y=0.6)+
  theme(panel.background =element_blank(),legend.position = c(0.15,0.15),
        axis.line=element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
dev.off()

cor.test(dist_inter$Freq,dist_inter$balance_value,method="spearman")
ggplot(dist_inter,aes(balance_value, Freq)) +
  geom_point() +
  geom_smooth(method='lm', formula= y~x) +  stat_cor(method = 'spearman', label.y=220)+theme_bw()

#========比较患者内和患者间====================
data_16S <- read.table("./final_16s_output/data_16S_genus1.txt",header=T)
amplicon_data_clr <- read.table("./final_16s_output/data_clr1.txt",header = T)
group_16s <- read.table("./final_16s_output/group_16s1.txt",header=T)

diversity<-vegan::diversity
shannon <- data.frame(Richness=specnumber(data_16S),
                      shannon=diversity(data_16S, index="shannon"),
                      simpson=diversity(data_16S, index="simpson"),
                      inverse_simpson=diversity(data_16S, index="invsimpson"),
                      Evenness=diversity(data_16S, index="shannon")/log(specnumber(data_16S)))
group_all <- merge(group_16s,shannon,by.x="sample",by.y="row.names")
group_all$batch <- as.factor(batch$batch[match(group_all$sample,batch$sample)])
group_all <- na.omit(group_all)
#group_all$Enterotype <- ifelse(group_all$Cluster=="4","ETE","ETnonE")
# summary(lm(shannon~Enterotype+batch+Age+Gender+BMI+Smoke+Drink,data=group_all))
# group_tumor$ad_shannon <- resid(lm(shannon~batch+Age+Gender+BMI+Smoke+Drink,data=group_tumor))

amplicon_single=data.frame(table(group_all$sampleID))
amplicon_single=as.character(amplicon_single$Var1[amplicon_single$Freq==1])
amplicon_phenotype=group_all[!group_all$sampleID %in% amplicon_single,]
amplicon_sample=unique(amplicon_phenotype$sampleID)
amplicon_all=unique(group_all$sampleID)

amplicon_diversity=vegdist((data_clr),method = "euclidean")
distmat <- as.matrix(amplicon_diversity)

i=226
dist_intra=foreach(i=1:length(amplicon_sample),.combine = rbind) %do% {
  cat(i,"\n")
  tmp_sample=amplicon_sample[i]
  tmp_dist=distmat[rownames(distmat) %like% tmp_sample,colnames(distmat) %like% tmp_sample]
  
  return.string=(data.frame(as.table(tmp_dist))[lower.tri(tmp_dist, diag = F), ])
}

i=208
dist_inter=foreach(i=1:length(amplicon_sample),.combine = rbind) %do% {
  cat(i,"\n")
  tmp_sample=amplicon_sample[i]
  tmp_dist=distmat[rownames(distmat) %like% tmp_sample,!colnames(distmat) %like% tmp_sample,drop=F]
  
  return.string=(setNames(melt(tmp_dist), c('Var1', 'Var2', 'Freq')))
}
dist_intra$Group="Intra"
dist_inter$Group="Inter"
dist_merge=rbind(dist_intra,dist_inter)
dist_merge=dist_merge[dist_merge$Var1 %like% "T" & dist_merge$Var2 %like% "T",]
table(dist_merge$Group)

tmp_comparisons <- list( c("Intra","Inter"))
pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig3/Heterogeneity250711.pdf",width = 4,height = 5)
ggboxplot(dist_merge, x = "Group", y = "Freq",
          color = "Group", palette = "jco")+
  stat_compare_means(comparisons = tmp_comparisons,method = "wilcox.test")+ylab("Heterogeneity")
dev.off()
save(dist_merge,file="/groups/ProHuShiX/home/liuyuyao/Figure_250105/16s_diversity/Heterogeneity.RData")

p_values1 <- compare_means(Freq ~ Group, data = dist_merge, 
                           group.by = NULL, # 若没有分组变量
                           comparisons = tmp_comparisons,
                           method = "wilcox.test") 

#====Fig3H&I: 16S_network.R====



