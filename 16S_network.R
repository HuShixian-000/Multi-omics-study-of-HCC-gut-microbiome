rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/")
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/functions.R")
library(vegan)
library(ggplot2)
library(foreach)
library(data.table)
library(cowplot)
library(ggsci)
library(dplyr)
library(tidyverse)
library(stringr)
library(ape)
library(multcomp)
library(extrafont)
library(patchwork)
library(ggtern)
library(reshape2)
library(ggpubr)
library(scales)
library(viridis)
library(rstatix)
library(readxl)
set.seed(12345)


#========16S数据===================
#=========每个患者随机抽取1个肿瘤样本并整理样本信息表=============
data_16S <- read.table("../final_16s_output/data_16S_genus1.txt",header=T)
# data_clr <- read.table("data/data_clr.txt",header = T)
group_16s <- read.table("../final_16s_output/group_16s1.txt",header=T)
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
group_tumor <- group_tumor %>%
  group_by(sampleID) %>%      # 按患者分组
  sample_n(1) %>%              # 随机抽取1个样本
  ungroup()                    # 取消分组

# group_tumor <- group_tumor[!duplicated(group_tumor$sampleID),]
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

#========E.coli:导出tsv计算fastspar网络===========
data_tumor_16S <- data_16S[row.names(data_16S)%in%group_tumor$sample,]
OTU <- data_tumor_16S[row.names(data_tumor_16S)%in%group_tumor[group_tumor$Cluster_new=="E.coli","sample"],]
OTU <- as.data.frame(t(OTU))
# OTU <- OTU*1000000
OTU <- OTU*2e8
mean(OTU == 0)
OTU <- mutate_all(OTU, function(x) (as.integer(x)))
OTU <- cbind(row.names(OTU),OTU)
colnames(OTU)[1] <- "#OTU_ID"
write_tsv(OTU,"../250711网络分析/fastp/E.coli_tumor/OTU.tsv")

#====Linux建网络====

corrlation <- read_tsv("./fastp/E.coli_tumor/median_correlation.tsv")
pvalue <- read_tsv("./fastp/E.coli_tumor/pvalues.tsv")

colnames(corrlation)[1] <- "OTU_1"
correlation <- corrlation %>%
  pivot_longer(
    cols = -OTU_1,                     # 排除id列
    names_to = "OTU_2",      # 原列名存入新列
    values_to = "correlation"        # 原数值存入新列
  )%>%
  filter(as.character(OTU_1) < as.character(OTU_2)) %>% # 只保留唯一组合
  arrange(OTU_1, OTU_2) %>%
  mutate(com =paste0(OTU_1,"_vs_",OTU_2))

colnames(pvalue)[1] <- "OTU_1"
pvalue<- pvalue %>%
  pivot_longer(
    cols = -OTU_1,                     # 排除id列
    names_to = "OTU_2",      # 原列名存入新列
    values_to = "pvalue"        # 原数值存入新列
  )%>%
  filter(as.character(OTU_1) < as.character(OTU_2)) %>% # 只保留唯一组合
  arrange(OTU_1, OTU_2)  %>%
  mutate(com =paste0(OTU_1,"_vs_",OTU_2))


E.coli_network <- merge(correlation,pvalue[,c("com","pvalue")],by="com")%>%
  filter(pvalue<0.05&abs(correlation)>0.4)

library(writexl)
write_xlsx(E.coli_network,"./fastp/E.coli_tumor/E.coli_network.xlsx")


#======E.coli计算fastspar网络属性===============
library(ggClusterNet)
library(phyloseq)
library(igraph)
corrlation <- read_tsv("./fastp/E.coli_tumor/median_correlation.tsv")
corrlation <- as.data.frame(corrlation)
row.names(corrlation) <- corrlation$`#OTU ID`
corrlation$`#OTU ID` <- NULL
##相关系数符合要求的才显示
corrlation[abs(corrlation)<0.4] <- 0
##排除自相关
for (i in 1:ncol(corrlation)){corrlation[i,i] <- 0}

##排除p>0.05
pvalue <- read_tsv("./fastp/E.coli_tumor/pvalues.tsv")
pvalue <- as.data.frame(pvalue)
row.names(pvalue) <- pvalue$`#OTU ID`
pvalue$`#OTU ID` <- NULL
# Step1:创建显著性掩码矩阵 (TRUE表示p<0.05)
sig_mask <- pvalue <0.05

# Step2:将不显著的相关系数设为NA (保留行列结构)
corrlation[!sig_mask] <- 0

corrlation <- as.matrix(corrlation)
result4 = nodeEdge(cor = corrlation)
#提取变文件
edge = result4[[1]]
#--提取节点文件
node = result4[[2]]

igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
net_pro <- net_properties.4(igraph,n.hub = T)%>%as.data.frame()
net_pro <-net_pro%>%mutate(propertie=row.names(net_pro))

nodepro = node_properties(igraph)%>%as.data.frame()
nodepro <- nodepro %>%mutate(node=row.names(nodepro))
head(nodepro)
library(writexl)
write_xlsx(nodepro,"./fastp/E.coli_tumor/E.coli_tumor_node_propertis.xlsx")
write_xlsx(net_pro,"./fastp/E.coli_tumor/E.coli_tumor_net_propertis.xlsx")

#======E.coli: 可视化fastspar网络===========
#--构建phyloseq对象
OTU <- data_tumor_16S[row.names(data_tumor_16S)%in%group_tumor[group_tumor$Cluster_new=="E.coli","sample"],]
OTU <- as.data.frame(t(OTU))
OTU_abundance <- as.data.frame(rowMeans(OTU))
OTU_abundance$OTU <- row.names(OTU_abundance)
colnames(OTU_abundance)[1] <- "mean_abundance"
#OTU <- OTU*100000
tax <- data.frame(Kingdom=paste0("K",1:nrow(OTU)),
                  Phylum=paste0("P",1:nrow(OTU)),Class=paste0("C",1:nrow(OTU)),Order=paste0("O",1:nrow(OTU)),
                  Family=paste0("F",1:nrow(OTU)),Genus=paste0("G",1:nrow(OTU)),Species=paste0("S",1:nrow(OTU)))
tax$Genus <- row.names(OTU)
row.names(tax) <- row.names(OTU)
OTU <- as.matrix(OTU)
TAX <- as.matrix(tax)
OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
ps = phyloseq(OTU, TAX)

#--提取相关矩阵
cor = corrlation
dim(cor)

result2 <-ggClusterNet::model_igraph(cor =  cor,
                                     method = "cluster_fast_greedy",
                                     seed = 123
)
node = result2[[1]]
head(node)

dat = result2[[2]]
head(dat)
tem = data.frame(mod = dat$model,col = dat$color) %>%  
  dplyr::distinct( mod, .keep_all = TRUE)  
col = tem$col
names(col) = tem$mod

#---node节点注释#-----------
otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)

tem2 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
head(tem2)

tem3 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                        manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
)
head(edge2)
col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% 
  select(color,manual)
col0 = col_edge$manual
names(col0) = col_edge$color

#如果不rm，后面开R会自动加载phyloseq包，就不能跑混合线性模型了
#rm(TAX,OTU,ps)

library(ggnewscale)
#a <- intersect(res_tumor[res_tumor$p.adj<0.05,]$genus,test[test$qval<0.05&test$statistic>0,]$feature)
#write.table(as.data.frame(a),"prevalance_abundance_overlap.txt",row.names = F,col.names = F,quote = F)

edge2$cor <- factor(edge2$cor,levels=c("+","-"))
p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color=cor),data = edge2, size = 0.5) +
  scale_color_manual(values=c("#FBEBB7","#BCDEE7"))
#scale_colour_manual(values = col0) 
library(ggpubr)
library(ggrepel)
# ggsave("./cs1.pdf",p1,width = 16,height = 14)

node_list <- lapply(unique(edge$OTU_1), function(x){ y <- edge[edge$OTU_1==x|edge$OTU_2==x,]
y$OTU <- x
return(y)}
)
names(node_list) <- unique(edge$OTU_1)
node_list <- do.call("rbind",node_list)

hub <- as.data.frame(table(node_list$OTU))
positive.hub <- as.data.frame(table(node_list[node_list$cor=="+",]$OTU))%>%arrange(-Freq)
negative.hub <- as.data.frame(table(node_list[node_list$cor=="-",]$OTU))%>%arrange(-Freq)

data_mgs_genus <- read.table("../final_16s_output/mgs.greengenes.txt",header=T)
common <- intersect(colnames(data_mgs_genus),colnames(data_16S))
dat$label <- ifelse(dat$OTU%in%common,dat$OTU,NA)
dat <- merge(dat,OTU_abundance,by="OTU")
dat$hub <- ifelse(dat$OTU%in%positive.hub[1:5,]$Var1,dat$OTU,ifelse(dat$OTU%in%negative.hub[1:5,]$Var1,dat$OTU,NA))
dat$direction <- ifelse(dat$OTU%in%positive.hub[1:5,]$Var1,"positive",ifelse(dat$OTU%in%negative.hub[1:5,]$Var1,"negative",NA))

table(dat$model)##判断展示的model

col <- c(col,"#C1C1C1")
names(col)[20] <- "mini_model"
# dat$model_new <- ifelse(dat$model%in%c("model_1","model_4"),"model_1",
#                         ifelse(dat$model%in%c("model_2","model_5","model_9"),"model_2",
#                                ifelse(dat$model=="model_3","model_3","mini_model")))

dat$model_new <- ifelse(dat$model%in%c("model_1","model_2","model_3","model_4","model_5","model_6",
                                       "model_7","model_8","model_9","model_10"),"model_1","mini_model")

col2 <- c("#DC4636","#3885B6","#FDA94D","#C1C1C1")
names(col2) <- c("model_1","model_2","model_3","mini_model")


#挑丰度前5的粪便也能检测到的菌来标记
dat2 <- dat[is.na(dat$label)==F,]%>%arrange(-mean_abundance)
dat2 <- dat2[1:5,]
dat$group <- ifelse(is.na(dat$label)==T,"gut absent","gut present")
dat$group <- factor(dat$group,levels = c("gut absent","gut present"))
E.coli_dat <- dat[is.na(dat$label)==F,]%>%mutate(Enterotype="ETE")


dat$group_new <- ifelse(dat$group=="gut present"&dat$model_new=="model_1","show_1",
                        ifelse(dat$group=="gut present"&dat$model_new=="mini_model","show_2",
                               ifelse(dat$group=="gut absent"&dat$model_new=="model_1","hide_1","hide_2")))
dat$group_new <- factor(dat$group_new,levels = c("show_1","show_2","hide_1","hide_2"))

p2 <- p1 +
  new_scale_color() +
  geom_point(aes(X1, X2,color =model_new,size=mean_abundance*100,fill=group_new), data = dat,shape=21) +
  scale_colour_manual(values = col2) +
  scale_fill_manual(values=c("#DC4636","#C1C1C1","white","white"))+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

pdf("./plot/E.coli_tumor_network.pdf",width=7,height=5)
p2+
  geom_text_repel(aes(X1, X2,label=label), data = dat2,min.segment.length = Inf, seed = 42,
                  box.padding = 0.2,fontface="bold")
dev.off()


write_xlsx(list(ETE_network_info1=dat,ETE_network_info2=E.coli_network),path="/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/supply_tables/ETE_network_info.xlsx")

#--随即取出任意比例节点-网络鲁棒性#---------
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/Network_stability_pre_encapsulation.R")
library(ggClusterNet)
library(phyloseq)

cor <- corrlation
#存在某些情况计算不出来相关系数，定义相关为0
cor[is.na(cor)]<-0
#-去除自相关的点
diag(cor)<-0  
#-查看网络边的数量
sum(abs(cor)>0)/2
#网络中节点的数量
sum(colSums(abs(cor))>0)  # node number: number of species with at least one linkage with others.
#去除没有任何相关的节点.
network.raw<-cor[colSums(abs(cor))>0,colSums(abs(cor))>0]
#对应的删除otu表格otu
#sp.ra2<-sp.ra[colSums(abs(cor))>0]
sp.ra2 <- OTU_abundance[OTU_abundance$OTU%in%rownames(network.raw),]
sp.ra2 <- sp.ra2[match(sp.ra2$OTU,rownames(network.raw)),]  #check if matched

## 鲁棒性评估robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained


#输入相关矩阵 OTU表格
Weighted.simu<-rmsimu(netRaw = network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance, 
                      abundance.weighted=T,nperm=100)
head(Weighted.simu)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance,
                        abundance.weighted=F,nperm=100)
head(Weighted.simu)
# tem = ps %>% sample_data() %>% .$Group %>% unique() %>% as.character()

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),
                 rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Group="E.coli")

head(dat1)

write.table(dat1,file="./fastp/E.coli_tumor/E.coli_tumor_robustness.txt",row.names = F,col.names = T,sep="\t",quote = F)

#==========保存每次抽样的数据画柱状图===============、
rmsimu_raw<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remains
  }))
}

Weighted.simu.raw<-rmsimu_raw(netRaw = network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance, 
                              abundance.weighted=T,nperm=100)

library(reshape2)
library(tidyverse)
Weighted.simu.raw<-as.data.frame(Weighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Weighted.simu.raw<-as.data.frame(Weighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Weighted.simu.raw<-reshape2::melt(Weighted.simu.raw,id.vars = "Proportion.removed")
Weighted.simu.raw <- Weighted.simu.raw%>%mutate(weighted="weighted",Group="E.coli")

Unweighted.simu.raw<-rmsimu_raw(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance,
                                abundance.weighted=F,nperm=100)

library(reshape2)
Unweighted.simu.raw<-as.data.frame(Unweighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Unweighted.simu.raw<-as.data.frame(Unweighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Unweighted.simu.raw<-reshape2::melt(Unweighted.simu.raw,id.vars = "Proportion.removed")
Unweighted.simu.raw <- Unweighted.simu.raw%>%mutate(weighted="unweighted",Group="E.coli")

ran_robust <- rbind(Weighted.simu.raw,Unweighted.simu.raw)
write.table(ran_robust,file="./fastp/E.coli_tumor/E.coli_tumor_ran_robustness_raw.txt",sep="\t",quote=F,col.names = T,row.names = F)

#====--去除关键节点-网络鲁棒性#------
library(ggClusterNet)
library(phyloseq)
library(igraph)

cor <- corrlation

#存在某些情况计算不出来相关系数，定义相关为0
cor[is.na(cor)]<-0
#-去除自相关的点
diag(cor)<-0  
#-查看网络边的数量
sum(abs(cor)>0)/2
#网络中节点的数量
sum(colSums(abs(cor))>0)  # node number: number of species with at least one linkage with others.
#去除没有任何相关的节点.
network.raw<-cor[colSums(abs(cor))>0,colSums(abs(cor))>0]

##read otu table

otutab<- ps %>% 
  scale_micro() %>%
  subset_taxa(
    row.names(tax_table(ps)) %in% row.names(network.raw)
  ) %>%
  vegan_otu() %>% 
  t() %>%
  as.data.frame()

#对应的删除otu表格otu
#sp.ra2<-sp.ra[colSums(abs(cor))>0]
sp.ra2 <- OTU_abundance[OTU_abundance$OTU%in%rownames(network.raw),]
sp.ra2 <- sp.ra2[match(sp.ra2$OTU,rownames(network.raw)),]  #check if matched
sp.ra <- sp.ra2$mean_abundance
names(sp.ra) <- sp.ra2$OTU
## robustness simulation 
#input network matrix, number of removed keystone species, keystonespecies list,  and ra of all species
#return the proportion of species remained

#get the keystone species list
igraph = make_igraph(cor)

degree = TRUE
zipi = FALSE
if (degree) {
  ret3 = node_properties(igraph) %>% 
    as.data.frame() %>%
    filter(!is.na(igraph.degree) ) %>% 
    arrange(desc(igraph.degree)) 
  head(ret3)
  tem = round(length(ret3$igraph.degree) * 0.10,0)
  module.hub = row.names(ret3)[1:tem]
}



if (zipi ) {
  res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
  p <- res[[1]]
  model = res[[2]] %>% filter(roles == "Module hubs")
  head(model)
  model$roles %>% unique()
  module.hub <- as.character(row.names(model))  
}

Weighted.simu<-rmsimu2(netRaw=network.raw,
                       rm.p.list=1:length(module.hub),
                       keystonelist=module.hub,
                       sp.ra=sp.ra,
                       abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu2(netRaw=network.raw, rm.p.list=1:length(module.hub),
                         keystonelist=module.hub,
                         sp.ra=sp.ra, abundance.weighted=F,nperm=100)

dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),
                 rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),
                              each=length(module.hub)),
                 Group= "E.coli")

head(dat1)
currentdat = dat1

write.table(dat1,file="./fastp/E.coli_tumor/E.coli_tumor_hub_robustness.txt",row.names = F,col.names = T,sep="\t",quote = F)

#======保存每次hub抽样的数据画柱状图============
rmsimu2_raw<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remains
  }))
}


Weighted.simu.raw<-rmsimu2_raw(netRaw=network.raw,
                               rm.p.list=1:length(module.hub),
                               keystonelist=module.hub,
                               sp.ra=sp.ra,
                               abundance.weighted=T,nperm=100)

library(reshape2)
library(tidyverse)
Weighted.simu.raw<-as.data.frame(Weighted.simu.raw)
Weighted.simu.raw<-Weighted.simu.raw%>%mutate(Number.hub.removed=row.names(Weighted.simu.raw))
Weighted.simu.raw<-melt(Weighted.simu.raw,id.vars = "Number.hub.removed")
Weighted.simu.raw <- Weighted.simu.raw%>%mutate(weighted="weighted",Group="E.coli")

Unweighted.simu.raw<-rmsimu2_raw(netRaw=network.raw,
                                 rm.p.list=1:length(module.hub),
                                 keystonelist=module.hub,
                                 sp.ra=sp.ra,
                                 abundance.weighted=F,nperm=100)

library(reshape2)
Unweighted.simu.raw<-as.data.frame(Unweighted.simu.raw)
Unweighted.simu.raw<-Unweighted.simu.raw%>%dplyr::mutate(Number.hub.removed=row.names(Unweighted.simu.raw))
Unweighted.simu.raw<-melt(Unweighted.simu.raw,id.vars = "Number.hub.removed")
Unweighted.simu.raw <- Unweighted.simu.raw%>%mutate(weighted="unweighted",Group="E.coli")

hub_robust <- rbind(Weighted.simu.raw,Unweighted.simu.raw)
write.table(hub_robust,file="./fastp/E.coli_tumor/E.coli_tumor_hub_robustness_raw.txt",sep="\t",quote=F,col.names = T,row.names = F)

#====加载易损性函数================
vulnerability = function(cor = cor){
  
  cor[abs(cor)>0]<-1 # adjacency matrix
  g = graph_from_adjacency_matrix(as.matrix(cor), 
                                  mode="undirected", 
                                  weighted = NULL, diag = FALSE,
                                  add.colnames = NULL) # note: this graph contains isolated nodes.
  # remove isolated nodes
  iso_node_id = which(igraph::degree(g)==0)
  g2 = igraph::delete_vertices(g, iso_node_id) # graph without isolated nodes
  
  #check node number and links
  length(V(g2));length(E(g2))   
  
  # calculate vulnerability of each node
  node.vul<-info.centrality.vertex(g2)
  return(max(node.vul))
}


#======计算网络易损性================      
i = 1
A = c()

cor = corrlation

cor[abs(cor)>0]<-1 # adjacency matrix
vulnerability(cor = cor)
##E.coli网络易损性为  0.0206265

#========non_E.coli:导出tsv计算fastspar网络===========
OTU <- data_tumor_16S[row.names(data_tumor_16S)%in%group_tumor[group_tumor$Cluster_new=="non_E.coli","sample"],]
OTU <- as.data.frame(t(OTU))
# OTU <- OTU*1000000
OTU <- OTU*2e8
mean(OTU == 0)
OTU <- mutate_all(OTU, function(x) (as.integer(x)))
OTU <- cbind(row.names(OTU),OTU)
colnames(OTU)[1] <- "#OTU_ID"
write_tsv(OTU,"./fastp/ETnonE_tumor/OTU.tsv")

#====Linux建网络====

corrlation <- read_tsv("./fastp/ETnonE_tumor/median_correlation.tsv")
pvalue <- read_tsv("./fastp/ETnonE_tumor/pvalues.tsv")

colnames(corrlation)[1] <- "OTU_1"
correlation <- corrlation %>%
  pivot_longer(
    cols = -OTU_1,                     # 排除id列
    names_to = "OTU_2",      # 原列名存入新列
    values_to = "correlation"        # 原数值存入新列
  )%>%
  filter(as.character(OTU_1) < as.character(OTU_2)) %>% # 只保留唯一组合
  arrange(OTU_1, OTU_2) %>%
  mutate(com =paste0(OTU_1,"_vs_",OTU_2))

colnames(pvalue)[1] <- "OTU_1"
pvalue<- pvalue %>%
  pivot_longer(
    cols = -OTU_1,                     # 排除id列
    names_to = "OTU_2",      # 原列名存入新列
    values_to = "pvalue"        # 原数值存入新列
  )%>%
  filter(as.character(OTU_1) < as.character(OTU_2)) %>% # 只保留唯一组合
  arrange(OTU_1, OTU_2)  %>%
  mutate(com =paste0(OTU_1,"_vs_",OTU_2))


nonE.coli_network <- merge(correlation,pvalue[,c("com","pvalue")],by="com")%>%
  filter(pvalue<0.05&abs(correlation)>0.4)
write_xlsx(nonE.coli_network,"./fastp/ETnonE_tumor/nonE.coli_network.xlsx")

#======non_E.coli:计算fastspar网络属性===============
corrlation <- read_tsv("./fastp/ETnonE_tumor/median_correlation.tsv")
corrlation <- as.data.frame(corrlation)
row.names(corrlation) <- corrlation$`#OTU ID`
corrlation$`#OTU ID` <- NULL
##相关系数符合要求的才显示
corrlation[abs(corrlation)<0.4] <- 0
##排除自相关
for (i in 1:ncol(corrlation)){corrlation[i,i] <- 0}

##排除p>0.05
pvalue <- read_tsv("./fastp/ETnonE_tumor/pvalues.tsv")
pvalue <- as.data.frame(pvalue)
row.names(pvalue) <- pvalue$`#OTU ID`
pvalue$`#OTU ID` <- NULL
# Step1:创建显著性掩码矩阵 (TRUE表示p<0.05)
sig_mask <- pvalue <0.05

# Step2:将不显著的相关系数设为NA (保留行列结构)
corrlation[!sig_mask] <- 0

corrlation <- as.matrix(corrlation)
result4 = nodeEdge(cor = corrlation)
#提取变文件
edge = result4[[1]]
#--提取节点文件
node = result4[[2]]

igraph  = igraph::graph_from_data_frame(edge, directed = FALSE, vertices = node)
net_pro <- net_properties.4(igraph,n.hub = T)%>%as.data.frame()
net_pro <-net_pro%>%mutate(propertie=row.names(net_pro))

nodepro = node_properties(igraph)%>%as.data.frame()
nodepro <- nodepro %>%mutate(node=row.names(nodepro))
head(nodepro)
library(writexl)
write_xlsx(nodepro,"./fastp/ETnonE_tumor/ETnonE_tumor_node_propertis.xlsx")
write_xlsx(net_pro,"./fastp/ETnonE_tumor/ETnonE_tumor_net_propertis.xlsx")

#=======non_E.coli:节点属性比较==============
library(readxl)
node_E_coli <- read_excel("./fastp/E.coli_tumor/E.coli_tumor_node_propertis.xlsx")%>%mutate(enterotype="E.coli")
node_nonE_coli <- read_excel("./fastp/ETnonE_tumor/ETnonE_tumor_node_propertis.xlsx")%>%mutate(enterotype="nonE.coli")

node_E_coli <- node_E_coli[node_E_coli$node%in%common,]
node_nonE_coli <- node_nonE_coli[node_nonE_coli$node%in%common,]

##补充在肠和肝都能测到的菌，但是没有成网络的部分
a <- common[!common%in%node_E_coli$node]
sup_E.coli <- data.frame(igraph.degree=rep(0,length(a)),
                         igraph.closeness=rep(0,length(a)),
                         igraph.betweenness=rep(0,length(a)),
                         igraph.cen.degree=rep(0,length(a)),
                         node=a,
                         enterotype="E.coli")
node_E_coli <- rbind(node_E_coli,sup_E.coli)


a <- common[!common%in%node_nonE_coli$node]
sup_nonE.coli <- data.frame(igraph.degree=rep(0,length(a)),
                            igraph.closeness=rep(0,length(a)),
                            igraph.betweenness=rep(0,length(a)),
                            igraph.cen.degree=rep(0,length(a)),
                            node=a,
                            enterotype="nonE.coli")
node_nonE_coli <- rbind(node_nonE_coli,sup_nonE.coli)

node_property<- rbind(node_E_coli,node_nonE_coli)

for (i in 1:4){
  p <- ggboxplot(node_property,x="enterotype",y=paste0(colnames(node_property)[i]),color="enterotype",palette = "jco")+
    stat_compare_means(comparisons = list(c("E.coli","nonE.coli")))+labs(title = gsub("igraph\\.","",colnames(node_property)[i]))
  print(p)
  ggsave(paste0("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/plot/网络属性图/",colnames(node_property)[i],".pdf"),width = 5,height = 8)
}


#======non_E.coli:可视化fastspar网络==========
#--构建phyloseq对象
OTU <- data_tumor_16S[row.names(data_tumor_16S)%in%group_tumor[group_tumor$Cluster_new=="non_E.coli","sample"],]
OTU <- as.data.frame(t(OTU))
OTU_abundance <- as.data.frame(rowMeans(OTU))
OTU_abundance$OTU <- row.names(OTU_abundance)
colnames(OTU_abundance)[1] <- "mean_abundance"
#--提取相关矩阵
cor = corrlation
head(cor)

result2 <- model_igraph(cor =  cor,
                        method = "cluster_fast_greedy",
                        seed = 123
)
node = result2[[1]]
head(node)

dat = result2[[2]]
head(dat)
tem = data.frame(mod = dat$model,col = dat$color) %>%  
  dplyr::distinct( mod, .keep_all = TRUE)  
col = tem$col
names(col) = tem$mod

#---node节点注释#-----------
otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = as.data.frame(vegan_tax(ps))
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)
#-----计算边#--------
edge = edgeBuild(cor = cor,node = node)
colnames(edge)[8] = "cor"
head(edge)

tem2 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
  dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
head(tem2)

tem3 = dat %>% 
  dplyr::select(OTU,model,color) %>%
  dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
  dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
head(tem3)

tem4 = tem2 %>%inner_join(tem3)
head(tem4)

edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                        manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
)
head(edge2)
col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% 
  select(color,manual)
col0 = col_edge$manual
names(col0) = col_edge$color

library(ggnewscale)
#a <- intersect(res_tumor[res_tumor$p.adj<0.05,]$genus,test[test$qval<0.05&test$statistic>0,]$feature)
#write.table(as.data.frame(a),"prevalance_abundance_overlap.txt",row.names = F,col.names = F,quote = F)

#p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
#                              data = edge2, size = 1) +
#  scale_colour_manual(values = col0) 


edge2$cor <- factor(edge2$cor,levels=c("+","-"))
p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color=cor),data = edge2, size = 0.5) +
  scale_color_manual(values=c("#FBEBB7","#BCDEE7"))


node_list <- lapply(unique(edge$OTU_1), function(x){ y <- edge[edge$OTU_1==x|edge$OTU_2==x,]
y$OTU <- x
return(y)}
)
names(node_list) <- unique(edge$OTU_1)
node_list <- do.call("rbind",node_list)

hub <- as.data.frame(table(node_list$OTU))
positive.hub <- as.data.frame(table(node_list[node_list$cor=="+",]$OTU))%>%arrange(-Freq)
negative.hub <- as.data.frame(table(node_list[node_list$cor=="-",]$OTU))%>%arrange(-Freq)

library(ggpubr)
library(ggrepel)
# ggsave("./cs1.pdf",p1,width = 16,height = 14)

dat$label <- ifelse(dat$OTU%in%common,dat$OTU,NA)
dat <- merge(dat,OTU_abundance,by="OTU")
dat$hub <- ifelse(dat$OTU%in%positive.hub[1:5,]$Var1,dat$OTU,ifelse(dat$OTU%in%negative.hub[1:5,]$Var1,dat$OTU,NA))
dat$direction <- ifelse(dat$OTU%in%positive.hub[1:5,]$Var1,"positive",ifelse(dat$OTU%in%negative.hub[1:5,]$Var1,"negative",NA))

table(dat$model)

col <- c(col,"#C1C1C1")
names(col)[22] <- "mini_model"
# dat$model_new <- ifelse(dat$model%in%c("model_1","model_5"),"model_1",
#                         ifelse(dat$model%in%c("model_2","model_4"),"model_2",
#                                ifelse(dat$model=="model_3","model_3","mini_model")))

dat$model_new <- ifelse(dat$model%in%c("model_1","model_2","model_3","model_4"),"model_1","mini_model")

col2 <- c("#DC4636","#3885B6","#FDA94D","#C1C1C1")
names(col2) <- c("model_1","model_2","model_3","mini_model")

p2 = p1 +
  new_scale_color() +
  geom_point(aes(X1, X2,color =model_new,size=mean_abundance*100), data = dat) +
  scale_colour_manual(values = col2) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

p2+geom_text_repel(aes(X1, X2,label=OTU,color=model_new), data = dat,min.segment.length = Inf, seed = 42,
                   box.padding = 0.2,fontface="bold")

#dat$direction <- factor(dat$direction,level=c("positive","negative"))
dat2 <- dat[is.na(dat$label)==F,]
p2+
  geom_text_repel(aes(X1, X2,label=label), data = dat2,min.segment.length = Inf, seed = 42,
                  box.padding = 0.2,fontface="bold")
# ggsave("./cs2.pdf",p2,width = 16,height = 14)


#挑丰度前5的粪便也能检测到的菌来标记
dat2 <- dat[is.na(dat$label)==F,]%>%arrange(-mean_abundance)
dat2 <- dat2[1:5,]
dat$group <- ifelse(is.na(dat$label)==T,"gut absent","gut present")
dat$group <- factor(dat$group,levels = c("gut absent","gut present"))
nonE.coli_dat <- dat[is.na(dat$label)==F,]%>%mutate(Enterotype="ETnonE")

dat$group_new <- ifelse(dat$group=="gut present"&dat$model_new=="model_1","show_1",
                        ifelse(dat$group=="gut present"&dat$model_new=="mini_model","show_2",
                               ifelse(dat$group=="gut absent"&dat$model_new=="model_1","hide_1","hide_2")))
dat$group_new <- factor(dat$group_new,levels = c("show_1","show_2","hide_1","hide_2"))

p2 <- p1 +
  new_scale_color() +
  geom_point(aes(X1, X2,color =model_new,size=mean_abundance*100,fill=group_new), data = dat,shape=21) +
  scale_colour_manual(values = col2) +
  scale_fill_manual(values=c("#DC4636","#C1C1C1","white","white"))+
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

pdf("./plot/ETnonE_network.pdf",width=7,height=5)
p2+
  geom_text_repel(aes(X1, X2,label=label), data = dat2,min.segment.length = Inf, seed = 42,
                  box.padding = 0.2,fontface="bold")
dev.off()

write_xlsx(list(ETnonE_network_info1=dat,ETnonE_network_info2=nonE.coli_network),path="/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/supply_tables/ETnonE_network_info.xlsx")

#========画参与组成网络的菌的占比比较图案================
dat_plot <- rbind(E.coli_dat[,c("OTU","group","Enterotype","model_new")],
                  nonE.coli_dat[,c("OTU","group","Enterotype","model_new")])

dat_plot$Enterotype <- factor(dat_plot$Enterotype,levels=c("ETnonE","ETE"))
pval1 <- chisq.test(table(dat_plot$Enterotype,dat_plot$model_new))[["p.value"]]
p1.1 <- ggplot(dat_plot,aes(x =Enterotype, fill = model_new)) +
  geom_bar(position = "fill")+scale_fill_manual(values = col2)+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "none")

p1.4 <- ggplot()+geom_text(aes(x =2,y = 0,
                               label =  paste0("p=",signif(pval1,2)),size =5.0,family = "sans",fontface = 1))+
  xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background=element_rect(I(0),linetype=0),
        plot.margin = margin(0,0,0,0, unit = "cm"),
        legend.position = "none")
p1 <- p1.4+p1.1+plot_layout(heights = c(0.1,1),widths = c(1,1),ncol = 1,nrow = 2)
p1
ggsave("plot/network_contribution.pdf",device = cairo_pdf,width =3, height =5,p1)

#--随即取出任意比例节点-网络鲁棒性#---------
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/Network_stability_pre_encapsulation.R")
library(ggClusterNet)
library(phyloseq)

#计算每个物种的平均丰度，使用测序深度标准化
#sp.ra<-colMeans(otutab)/mean(rowSums(otutab))   #relative abundance of each species
#这个我们算了，是OTU_abundance 

cor = corrlation
head(cor)

#存在某些情况计算不出来相关系数，定义相关为0
cor[is.na(cor)]<-0
#-去除自相关的点
diag(cor)<-0  
#-查看网络边的数量
sum(abs(cor)>0)/2
#网络中节点的数量
sum(colSums(abs(cor))>0)  # node number: number of species with at least one linkage with others.
#去除没有任何相关的节点.
network.raw<-cor[colSums(abs(cor))>0,colSums(abs(cor))>0]
#对应的删除otu表格otu
#sp.ra2<-sp.ra[colSums(abs(cor))>0]
sp.ra2 <- OTU_abundance[OTU_abundance$OTU%in%rownames(network.raw),]
sp.ra2 <- sp.ra2[match(sp.ra2$OTU,rownames(network.raw)),]  #check if matched

## 鲁棒性评估robustness simulation 
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained


#输入相关矩阵 OTU表格
Weighted.simu<-rmsimu(netRaw = network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance, 
                      abundance.weighted=T,nperm=100)
head(Weighted.simu)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance,
                        abundance.weighted=F,nperm=100)
head(Weighted.simu)
# tem = ps %>% sample_data() %>% .$Group %>% unique() %>% as.character()

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),
                 rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 Group="nonE.coli")

head(dat1)


write.table(dat1,file="./fastp/ETnonE_tumor/nonE.coli_tumor_robustness.txt",row.names = F,col.names = T,sep="\t",quote = F)



#==========保存每次抽样的数据画柱状图===============、
rmsimu_raw<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remains
  }))
}

Weighted.simu.raw<-rmsimu_raw(netRaw = network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance, 
                              abundance.weighted=T,nperm=100)

library(reshape2)
library(tidyverse)
Weighted.simu.raw<-as.data.frame(Weighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Weighted.simu.raw<-as.data.frame(Weighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Weighted.simu.raw<-melt(Weighted.simu.raw,id.vars = "Proportion.removed")
Weighted.simu.raw <- Weighted.simu.raw%>%mutate(weighted="weighted",Group="nonE.coli")

Unweighted.simu.raw<-rmsimu_raw(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2$mean_abundance,
                                abundance.weighted=F,nperm=100)

library(reshape2)
Unweighted.simu.raw<-as.data.frame(Unweighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Unweighted.simu.raw<-as.data.frame(Unweighted.simu.raw)%>%dplyr::mutate(Proportion.removed=seq(0.05,1,by=0.05))
Unweighted.simu.raw<-melt(Unweighted.simu.raw,id.vars = "Proportion.removed")
Unweighted.simu.raw <- Unweighted.simu.raw%>%mutate(weighted="unweighted",Group="nonE.coli")

ran_robust <- rbind(Weighted.simu.raw,Unweighted.simu.raw)
write.table(ran_robust,file="./fastp/ETnonE_tumor/nonE.coli_tumor_ran_robustness_raw.txt",sep="\t",quote=F,col.names = T,row.names = F)

#====--去除关键节点-网络鲁棒性#------
library(ggClusterNet)
library(phyloseq)
library(igraph)

cor = corrlation
head(cor)

#存在某些情况计算不出来相关系数，定义相关为0
cor[is.na(cor)]<-0
#-去除自相关的点
diag(cor)<-0  
#-查看网络边的数量
sum(abs(cor)>0)/2
#网络中节点的数量
sum(colSums(abs(cor))>0)  # node number: number of species with at least one linkage with others.
#去除没有任何相关的节点.
network.raw<-cor[colSums(abs(cor))>0,colSums(abs(cor))>0]

##read otu table

otutab<- ps %>% 
  scale_micro() %>%
  subset_taxa(
    row.names(tax_table(ps)) %in% row.names(network.raw)
  ) %>%
  vegan_otu() %>% 
  t() %>%
  as.data.frame()



#对应的删除otu表格otu
#sp.ra2<-sp.ra[colSums(abs(cor))>0]
sp.ra2 <- OTU_abundance[OTU_abundance$OTU%in%rownames(network.raw),]
sp.ra2 <- sp.ra2[match(sp.ra2$OTU,rownames(network.raw)),]  #check if matched
sp.ra <- sp.ra2$mean_abundance
names(sp.ra) <- sp.ra2$OTU
## robustness simulation 
#input network matrix, number of removed keystone species, keystonespecies list,  and ra of all species
#return the proportion of species remained

#get the keystone species list
igraph = make_igraph(cor)

degree = TRUE
zipi = FALSE
if (degree) {
  ret3 = node_properties(igraph) %>% 
    as.data.frame() %>%
    filter(!is.na(igraph.degree) ) %>% 
    arrange(desc(igraph.degree)) 
  head(ret3)
  tem = round(length(ret3$igraph.degree) * 0.20,0)    ##这里设置对于hub菌的定义，为了画图好看，改成20%(原本是10%)
  module.hub = row.names(ret3)[1:tem]
}



if (zipi ) {
  res = ZiPiPlot(igraph = igraph,method = "cluster_fast_greedy")
  p <- res[[1]]
  model = res[[2]] %>% filter(roles == "Module hubs")
  head(model)
  model$roles %>% unique()
  module.hub <- as.character(row.names(model))  
}

Weighted.simu<-rmsimu2(netRaw=network.raw,
                       rm.p.list=1:length(module.hub),
                       keystonelist=module.hub,
                       sp.ra=sp.ra,
                       abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu2(netRaw=network.raw, rm.p.list=1:length(module.hub),
                         keystonelist=module.hub,
                         sp.ra=sp.ra, abundance.weighted=F,nperm=100)

dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),
                 rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),
                              each=length(module.hub)),
                 Group= "nonE.coli")

head(dat1)
currentdat = dat1

write.table(dat1,file="./fastp/ETnonE_tumor/nonE.coli_tumor_hub_robustness.txt",row.names = F,col.names = T,sep="\t",quote = F)


#======保存每次hub抽样的数据画柱状图============
rmsimu2_raw<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remains
  }))
}


Weighted.simu.raw<-rmsimu2_raw(netRaw=network.raw,
                               rm.p.list=1:length(module.hub),
                               keystonelist=module.hub,
                               sp.ra=sp.ra,
                               abundance.weighted=T,nperm=100)

library(reshape2)
library(tidyverse)
Weighted.simu.raw<-as.data.frame(Weighted.simu.raw)
Weighted.simu.raw<-Weighted.simu.raw%>%dplyr::mutate(Number.hub.removed=row.names(Weighted.simu.raw))
Weighted.simu.raw<-melt(Weighted.simu.raw,id.vars = "Number.hub.removed")
Weighted.simu.raw <- Weighted.simu.raw%>%mutate(weighted="weighted",Group="nonE.coli")

Unweighted.simu.raw<-rmsimu2_raw(netRaw=network.raw,
                                 rm.p.list=1:length(module.hub),
                                 keystonelist=module.hub,
                                 sp.ra=sp.ra,
                                 abundance.weighted=F,nperm=100)

library(reshape2)
Unweighted.simu.raw<-as.data.frame(Unweighted.simu.raw)
Unweighted.simu.raw<-Unweighted.simu.raw%>%dplyr::mutate(Number.hub.removed=row.names(Unweighted.simu.raw))
Unweighted.simu.raw<-melt(Unweighted.simu.raw,id.vars = "Number.hub.removed")
Unweighted.simu.raw <- Unweighted.simu.raw%>%mutate(weighted="unweighted",Group="nonE.coli")

hub_robust <- rbind(Weighted.simu.raw,Unweighted.simu.raw)
write.table(hub_robust,file="./fastp/ETnonE_tumor/nonE.coli_tumor_hub_robustness_raw.txt",sep="\t",quote=F,col.names = T,row.names = F)

#====加载易损性函数================
vulnerability = function(cor = cor){
  
  cor[abs(cor)>0]<-1 # adjacency matrix
  g = graph_from_adjacency_matrix(as.matrix(cor), 
                                  mode="undirected", 
                                  weighted = NULL, diag = FALSE,
                                  add.colnames = NULL) # note: this graph contains isolated nodes.
  # remove isolated nodes
  iso_node_id = which(igraph::degree(g)==0)
  g2 = igraph::delete_vertices(g, iso_node_id) # graph without isolated nodes
  
  #check node number and links
  length(V(g2));length(E(g2))   
  
  # calculate vulnerability of each node
  node.vul<-info.centrality.vertex(g2)
  return(max(node.vul))
}

#======计算网络易损性================      
i = 1
A = c()

cor = corrlation

cor[abs(cor)>0]<-1 # adjacency matrix
vulnerability(cor = cor)
##nonE.coli网络易损性为 0.07617797

#=====网络属性图从此画====
#======colorset=============
col1 <- c("#2482BC","#f99170","#e73133")
col2 <- c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
col3 <- c("#fcd364","#b66d32")
col4 <- c("#b39658","#b66d32")

#=======综合画两个随机网络稳定性的图=================
robust_E.coli <- read.table("./fastp/E.coli_tumor/E.coli_tumor_robustness.txt",header=T)
robust_nonE.coli <- read.table("./fastp/ETnonE_tumor/nonE.coli_tumor_robustness.txt",header=T)

robust_weight <- rbind(robust_E.coli[robust_E.coli$weighted=="weighted",],robust_nonE.coli[robust_nonE.coli$weighted=="weighted",])

ggplot(data=robust_weight, 
       aes(x=Proportion.removed, y=remain.mean, group=Group, color=Group)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="weighted")+
  theme_light()

robust_unweight <- rbind(robust_E.coli[robust_E.coli$weighted=="unweighted",],robust_nonE.coli[robust_nonE.coli$weighted=="unweighted",])

ggplot(data=robust_unweight, 
       aes(x=Proportion.removed, y=remain.mean, group=Group, color=Group)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="unweighted")+
  theme_light()

ggsave("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/plot/网络属性图/ETE_unweighted.pdf",width = 8,height = 6)

#=======综合画两个随机网络稳定性的柱状图=================
robust_E.coli <- read.table("./fastp/E.coli_tumor/E.coli_tumor_robustness.txt",header=T)
robust_nonE.coli <- read.table("./fastp/ETnonE_tumor/nonE.coli_tumor_robustness.txt",header=T)
robust_E.coli_raw <- read.table("./fastp/E.coli_tumor/E.coli_tumor_ran_robustness_raw.txt",header=T)
robust_nonE.coli_raw <- read.table("./fastp/ETnonE_tumor/nonE.coli_tumor_ran_robustness_raw.txt",header=T)

#########加权
robust_weight <- rbind(robust_E.coli[robust_E.coli$weighted=="weighted",],robust_nonE.coli[robust_nonE.coli$weighted=="weighted",])
robust_weight_raw <- rbind(robust_E.coli_raw[robust_E.coli_raw$weighted=="weighted",],robust_nonE.coli_raw[robust_nonE.coli_raw$weighted=="weighted",])

robust_weight$pvalue <- NA

for (i in unique(robust_weight$Proportion.removed)){
  pvalue <- wilcox.test(value~Group,data=robust_weight_raw[robust_weight_raw$Proportion.removed==i,])[["p.value"]]
  robust_weight[robust_weight$Proportion.removed==i,"pvalue"] <- pvalue
}
###注意要根据robust_weight的pvalue后续修改图中的p值

robust_weight$Proportion.removed<- as.factor(robust_weight$Proportion.removed)
robust_weight$group <- ifelse(robust_weight$Group=="nonE.coli","ETnonE","ETE")
robust_weight$group <- factor(robust_weight$group,levels=c("ETnonE","ETE"))
###根据robust_weight的pvalue后续修改图中的p值
robust_label <- robust_weight[robust_weight$Group=="E.coli",] 
robust_label$star <- cut(robust_label$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
ggplot()+
  geom_errorbar(data=robust_weight,
                aes(x=Proportion.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_weight,
           aes(x=Proportion.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="weighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  # coord_cartesian(ylim = c(0.8, 1))+
  theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Proportion.removed,label=star),position=position_dodge(0.9),size=3,fontface="bold")


##下面这段代码画的图是标具体P值的
robust_label <- robust_weight[robust_weight$Group=="E.coli",]  
ggplot()+
  geom_errorbar(data=robust_weight,
                aes(x=Proportion.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_weight,
           aes(x=Proportion.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="weighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  coord_cartesian(ylim = c(0, 1))+theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Proportion.removed,label=signif(pvalue,digits = 3)),position=position_dodge(0.9),size=3,fontface="bold")


#########未加权
robust_unweight <- rbind(robust_E.coli[robust_E.coli$weighted=="unweighted",],robust_nonE.coli[robust_nonE.coli$weighted=="unweighted",])
robust_unweight_raw <- rbind(robust_E.coli_raw[robust_E.coli_raw$weighted=="unweighted",],robust_nonE.coli_raw[robust_nonE.coli_raw$weighted=="unweighted",])

robust_unweight$pvalue <- NA

for (i in unique(robust_unweight$Proportion.removed)){
  pvalue <- wilcox.test(value~Group,data=robust_unweight_raw[robust_unweight_raw$Proportion.removed==i,])[["p.value"]]
  robust_unweight[robust_unweight$Proportion.removed==i,"pvalue"] <- pvalue
}

robust_unweight$Proportion.removed<- as.factor(robust_unweight$Proportion.removed)
robust_unweight$group <- ifelse(robust_unweight$Group=="nonE.coli","ETnonE","ETE")
robust_unweight$group <- factor(robust_unweight$group,levels=c("ETnonE","ETE"))
###根据robust_unweight的pvalue后续修改图中的p值
robust_label <- robust_unweight[robust_unweight$Group=="E.coli",] 
robust_label$star <- cut(robust_label$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns"))
p1 <- ggplot()+
  geom_errorbar(data=robust_unweight,
                aes(x=Proportion.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_unweight,
           aes(x=Proportion.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="unweighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  #coord_cartesian(ylim = c(0.8, 1))+
  theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Proportion.removed,label=star),position=position_dodge(0.9),size=3,fontface="bold")


###下面是标具体P值的
p2 <- ggplot()+
  geom_errorbar(data=robust_unweight,
                aes(x=Proportion.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_unweight,
           aes(x=Proportion.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="unweighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  coord_cartesian(ylim = c(0, 1))+theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Proportion.removed,label=signif(pvalue,digits = 3)),position=position_dodge(0.9),size=3,fontface="bold")

library(gridExtra)
pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/plot/网络属性图/柱状图(两张一样).pdf",width = 18,height = 5)
grid.arrange(p1, p2, ncol = 2)
dev.off()

#=======综合画两个hub网络稳定性的图=================
##只画6个菌的
robust_E.coli <- read.table("./fastp/E.coli_tumor/E.coli_tumor_hub_robustness.txt",header=T)
robust_E.coli <- robust_E.coli[robust_E.coli$Number.hub.removed%in%c(1:6),]
robust_nonE.coli <- read.table("./fastp/ETnonE_tumor/nonE.coli_tumor_hub_robustness.txt",header=T)

robust_weight <- rbind(robust_E.coli[robust_E.coli$weighted=="weighted",],robust_nonE.coli[robust_nonE.coli$weighted=="weighted",])

ggplot(data=robust_weight, 
       aes(x=Number.hub.removed, y=remain.mean, group=Group, color=Group)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="weighted")+
  theme_light()


robust_unweight <- rbind(robust_E.coli[robust_E.coli$weighted=="unweighted",],robust_nonE.coli[robust_nonE.coli$weighted=="unweighted",])

ggplot(data=robust_unweight, 
       aes(x=Number.hub.removed, y=remain.mean, group=Group, color=Group)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue","red"))+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+labs(title="unweighted")+
  theme_light()

ggsave("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/plot/网络属性图/hub_bac.pdf",width = 8,height = 5)

#=======综合画两个hub网络稳定性的柱状图=================
robust_E.coli <- read.table("./fastp/E.coli_tumor/E.coli_tumor_hub_robustness.txt",header=T)
robust_nonE.coli <- read.table("./fastp/ETnonE_tumor/nonE.coli_tumor_hub_robustness.txt",header=T)
robust_E.coli_raw <- read.table("./fastp/E.coli_tumor/E.coli_tumor_hub_robustness_raw.txt",header=T)
robust_nonE.coli_raw <- read.table("./fastp/ETnonE_tumor/nonE.coli_tumor_hub_robustness_raw.txt",header=T)

#########加权
robust_weight <- rbind(robust_E.coli[robust_E.coli$weighted=="weighted",],robust_nonE.coli[robust_nonE.coli$weighted=="weighted",])
robust_weight_raw <- rbind(robust_E.coli_raw[robust_E.coli_raw$weighted=="weighted",],robust_nonE.coli_raw[robust_nonE.coli_raw$weighted=="weighted",])

table(robust_weight_raw$Number.hub.removed)
##根据两个型中都有的踢菌的数目来设置循环
robust_weight <-robust_weight[robust_weight$Number.hub.removed%in%c(1:6),] 
robust_weight_raw <-robust_weight_raw[robust_weight_raw$Number.hub.removed%in%c(1:6),] 

robust_weight$pvalue <- NA

for (i in 1:6){
  pvalue <- wilcox.test(value~Group,data=robust_weight_raw[robust_weight_raw$Number.hub.removed==i,])[["p.value"]]
  robust_weight[robust_weight$Number.hub.removed==i,"pvalue"] <- pvalue
  cat(i,"\n")
}

robust_weight$Number.hub.removed <- as.factor(robust_weight$Number.hub.removed)
robust_weight$group <- ifelse(robust_weight$Group=="nonE.coli","ETnonE","ETE")
robust_weight$group <- factor(robust_weight$group,levels=c("ETnonE","ETE"))
###根据robust_weight的pvalue后续修改图中的p值
robust_label <- robust_weight[robust_weight$Group=="E.coli",] 
robust_label$star <- cut(robust_label$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "ns"))
robust_weight$remain.mean <- as.numeric(robust_weight$remain.mean)

p1 <- ggplot()+
  geom_errorbar(data=robust_weight,
                aes(x=Number.hub.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_weight,
           aes(x=Number.hub.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Number of hub species removed")+
  ylab("Proportion of species remained")+labs(title="weighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  coord_cartesian(ylim = c(0, 1))+theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Number.hub.removed,label=star),position=position_dodge(0.9),size=3,fontface="bold")


##这是标P值的
p2 <- ggplot()+
  geom_errorbar(data=robust_weight,
                aes(x=Number.hub.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_weight,
           aes(x=Number.hub.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Number of hub species removed")+
  ylab("Proportion of species remained")+labs(title="weighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  coord_cartesian(ylim = c(0, 1))+theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Number.hub.removed,label=signif(pvalue,digits = 3)),position=position_dodge(0.9),size=3,fontface="bold")


#########未加权
robust_unweight <- rbind(robust_E.coli[robust_E.coli$weighted=="unweighted",],robust_nonE.coli[robust_nonE.coli$weighted=="unweighted",])
robust_unweight_raw <- rbind(robust_E.coli_raw[robust_E.coli_raw$weighted=="unweighted",],robust_nonE.coli_raw[robust_nonE.coli_raw$weighted=="unweighted",])

robust_unweight <-robust_unweight[robust_unweight$Number.hub.removed%in%c(1:6),] 
robust_unweight_raw <-robust_unweight_raw[robust_unweight_raw$Number.hub.removed%in%c(1:6),] 

robust_unweight$pvalue <- NA

for (i in 1:6){
  pvalue <- wilcox.test(value~Group,data=robust_unweight_raw[robust_unweight_raw$Number.hub.removed==i,])[["p.value"]]
  robust_unweight[robust_unweight$Number.hub.removed==i,"pvalue"] <- pvalue
}


robust_unweight$Number.hub.removed <- as.factor(robust_unweight$Number.hub.removed)
robust_unweight$group <- ifelse(robust_unweight$Group=="nonE.coli","ETnonE","ETE")
robust_unweight$group <- factor(robust_unweight$group,levels=c("ETnonE","ETE"))
###根据robust_unweight的pvalue后续修改图中的p值
robust_label <- robust_unweight[robust_unweight$Group=="E.coli",] 
robust_label$star <- cut(robust_label$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))

p1 <- ggplot()+
  geom_errorbar(data=robust_unweight,
                aes(x=Number.hub.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_unweight,
           aes(x=Number.hub.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Number of hub species removed")+
  ylab("Proportion of species remained")+labs(title="unweighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  coord_cartesian(ylim = c(0, 1))+theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Number.hub.removed,label=star),position=position_dodge(0.9),size=3,fontface="bold")


##这是标P值的
p2 <- ggplot()+
  geom_errorbar(data=robust_unweight,
                aes(x=Number.hub.removed,
                    ymin=remain.mean-remain.sd,
                    ymax=remain.mean+remain.sd,color=group),
                position = position_dodge(0.8), width = 0.3)+
  geom_bar(data=robust_unweight,
           aes(x=Number.hub.removed,y=remain.mean,fill=group),stat = "identity",position = "dodge",width = 0.8
  )+
  scale_color_manual(values=col3)+
  scale_fill_manual(values=col3)+
  xlab("Number of hub species removed")+
  ylab("Proportion of species remained")+labs(title="unweighted")+
  # scale_y_continuous(limits = c(0.5, 1))+
  coord_cartesian(ylim = c(0, 1))+theme_light()+
  geom_text(data=robust_label,aes(y=1,x=Number.hub.removed,label=signif(pvalue,digits = 3)),position=position_dodge(0.9),size=3,fontface="bold")

pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/16S_250509/250711网络分析/plot/网络属性图/hub_bac柱状图(2张一样).pdf",height = 5,width = 15)
p1+p2
#grid.arrange(p1, p2, ncol = 2)
dev.off()
