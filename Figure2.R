rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/")
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/functions.R")
#==================
#Fig2-解释肠型
#===========================
library(ggpubr)
library(vegan)
library(ggplot2)
library(ggsci)
library(data.table)
library(foreach)
library(dplyr)
library(crayon)
library(nlme)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggprism)
library(cowplot)
library(survminer)
library(survival)
library(mRMRe)
set.seed(12345)

#======colorset=============
col1 <- c("#2482BC","#f99170","#e73133") ##健康人，肝硬化，HCC
col2 <- c("#5ba787","#fa9fcb","#b35c31","#6e7ca5")  ##4种肠型
col3 <- c("#fcd364","#b35c31")  ##ETnonE,ETE
col4 <- c("#BD9E4B","#BF6C0F")  ##画小提琴图加粗的线的颜色

#======Fig2A_balance feature======================
mi <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Enterotype/mrmr_MI.txt")
taxa_up <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/Score.feature.up1.txt")
taxa_down <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/Score.feature.down1.txt")

taxa_up$effect <- mi$X.results.mi_matrix.1....[match(taxa_up$x,rownames(mi))]
taxa_down$effect <- mi$X.results.mi_matrix.1....[match(taxa_down$x,rownames(mi))]

taxa <- rbind(taxa_up,taxa_down)

taxa <- taxa%>%mutate(EAE_DA= ifelse(effect > 0, "ETnonE",'ETE')) 

taxa= taxa[order(taxa$effect,decreasing = T),]
colnames(taxa)[1] <- "species"
taxa$species= factor(taxa$species, levels= taxa$species)
##
table(taxa$EAE_DA)
taxa$EAE_DA <- factor(taxa$EAE_DA,levels=c("ETnonE","ETE"))
##
p <- ggplot(taxa, aes(y = species, x = effect, label = species)) +
  geom_col(aes(fill = EAE_DA), width=0.5) +
  coord_flip()+
  #geom_point() +
  geom_point(aes(color=EAE_DA))+
  # scale_fill_manual(values = c(col1[1],col1[3]))+ 
  #scale_color_manual(values=c(col1[1],col1[3]))+
  scale_fill_manual(values = col3)+ 
  scale_color_manual(values=col3)+
  #coord_flip()+ 
  theme_classic()+
  # ylim(-2, 2) +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  ) +xlim(-0.8, 0.8) +
  geom_text_repel(
    seed = 233,
    size = 3.5,
    color = 'black',
    min.segment.length = 0,
    force = 2,
    force_pull = 2,
    box.padding = 0.1,
    max.overlaps = Inf,
    segment.linetype = 2, #线段类型,1为实线,2-6为不同类型虚线
    segment.color = 'black', #线段颜色
    segment.alpha = 0.7, #线段不透明度
    #nudge_x = 1 - abs(dadraw$cohens_d), #标签x轴起始位置调整
    nudge_x = 0.3, #标签x轴起始位置调整
    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 0.5
  ) +
  ggtitle("ETnonE VS ETE")+
  theme(legend.position = "none",
        aspect.ratio = 0.3,
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #axis.text.x = element_blank(),
        plot.title = element_text(size=8,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.1, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5)
  )

pdf("./Fig2A_feature_taxa.pdf",height = 9,width = 12)
cowplot::plot_grid(p)
dev.off()
save(taxa,file="./Fig2A.RData")

#===========Fig2B_humann_PCA==============
#========预处理humann================
library(readr)
abu_table <- read_tsv("/groups/ProHuShiX/home/share/BIGMetaG_result/humann/pathabundance/pathabundance_all.tsv")
colnames(abu_table) <- gsub("_.*$","",colnames(abu_table))
colnames(abu_table) <- gsub("-",".",colnames(abu_table))
abu_table <- abu_table[1:34305,]##最后面这几行行名不知道为什么是脚本命令，把它删掉

score <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt",header=T)
mgs_path <- abu_table[,c("# Pathway",score$Row.names)]
mgs_path <- as.data.frame(mgs_path)
row.names(mgs_path) <- mgs_path$`# Pathway`
mgs_path$`# Pathway` <- NULL
mgs_path=mgs_path[!rownames(mgs_path) %like% "\\|",]
mgs_path=mgs_path[rowSums(mgs_path>0)>(ncol(mgs_path)*0.1),]
mgs_path=mgs_path[!rownames(mgs_path) %like% c("UNMAPPED"),]
mgs_path=mgs_path[!rownames(mgs_path) %like% c("UNINTEGRATED"),]


mgs_path=apply(mgs_path,2,function(x){
  x=x/sum(x)
  return(x)
})
colSums(mgs_path)[1:10]


#先修改一下名字，不然不能clr
mgs_path_name <- as.data.frame(rownames(mgs_path))
colnames(mgs_path_name) <- c("name")
mgs_path_name$num <- paste("H",1:nrow(mgs_path_name),sep = "")
row.names(mgs_path) <- mgs_path_name$num

#对丰度进行clr转化
mgs_path_clr <- mgs_path
mgs_path_clr=transform_and_filter_taxa(mgs_path_clr,samples_row = F,method = "clr",missing_filter = 0)

humann_f <- mgs_path

beta_diversity=vegdist(humann_f,method = "bray")
pcoa_analysis=cmdscale(beta_diversity,k=4,eig=TRUE)
pc<- round(pcoa_analysis$eig/sum(pcoa_analysis$eig)*100,digits=2)#这一步是算解释度的，pc的第一和第二位分别代表X轴和y轴的解释度
pcoa_analysis <- as.data.frame(pcoa_analysis[["points"]])
pcoa_analysis <- merge(pcoa_analysis,cluster,by="row.names",all=F)
pcoa_analysis $Enterotype <- ifelse(pcoa_analysis$Cluster=="4","ETE","ETnonE")
pcoa_analysis $Enterotype<- factor(pcoa_analysis $Enterotype,levels = c("ETnonE","ETE"))

cov.data <- pcoa_analysis[match(rownames(humann_f),pcoa_analysis$Row.names),]
cov.data$Enterotype <- as.factor(cov.data$Enterotype)
otu.adonis=adonis2(humann_f ~ cov.data$Enterotype ,  permutations = 999, method = "bray")

bdiversity1 <- ggplot(pcoa_analysis,aes(V1, V2, color = Enterotype))+
  geom_point(size = 2)+#ggtitle("France")+
  stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)+
  xlab(paste0("PCoA1(",pc[1],"%)"))+ylab(paste0("PCoA2(",pc[2],"%)")) +
  scale_color_manual(values =c(col3))+
  #guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position =c(.85, .88),
        legend.title = element_text(size =12,family ="bold",color="black"),
        legend.text = element_text(size =9,family ="sans",color="black"),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))+
  annotate("text",x=-0.25,y=0.13,label = paste("PERMANOVA:\np=",otu.adonis[1,5],
                                               "\nR2 = ", round(otu.adonis[1,3],4),sep = ""
  ))
bdiversity1

bdiversity1 <- ggExtra::ggMarginal(bdiversity1, type = "density", groupColour = F, groupFill = T,
                                   xparams = list(alpha = 0.5,color=NA),
                                   yparams = list(alpha = 0.5,color=NA))


ggsave("./Fig2B_humann_PCA.pdf",device = cairo_pdf,width =5.5, height =5,bdiversity1)
save(pcoa_analysis,otu.adonis,file="./Fig2B.RData")

#========Fig2S_humann_score_PCA=============
otu.adonis=adonis2(humann_f ~ cov.data$balance_value,permutations = 999, method = "bray")

pdf("./plots/Fig2S_humann_score_PCA250913.pdf",width = 4.8,height = 5)
ggplot(pcoa_analysis,aes(V1, V2, color = balance_value))+
  geom_point(size = 2)+#ggtitle("France")+
  #stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)
  xlab(paste0("PCoA1(",pc[1],"%)"))+ylab(paste0("PCoA2(",pc[2],"%)"))+
  labs(color="balance_score") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.6,9.6))+
  theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.title = element_text(hjust = 0.5),
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
  annotate("text",x=0.23,y=0.13,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis[1,5],2),
                                            "\nR2=",signif(otu.adonis[1,3],2)))+ggtitle("FAH corhort1")
dev.off()
save(pcoa_analysis,otu.adonis,file="./Fig2S_human_score_PCA.RData")

#=========Fig2C_humann_feature============
#==========差异分析=============
mgs_path_clr <- mgs_path_clr[match(score$Row.names,row.names(mgs_path_clr)),]
colnames(score) <- gsub("Cluster","cluster",colnames(score))
mgs_path <- as.data.frame(t(mgs_path))
mgs_path <- mgs_path[match(row.names(mgs_path_clr),row.names(mgs_path)),match(colnames(mgs_path_clr),colnames(mgs_path))]
identical(score$Row.names,row.names(mgs_path_clr))

res <- list()
i=1
for (i in 1:ncol(mgs_path_clr)){
  tmp.data <- as.data.frame(cbind(mgs_path_clr[,i],mgs_path[i],score$cluster))
  colnames(tmp.data) <- c("value","value_raw","cluster")
  tmp.data$cluster <- ifelse(tmp.data$cluster==4,"ETE","ETnonE")
  p.value<- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  a <- tmp.data %>% filter(cluster == "ETE") %>% dplyr::select(value_raw)
  b <- tmp.data%>% filter(cluster == "ETnonE") %>% dplyr::select(value_raw)
  # a[a==0] <- NA
  #b[b==0] <- NA
  fold <- mean(na.omit(b[,1]))/mean(na.omit(a[,1]))
  res[[i]] <- data.frame(mgs_path_num=colnames(mgs_path_clr)[i],FC=fold,pvalue=p.value,
                         mean_ETE=mean(a[,1]),mean_ETnonE=mean(b[,1]),mean_all=mean(tmp.data$value_raw))
}
res <- do.call("rbind",res)
res$padj <- p.adjust(res$pvalue,"BH")
res <- res[res$padj<0.05,]
res <- merge(res,mgs_path_name,by.x="mgs_path_num",by.y="num",all.x=T)
res$short <- gsub("^.*:","",res$name)

# supplymentary放所有humann pathway的diff结果
write_xlsx(res,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/Fig2C_diff_humann_res.xlsx")


up<- res[res$FC>1,]
down<- res[res$FC<1,]

taxa_up <- up
taxa_down <- down
taxa_up <-taxa_up%>%arrange(padj)%>%mutate(feature="ETnonE")
taxa_down <- taxa_down%>%arrange(padj)%>%mutate(feature="ETE")

taxa_down <- taxa_down[1:5,]%>%arrange(padj)
forplot <- rbind(taxa_up[1:5,],taxa_down)
forplot$name <- gsub("^.*:","",forplot$name)
forplot <- forplot%>%arrange(-FC)
forplot$name <- factor(forplot$name,levels = rev(forplot$name))
forplot$feature <- factor(forplot$feature,levels = c("ETnonE","ETE"))
forplot$effect <- log2(forplot$FC)
forplot$stars <- cut(forplot$padj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
col3 <- c("#fcd364","#b35c31") 

p <- ggplot(forplot, aes(x = name, y = effect))+ 
  geom_segment( aes(x=name, xend=name, y=0, yend=effect,color=feature),size=0.8,linetype=2)+#使用reorder()排序变量
  geom_point(aes(color=feature),size =6) +
  geom_text(aes(label=stars), fontface="bold", color="white",size=3 ,nudge_y=0.005,nudge_x = 0)+
  #geom_col(aes(fill=feature))+
  # scale_y_reverse()+ #控制线段的参数，见下
  scale_color_manual(values = c(col3))+
  theme(panel.border = element_rect(size=1,fill=NA),
        axis.text.y = element_text(face = "bold",size = 10),
        legend.position = c(.75,.15),
        panel.background =element_blank(),
        #axis.line=element_line(),
        panel.grid.major.y = element_blank(),   ##设置x,y轴的主次提示线都为空白
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  geom_hline(yintercept=c(0),lty=3,col="darkgrey",lwd=0.5)+
  ylab("Log2(Fold Change)")+xlab("")+
  #ylim(-0.6,0.6)+
  coord_flip()
p

pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2C_feature_pathway250702.pdf",width = 8,height = 8)
p
dev.off()

#====correlation with BS/ETE====
load("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/442MGS_clinical.RData")
identical(MGS_meta$sample,rownames(mgs_path_clr))

MGS_meta$Age <- as.numeric(MGS_meta$Age);MGS_meta$Gender<-as.factor(MGS_meta$Gender);MGS_meta$BMI<-as.numeric(MGS_meta$BMI);MGS_meta$Smoke<-as.factor(MGS_meta$Smoke);MGS_meta$Drink<-as.factor(MGS_meta$Drink)
res1 <- list()
i=1
for (i in 1:ncol(mgs_path_clr)){
  tmp.data <- data.frame(value=mgs_path_clr[,i],value_raw=mgs_path[,i],
                         cluster=score$cluster,BS=score$balance_value,
                         age=MGS_meta$Age,gender=MGS_meta$Gender,BMI=MGS_meta$BMI,
                         smoke=MGS_meta$Drink,drink=MGS_meta$Smoke)
  tmp.data$BS <- as.numeric(tmp.data$BS)
  tmp.data$cluster <- ifelse(tmp.data$cluster==4,"ETE","ETnonE")
  p.value1<- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  p.value2<- summary(lm(value~BS+gender+age+BMI+smoke+drink, data=tmp.data))$coefficients["BS", "Pr(>|t|)"]
  estimate<- summary(lm(value~BS+gender+age+BMI+smoke+drink, data=tmp.data))$coefficients["BS", "Estimate"]
  a <- tmp.data %>% filter(cluster == "ETE") %>% select(value_raw)
  b <- tmp.data%>% filter(cluster == "ETnonE") %>% select(value_raw)
  fold <- mean(na.omit(b[,1]))/mean(na.omit(a[,1]))
  res1[[i]] <- data.frame(mgs_path_num=colnames(mgs_path_clr)[i],FC=fold,pvalue_cluster=p.value1,estimate=estimate,pvalue_BS=p.value2,
                          mean_ETE=mean(a[,1]),mean_ETnonE=mean(b[,1]),mean_all=mean(tmp.data$value_raw))
}
res1 <- do.call("rbind",res1)
res1$padj_BS <- p.adjust(res1$pvalue_BS,"BH")
res1$padj_cluster <- p.adjust(res1$pvalue_cluster,"BH")
#res1 <- res1[res1$padj<0.05,]
res1 <- merge(res1,mgs_path_name,by.x="mgs_path_num",by.y="num",all.x=T)

writexl::write_xlsx(res1,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/FigS2C_humann_BS.xlsx")



#========Fig2D_VF和ARG总量比较================
VF <- read.table("/groups/ProHuShiX/home/xiashufen/Metagenomics/VF/VF_merged_new.txt",row.names = 1,header=T,sep="\t")
VF <- VF[,colnames(VF)%in%balance$Row.names]
VF_sum <- colSums(VF)%>%as.data.frame()
colnames(VF_sum) <- "VF"
VF_sum <- merge(VF_sum,cluster,by="row.names")
VF_sum$Enterotype <- ifelse(VF_sum$Cluster=="4","ETE","ETnonE")
VF_sum$Enterotype <- factor(VF_sum$Enterotype,levels = c("ETnonE","ETE"))

ARG <- read.table("/groups/ProHuShiX/home/share/BIGMetaG_result/ARGs/ARGs_merged.txt",header=T,check.names = F,row.names = 1)
colnames(ARG) <- gsub("-",".",colnames(ARG))
ARG <- ARG[,colnames(ARG)%in%balance$Row.names]
ARG_sum <- colSums(ARG)%>%as.data.frame()
colnames(ARG_sum) <- "ARG"
ARG_sum <- merge(ARG_sum,cluster,by="row.names")
ARG_sum$Enterotype <- ifelse(ARG_sum$Cluster=="4","ETE","ETnonE")
ARG_sum$Enterotype <- factor(ARG_sum$Enterotype,levels = c("ETnonE","ETE"))

sum1 <- ggplot(VF_sum, aes(x=Enterotype, y=VF, color=Enterotype)) + 
  geom_violin(trim=FALSE,aes(fill=Enterotype,alpha=0.8))+
  geom_boxplot(width=0.1, fill="white")+
  theme_classic()+theme(plot.title = element_text(hjust=0.5),legend.position = "none")+
  stat_compare_means(comparisons =list(c("ETE","ETnonE")))+
  scale_color_manual(values = col4)+
  scale_fill_manual(values=col3)+labs(title = "VF")

sum2 <- ggplot(ARG_sum, aes(x=Enterotype, y=ARG, color=Enterotype)) + 
  geom_violin(trim=FALSE,aes(fill=Enterotype,alpha=0.8))+
  geom_boxplot(width=0.1, fill="white")+
  theme_classic()+theme(plot.title = element_text(hjust=0.5),legend.position = "none")+
  stat_compare_means(comparisons =list(c("ETE","ETnonE")))+
  scale_color_manual(values = col4)+
  scale_fill_manual(values=col3)+labs(title = "ARG")

library(patchwork)
pdf("./Fig2D_sum.pdf",width = 3.5,height = 10)
sum1+sum2+plot_layout(heights = c(1,1),widths = c(1,1),ncol = 1,nrow = 2)
dev.off()
save(VF_sum,ARG_sum,file="./Fig2D.RData")


#======Fig2S_VF_ARG_score_lm==============
lm1 <- ggplot(data =VF_sum,aes(x=balance_value,y=VF))+
  geom_point(aes(color=Enterotype))+geom_smooth(method = 'lm',color=col4[2],fill=col3[1])+
  scale_color_manual(values = col3)+
  stat_cor(method = 'spearman', label.y=30000)+
  theme(panel.background =element_blank(),
        axis.line=element_line(),
        plot.title = element_text(hjust=0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position =c(.85, .88)
  )+ggtitle("VF")

lm2 <- ggplot(data =ARG_sum,aes(x=balance_value,y=ARG))+
  geom_point(aes(color=Enterotype))+geom_smooth(method = 'lm',color=col4[2],fill=col3[1])+
  scale_color_manual(values = col3)+
  stat_cor(method = 'spearman', label.y=30000)+
  theme(panel.background =element_blank(),
        axis.line=element_line(),
        plot.title = element_text(hjust=0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position =c(.85, .88)
  )+ggtitle("ARG")

pdf("./Fig2S_sum_lm.pdf",width = 10,height = 5)
lm1+lm2+plot_layout(heights = c(1,1),widths = c(1,1),ncol = 2,nrow = 1)
dev.off()
save(VF_sum,ARG_sum,file="./Fig2S_sum_lm.RData")

#========Fig2E_VF_ARG_PCA============
VF=VF[rowSums(VF>0)>(ncol(VF)*0.1),]
VF=as.data.frame(apply(VF,2,function (x){
  x=x/sum(x)
  return(x)
}))
VF <- as.data.frame(t(VF))
save(VF,file="./trim_VF_new.RData")

Cluster <- balance
cluster <-  Cluster[,c("Row.names","Cluster","balance_value")]
row.names(cluster) <- cluster$Row.names
cluster$Row.names <- NULL

beta_diversity=vegdist(VF,method = "bray")
pcoa_analysis1=cmdscale(beta_diversity,k=4,eig=TRUE)
pc1<- round(pcoa_analysis1$eig/sum(pcoa_analysis1$eig)*100,digits=2)#这一步是算解释度的，pc的第一和第二位分别代表X轴和y轴的解释度
pcoa_analysis1 <- as.data.frame(pcoa_analysis1[["points"]])

pcoa_analysis1 <- merge(pcoa_analysis1,cluster,by="row.names",all=F)
pcoa_analysis1 $Enterotype <- ifelse(pcoa_analysis1$Cluster=="4","ETE","ETnonE")
pcoa_analysis1 $Enterotype<- factor(pcoa_analysis1 $Enterotype,levels = c("ETnonE","ETE"))

cov.data <- pcoa_analysis1[match(row.names(VF),pcoa_analysis1$Row.names),]
cov.data$Enterotype <- as.factor(cov.data$Enterotype)
otu.adonis1=adonis2(VF ~ cov.data$Enterotype ,  permutations = 999, method = "bray")

bdiversity1 <- ggplot(pcoa_analysis1,aes(V1, V2, color = Enterotype))+
  geom_point(size = 2)+#ggtitle("France")+
  stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)+
  xlab(paste0("PCoA1(",pc1[1],"%)"))+ylab(paste0("PCoA2(",pc1[2],"%)")) +
  scale_color_manual(values =c(col3))+
  #guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.title = element_text(hjust=0.5),
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position =c(.87, .15),
        legend.title = element_text(size =12,family ="bold",color="black"),
        legend.text = element_text(size =9,family ="sans",color="black"),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA,color="black")
  )+
  annotate("text",x=-0.45,y=0.28,label = paste("PERMANOVA:\np=",otu.adonis1[1,5],
                                               "\nR2 = ", round(otu.adonis1[1,3],4),sep = ""
  ))+ggtitle("VF")
bdiversity1


ARG=ARG[rowSums(ARG>0)>(ncol(ARG)*0.1),]
ARG=as.data.frame(apply(ARG,2,function (x){
  x=x/sum(x)
  return(x)
}))
ARG <- as.data.frame(t(ARG))

beta_diversity=vegdist(ARG,method = "bray")
pcoa_analysis2=cmdscale(beta_diversity,k=4,eig=TRUE)
pc2<- round(pcoa_analysis2$eig/sum(pcoa_analysis2$eig)*100,digits=2)#这一步是算解释度的，pc的第一和第二位分别代表X轴和y轴的解释度
pcoa_analysis2 <- as.data.frame(pcoa_analysis2[["points"]])

pcoa_analysis2 <- merge(pcoa_analysis2,cluster,by="row.names",all=F)
pcoa_analysis2 $Enterotype <- ifelse(pcoa_analysis2$Cluster=="4","ETE","ETnonE")
pcoa_analysis2 $Enterotype<- factor(pcoa_analysis2 $Enterotype,levels = c("ETnonE","ETE"))

cov.data <- pcoa_analysis2[match(row.names(ARG),pcoa_analysis2$Row.names),]
cov.data$Enterotype <- as.factor(cov.data$Enterotype)
otu.adonis2=adonis2(ARG ~ cov.data$Enterotype ,  permutations = 999, method = "bray")

bdiversity2 <- ggplot(pcoa_analysis2,aes(V1, V2, color = Enterotype))+
  geom_point(size = 2)+#ggtitle("France")+
  stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)+
  xlab(paste0("PCoA1(",pc2[1],"%)"))+ylab(paste0("PCoA2(",pc2[2],"%)")) +
  scale_color_manual(values =c(col3))+
  #guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.title = element_text(hjust=0.5),
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position =c(.87, .15),
        legend.title = element_text(size =12,family ="bold",color="black"),
        legend.text = element_text(size =9,family ="sans",color="black"),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA,color="black")
  )+
  annotate("text",x=-0.55,y=0.28,label = paste("PERMANOVA:\np=",otu.adonis2[1,5],
                                               "\nR2 = ", round(otu.adonis2[1,3],4),sep = ""
  ))+ggtitle("ARG")
bdiversity2

bdiversity_all <- ggarrange(bdiversity1,bdiversity2, ncol = 1, nrow = 2)
ggsave("./Fig2E_VF_ARG_PCA.pdf",device = cairo_pdf,width =5, height =10,bdiversity_all)
save(pcoa_analysis1,pcoa_analysis2,otu.adonis1,otu.adonis2,pc1,pc2,file="./Fig2E.RData")


#========Fig2S_VF_ARG_score_PCA===============
cov.data <- pcoa_analysis1[match(row.names(VF),pcoa_analysis1$Row.names),]
cov.data$Enterotype <- as.factor(cov.data$Enterotype)
otu.adonis1=adonis2(VF ~ cov.data$balance_value,permutations = 999, method = "bray")

p1 <- ggplot(pcoa_analysis1,aes(V1, V2, color = balance_value))+
  geom_point(size = 2)+#ggtitle("France")+
  #stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)
  xlab(paste0("PCoA1(",pc1[1],"%)"))+ylab(paste0("PCoA2(",pc1[2],"%)"))+
  labs(color="balance_score") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.6,9.6))+
  theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.title = element_text(hjust=0.5),
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
  annotate("text",x=-0.15,y=-0.4,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis1[1,5],2),
                                             "\nR2=",signif(otu.adonis1[1,3],2)))+
  ggtitle("VF")

cov.data <- pcoa_analysis2[match(row.names(ARG),pcoa_analysis2$Row.names),]
cov.data$Enterotype <- as.factor(cov.data$Enterotype)
otu.adonis2=adonis2(ARG ~ cov.data$balance_value,  permutations = 999, method = "bray")
range(cov.data$balance_value)
p2 <- ggplot(pcoa_analysis2,aes(V1, V2, color = balance_value))+
  geom_point(size = 2)+#ggtitle("France")+
  #stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)
  xlab(paste0("PCoA1(",pc2[1],"%)"))+ylab(paste0("PCoA2(",pc2[2],"%)"))+
  labs(color="balance_score") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.6,9.6))+
  theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.title = element_text(hjust=0.5),
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
  annotate("text",x=0.5,y=-0.4,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis2[1,5],2),
                                           "\nR2=",signif(otu.adonis2[1,3],2)))+
  ggtitle("ARG")


pdf("./Fig2S_VF_ARG_score_PCA250708.pdf",width = 10.5,height = 5.5)
p1+p2+plot_layout(heights = c(1,1),widths = c(1,1),ncol = 2,nrow = 1)
dev.off()
save(pcoa_analysis1,otu.adonis1,pcoa_analysis2,otu.adonis2,pc1,pc2,file="./Fig2S_VF_ARG_score_PCA.RData")


#=========Fig2F_VF_ARG_diff===========
#========预处理VF和ARG================
##挑样本
VF <- read.table("/groups/ProHuShiX/home/share/BIGMetaG_result/VF_2022/HCC+health/VF_merged_new.txt",header=T)
#colnames(VF) <- gsub("\\.","-",colnames(VF))
score <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt",header=T)
#score$sample <- gsub("\\.","-",score$Row.names)
VF <- VF[,colnames(VF)%in%score$Row.names]

ARG <- read.table("/groups/ProHuShiX/home/share/BIGMetaG_result/ARGs/ARGs_merged.txt",header=T,check.names = F,row.names = 1)
colnames(ARG) <- gsub("-",".",colnames(ARG))
ARG <- ARG[,colnames(ARG)%in%score$Row.names]

##剔除检出率
VF=VF[rowSums(VF>0)>(ncol(VF)*0.1),]
VF=as.data.frame(apply(VF,2,function (x){
  x=x/sum(x)
  return(x)
}))
VF <- as.data.frame(t(VF))
save(VF,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/trim_VF_new.RData")
VF_clr=transform_and_filter_taxa(VF,samples_row = T,method = "clr",missing_filter = 0)
save(VF_clr,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/VF_clr.RData")

ARG=ARG[rowSums(ARG>0)>(ncol(ARG)*0.1),]
ARG=as.data.frame(apply(ARG,2,function (x){
  x=x/sum(x)
  return(x)
}))
ARG <- as.data.frame(t(ARG))
save(ARG,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/trim_ARG_new.RData")
ARG_clr=transform_and_filter_taxa(ARG,samples_row = T,method = "clr",missing_filter = 0)
save(ARG_clr,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/ARG_clr.RData")

#=========VF===========
VF_clr <- VF_clr[match(score$Row.names,row.names(VF_clr)),]
colnames(score) <- gsub("Cluster","cluster",colnames(score))
identical(colnames(VF),colnames(VF_clr))

res <- list()
for (i in 1:ncol(VF_clr)){
  tmp.data <- as.data.frame(cbind(VF_clr[,i],VF[,i],score$cluster))
  colnames(tmp.data) <- c("value","value_raw","cluster")
  tmp.data$cluster <- ifelse(tmp.data$cluster==4,"ETE","ETnonE")
  p.value<- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  a <- tmp.data %>% filter(cluster == "ETE") %>% select(value_raw)
  b <- tmp.data%>% filter(cluster == "ETnonE") %>% select(value_raw)
  # a[a==0] <- NA
  #b[b==0] <- NA
  fold <- mean(na.omit(a[,1]))/mean(na.omit(b[,1]))
  res[[i]] <- data.frame(VF_clr=colnames(VF_clr)[i],FC=fold,pvalue=p.value,
                         mean_ETE=mean(a[,1]),mean_ETnonE=mean(b[,1]),mean_all=mean(tmp.data$value_raw))
}
res <- do.call("rbind",res)
res$padj <- p.adjust(res$pvalue,"BH")
res <- res[res$pvalue<0.05,]
res <- res[res$FC>1,]
res$ETE_fc <- res$mean_ETE/res$mean_all
res$ETnonE_fc <- res$mean_ETnonE/res$mean_all
res$ETE_log2fc <- log2(res$ETE_fc)
res$ETnonE_log2fc <- log2(res$ETnonE_fc)

library(readxl)
VF_info <- read_excel("/groups/ProHuShiX/home/share/VFARDB/VFDB/VFDB_2022_header.xlsx")

res <- merge(res,VF_info[,c("VFG","Des","Description","VF_class","VFC_class")],by.x="VF_clr",by.y="VFG",all.x=T)
forplot <- res%>%arrange(-FC)

# 注意这里导出的VF表是全部所有的VF，不只是挑出来用于画图的部分（不筛ｐ值）
write_xlsx(forplot,"/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/select_VF250616.xlsx")

# forplot <- forplot[1:10,]
forplot <- forplot[!forplot$VFC_class%in%c("Immune modulation","Adherence"),]
table(forplot$VFC_class)

#========批量画VF图=====================
VFC_class <- c("Nutritional/Metabolic factor","Effector delivery system","Invasion")

plot_function <- function(i){
  tmp.forplot <- forplot[forplot$VFC_class==VFC_class[i],]%>%arrange(-FC)
  colnames(tmp.forplot) <- gsub("VF_clr","VF",colnames(tmp.forplot))
  
  forplot1 <- tmp.forplot[,c("VF","FC","ETE_log2fc","Des")]
  colnames(forplot1) <- c("VF","FC","log2fc","Des")
  forplot1$Enterotype <- "ETE"
  forplot2 <- tmp.forplot[,c("VF","FC","ETnonE_log2fc","Des")]
  colnames(forplot2) <- c("VF","FC","log2fc","Des")
  forplot2$Enterotype <- "ETnonE"
  
  fordot1 <- rbind(forplot1,forplot2)
  fordot1$Des <- factor(fordot1$Des,levels = tmp.forplot$Des)
  fordot1$Enterotype <- factor(fordot1$Enterotype,levels=c("ETE","ETnonE"))
  
  p <- ggplot(fordot1,aes(y=Des,x=Enterotype))+ 
    geom_point(aes(color=log2fc,size=5))+ ##设置框线周边的线为灰色，宽度为1
    # scale_fill_gradient2(low="#0E5B5B",mid = "#FFFFFF",high= "#862627")
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1,face="bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom")+labs(x="",y="")+coord_flip()+
    scale_color_gradientn(
      colors = c("#2739FF", "white","#d73027"),
      values = scales::rescale(c(-0.15, 0, 0.8)),  # 非对称设置：让-0.08~0之间用一半梯度，0~0.8用另一半
      limits = c(-0.08, 0.8),                     # 可根据数据自动设
      oob = scales::squish                        # 防止极值报错
    ) +ggtitle(VFC_class[i])
  return(p)
}

library(patchwork)
plist <- lapply(1:3, plot_function)
p1 <- plist[[1]]
p2 <- plist[[2]]
p3 <- plist[[3]]

VF_p <- p1+p2+p3+plot_layout(heights = c(1,1,1),widths = c(15,6,4),ncol = 3,nrow = 1)
VF_p

#====VF correlation with BS====
load("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/442MGS_clinical.RData")
identical(MGS_meta$sample,rownames(VF_clr))

MGS_meta$Age <- as.numeric(MGS_meta$Age);MGS_meta$Gender<-as.factor(MGS_meta$Gender);MGS_meta$BMI<-as.numeric(MGS_meta$BMI);MGS_meta$Smoke<-as.factor(MGS_meta$Smoke);MGS_meta$Drink<-as.factor(MGS_meta$Drink)
res1 <- list()
i=1
for (i in 1:ncol(VF_clr)){
  tmp.data <- data.frame(value=VF_clr[,i],value_raw=VF[,i],
                         cluster=score$cluster,BS=score$balance_value,
                         age=MGS_meta$Age,gender=MGS_meta$Gender,BMI=MGS_meta$BMI,
                         smoke=MGS_meta$Drink,drink=MGS_meta$Smoke)
  tmp.data$BS <- as.numeric(tmp.data$BS)
  tmp.data$cluster <- ifelse(tmp.data$cluster==4,"ETE","ETnonE")
  p.value1<- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  p.value2<- summary(lm(value~BS+gender+age+BMI+smoke+drink, data=tmp.data))$coefficients["BS", "Pr(>|t|)"]
  estimate<- summary(lm(value~BS+gender+age+BMI+smoke+drink, data=tmp.data))$coefficients["BS", "Estimate"]
  a <- tmp.data %>% filter(cluster == "ETE") %>% select(value_raw)
  b <- tmp.data%>% filter(cluster == "ETnonE") %>% select(value_raw)
  fold <- mean(na.omit(a[,1]))/mean(na.omit(b[,1]))
  res1[[i]] <- data.frame(VF_clr=colnames(VF_clr)[i],FC=fold,pvalue_cluster=p.value1,estimate=estimate,pvalue_BS=p.value2,
                          mean_ETE=mean(a[,1]),mean_ETnonE=mean(b[,1]),mean_all=mean(tmp.data$value_raw))
}
res1 <- do.call("rbind",res1)
res1$padj_BS <- p.adjust(res1$pvalue_BS,"BH")
res1$padj_cluster <- p.adjust(res1$pvalue_cluster,"BH")

writexl::write_xlsx(res1,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/FigS2F_VF_BS.xlsx")

# 看一下挑来画图的VF是否都与BS负相关
res1_f <- res1[res1$VF_clr%in%forplot$VF_clr,]
writexl::write_xlsx(res1_f,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/select_VF_BS250709.xlsx")

#========挑ARG================
ARG <- ARG[match(score$Row.names,row.names(ARG)),]
ARG_clr <- ARG_clr[match(score$Row.names,row.names(ARG_clr)),]
identical(row.names(ARG),row.names(ARG_clr))

res <- list()
i=1
for (i in 1:ncol(ARG)){
  tmp.data <- as.data.frame(cbind(ARG_clr[,i],ARG[,i],score$cluster))
  colnames(tmp.data) <- c("value","value_raw","cluster")
  tmp.data$cluster <- ifelse(tmp.data$cluster==4,"ETE","ETnonE")
  p.value<- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  a <- tmp.data %>% filter(cluster == "ETE") %>% select(value_raw)
  b <- tmp.data%>% filter(cluster == "ETnonE") %>% select(value_raw)
  # a[a==0] <- NA
  #b[b==0] <- NA
  fold <- mean(na.omit(a[,1]))/mean(na.omit(b[,1]))
  res[[i]] <- data.frame(ARG=colnames(ARG)[i],FC=fold,pvalue=p.value,
                         mean_ETE=mean(a[,1]),mean_ETnonE=mean(b[,1]),mean_all=mean(tmp.data$value_raw))
}
res <- do.call("rbind",res)
res$padj <- p.adjust(res$pvalue,"BH")
res <- res[res$pvalue<0.05,]
res <- res[res$FC>1,]
res$ETE_fc <- res$mean_ETE/res$mean_all
res$ETnonE_fc <- res$mean_ETnonE/res$mean_all
res$ETE_log2fc <- log2(res$ETE_fc)
res$ETnonE_log2fc <- log2(res$ETnonE_fc)


library(readxl)
# 只跑一遍即可
ARG_name <- as.data.frame(colnames(ARG))
colnames(ARG_name) <- c("name")
ARG_name$num <- paste("ARG",1:nrow(ARG_name),sep = "")
colnames(ARG) <- paste("ARG",1:ncol(ARG),sep = "")
ARG_name$Des <- gsub("^.*\\|","",ARG_name$name)
ARG_name$Des <- gsub("_.*$","",ARG_name$Des)
ARG_name$Species <- NA

res <- merge(res,ARG_name[,c("name","num","Des")],by.x="ARG",by.y="num",all.x=T)
forplot <- res%>%arrange(-FC)

# 注意这里导出的ARG也是全部的ARG（不筛ｐ值）
# library(writexl)
# write_xlsx(forplot,"/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/select_ARG_all250709.xlsx")

library(readxl)
forplot <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/select_ARG_plot250709.xlsx",sheet = "Sheet1")
table(forplot$`Resistance Mechanism`)
# forplot <- forplot[1:10,]
forplot <- forplot %>%
  arrange(desc(FC)) %>%                                  # 按 value 从大到小排序
  slice(1:30) # 挑前30个
forplot <- forplot[!forplot$`Resistance Mechanism`%in%c("antibiotic target protection",
                                                        "antibiotic target replacement","reduced permeability to antibiotic, antibiotic efflux"),]

table(forplot$`Resistance Mechanism`)

forplot$ETE_fc <- forplot$mean_ETE/forplot$mean_all
forplot$ETnonE_fc <- forplot$mean_ETnonE/forplot$mean_all
forplot$ETE_log2fc <- log2(forplot$ETE_fc)
forplot$ETnonE_log2fc <- log2(forplot$ETnonE_fc)

#========批量画ARG图=====================
ARG_class <- c("antibiotic efflux","antibiotic inactivation","antibiotic target alteration")
i=1

min(forplot$ETE_log2fc,forplot$ETnonE_log2fc);max(forplot$ETE_log2fc,forplot$ETnonE_log2fc)

plot_function <- function(i){
  tmp.forplot <- forplot[forplot$`Resistance Mechanism`==ARG_class[i],]%>%arrange(-FC)
  colnames(tmp.forplot) <- gsub("Resistance Mechanism","ARG_class",colnames(tmp.forplot))
  
  forplot1 <- tmp.forplot[,c("ARG","FC","ETE_log2fc","Des")]
  colnames(forplot1) <- c("ARG","FC","log2fc","Des")
  forplot1$Enterotype <- "ETE"
  forplot2 <- tmp.forplot[,c("ARG","FC","ETnonE_log2fc","Des")]
  colnames(forplot2) <- c("ARG","FC","log2fc","Des")
  forplot2$Enterotype <- "ETnonE"
  
  fordot2 <- rbind(forplot1,forplot2)
  fordot2$Des <- factor(fordot2$Des,levels = tmp.forplot$Des)
  fordot2$Enterotype <- factor(fordot2$Enterotype,levels=c("ETE","ETnonE"))
  
  p <- ggplot(fordot2,aes(y=Des,x=Enterotype))+ 
    geom_point(aes(color=log2fc,size=5))+ ##设置框线周边的线为灰色，宽度为1
    # scale_fill_gradient2(low="#0E5B5B",mid = "#FFFFFF",high= "#862627")
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1,face="bold"),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom")+labs(x="",y="")+coord_flip()+
    scale_color_gradientn(
      colors = c("#2739FF", "white","#d73027"),
      values = scales::rescale(c(-0.3, 0, 1.2)),  # 非对称设置：让-0.08~0之间用一半梯度，0~0.8用另一半
      limits = c(-0.3,1.2),                     # 可根据数据自动设
      oob = scales::squish                        # 防止极值报错
    ) +
    ggtitle(ARG_class[i])
  return(p)
}


library(patchwork)
plist <- lapply(1:3, plot_function)
p1 <- plist[[1]]
p2 <- plist[[2]]
p3 <- plist[[3]]

ARG_p <- p1+p2+p3+plot_layout(heights = c(1,1,1),widths = c(16,6,3),ncol = 3,nrow = 1)
ARG_p



pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/diff_VF_ARG.pdf",width=12,height = 3)
VF_p
ARG_p
dev.off()

pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/diff_ARG250709.pdf",width=12,height = 3)
ARG_p
dev.off()

#=====ARG correlation with BS====
identical(MGS_meta$sample,rownames(ARG_clr))

MGS_meta$Age <- as.numeric(MGS_meta$Age);MGS_meta$Gender<-as.factor(MGS_meta$Gender);MGS_meta$BMI<-as.numeric(MGS_meta$BMI);MGS_meta$Smoke<-as.factor(MGS_meta$Smoke);MGS_meta$Drink<-as.factor(MGS_meta$Drink)
res1 <- list()
i=1
for (i in 1:ncol(ARG_clr)){
  tmp.data <- data.frame(value=ARG_clr[,i],value_raw=ARG[,i],
                         cluster=score$cluster,BS=score$balance_value,
                         age=MGS_meta$Age,gender=MGS_meta$Gender,BMI=MGS_meta$BMI,
                         smoke=MGS_meta$Drink,drink=MGS_meta$Smoke)
  tmp.data$BS <- as.numeric(tmp.data$BS)
  tmp.data$cluster <- ifelse(tmp.data$cluster==4,"ETE","ETnonE")
  p.value1<- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  p.value2<- summary(lm(value~BS+gender+age+BMI+smoke+drink, data=tmp.data))$coefficients["BS", "Pr(>|t|)"]
  estimate<- summary(lm(value~BS+gender+age+BMI+smoke+drink, data=tmp.data))$coefficients["BS", "Estimate"]
  a <- tmp.data %>% filter(cluster == "ETE") %>% select(value_raw)
  b <- tmp.data%>% filter(cluster == "ETnonE") %>% select(value_raw)
  fold <- mean(na.omit(a[,1]))/mean(na.omit(b[,1]))
  res1[[i]] <- data.frame(ARG=colnames(ARG)[i],FC=fold,pvalue_cluster=p.value1,estimate=estimate,pvalue_BS=p.value2,
                          mean_ETE=mean(a[,1]),mean_ETnonE=mean(b[,1]),mean_all=mean(tmp.data$value_raw))
}
res1 <- do.call("rbind",res1)
res1$padj_BS <- p.adjust(res1$pvalue_BS,"BH")
res1$padj_cluster <- p.adjust(res1$pvalue_cluster,"BH")
res1 <- merge(res1,ARG_name[,c("name","num","Des")],by.x="ARG",by.y="num",all.x=T)
writexl::write_xlsx(res1,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/FigS2F_ARG_BS.xlsx")

res1_f <- res1[res1$name%in%res$ARG,]

# select画图用
write_xlsx(res1_f,"/groups/ProHuShiX/home/xiashufen/bigmeta/Fig2/Fig2F/select_ARG_plot250709.xlsx")


