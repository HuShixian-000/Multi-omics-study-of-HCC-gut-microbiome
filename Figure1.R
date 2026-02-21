rm(list=ls())
library(microbiome)
library(crayon)
library(nlme)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggprism)
library(vegan)
library(ggplot2)
library(ggsci)
library(data.table)
library(foreach)
library(dplyr)
library(gghalves)
library(ggdist)
library(purrr)
library(tidyverse)
library(magrittr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

#=========Fig1A big Adonis==============
# 数据准备见R_adonis.R
suffix <- "_ad1.txt"
file_name <- list.files('/groups/ProHuShiX/home/xiashufen/bigmeta/big_adonis',pattern = suffix)
file_name2 <- map_chr(file_name,function(x){stringr::str_remove_all(x,suffix)})
all_files <- lapply(file_name,function(x){y <- read.table(x,header = T)})
names(all_files) <- file_name2
all_files <- do.call("rbind",all_files)


df <- all_files
df$P <- p.adjust(df$P, method="BH")

df$varExpPct <- sprintf("%.1f%%", 100*df$R2)
df$varExpPct[is.na(df$R2)] <- ""
df$NAtext <- ""
df$NAtext[is.na(df$R2)] <- "N/A"
df$stars <- cut(df$Pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
df$stars[is.na(df$R2)] <- ""

# Try to make a reasonable color scheme that has contrast where needed
colors <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100)
#R2Q <- quantile(df$R2, c(0.25, 0.75), na.rm=T)

#=====β分布中的颜色映射====
#这两个参数控制颜色映射分布（β分布）
#这里大多数R2分布在0-0.013，但sample都为1
#所以想要颜色分布的效果是：在低值段的颜色变化更加细腻
#beta取小、alpha取大即可
alpha <- 2
beta <- 0.08

# 自动拟合 or 手动拟合
if (is.na(alpha) || is.na(beta)) {
  labhat <- optim(par=c(0, 0), method="Nelder-Mead",
                  fn=function(lab) sum((pbeta(c(0.25, 0.75), exp(lab[1]), exp(lab[2])) - R2Q)^2))
  abhat <- exp(labhat$par)
  colorvalues <- pbeta(seq(0, 1, length=length(colors)), abhat[1], abhat[2])
  
  # Label colors flip to white when the color is too dark
  df$lblcolor <- ifelse(qbeta(df$R2, abhat[1], abhat[2]) < 0.8, "black", "white")
  #cat(sprintf("Best-fit alpha = %g, beta = %g\n", abhat[1], abhat[2]))
} else {
  colorvalues <- pbeta(seq(0, 1, length=length(colors)), alpha, beta)
  
  # Label colors flip to white when the color is too dark
  df$lblcolor <- ifelse(qbeta(df$R2, alpha, beta) < 0.8, "black", "white")
}

df$Dataset <- factor(df$Dataset,levels = c("MGS","MetaCyc","VF","ARG"))
df$Cov <- factor(df$Cov,levels = rev(c("sample" ,"group" , "Age", "Gender" , "BMI" , "Smoke" ,"Drink",
                                       "Antibiotics" , "HBsAg" ,"Fatty_liver", "Cirrhosis","ALT" , "AST",
                                       "AFP" , "BCLC")))  

fontsize=15
source("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/theme_nature.r")
p <- ggplot(data=df, aes(x=Dataset, y=Cov)) +
  geom_tile(aes(fill=R2),width=0.97,height=0.95) +
  geom_text(aes(label=varExpPct, color=lblcolor), size=fontsize/(20/5), nudge_y=-0.15) +
  geom_text(aes(label=NAtext), color="grey", size=fontsize/(20/5), nudge_y=-0.15) +
  geom_text(aes(label=stars, color=lblcolor), fontface="bold", size=1.25*fontsize/(14/5), nudge_y=0.12) +
  scale_fill_gradientn(colors=colors, values=colorvalues, limits=c(0, 1), na.value="white") +
  scale_color_manual(values=c(black="black", white="white")) +
  scale_x_discrete(expand=c(0,0)) + xlab(NULL) +
  scale_y_discrete(expand=c(0,0), position = "right", limits = rev(levels(df$Feature))) + ylab(NULL) +
  guides(color="none",
         fill=guide_colourbar(title=NULL, barheight=unit(40,"mm"), 
                              label.position = "left")) +
  theme_nature() +
  theme(axis.text.x = element_text(angle=-17, hjust=0),
        panel.border=element_rect(fill=NA),
        legend.position = "left", axis.ticks.y = element_blank())

pdf("Fig1b.adonis250901.pdf",width = 5,height = 7)
print(p)
dev.off()



#======colorset=============
col1 <- c("#2482BC","#f99170","#e73133")
col2 <- c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
col3 <- c("#fcd364","#fe9014")
col4 <- c("#2482BC","#FFCC33","#f99170","#e73133")

data <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/species_raw_filter.txt")
clinical <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/clinical_filter.txt")

# 90% prevelance + relative abundance
data <- data[rowSums(data!=0)>ncol(data)*0.1,]
data <- apply(data,2,function(x){
  y <- x/sum(x)
  return(y)
})
data <- as.data.frame(t(data))

#====α diversity====
shannon <- as.data.frame(vegan::diversity((data), index="shannon"))
shannon <- shannon %>% mutate(group=clinical$group[match(rownames(shannon),clinical$Row.names)])
shannon$group <- factor(shannon$group,levels=c("Health","HCC"))
colnames(shannon)[1] <- "shannon"
ggplot(shannon, aes(x=group, y=shannon, fill=group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=c("#31AEEF","#FF9A9C"))+
  xlab("")+ylim(0,5)+
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')+
  geom_signif(# 添加显著性标签
    comparisons=list(c("HCC","Health")),
    step_increase = 0.1,
    test="wilcox.test",   # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=F,   # 标签样式F为数字，T为*号
    size=1,textsize = 7,y_position = 4.5
  ) 

adiversity <- ggplot(shannon,aes(x= group, y= shannon,fill=group,color=group))+
  geom_half_violin(position=position_nudge(x=0.1,y=0), side='R',adjust=1.2,trim=F,color=NA,alpha=0.79)+ 
  geom_half_point(position=position_nudge(x=-0.35,y=0),size =1, shape =19,range_scale = 0.5,alpha=0.9)+
  geom_boxplot(outlier.shape = NA, #隐藏离群点； 
               width =0.1, alpha=0.9)+ scale_fill_manual(values = c(col1[1], col1[3]))+
  scale_colour_manual(values =c(col1[1], col1[3]))+
  stat_compare_means(comparisons = list(c('HCC','Health')),label.y =4.8,textsize=6,fontface='bold')+
  stat_summary(fun="median",geom= "crossbar", width = 0.09,color="white")+
  theme_bw()+theme(panel.grid=element_blank(),legend.position = "none",
                   axis.text = element_text(colour = "black",size=12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(colour = "black",size=13))
pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/FigS1A_a_diversity.pdf",width = 3.5,height = 5)
print(adiversity)
dev.off()

#====重复抽样10次====
set.seed(666)
resampled_data <- map_df(1:10, function(i) {
  sampled_HCC <- sample_n(shannon %>% filter(group == "HCC"), 54)
  sampled_Health <- shannon %>% filter(group == "Health") 
  sampled_data <- bind_rows(sampled_Health, sampled_HCC) 
  sampled_data$iteration <- i # 添加标记，表示第几次抽样
  return(sampled_data)
})
resample_adiversity <- ggplot(resampled_data,aes(x= group, y= shannon,fill=group,color=group))+
  geom_half_violin(position=position_nudge(x=0.1,y=0), side='R',adjust=1.2,trim=F,color=NA,alpha=0.79)+ 
  geom_half_point(position=position_nudge(x=-0.35,y=0),size =1, shape =19,range_scale = 0.5,alpha=0.9)+
  geom_boxplot(outlier.shape = NA, #隐藏离群点； 
               width =0.1, alpha=0.9)+ scale_fill_manual(values = c(col1[1], col1[3]))+
  scale_colour_manual(values =c(col1[1], col1[3]))+
  stat_compare_means(comparisons = list(c('HCC','Health')),label.y =4.5,textsize=6,fontface='bold')+
  stat_summary(fun="median",geom= "crossbar", width = 0.09,color="white")+
  theme_bw()+theme(panel.grid=element_blank(),legend.position = "none",
                   axis.text = element_text(colour = "black",size=12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(colour = "black",size=13))+
  facet_wrap(~iteration,ncol=5) # 按照抽样批次分别绘制
pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/resample_a_diversity.pdf",width = 10,height = 4)
print(resample_adiversity)
dev.off()


#====β diversity====
library(ape)
library(multcomp)
library(extrafont)
library(patchwork)
library(vegan)
library(ggtern)
library(reshape2)
library(ggsci)
library(scales)
library(patchwork)

set.seed(666)
dist_BC <- vegdist(data,method = "bray")
pcoa <- pcoa(dist_BC,correction = "none", rn = NULL) 
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]

groups<-clinical$group[match(rownames(pcoa$vectors),clinical$Row.names)]
pcoadata <- data.frame(rownames(pcoa$vectors),PC1,PC2,groups)
colnames(pcoadata) <- c("sample","PC1","PC2","group")
pcoadata$group <- factor(pcoadata$group,levels =c("Health","HCC"))

dist_BC <- as.matrix(dist_BC)
#Adonis test
otu.adonis <- adonis2(dist_BC~group,data = pcoadata,permutations = 999)
#write.table(otu.adonis$aov.tab,'beta_bray-curtis_adonis.tsv',sep = '\t',quote = F)

otu_pvalue <- as.numeric(otu.adonis$`Pr(>F)`[1])
otu_r2 <- as.numeric(otu.adonis$R2[1])

bdiversity <- ggplot(pcoadata, aes(PC1, PC2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(PC1, PC2), level = 0.95, linetype = 2) +
  xlim(-0.6, 0.6) + ylim(-0.5, 0.5) +
  labs(
    x = paste0("PCoA1 ( ", floor(pcoa$values$Relative_eig[1] * 100), "% )"),
    y = paste0("PCoA2 ( ", floor(pcoa$values$Relative_eig[2] * 100), "% )"),
    color = "Group"
  ) +
  scale_color_manual(values = c(col1[1], col1[3])) +
  theme(
    plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.position = c(.15, .15),  # 修改 ggplot2 3.5.0 兼容性问题
    legend.title = element_text(size = 12, family = "bold", color = "black"),
    legend.text = element_text(size = 9, family = "sans", color = "black"),
    legend.key = element_blank(), 
    axis.text.x = element_text(size = 10, family = "sans", color = "black"),
    axis.text.y = element_text(size = 10, family = "sans", color = "black"),
    panel.grid.major = element_line(colour = NA),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill = NA)
  ) +
  annotate("text", x = -0.49, y = 0.45, 
           label = paste("PERMANOVA:\np-value =", signif(otu_pvalue, 2),
                         "\nR² =", signif(otu_r2, 2)))

bdiversity <- ggExtra::ggMarginal(bdiversity, type = "histogram", groupColour = F, groupFill = T,
                                  xparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'),
                                  yparams = list(bins = 60, alpha = 0.8,position = 'identity', color = 'white'))

bdiversity <- ggExtra::ggMarginal(bdiversity, type = "density", groupColour = F, groupFill = T,
                                  xparams = list(alpha = 0.6,color=NA),
                                  yparams = list(alpha = 0.6,color=NA))


ggsave("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/FigS1B_betadiversity0.pdf",device = cairo_pdf,width =5.5, height =5,bdiversity)

rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/Enterotype/")

##====Enterotype====
#======colorset=============
col1 <- c("#2482BC","#f99170","#e73133")
col2 <- c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
col3 <- c("#fcd364","#b66d32")
col4 <- c("#b39658","#b66d32")

source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/functions.R")

# MGS
data_taxa=read.table("../Final_input/species_raw_filter.txt",sep = " ",header = T,row.names = 1)
data_phenotype=read.table("../Final_input/clinical_filter.txt",sep = " ",header = T,stringsAsFactors = F,row.names = 1)
data_phenotype$Cluster=NULL;data_phenotype$alpha_diversity=NULL;data_phenotype$balance_value=NULL;data_phenotype$new_group=NULL

data_taxa=data_taxa[rowSums(data_taxa>0)>(ncol(data_taxa)*0.1),]
data_taxa=as.data.frame(apply(data_taxa,2,function (x){
  x=x/sum(x)
  return(x)
}))
data_taxa=as.data.frame(t(data_taxa))
#write.table(data_taxa,"/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/species_0.1prevelance.txt")
data_taxa_f=data_taxa[,colMeans(data_taxa)>0.001] 
data_taxa_clr=transform_and_filter_taxa(data_taxa,samples_row = T,method = "clr",missing_filter = 0)
#write.table(data_taxa_clr,"/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/data_mgs_clr.txt")

# shannon
data_shannon <- as.data.frame(vegan::diversity((data_taxa), index="shannon"))
colnames(data_shannon)[1] <- "shannon"
data_shannon=merge(data_shannon, data_phenotype,by.x="row.names",by.y="Row.names",all=F)


#rxx
data_taxa=read.table("/groups/ProHuShiX/home/xiashufen/for_shufen/rxx_raw_species.txt",sep = "\t",header = T,row.names = 1,check.names = F)
data_taxa_f_rxx <- data_taxa[rownames(data_taxa)%in%colnames(data_taxa_f),]
HCC114_Health37 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
data_taxa_f_rxx <- data_taxa_f_rxx %>% 
  select(HCC114_Health37$sample)
data_taxa_f_rxx=as.data.frame(apply(data_taxa_f_rxx,2,function (x){
  x=x/sum(x)
  return(x)
}))
data_taxa_f_rxx=as.data.frame(t(data_taxa_f_rxx))
#data_taxa_f <- data_taxa_f_rxx

# enterotype
#====DMM====
set.seed(12345)
min(data_taxa_f[data_taxa_f>0]) # check 10^n
data_taxa_reads=as.data.frame((data_taxa_f*10000000))

fit1 <- mclapply(1:6, dmn, count = as.matrix(data_taxa_reads), verbose=TRUE,mc.cores=1)
lplc <- sapply(fit1, laplace)
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")

best <- fit1[[which.min(lplc)]]
cluster <- data.frame(Cluster=apply(mixture(best), 1, which.max))
k=4
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("Species", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    arrange(value) %>%
    mutate(Species = factor(Species, levels = unique(Species))) %>%
    filter(abs(value) > quantile(abs(value), 0.8))
  
  p <- ggplot(d, aes(x = Species, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers1: community type", k))
  print(p)
  #ggsave(paste(k,"lplc.core.species1.pdf",sep="_"))
}
heatmapdmn(as.matrix(data_taxa_reads), fit[[1]], best, 8)

k=1
for (k in seq(ncol(fitted(best)))) {
  # 转长格式并命名
  d <- reshape2::melt(fitted(best))
  colnames(d) <- c("Species", "cluster", "value")
  
  # 选定当前 cluster 的数据
  d <- subset(d, cluster == k) %>%
    arrange(desc(value)) %>%                                  # 按 value 从大到小排序
    slice(1:20) %>%                                           # 取前20个物种
    mutate(Species = factor(Species, levels = rev(Species))) # 保持排序顺序画图
  
  # 画图
  p <- ggplot(d, aes(x = Species, y = value)) +
    geom_bar(stat = "identity", fill = col2[k]) +             # 使用 col2[k] 作为填充色
    coord_flip() +
    labs(
      title = paste("Top 20 drivers: Community Type", k),
      y = "value",
      x = NULL
    ) +
    theme_minimal(base_size = 14) +                           # 简洁背景+大字号
    theme(
      axis.text.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  print(p)
  #ggsave(paste0(k, "_lplc.core.species250630.pdf"), plot = p, width = 6, height = 5)
}

# check associations
tmp_pheno=merge(data_shannon,cluster,by.x="Row.names",by.y="row.names",all=F)
tmp_pheno$Cluster=as.factor(tmp_pheno$Cluster)
tmp_pheno$height=tmp_pheno$height/100
tmp_pheno$weight=tmp_pheno$weight/(tmp_pheno$height^2)
tmp_pheno$weight[tmp_pheno$weight<0 | tmp_pheno$weight>35]=NA

ggplot(tmp_pheno) +
  aes(x = group, fill = factor(Cluster)) +
  geom_bar(position = "fill")
table(tmp_pheno$group,tmp_pheno$Cluster)
ggboxplot(tmp_pheno, x = "Cluster", y = "shannon",
          color = "Cluster", palette = "jco")

# set top selected bugs
#tmp_pheno=merge(data_phenotype,cluster,by.x="Row.names",by.y="row.names",all=F)
tmp_pheno$new_group[tmp_pheno$Cluster==1]=1
tmp_pheno$new_group[tmp_pheno$Cluster==2]=1
tmp_pheno$new_group[tmp_pheno$Cluster==3]=1
tmp_pheno$new_group[tmp_pheno$Cluster==4]=0
tmp.data=merge(tmp_pheno[,c("Row.names","new_group")],data_taxa_clr,by.x="Row.names",by.y="row.names",all=F)
tmp.data$sample=NULL
tmp.data$Row.names=NULL
tmp.data[[1]] <- as.numeric(tmp.data[[1]])

#set.seed(666)
#n=50
f_data <- mRMR.data(data = data.frame(tmp.data))
results <- mRMR.classic("mRMRe.Filter", data = f_data, target_indices = 1,
                        feature_count = 50) # select top n predictors
tmp.taxa=colnames(tmp.data)[(solutions(results)$`1`)]
mm=as.data.frame((results@mi_matrix[1,]))
# colnames(mm)[1] <- "MI"
#write.table(mm,"/groups/ProHuShiX/home/xiashufen/bigmeta/Enterotype/mrmr_MI.txt")

taxa_up=tmp.taxa[tmp.taxa %in% rownames(mm)[mm$`(results@mi_matrix[1, ])`>0]]
taxa_down=tmp.taxa[tmp.taxa %in% rownames(mm)[mm$`(results@mi_matrix[1, ])`<0]]
taxa_up
taxa_down
taxa_up_index <- match(taxa_up, colnames(data_taxa))
taxa_down_index <- match(taxa_down, colnames(data_taxa))

write.table(taxa_up,file = "/groups/ProHuShiX/home/xiashufen/bigmeta/BS/Score.feature.up1.txt",quote = F)
write.table(taxa_down,file = "/groups/ProHuShiX/home/xiashufen/bigmeta/BS/Score.feature.down1.txt",quote = F)

# Compute balance
balance_df=data.frame(sample=rownames(data_taxa),balance_value=NA)
data_taxa[data_taxa==0]=min(data_taxa[data_taxa>0])
for(sample in rownames(data_taxa)){
  balance_df$balance_value[balance_df$sample==sample] = log(exp(mean(as.numeric(log(data_taxa[sample,taxa_up_index])),na.rm = T))) - log(exp(mean(as.numeric(log(data_taxa[sample,taxa_down_index])),na.rm = T)))  
}
rownames(balance_df)=balance_df$sample
balance_df$sample=NULL

# link to phenotype
tmp_pheno=merge(tmp_pheno,balance_df,by="Row.names",by.y="row.names",all=F)
tmp_pheno$Cluster=as.factor(tmp_pheno$Cluster)
ggboxplot(tmp_pheno, x = "Cluster", y = "shannon",
          color = "Cluster", palette = "jco")
# BS明显比shannon更能区分4种肠型
ggboxplot(tmp_pheno, x = "Cluster", y = "balance_value",
          color = "Cluster", palette = "jco")
ggboxplot(tmp_pheno, x = "BCLC", y = "balance_value",
          color = "BCLC", palette = "jco")
ggboxplot(tmp_pheno, x = "new_group", y = "balance_value",
          color = "new_group", palette = "jco")
tmp_comparisons <- list( c("HCC","Health"))
ggboxplot(tmp_pheno, x = "group", y = "balance_value",
          color = "group", palette = "jco")+
  stat_compare_means(comparisons = tmp_comparisons)
ggsave("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/BS_HCC_Health.pdf",width = 5,height = 6)
write.table(tmp_pheno,file="../Final_output/Balance.score1.txt",sep = "\t",quote = F)

#====enterotype&phenotype
data_phnotype=read.table("../Final_output//Balance.score1.txt",sep = "\t",stringsAsFactors = F,header = T)
data_survive=read.table("../Final_input/MGSinfo_0326.txt",sep=" ",header = T,stringsAsFactors = F)
data_survive$OS=as.numeric(ifelse(data_survive$OS=="生存",1,0))
data_survive$RFS0=as.numeric(data_survive$RFS0)
data_survive$RFSTime0=as.numeric(data_survive$RFSTime0)

data_phnotype=merge(data_phnotype,data_survive,by.x="Row.names",by.y="sample",all=F)
data_phnotype$Cluster=as.factor(data_phnotype$Cluster)

ggplot(data_phnotype[!is.na(data_phnotype$MVI),]) +
  aes(x = Cluster, fill = factor(MVI)) +
  geom_bar(position = "fill")
table(data_phnotype$MVI,data_phnotype$Cluster)

tmp_comparisons <- list( c("0", "1") )
ggboxplot(data_phnotype, x = "MVI", y = "balance_value",
          color = "MVI")+
  stat_compare_means(comparisons = tmp_comparisons)
summary(lm(Cluster~smoking,data_phnotype))

data_phnotype <- data_phnotype[!is.na(data_phnotype$Antibiotics),]
data_phnotype$Antibiotics <- as.factor(data_phnotype$Antibiotics)
summary(lm(Cluster~Antibiotics,data_phnotype)) # 0.521
p <- chisq.test(table(data_phnotype$Antibiotics,data_phnotype$Cluster))$p.value # 0.516

ggplot(data_phnotype,aes(x = Antibiotics, fill = Cluster)) +
  geom_bar(position = "fill")+scale_fill_manual(values =  col2)+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "none")

summary(lm(Cluster~drinking,data_phnotype))
summary(lm(Cluster~Gender,data_phnotype))
summary(lm(Cluster~Age,data_phnotype))

data_phnotype <- na.omit(data_phnotype)
data_phnotype$OSTime
cox <- coxph(Surv(OSTime,OS) ~ MVI + age + male + height+ weight+antibiotics+drinking + smoking, data = data_phnotype)   
summary(cox)

# continues analysis
model <- coxph(Surv(OSTime,OS) ~ Cluster + age + male + height+ weight+antibiotics+drinking + smoking, data=data_phnotype, x=TRUE)
model
model <- coxph(Surv(OSTime,OS) ~ balance_value + age + male + height+ weight+antibiotics+drinking + smoking, data=data_phnotype, x=TRUE)
model

# plot survival area
plot_surv_area(time="OSTime",
               status="OS",
               variable="balance_value",
               data=data_phnotype,start_color = "#BB5566",end_color="#44AA99",
               model=model,discrete=TRUE,bins = 5)

#====check PAM====
data_taxa_t=as.data.frame(t(data_taxa_f))
# 3种dist都试一下
data_taxa.dist=dist.JSD(data_taxa_t) # 行为样本
data_taxa.dist_j=vegdist(data_taxa_f,method = "jaccard") # 列为样本
data_taxa.dist_b=vegdist(data_taxa_f,method = "bray")

data_taxa.dist <- data_taxa.dist_j

# perform clustering
nclusters=NULL
# > calculate CH index to identify optimal number of clusters
for (k in 1:10) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    data_taxa.cluster_temp=pam.clustering(data_taxa.dist, k)
    nclusters[k]=index.G1(data_taxa_f,data_taxa.cluster_temp,  d = data_taxa.dist, # 列为样本
                          centrotypes = "medoids")# "medoids" PAM/"centroids" K means
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

# 最佳k为4
data_taxa.cluster=pam.clustering(data_taxa.dist, k=4)
cluster=data.frame(row.names = colnames(data_taxa_t),Cluster=data_taxa.cluster)
tmp_pheno=merge(data_shannon,cluster,by.x="Row.names",by.y="row.names",all=F)
#write.table(tmp_pheno,file="../Final_output/Balance.score.pam_jaccard.txt",sep = "\t",quote = F)

# top drivers
# 3型为E.coli
# set top selected bugs
tmp_pheno$new_group[tmp_pheno$Cluster==1]=1
tmp_pheno$new_group[tmp_pheno$Cluster==2]=1
tmp_pheno$new_group[tmp_pheno$Cluster==3]=0
tmp_pheno$new_group[tmp_pheno$Cluster==4]=1
tmp.data=merge(tmp_pheno[,c("Row.names","new_group")],data_taxa_clr,by.x="Row.names",by.y="row.names",all=F)
tmp.data$sample=NULL
tmp.data$Row.names=NULL
tmp.data[[1]] <- as.numeric(tmp.data[[1]])
f_data <- mRMR.data(data = data.frame(tmp.data))
results <- mRMR.classic("mRMRe.Filter", data = f_data, target_indices = 1,
                        feature_count = 50) # select top n predictors
tmp.taxa=colnames(tmp.data)[(solutions(results)$`1`)]
mm=as.data.frame((results@mi_matrix[1,]))

#n=50
taxa_up=tmp.taxa[tmp.taxa %in% rownames(mm)[mm$`(results@mi_matrix[1, ])`>0]] # 18
taxa_down=tmp.taxa[tmp.taxa %in% rownames(mm)[mm$`(results@mi_matrix[1, ])`<0]] # 32
taxa_up
taxa_down
taxa_up_index <- match(taxa_up, colnames(data_taxa))
taxa_down_index <- match(taxa_down, colnames(data_taxa))

write.table(taxa_up,file = "/groups/ProHuShiX/home/xiashufen/bigmeta/BS/PAM_jaccard.Score.feature.up.txt",quote = F)
write.table(taxa_down,file = "/groups/ProHuShiX/home/xiashufen/bigmeta/BS/PAM_jaccard.Score.feature.down.txt",quote = F)

# Compute balance
balance_df=data.frame(sample=rownames(data_taxa),balance_value=NA)
data_taxa[data_taxa==0]=min(data_taxa[data_taxa>0])
for(sample in rownames(data_taxa)){
  balance_df$balance_value[balance_df$sample==sample] = log(exp(mean(as.numeric(log(data_taxa[sample,taxa_up_index])),na.rm = T))) - log(exp(mean(as.numeric(log(data_taxa[sample,taxa_down_index])),na.rm = T)))  
}
rownames(balance_df)=balance_df$sample
balance_df$sample=NULL

tmp_pheno=merge(tmp_pheno,balance_df,by="Row.names",by.y="row.names",all=F)
tmp_pheno$Cluster=as.factor(tmp_pheno$Cluster)
tmp_comparisons <- list( c("HCC","Health"))
ggboxplot(tmp_pheno, x = "Cluster", y = "shannon",
          color = "Cluster", palette = "jco")
# BS明显比shannon更能区分4种肠型
ggboxplot(tmp_pheno, x = "Cluster", y = "balance_value",
          color = "Cluster", palette = "jco")
ggboxplot(tmp_pheno, x = "group", y = "balance_value",
          color = "group", palette = "jco")+
  stat_compare_means(comparisons = tmp_comparisons)
ggsave("/groups/ProHuShiX/home/xiashufen/bigmeta/Enterotype/pamscore.BS_HCC_Health.pdf",width = 5,height = 6)
write.table(tmp_pheno,file="../Final_output/pamjaccard.Balance.score.txt",sep = "\t",quote = F)

##====enterotpye diversity====
#======colorset=============
col1 <- c("#2482BC","#f99170","#e73133") ##健康人，肝硬化，HCC
col2 <- c("#5ba787","#fa9fcb","#fe9014","#6e7ca5") ##4种肠型
col3 <- c("#fcd364","#b35c31")  ##ETnonE,ETE
col4 <- c("#2482BC","#FFCC33","#f99170","#e73133")
# col4 <- c("#BD9E4B","#BF6C0F")  ##画小提琴图加粗的线的颜色
#=======Fig1C 4种肠型PCA图======================
data <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/species_raw_filter.txt",header = T)
#Cluster <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score.txt",header=T)
Cluster <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt",header=T)
# 90% prevelance + relative abundance
data <- data[rowSums(data!=0)>ncol(data)*0.1,]
data <- apply(data,2,function(x){
  y <- x/sum(x)
  return(y)
})
data <- as.data.frame(t(data))
dist_BC <- vegdist((data),method = "bray")
#dist_BC <- cmdscale(dist_BC,k=2,eig=TRUE)
pcoa <- pcoa(dist_BC,correction = "none", rn = NULL) 
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]
Enterotype<-Cluster$Cluster[match(rownames(pcoa$vectors),Cluster$Row.names)]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,Enterotype)
colnames(pcoadata) <- c("sample","PC1","PC2","Cluster")

pcoadata$group <- factor(pcoadata$Cluster,levels =c("1","2","3","4"))

#Adonis test
otu.adonis <- adonis2(dist_BC~Cluster,data = pcoadata,permutations = 999)
#write.table(otu.adonis$aov.tab,'beta_bray-curtis_adonis.tsv',sep = '\t',quote = F)

pcoadata$Cluster <- as.factor(pcoadata$Cluster)

p <- ggplot(pcoadata,aes(PC1, PC2, color = Cluster))+
  geom_point(size = 2)+#ggtitle("France")+
  stat_ellipse(aes(PC1,PC2) ,level=0.95,linetype = 2)+
  xlim(-0.6,0.6)+ylim(-0.5,0.5)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )"),
       color="Group") +
  scale_color_manual(values =col2,
                     labels=c("ETF(n=147)","ETR(n=130)","ETP(n=110)","ETE(n=55)"))+
  labs(color="Enterotype")+theme_classic()+
  #guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        legend.position =c(.14, .17),
        legend.title = element_text(size =12,face  ="bold",color="black"),
        legend.text = element_text(size =9,face ="bold",color="black"),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10,face ="bold",color="black"),
        axis.text.y = element_text(size = 10,face ="bold",color="black"))+
  annotate("text",x=-0.49,y=0.45,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis$`Pr(>F)`[1],2),
                                             "\nR2=",signif(otu.adonis$R2[1],2)))
p
ggsave("Fig1C_ent_betadiversity.pdf",device = cairo_pdf,width =5.5, height =5,p)
save(pcoadata,otu.adonis,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/Fig1C_ent_bdiv.RData")

#======FigS1 ent_alpha_diversity=============
shannon <- as.data.frame(vegan::diversity((data), index="shannon"))
shannon$Cluster <- Cluster$Cluster[match(rownames(shannon),Cluster$Row.names)]

shannon$Cluster <- factor(shannon$Cluster,levels = c("1","2","3","4"))
shannon$Enterotype <- ifelse(shannon$Cluster=="1","ETF",ifelse(shannon$Cluster=="2","ETR",
                                                               ifelse(shannon$Cluster=="3","ETP","ETE")))
shannon$Enterotype <- factor(shannon$Enterotype,levels = c("ETF","ETR","ETP","ETE"))
colnames(shannon)[1] <- "shannon"
p <- ggplot(shannon, aes(x=Enterotype, y=shannon, color=Enterotype)) + 
  geom_violin(trim=FALSE,fill="white")+
  geom_boxplot(width=0.15, fill="white")+
  #scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  scale_color_manual(values=col2)+
  xlab("")+labs(y="Shannon index")+
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'none',
        axis.title.y = element_text(size=15))+
  geom_signif(                         # 添加显著性标签
    comparisons=list(c("ETF","ETR"),c("ETF","ETP"),c("ETF","ETE"),c("ETR","ETP"),c("ETR","ETE"),c("ETP","ETE")),
    step_increase = 0.1,
    test="wilcox.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=F,   # 标签样式F为数字，T为*号
    size=0.6,textsize = 5,color="black"#,y_position = c(1)
  ) 
p

# 获取精确p值
a <- wilcox.test(shannon ~ Cluster, 
                 data = subset(shannon, Cluster %in% c("1", "4")),
                 exact = FALSE)
a$p.value



ggsave("FigS1_ent_alphadiversity.pdf",device = cairo_pdf,width =4.5, height =5,p)
save(shannon,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/FigS1_ent_alphadiversity.RData")

#========Fig1D_ent_distrubution==============
Cluster$group <- factor(Cluster$group,levels = c("Health","HCC"))
Cluster$Enterotype <- ifelse(Cluster$Cluster=="1","ETF1",ifelse(Cluster$Cluster=="2","ETF2",
                                                                ifelse(Cluster$Cluster=="3","ETP","ETE")))
Cluster$Enterotype <- factor(Cluster$Enterotype,levels = c("ETF1","ETF2","ETP","ETE"))
# 卡方检验用于检验肠型在HCC与Health之间的分布/组成是否有差异
# clinical <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/clinical_filter.txt",header = T)
# Cluster$Group <- clinical$group[match(rownames(Cluster), clinical$Row.names)]
# Cluster$shannon <- shannon$shannon[match(rownames(Cluster),rownames(shannon))]
# #write.table(Cluster,"/groups/ProHuShiX/home/xiashufen/bigmeta/Enterotype/a_diver_in_4enteritypes.txt")

col2 <- c("#5ba787", "#fcd364", "#f56a5f", "#2b6cbe")
col2 <- c("#F7CCC1", "#BCD1EC", "#ECA579", "#BBD5A6")
col2 <- c("#99c48a","#88adcf","#f4b584","#e3867f","#a18cbf")
col2 <- c("#A8D5BA", "#9DC3E6", "#F9CB9C", "#F4A6A6")
col2 <- c("#B5EAD7", "#C7CEEA", "#FFDAC1", "#FF9AA2")
col2 <- c( "#66C2A5","#4D7EA8", "#EFC000","#D33F49")
col2 <- c("#A8D5BA",  "#90B4E0",  "#F7C59F", "#E28F8F")
pval1 <- round(chisq.test(table(Cluster$Enterotype,Cluster$group))[["p.value"]],4)
p1.1 <- ggplot(Cluster,aes(x =group, fill = Enterotype)) +
  geom_bar(position = "fill")+scale_fill_manual(values = col2)+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "right")

p1.4 <- ggplot()+geom_text(aes(x =2,y = 0,
                               label =  paste0("p=",signif(pval1,4)),size =5.0,family = "sans",fontface = 1))+
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
ggsave("Fig1D_ent_distrubution.pdf",device = cairo_pdf,width =3, height =5,p1)
save(Cluster,pval1,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/Fig1D_ent_distrubution.RData")

#====ent_distribution_ETnonE/ETE====
Cluster$enterotype <- as.factor(ifelse(Cluster$new_group==1,"ETE","ETnonE"))
Cluster$enterotype <- factor(Cluster$enterotype,levels=c("ETnonE","ETE"))
pval3 <- chisq.test(table(Cluster$enterotype,Cluster$group))$p.value
c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
c( "#F7CCC1", "#BCD1EC")
p3.1 <- ggplot(Cluster,aes(x =group, fill = enterotype)) +
  geom_bar(position = "fill")+scale_fill_manual(values =  c("#C1C1C1","#fe9014"))+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "right")

p3.4 <- ggplot()+geom_text(aes(x =2,y = 0,
                               label =  paste0("p=",signif(pval3,4)),size =5.0,family = "sans",fontface = 1))+
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
p3 <- p3.4+p3.1+plot_layout(heights = c(0.1,1),widths = c(1,1),ncol = 1,nrow = 2)
p3
ggsave("./ent_distrubution_ETEETnonE.pdf",device = cairo_pdf,width =3, height =5,p3)

#====ent_distribution_ETnonF1/ETF1====
c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
Cluster$enterotype1 <- as.factor(ifelse(Cluster$Cluster==1,"ETF1","ETnonF1"))
Cluster$enterotype1 <- factor(Cluster$enterotype1,levels=c("ETnonF1","ETF1"))
pval4 <- chisq.test(table(Cluster$enterotype1,Cluster$group))$p.value

p4.1 <- ggplot(Cluster,aes(x =group, fill = enterotype1)) +
  geom_bar(position = "fill")+scale_fill_manual(values =  c("#C1C1C1","#5ba787"))+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "right")

p4.4 <- ggplot()+geom_text(aes(x =2,y = 0,
                               label =  paste0("p=",signif(pval4,4)),size =5.0,family = "sans",fontface = 1))+
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
p4 <- p4.4+p4.1+plot_layout(heights = c(0.1,1),widths = c(1,1),ncol = 1,nrow = 2)
p4
ggsave("./ent_distrubution_ETF1ETnonF1.pdf",device = cairo_pdf,width =3, height =5,p4)

#====ent_distribution_ETnonF2/ETF2====
c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
Cluster$enterotype2 <- as.factor(ifelse(Cluster$Cluster==2,"ETF2","ETnonF2"))
Cluster$enterotype2 <- factor(Cluster$enterotype2,levels=c("ETnonF2","ETF2"))
pval5 <- chisq.test(table(Cluster$enterotype2,Cluster$group))$p.value

p5.1 <- ggplot(Cluster,aes(x =group, fill = enterotype2)) +
  geom_bar(position = "fill")+scale_fill_manual(values =  c("#C1C1C1","#fa9fcb"))+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "right")

p5.4 <- ggplot()+geom_text(aes(x =2,y = 0,
                               label =  paste0("p=",signif(pval5,4)),size =5.0,family = "sans",fontface = 1))+
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
p5 <- p5.4+p5.1+plot_layout(heights = c(0.1,1),widths = c(1,1),ncol = 1,nrow = 2)
p5
ggsave("./ent_distrubution_ETF2ETnonF2.pdf",device = cairo_pdf,width =3, height =5,p5)

#====ent_distribution_ETnonP/ETP====
c("#5ba787","#fa9fcb","#fe9014","#6e7ca5")
Cluster$enterotype3 <- as.factor(ifelse(Cluster$Cluster==3,"ETP","ETnonP"))
Cluster$enterotype3 <- factor(Cluster$enterotype3,levels=c("ETnonP","ETP"))
pval6 <- chisq.test(table(Cluster$enterotype3,Cluster$group))$p.value

p6.1 <- ggplot(Cluster,aes(x =group, fill = enterotype3)) +
  geom_bar(position = "fill")+scale_fill_manual(values =  c("#C1C1C1","#fe9014"))+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "right")

p6.4 <- ggplot()+geom_text(aes(x =2,y = 0,
                               label =  paste0("p=",signif(pval6,4)),size =5.0,family = "sans",fontface = 1))+
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
p6 <- p6.4+p6.1+plot_layout(heights = c(0.1,1),widths = c(1,1),ncol = 1,nrow = 2)
p6
ggsave("./ent_distrubution_ETPETnonP.pdf",device = cairo_pdf,width =3, height =5,p6)

#====ent_distribution_Gender/Age/BMI/Smoke/Drink====
Cluster <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt",header=T)
Cluster$Enterotype <- ifelse(Cluster$Cluster=="1","ETF1",ifelse(Cluster$Cluster=="2","ETF2",
                                                                ifelse(Cluster$Cluster=="3","ETP","ETE")))
Cluster$Enterotype <- factor(Cluster$Enterotype,levels = c("ETF1","ETF2","ETP","ETE"))
Cluster$male <- factor(ifelse(Cluster$male==1,"Male","Female"),levels = c("Female","Male"))
Cluster$age <- as.numeric(Cluster$age)
Cluster$weight <- as.numeric(Cluster$weight)
Cluster$smoking <- factor(ifelse(Cluster$smoking==1,"Yes","No"),levels = c("No","Yes"))
Cluster$drinking <- factor(ifelse(Cluster$drinking==1,"Yes","No"),levels = c("No","Yes"))

tmp.pheno1 <- c("male","smoking","drinking")
tmp.pheno2 <- c("age","weight")

setwd("./ent_distribution")

i="drinking"
for (i in tmp.pheno1){
  tmp_df <- Cluster[, c("Enterotype", i)]
  colnames(tmp_df) <- c("Enterotype", "phenotype")  # 重命名便于通用使用
  tmp_df <- tmp_df[!is.na(tmp_df$phenotype), ]      # 去除 NA
  
  pval <- chisq.test(table(tmp_df$Enterotype,tmp_df$phenotype))$p.value
  p6.1 <- ggplot(tmp_df,aes(x = phenotype, fill = Enterotype)) +
    geom_bar(position = "fill")+scale_fill_manual(values = col2)+theme_bw()+
    ylab('Pecantage')+xlab("")+theme(legend.position = "right")
  p6.4 <- ggplot()+geom_text(aes(x =2,y = 0,
                                 label =  paste0(i," p=",signif(pval,4)),size =5.0,family = "sans",fontface = 1))+
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
  p6 <- p6.4+p6.1+plot_layout(heights = c(0.1,1),widths = c(1,1),ncol = 1,nrow = 2)
  print(p6)
  ggsave(paste0(i,"_ent_pheno_distribution250705.pdf"),width = 8,height = 6,p6)
}
for (i in tmp.pheno2){
  tmp_df <- Cluster[, c("Enterotype", i)]
  colnames(tmp_df) <- c("Enterotype", "phenotype")  # 重命名便于通用使用
  tmp_df <- tmp_df[!is.na(tmp_df$phenotype), ]      # 去除 NA
  
  p <- ggplot(tmp_df, aes(x=Enterotype, y=phenotype, color=Enterotype)) + 
    geom_violin(trim=FALSE,fill="white")+
    geom_boxplot(width=0.15, fill="white")+
    scale_color_manual(values=col2)+
    xlab("")+ylab(i)+
    theme_classic() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),legend.position = 'none',
          axis.title.y = element_text(size=15))+
    geom_signif(                         # 添加显著性标签
      comparisons=list(c("ETF1","ETF2"),c("ETF1","ETP"),c("ETF1","ETE"),c("ETF2","ETP"),c("ETF2","ETE"),c("ETP","ETE")),
      step_increase = 0.1,
      test="wilcox.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
      map_signif_level=F,   # 标签样式F为数字，T为*号
      size=0.6,textsize = 5,color="black"#,y_position = c(1)
    ) 
  print(p)
  ggsave(paste0(i,"_ent_distribution250705.pdf"),width = 8,height = 6,p)
}


#====ent_distribution_curetype====
sur <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/MGSinfo_0326.txt")
# sur$患者当次治疗类型 <- as.factor(ifelse(sur$患者当次治疗类型=="原发初治","Primary treatment",
#                                  ifelse(sur$患者当次治疗类型=="原发再次治疗","Primary retreatment",
#                                         ifelse(sur$患者当次治疗类型=="复发首次治疗","Relapse initial treatment",
#                                                "Relapse retreatment"))))
# sur$患者当次治疗类型 <- factor(sur$患者当次治疗类型,levels = c("Primary treatment","Primary retreatment","Relapse initial treatment","Relapse retreatment"))
sur$患者当次治疗类型 <- as.factor(ifelse(sur$患者当次治疗类型%in%c("原发初治","原发再次治疗"),"Primary treatment","Recurrence retreatment"))
sur$患者当次治疗类型 <- factor(sur$患者当次治疗类型,levels = c("Primary treatment","Recurrence retreatment"))

sur$Enterotype <- as.factor(Cluster$Cluster[match(sur$sample,Cluster$Row.names)])
pval2 <- chisq.test(table(sur$患者当次治疗类型,sur$Enterotype))[["p.value"]]
p2.1 <- ggplot(sur,aes(x =`患者当次治疗类型`, fill = Enterotype)) +
  geom_bar(position = "fill")+scale_fill_manual(values = col2)+theme_bw()+
  ylab('Pecantage')+xlab("")+theme(plot.title = element_text(hjust=0.5),legend.position = "none")+
  ggtitle(label=paste0("p=",signif(pval2,2)))#,size =5.0,family = "sans",fontface = 1)

p2.1
ggsave("ent_distrubution_curetype250623.pdf",device = cairo_pdf,width =8, height =5,p2.1)
#save(Cluster,pval1,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/Fig1D_ent_distrubution.RData")

#========FigS1_ent_survival===============
library(readxl)
library(survival)
library(survminer)
sur <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/MGSinfo_0326.txt")
sur$status <- ifelse(sur$OS=="死亡",1,0)
sur$Enterotype <- Cluster$Cluster[match(sur$sample,Cluster$Row.names)]
pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/FigS1_ent_survival.pdf",width =4, height =5)
ggsurv <- ggsurvplot(survfit(Surv(OSTime,status) ~ Enterotype, # 创建的拟合对象
                             data = sur),  # 指定变量数据来源
                     conf.int = F, # 显示置信区间
                     pval = TRUE, # 添加P值
                     risk.table = T, # 添加风险表
                     #surv.median.line = "hv",  # 添加中位生存时间线
                     xlab = "Follow up time(m)",
                     #title = "5-year-follow",
                     palette =col2,
                     break.x.by = 12)
ggsurv$plot <- ggsurv$plot + 
  guides(color = guide_legend(title="Enterotype",labels= c("ETF","ETR","ETP","ETE"))) 
ggsurv
dev.off()
save(sur,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/FigS1_ent_survival.RData")

#========Fig1E_ent_survival===============
library(survival)
library(survminer)
#====enterotype1 vs enterotype 4====
sur <- sur[sur$Enterotype%in%c("1","4"),]
sur$status <- ifelse(sur$OS=="死亡",1,0)
pdf("Fig1E_ent_survival.pdf",width =4, height =5)
ggsurv <- ggsurvplot(survfit(Surv(OSTime,status) ~ Enterotype, # 创建的拟合对象
                             data = sur),  # 指定变量数据来源
                     conf.int = F, # 显示置信区间
                     pval = TRUE, # 添加P值
                     risk.table = T, # 添加风险表
                     #surv.median.line = "hv",  # 添加中位生存时间线
                     xlab = "Follow up time(m)",
                     #title = "5-year-follow",
                     palette =c(col2[1],col2[3]),
                     break.x.by = 12)
ggsurv$plot <- ggsurv$plot + 
  guides(color = guide_legend(title="Enterotype",labels= c("ETF","ETR","ETP","ETE"))) 
ggsurv
dev.off()
save(sur,file="./Fig1E_ent_survival.RData")

#====Fig 1F 3cohort GBS====
#=======FigS1_socre_group=============
#----bigmeta----
Cluster <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt")
Cluster$group <- factor(Cluster$group,levels = c("Health","HCC"))
p1 <- ggplot(Cluster, aes(x=group, y=balance_value, color=group)) + 
  geom_violin(trim=FALSE,aes(fill=group,alpha=0.8))+
  geom_boxplot(width=0.15, fill="white")+
  #scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  scale_fill_manual(values=c(col1[1],col1[3]))+
  scale_color_manual(values=c(col1[1],col1[3]))+
  xlab("")+labs(y="balance score")+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.line = element_line(colour = "black"),
        legend.position = 'none',
        axis.title.y = element_text(size=15))+
  geom_signif(                         # 添加显著性标签
    #comparisons=list(c("ETM","ETF"),c("ETM","ETE"),c("ETM","ETB"),c("ETF","ETE"),c("ETF","ETB"),c("ETE","ETB")),
    comparisons=list(c("Health","HCC")),
    step_increase = 0.1,
    test="wilcox.test",   # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=F,   # 标签样式F为数字，T为*号
    size=0.6,textsize = 5,color="black"#,y_position = c(1)
  ) +ggtitle("Cohort1")

save(Cluster,p1,file="./FigS1_bigmeta_score_group.RData")

#----renxx----
#balance_rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon1.txt")
balance_rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon_0.1pre.txt")

HCC114_Health37 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
balance_rxx <- balance_rxx %>% 
  filter(sample %in% HCC114_Health37$sample)
balance_rxx$group <- as.factor(ifelse(grepl("H",balance_rxx$sample),"Health","HCC"))
balance_rxx$group <- factor(balance_rxx$group,levels = c("Health","HCC"))
p2 <- ggplot(balance_rxx, aes(x=group, y=balance_value, color=group)) + 
  geom_violin(trim=FALSE,aes(fill=group,alpha=0.8))+
  geom_boxplot(width=0.15, fill="white")+
  #scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  scale_fill_manual(values=c(col1[1],col1[3]))+
  scale_color_manual(values=c(col1[1],col1[3]))+
  xlab("")+labs(y="balance score")+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.line = element_line(colour = "black"),
        legend.position = 'none',
        axis.title.y = element_text(size=15))+
  geom_signif(                         # 添加显著性标签
    #comparisons=list(c("ETM","ETF"),c("ETM","ETE"),c("ETM","ETB"),c("ETF","ETE"),c("ETF","ETB"),c("ETE","ETB")),
    comparisons=list(c("Health","HCC")),
    step_increase = 0.1,
    test="wilcox.test",                     # "t 检验，比较两组（参数）" = "t.test","Wilcoxon 符号秩检验，比较两组（非参数）" = "wilcox.test"
    map_signif_level=F,   # 标签样式F为数字，T为*号
    size=0.6,textsize = 5,color="black"#,y_position = c(1)
  ) +ggtitle("Cohort2")
save(balance_rxx,p2,file="./FigS1_rxx_score_group0.1.RData")
p2
ggsave("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/plots/Fig1F_FAHcohort2_BS_250915.pdf",width = 5,height = 6)

load("./RData/FigS1_bigmeta_score_group.RData")
load("./RData/FigS1_rxx_score_group0.1.RData")

# data3 <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/newMGS_BS&shanno")
data3 <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/newMGS_BS&shannon_0.1pre.txt")
clinical3 <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/分析用临床信息_250302.xlsx")
data3$group <- factor(clinical3$group[match(data3$sample,clinical3$sample)],levels=c("Health","Hepatitis","Cirrhosis","HCC"))

p3 <- ggplot(data3, aes(x=group, y=balance_value, color=group)) + 
  geom_violin(trim=FALSE,aes(fill=group,alpha=0.8))+
  geom_boxplot(width=0.15, fill="white")+
  #scale_fill_manual(values=c("#44BB99","#EE8866","#AA4499"))+
  scale_fill_manual(values=c(col4))+
  scale_color_manual(values=c(col4))+
  xlab("")+labs(y="balance score")+
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.line = element_line(colour = "black"),
        legend.position = 'none',
        axis.title.y = element_blank())+
  stat_compare_means(comparisons = list(c("Cirrhosis","HCC"),c("Hepatitis","HCC"),c('Health','HCC')))+
  ggtitle("Cohort3")

pdf("./plots/Fig1F_3cohort_score_group250915.pdf",width = 7,height = 5)
p1+p2+p3+plot_layout(heights = c(1,1,1),widths = c(1,1,1.6),ncol = 3,nrow = 1)
dev.off()


#====ETE vs ETnonE====
sur <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/MGSinfo_0514.xlsx")
sur <- sur %>% mutate(Cluster=Cluster$Cluster[match(sur$sample,Cluster$Row.names)])
sur$Ent_new <- ifelse(sur$Cluster=="4","ETE","ETnonE")
sur$Ent_new <- factor(sur$Ent_new,levels=c("ETnonE","ETE"))
sur$status <- ifelse(sur$OS=="死亡",1,0)
pdf("Fig1E_ent_survival_V2.pdf",width =4, height =5)
ggsurv <- ggsurvplot(survfit(Surv(OSTime,status) ~ Ent_new, # 创建的拟合对象
                             data = sur),  # 指定变量数据来源
                     conf.int = F, # 显示置信区间
                     pval = TRUE, # 添加P值
                     risk.table = T, # 添加风险表
                     #surv.median.line = "hv",  # 添加中位生存时间线
                     xlab = "Follow up time(m)",
                     #title = "5-year-follow",
                     palette =col3,
                     break.x.by = 12)
ggsurv$plot <- ggsurv$plot + 
  guides(color = guide_legend(title="Enterotype")) 
ggsurv
dev.off()
save(sur,file="/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/Fig1E_ent_survival_V2.RData")

#=========Fig1H_socre_survival===============
library(contsurvplot)
library(ggplot2)
library(survival)
library(riskRegression)
library(pammtools)

sur$BMI[is.na(sur$BMI)]=median(sur$BMI[!is.na(sur$BMI)])
sur$Gender <- as.factor(sur$Gender)
sur$Smoke <- as.factor(sur$Smoke)
sur$Drink <- as.factor(sur$Drink)
sur$status <- ifelse(sur$OS=="死亡",1,0)
sur$balance_value <- as.numeric(Cluster$balance_value[match(sur$sample,Cluster$Row.names)])
model <- coxph(Surv(OSTime,status) ~ balance_value + Age+ Gender+BMI+Smoke+Drink , data=sur, x=TRUE)
model
a <- as.data.frame(summary(model)$coefficients)
# 提取 HR 和 95% CI
hr_df <- as.data.frame(summary(model)$conf.int)

#这个是为了得到risk.table的
ggsurvplot(survfit(Surv(OSTime,status) ~ 1, # 创建的拟合对象
                   data = sur),  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = T, # 添加风险表
           #surv.median.line = "hv",  # 添加中位生存时间线
           xlab = "Follow up time(m)",
           title = "5-year-follow",
           break.x.by = 12)

p1 <- plot_surv_area(time="OSTime",
                     status="status",
                     variable="balance_value",
                     data=sur,start_color ="#92B096",end_color="#5D7160",
                     model=model,discrete=TRUE,bins = 3)+
  guides(color = guide_legend(title = "balance score", 
                              labels = c("Low(-5.45--1.10)", "Medium(-1.10-3.23)", "High(3.23-7.57)")))

p1 <- p1 +
  annotate("text", x = 50, y = 1,size=6,label = paste0("pvalue=",signif(a[1,5],2)))+
  theme(legend.position =c(0.2, 0.2))

p1
ggplot_build(p1)$data

table(Cluster$group)
Cluster$group <- as.factor(Cluster$group)
ggboxplot(Cluster, x = "group", y = "balance_value",
          color = "group", palette = "jco") +
  stat_compare_means(comparisons = list(c("HCC","Health")))
ggsave("./BS_HCC_Health.pdf",width = 5,height = 6)
#====Renxx data
# meta_all <- meta
# load("/groups/ProHuShiX/home/liuyuyao/renxx_cohort/metaG/meta_sur.RData")
# 114HCC 37Health
sur_rxx <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
# 114HCC 53Health 存在一些不那么健康的健康人
balance_rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon_0.1pre.txt")
# 131HCC 37Health
Group <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxCohortInfo_131HCC_37Health_250310.xlsx")

balance_rxx <- balance_rxx %>% mutate(group=as.factor(ifelse(grepl("H",balance_rxx$sample),"Health","HCC")))
# balance_rxx <- balance_rxx %>% mutate(median_group = as.factor(ifelse(balance_rxx$balance_value > median(balance_rxx$balance_value),"High","Low")))

# 按照sur_rxx的样本
sur_rxx <- merge(sur_rxx,balance_rxx[,c("sample","balance_value")],by="sample")
# 相当于剔除Health
sur_rxx <- sur_rxx[is.na(sur_rxx$OS)==F,]
#sur_rxx$OS[which(sur_rxx$sample%in%c("1815","1068"))] <- "死亡"
sur_rxx$status <- ifelse(sur_rxx$OS=="死亡",1,0)
# 剔除6个HCC后 HCC=108 Health=37
data_survival_health <- sur_rxx %>% filter(group=="Health")
data_survival<- sur_rxx#[!sur_rxx$sample%in%c("1815","1068","2140","1302","1769","2155","1805","1999","1401","1650","1508","1379","2014","2157","2080"),]#,,

data_survival$BMI[is.na(data_survival$BMI)]=median(data_survival$BMI[!is.na(data_survival$BMI)])
data_survival$Gender <- as.factor(data_survival$Gender)
data_survival$Smoke <- as.factor(data_survival$Smoke)
data_survival$Alcohol <- as.factor(data_survival$Alcohol)

model <- coxph(Surv(OSTime,status) ~ balance_value + Age+ Gender+BMI+Smoke+Alcohol , data=data_survival, x=TRUE)
model
a1 <- as.data.frame(summary(model)$coefficients)

# 提取 HR 和 95% CI
hr_df1 <- as.data.frame(summary(model)$conf.int)


pdf("./plots/Fig1F_risktable250915.pdf",width = 4,height = 5)
#这个是为了得到risk.table的
ggsurvplot(survfit(Surv(OSTime,status) ~ 1, # 创建的拟合对象
                   data = sur),  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = T, # 添加风险表
           #surv.median.line = "hv",  # 添加中位生存时间线
           xlab = "Follow up time(m)",
           title = "5-year-follow",
           break.x.by = 20)

#这个是为了得到risk.table的
ggsurvplot(survfit(Surv(OSTime,status) ~ 1, # 创建的拟合对象
                   data = data_survival),  # 指定变量数据来源
           conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
           risk.table = T, # 添加风险表
           #surv.median.line = "hv",  # 添加中位生存时间线
           xlab = "Follow up time(m)",
           title = "5-year-follow",
           break.x.by = 20)
dev.off()

p2 <- plot_surv_area(time="OSTime",
                     status="status",
                     variable="balance_value",
                     data=data_survival,start_color ="#92B096",end_color="#5D7160",
                     model=model,discrete=TRUE,bins = 3)+
  annotate("text", x = 70, y = 1, size=6,label =paste0("pvalue=",signif(a1[1,5],3)))
p2 <- p2+theme(legend.position =c(0.2, 0.2))+
  guides(color = guide_legend(title = "balance score", 
                              labels = c("Low(-3.80-0.35)", "Medium(0.35-4.50)", "High(4.50-8.65)")))

p2
pdf("Fig1F_score_survivial.pdf",width = 8,height = 4)
(p1+ggtitle("Cohort1"))+
  (p2+ggtitle("Cohort2"))+plot_layout(heights = c(1,1),widths = c(1,1),ncol = 2,nrow = 1)
dev.off()
save(data_survival,sur,file="./Fig1F.Rdata")

pdf("./plots/Fig1H_rxx_score_survivial250915.pdf",width = 4,height = 4)
(p2+ggtitle("FAH cohort2"))
dev.off()

HCC114_Health37 <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
balance_rxx <- balance_rxx %>% 
  filter(sample %in% HCC114_Health37$sample)
#balance_rxx_f <- balance_rxx[!balance_rxx$sample%in%c("1815","1068","2140","1302","1769","2155","1805","1999","1401","1650","1508","1379","2014","2157","2080"),]
table(balance_rxx$group)
ggboxplot(balance_rxx, x = "group", y = "balance_value",
          color = "group", palette = "jco") +
  stat_compare_means(comparisons = list(c("HCC","Health"))) 
ggsave("./BS_HCC_Health_rxx.pdf",width = 5,height = 6)


#=====Fig 1G====
#=======整理bigmeta数据==============
MGS_meta <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt")
colnames(MGS_meta)[1] <- "sample"
#MGS_meta <- subset(MGS_meta[,c("sample","balance_value","group"),drop=F])
##shannon
load("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/RData/FigS1_bigmeta_adiv.RData")
MGS_meta$shannon <- shannon$shannon[match(MGS_meta$sample,rownames(shannon))]

##HCC
MGSinfo <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/MGSinfo_0514.xlsx")
bigmeta_HCC <- MGSinfo[,c("sample","Age","Gender","BMI")]
bigmeta_HCC <- merge(bigmeta_HCC,MGS_meta[MGS_meta$group=="HCC",],by="sample")
bigmeta_HCC <- bigmeta_HCC[,c("sample","Age","Gender","BMI","balance_value","shannon","group")]
##Health
Heainfo <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/健康人临床信息.xlsx",sheet = "整理2")
colnames(Heainfo)[3] <- "Age"
Heainfo$sample <- MGS_meta$sample[match(Heainfo$Name,MGS_meta$name)]
Heainfo$shannon <- shannon$shannon[match(Heainfo$sample,rownames(shannon))]
Heainfo$balance_value <- MGS_meta$balance_value[match(Heainfo$sample,rownames(shannon))]
bigmeta_Hea <- Heainfo[,c("sample","Age","Gender","BMI","balance_value","shannon")]
bigmeta_Hea$group <- "Health"

bigmeta <- rbind(bigmeta_HCC,bigmeta_Hea)

bigmeta$BMI[is.na(bigmeta$BMI)==T] <- median(bigmeta$BMI,na.rm = T)
bigmeta$outcome <- ifelse(bigmeta$group=="HCC",1,0)
bigmeta$Age <- as.numeric(bigmeta$Age)
colnames(bigmeta) <- gsub("balance_value","balance_score",colnames(bigmeta))
bigmeta$Gender <- ifelse(bigmeta$Gender=="male",1,0)

#=======整理Renxx队列信息=============
RenxxInfo <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
load("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/RData/FigS1_rxx_score_group.RData")
old <- balance_rxx
balance_rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon_0.1pre.txt")
colnames(balance_rxx)[2] <- "balance_score"
balance_rxx$group <- old$group[match(balance_rxx$sample,old$sample)]
##balance_rxx
#meta <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
##meta
#load("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig1/RData/FigS1_rxx_adiv.RData")
##data2
renxx <- merge(balance_rxx[,c("sample","balance_score","shannon")],RenxxInfo[,c("sample","group","Gender","Age","BMI")],by="sample")
#renxx <- merge(renxx,data2[,c("sample","shannon")],by="sample")
renxx$Age[is.na(renxx$Age)] <- median(renxx$Age,na.rm = T)
renxx$BMI[is.na(renxx$BMI)==T] <- median(renxx$BMI,na.rm = T)
renxx$outcome <- ifelse(renxx$group=="HCC",1,0)
#renxx$Gender <- ifelse(renxx$Gender=="male",1,0)

#=======整理新收HCC队列信息============
new_HCC_balance <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/newMGS_BS&shannon_0.1pre.txt",header=T)
colnames(new_HCC_balance)[2] <- "balance_score"

# 加上漏掉的被分到Health的8个Hep
hep <- c("HC128","HC129","HC131_1","HC287","HC405","HC525","HC608","HC616")

new_HCC_clinical1 <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/分析用临床信息_250414.xlsx") %>%
  filter(group%in%c("Health","HCC")|sample%in%hep) %>%
  filter(sample%in%new_HCC_balance$sample)
new_HCC_clinical1$BMI[new_HCC_clinical1$sample=="HC407"] <- as.numeric(new_HCC_clinical1$Weight[new_HCC_clinical1$sample=="HC407"])/as.numeric(new_HCC_clinical1$Height[new_HCC_clinical1$sample=="HC407"])^2
new_HCC_clinical1$BMI[new_HCC_clinical1$sample=="4290"] <- 26.67

new_HCC_clinical <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/肝炎肝硬化整理250904.xlsx") %>% 
  filter(sample%in%new_HCC_balance$sample)


new_HCC_clinical$Height <- gsub("无可录信息",NA,new_HCC_clinical$Height);new_HCC_clinical$Height <- as.numeric(new_HCC_clinical$Height)
new_HCC_clinical$Weight <- as.numeric(new_HCC_clinical$Weight)
new_HCC_clinical$BMI <- new_HCC_clinical$Weight/(new_HCC_clinical$Height/100)^2

new_HCC_clinical <- new_HCC_clinical[,c("sample","group","Gender","Age","BMI")]
new_HCC_clinical1 <- new_HCC_clinical1[,c("sample","group","Gender","Age","BMI")]

new_HCC_clinical <- rbind(new_HCC_clinical,new_HCC_clinical1)

new_HCC <- new_HCC_clinical[new_HCC_clinical$group%in%c("HCC","Health"),c("sample","group","Age","Gender","BMI")]
new_HCC <- merge(new_HCC,new_HCC_balance[,c("sample","shannon","balance_score")],by="sample")
new_HCC$BMI[new_HCC$BMI=="."] <- NA
new_HCC$BMI <- as.numeric(new_HCC$BMI)
new_HCC$BMI[is.na(new_HCC$BMI)==T] <- median(new_HCC$BMI,na.rm = T)
new_HCC$Age <- as.numeric(new_HCC$Age)
new_HCC$outcome <- ifelse(new_HCC$group=="HCC",1,0)
new_HCC$Gender <- ifelse(new_HCC$Gender%in%c("male","男"),1,0)

#======整理肝炎队列信息===========
new_MGS <- new_HCC_clinical[new_HCC_clinical$group%in%c("HCC","Hepatitis","Cirrhosis"),c("sample","group","Age","Gender","BMI")]
new_MGS <- merge(new_MGS,new_HCC_balance[,c("sample","shannon","balance_score")],by="sample")
new_MGS$BMI[new_MGS$BMI=="."] <- NA
new_MGS$BMI <- as.numeric(new_MGS$BMI)
new_MGS$BMI[new_MGS$sample=="LC017_MJP"] <- NA
#new_MGS$BMI[is.na(new_MGS$BMI)==T] <- median(new_MGS$BMI,na.rm = T)
new_MGS$Age <- as.numeric(new_MGS$Age)
new_MGS$outcome <- ifelse(new_MGS$group=="HCC",1,0)
new_MGS$Gender <- ifelse(new_MGS$Gender%in%c("male","男"),1,0)

Hep <- new_MGS[new_MGS$group%in%c("HCC","Hepatitis"),]
Hep <- Hep[!(Hep$group=="Hepatitis"&is.na(Hep$BMI)==T),]
#HCC Hepatitis 
#134      86 

# 只有Hep中BMI含NA，删去则无
#Hep$BMI[is.na(Hep$BMI)==T] <- median(Hep$BMI,na.rm=T)

#=======整理肝硬化患者信息===============
Cir <- new_MGS[new_MGS$group%in%c("HCC","Cirrhosis"),]
Cir$outcome <- ifelse(Cir$group=="HCC",1,0) 
Cir <- Cir[!(Cir$group=="Cirrhosis"&is.na(Cir$BMI)==T),]
#Cirrhosis       HCC 
# 63      134  

# 只有Cir中BMI含NA，删去则无
#Cir$BMI[is.na(Cir$BMI)==T] <- median(Cir$BMI,na.rm=T)

#=========整理CLD信息============
CLD <- new_MGS
CLD <- CLD[!(CLD$group=="Cirrhosis"&is.na(CLD$BMI)==T),]
CLD <- CLD[!(CLD$group=="Hepatitis"&is.na(CLD$BMI)==T),]
CLD$BMI[is.na(CLD$BMI)==T] <- median(CLD$BMI,na.rm=T)
table(CLD$outcome)

#======重跑开始处====
#bigmeta_new <- bigmeta
save(bigmeta,renxx,new_MGS,new_HCC,Hep,Cir,CLD,file="/groups/ProHuShiX/home/xiashufen/bigmeta/xgboost/data_prepare250916.RData")

load("D:/Desktop/HCC_bigmeta/xgboost/data_prepare250916.RData")
outpath <- "D:/Desktop/HCC_bigmeta/xgboost/test250916/"

#=========2:8================
set.seed(5364)
train_index <- createDataPartition(bigmeta$outcome, p = 0.8, list = FALSE)
train_data <- bigmeta[train_index, ]  # 80%训练集
valid_data <- bigmeta[-train_index, ] # 20%验证集

train_outcome=train_data$outcome
train_feature=as.matrix(train_data[,c("Age","Gender","BMI","balance_score"),drop=F])

set.seed(85806)
# 定义训练控制（10折交叉验证）
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  # p=0.9,
  classProbs = TRUE,          # 输出概率（用于计算AUC）
  summaryFunction = twoClassSummary,  # 使用AUC评估
  allowParallel = TRUE,      # 启用并行加速
  search = "grid"
)

grid <- expand.grid(
  nrounds = 50,
  max_depth = seq(1,3,1),
  eta = seq(0.1,0.8,0.1),
  gamma = 0,
  colsample_bytree =  seq(0.4,0.8,0.1),
  min_child_weight = 1,
  subsample = seq(0.5,0.8,0.1)
)
#     
set.seed(85806)
xgb_model <- train(
  x = train_data[,c("Age","Gender","BMI","balance_score"),drop=F],
  y = factor(train_data$outcome,levels=c("1","0"),labels=c("HCC","Health")),
  method = "xgbTree",
  trControl = ctrl,
  # tuneGrid = grid,
  metric = "ROC"    # 以AUC作为优化指标
)

best_tune <- xgb_model$bestTune

# 用最佳组合的参数来计算最佳迭代次数
params <- list(
  booster="gbtree",
  objective = "binary:logistic",
  max_depth = best_tune$max_depth,  ##树的深度，数值大容易过拟合
  eta = best_tune$eta,  ##learing rate,更新中减少的步长来防止过拟合
  subsample = best_tune$subsample, ##对于每棵树，随机采样的比例
  gamma = best_tune$gamma,
  min_child_weight=best_tune$min_child_weight,
  colsample_bytree=best_tune$colsample_bytree,
  #nrounds = 50  ##最大迭代次数
  #lambda=cvfit$lambda.min, ##默认为1
  alpha =0 , ##默认为0
  random_state = 555 , # 设置XGBoost内部种子
  nthread =1
)

dtrain <- xgb.DMatrix(data = train_feature, label = train_outcome)
#     
set.seed(85806)

cv_results <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 100,
  nfold = 10,
  metrics = "logloss",
  early_stopping_rounds = 10,
  verbose = 0
)

# 输出最佳迭代次数
best_iteration <- cv_results$best_iteration
print(best_iteration)
#     
set.seed(85806)

model <- xgboost(
  data = train_feature,
  label = train_outcome,
  params = params,
  nrounds = best_iteration,
  verbose = 0
)


# 快速评估
y_pred_prob <- predict(model, train_feature,type="response")
train_auc <- roc(train_outcome, y_pred_prob)$auc  # 预期值约0.6-0.8（取决于数据分布）

test_outcome=valid_data$outcome
test_feature=as.matrix(valid_data[,c("Age","Gender","BMI","balance_score"),drop=F])
# 快速评估
y_pred_prob <- predict(model, test_feature,type="response")
test_auc <- roc(test_outcome, y_pred_prob)$auc  # 预期值约0.6-0.8（取决于数据分布）

##旭鑫师兄验证
rxx_outcome=renxx$outcome
rxx_feature=as.matrix(renxx[,c("Age","Gender","BMI","balance_score"),drop=F])
# 快速评估
y_pred_prob <- predict(model, rxx_feature,type="response")
rxx_auc <- roc(rxx_outcome, y_pred_prob)$auc  # 预期值约0.6-0.8（取决于数据分布）

# 快速评估
new_outcome=new_HCC$outcome
new_feature=as.matrix(new_HCC[,c("Age","Gender","BMI","balance_score"),drop=F])
y_pred_prob <- predict(model, new_feature,type="response")
new_auc <- roc(new_outcome, y_pred_prob)$auc # 预期值约0.6-0.8（取决于数据分布）

# 肝炎验证
Hep_outcome=Hep$outcome
Hep_feature=as.matrix(Hep[,c("Age","Gender","BMI","balance_score"),drop=F])
y_pred_prob <- predict(model, Hep_feature,type="response")
Hep_auc <- roc(Hep_outcome, y_pred_prob)$auc # 预期值约0.6-0.8（取决于数据分布）
Hep_auc

# 肝硬化验证
Cir_outcome=Cir$outcome
Cir_feature=as.matrix(Cir[,c("Age","Gender","BMI","balance_score"),drop=F])
y_pred_prob <- predict(model, Cir_feature,type="response")
Cir_auc <- roc(Cir_outcome, y_pred_prob)$auc # 预期值约0.6-0.8（取决于数据分布）
Cir_auc

# 肝病验证
CLD_outcome=CLD$outcome
CLD_feature=as.matrix(CLD[,c("Age","Gender","BMI","balance_score"),drop=F])
y_pred_prob <- predict(model, CLD_feature,type="response")
CLD_auc <- roc(CLD_outcome, y_pred_prob)$auc # 预期值约0.6-0.8（取决于数据分布）
CLD_auc

res <- data.frame(ratio="8:2",#seed1=random_numbers[j],seed2=random_numbers2[i],
                  train_auc=train_auc,test_auc=test_auc,
                  rxx_auc=rxx_auc,new_auc=new_auc,
                  Hep_auc=Hep_auc,Cir_auc=Cir_auc,CLD_auc=CLD_auc,
                  best_iteration=best_iteration)

#======整理画图代码=============
# 计算每个模型的ROC对象
roc1 <-  roc(test_outcome,  predict(model, test_feature,type="response"),ci=T)
roc2 <-roc(rxx_outcome, predict(model, rxx_feature,type="response"),ci=T)
roc3 <- roc(new_outcome,  predict(model, new_feature,type="response"), ci = TRUE)

# 将ROC对象存入列表
roc_list <- list(`FAH cohort1 test`= roc1, 
                 `FAH cohort2` = roc2, 
                 `TAH cohort` = roc3)

test_auc <- round(auc(roc1),4) 
rxx_auc <- round(auc(roc2),4)
new_auc<- round(auc(roc3),4)

lab1<- paste0( "AUC = ",test_auc )
lab2<- paste0( "AUC = ",rxx_auc)
lab3<- paste0( "AUC = ",new_auc)

library(cols4all)
library(ggplot2)
library(pROC)
mycol3 <- c4a('set1',3)

g2 <- ggroc(roc_list,
            legacy.axes = TRUE,  
            size = 1)+
  theme_bw()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(legend.position = c(.5,.148),
        legend.key.width = unit(0.7,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(color=guide_legend(title = ""))+
  ggtitle("HCC vs HC") + 
  theme(title = element_text(size = 15,face = "bold",color="black"),
        plot.title = element_text(hjust=0.5))+
  scale_x_continuous(expand = c(0.02,0))+
  scale_y_continuous(expand = c(0.02,0))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.text.x=element_text(size = 12,color="black",hjust=0.8))+
  theme(axis.text.y=element_text(size = 13,color="black"))+
  theme(axis.title.y=element_text(size = 15))+ 
  theme(legend.text = element_text(size = 15))+
  theme(panel.border = element_rect(fill=NA,color=NA),
        panel.grid.major.y = element_line(color="lightgrey",linetype=2,size=0.5),   ##设置x,y轴的主次提示线都为空白
        panel.grid.major.x =  element_line(color="lightgrey",linetype=2,size=0.5), 
        panel.grid.minor.x = element_line(color="lightgrey",linetype=3,size=0.5), 
        panel.grid.minor.y =  element_line(color="lightgrey",linetype=3,size=0.5))+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", size=0.6,
               linetype=2)+
  scale_colour_manual(values = mycol3)+ 
  annotate( "text",x = 0.87,y = 0.15, label=lab1,size= 5.5)+
  annotate( "text",x = 0.87,y = 0.10, label=lab2,size= 5.5)+
  annotate( "text",x = 0.87,y = 0.05, label=lab3,size= 5.5)
g2

ci.list <- lapply(roc_list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

color1=list("#E41A1C" ,"#377EB8", "#4DAF4A")
p1 <- g2 + 
  geom_ribbon(
    data = dat.ci.list[[1]],
    aes(x = 1-x, ymin = lower, ymax = upper),
    fill = color1[1],
    alpha = 0.1,
    inherit.aes = F) +
  geom_ribbon(
    data = dat.ci.list[[2]],
    aes(x = 1-x, ymin = lower, ymax = upper),
    fill = color1[2],
    alpha = 0.1,
    inherit.aes = F) +
  geom_ribbon(
    data = dat.ci.list[[3]],
    aes(x = 1-x, ymin = lower, ymax = upper),
    fill = color1[3],
    alpha = 0.2,
    inherit.aes = F)
p1


############肝炎肝硬化
# 计算每个模型的ROC对象
roc1 <-  roc(Hep_outcome,  predict(model, Hep_feature,type="response"),ci=T)
roc2 <-roc(Cir_outcome, predict(model, Cir_feature,type="response"),ci=T)

# 将ROC对象存入列表
roc_list <- list(`Hepatitis vs HCC`= roc1, 
                 `Cirrhosis vs HCC` = roc2)

Hep_auc <- round(auc(roc1),4) 
Cir_auc <- round(auc(roc2),4)

lab1<- paste0( "AUC = ",Hep_auc )
lab2<- paste0( "AUC = ",Cir_auc)

library(cols4all)
library(ggplot2)
library(pROC)
#mycol3 <-c("#BB6449", "#9A68A4")
mycol3 <-c("#DD9123", "#9A68A4")

g2 <- ggroc(roc_list,
            legacy.axes = TRUE,  
            size = 1)+
  theme_bw()+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  theme(legend.position = c(.5,.148),
        legend.key.width = unit(0.7,"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(color=guide_legend(title = ""))+
  ggtitle("HCC vs Hepatitis&Cirrhosis ") + 
  theme(title = element_text(size = 15,face = "bold",color="black"),
        plot.title = element_text(hjust=0.5))+
  scale_x_continuous(expand = c(0.02,0))+
  scale_y_continuous(expand = c(0.02,0))+
  theme(axis.title.x = element_text(size = 14))+
  theme(axis.text.x=element_text(size = 12,color="black",hjust=0.8))+
  theme(axis.text.y=element_text(size = 13,color="black"))+
  theme(axis.title.y=element_text(size = 15))+ 
  theme(legend.text = element_text(size = 15))+
  theme(panel.border = element_rect(fill=NA,color=NA),
        panel.grid.major.y = element_line(color="lightgrey",linetype=2,size=0.5),   ##设置x,y轴的主次提示线都为空白
        panel.grid.major.x =  element_line(color="lightgrey",linetype=2,size=0.5), 
        panel.grid.minor.x = element_line(color="lightgrey",linetype=3,size=0.5), 
        panel.grid.minor.y =  element_line(color="lightgrey",linetype=3,size=0.5))+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", size=0.6,
               linetype=2)+
  scale_colour_manual(values = mycol3)+ 
  annotate( "text",x = 0.87,y = 0.15, label=lab1,size= 5.5)+
  annotate( "text",x = 0.87,y = 0.10, label=lab2,size= 5.5)

g2

ci.list <- lapply(roc_list, ci.se, specificities = seq(0, 1, l = 25))

dat.ci.list <- lapply(ci.list, function(ciobj) 
  data.frame(x = as.numeric(rownames(ciobj)),
             lower = ciobj[, 1],
             upper = ciobj[, 3]))

color1=list("#DD9123", "#9A68A4")
p2 <- g2 + 
  geom_ribbon(
    data = dat.ci.list[[1]],
    aes(x = 1-x, ymin = lower, ymax = upper),
    fill = color1[1],
    alpha = 0.2,
    inherit.aes = F) +
  geom_ribbon(
    data = dat.ci.list[[2]],
    aes(x = 1-x, ymin = lower, ymax = upper),
    fill = color1[2],
    alpha = 0.1,
    inherit.aes = F) 
p2


library(patchwork)
pdf("./test250916/model_sum_251009.pdf",width=10,height=5)
p1+p2+plot_layout(ncol=2,nrow=1,widths = c(1,1),heights = c(1,1))
dev.off()

#=========Fig1G_score_pheno===========
#=====大样本队列表型图==================
sur <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/MGSinfo_250310.xlsx")

sur$BMI[is.na(sur$BMI)]=median(sur$BMI[!is.na(sur$BMI)])
sur$log2AFP <- log2(sur$AFP)
sur$BCLC <- ifelse(sur$BCLC=="0",0,ifelse(sur$BCLC=="A",1,ifelse(sur$BCLC=="B",2,ifelse(sur$BCLC=="C",3,NA))))
sur$TumorNumber<- ifelse(sur$TumorNumber==">3",4,sur$TumorNumber)
sur$Edmondson <-  ifelse(sur$病理简表Edmondson分级=="1级",1,ifelse(sur$病理简表Edmondson分级=="1-2级",2,ifelse(sur$病理简表Edmondson分级=="2级",3,
                                                                                                 ifelse(sur$病理简表Edmondson分级=="2-3级",4,ifelse(sur$病理简表Edmondson分级=="3级",5,
                                                                                                                                             ifelse(sur$病理简表Edmondson分级=="3-4级",6,ifelse(sur$病理简表Edmondson分级=="4级",7,NA)))))))

sur$`Tumor thromb` <- ifelse(sur$病理简表肉眼癌栓=="是",1,ifelse(sur$病理简表肉眼癌栓=="否",0,NA))
sur$`Satellite leision` <-  ifelse(sur$病理简表卫星病灶=="是",1,ifelse(sur$病理简表卫星病灶=="否",0,NA))
sur$sample <- gsub("-",".",sur$sample)
balance <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt")
sur <- sur[sur$sample%in%balance$Row.names,]
sur$balance_value <- as.numeric(balance$balance_value[match(sur$sample,balance$Row.names)])
sur$cluster <- as.numeric(balance$Cluster[match(sur$sample,balance$Row.names)])
##########连续型变量
pheno <- c("TBil", "ALT", "AST","log2AFP","TumorSize","NLR","BCLC","TumorNumber","Edmondson")

sur$TBil <- as.numeric(sur$TBil)
sur$BCLC <- as.numeric(sur$BCLC)
sur$TumorNumber <- as.numeric(sur$TumorNumber )
sur$Edmondson <- as.numeric(sur$Edmondson)
colnames(sur)[colnames(sur)=="Tumor size"] <- "TumorSize"
colnames(sur)[colnames(sur)=="Tumor number"] <- "TumorNumber"

res_cor <- list()
res_lm <- list()
res_lm_ad <- list()
# 用lm时必矫正性别年龄吸烟喝酒抗生素等
i=1
for (i in 1:length(pheno)){
  tmp.data <- sur[,c(pheno[i],"balance_value")]
  tmp.res <-cor.test(tmp.data[[1]],tmp.data$balance_value,method = "spearman") 
  res_cor[[i]] <- data.frame(rho=tmp.res[["estimate"]],pvalue=tmp.res[["p.value"]],pheno=pheno[i])
  mod1 <- lm(as.formula(paste0(pheno[i],"~balance_value")),data=sur)
  a <- as.data.frame(summary(mod1)[["coefficients"]])
  a<-  as.data.frame(a[-1,])
  a$pheno <- pheno[i]
  res_lm[[i]] <-a
  mod2<-lm(as.formula(paste0(pheno[i],"~balance_value+Age+Gender+BMI+Smoke+Drink")),data = sur)
  a <- as.data.frame(summary(mod2)[["coefficients"]])
  a<-  as.data.frame(a[2,])
  a$pheno <- pheno[i]
  res_lm_ad[[i]] <-a
}
res_cor <- do.call("rbind",res_cor)
res_lm <- do.call("rbind",res_lm)
res_lm_ad <- do.call("rbind",res_lm_ad)
res_lm$adj.p <- p.adjust(res_lm$`Pr(>|t|)`,"fdr")
res_lm_ad$adj.p <- p.adjust(res_lm_ad$`Pr(>|t|)`,"fdr")

###########分类型变量
sur$Child_Pugh <- ifelse(sur$Child_Pugh=="A",0,ifelse(sur$Child_Pugh=="B",1,NA))
colnames(sur)[colnames(sur)=="Tumor thromb"] <- "Tumor_thrombus"
colnames(sur)[colnames(sur)=="Satellite leision"] <- "Satellite_lesions"
colnames(sur)[colnames(sur)=="Asites"] <- "Ascites"
colnames(sur)[colnames(sur)=="pathology_MVI"] <- "MVI"

write.xlsx(sur,"../Final_input/MGSinfo_0514.xlsx")

pheno <- c( "PVTT", "HVTT" , "Ascites", "Child_Pugh","MVI" ,"Tumor_thrombus" , "Satellite_lesions" )


sur$HVTT <- factor(sur$HVTT,levels=c(0,1))
sur$PVTT <- factor(sur$PVTT,levels = c(0,1))
sur$Ascites <- factor(sur$Ascites,levels = c(0,1))
sur$Child_Pugh <- factor(sur$Child_Pugh,levels = c(0,1))
sur$MVI <- factor(sur$MVI,levels = c(0,1))
sur$Tumor_thrombus <- factor(sur$Tumor_thrombus,levels = c(0,1))
sur$Satellite_lesions <- factor(sur$Satellite_lesions,levels = c(0,1))

res_glm <- list()
res_glm_ad <- list()
i = 2
for (i in 1:length(pheno)){
  # glm.md <- glm(as.formula(paste0(pheno[i],"~balance_value")),data=sur, family = "binomial")
  # a <- summary(glm.md)$coefficients %>% as.data.frame()
  # a<-  as.data.frame(a[-1,])
  # a$pheno <- pheno[i]
  # res_glm[[i]] <-a
  mod2<-glm(as.formula(paste0(pheno[i],"~balance_value+Age+Gender+BMI+Drink+Smoke")),data = sur,family = "binomial")
  a <- summary(mod2)$coefficients %>% as.data.frame()
  a<-  as.data.frame(a[2,])
  a$pheno <- pheno[i]
  res_glm_ad[[i]] <-a
}
res_glm <- do.call("rbind",res_glm)
res_glm_ad <- do.call("rbind",res_glm_ad)
res_glm$adj.p <- p.adjust(res_glm$`Pr(>|z|)`,"BH")
res_glm_ad$adj.p <- p.adjust(res_glm_ad$`Pr(>|z|)`,"BH")

#save(sur,file="bigmeta_sur_250107.RData")
#==========画大队列表型棒棒糖图====================
colnames(res_lm_ad)[4] <- "pvalue"
colnames(res_glm_ad)[4] <- "pvalue"
forplot1 <- rbind(res_lm_ad[,c("Estimate","pvalue", "pheno", "adj.p")],res_glm_ad[,c("Estimate","pvalue", "pheno", "adj.p")])

forplot1 <- forplot1[!forplot1$Estimate>0,]
forplot1 <- forplot1%>%arrange(Estimate)
forplot1$pheno <- factor(forplot1$pheno,levels = rev(forplot1$pheno))
forplot1$stars <- cut(forplot1$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))


cl <- c("#F39B7F","#F7DB70","#71D0F5FF" ,"#D2AF81FF","#D5E4A2FF" ,"#197EC0FF","#709AE1FF"  ,
        "#F05C3BFF" , "#46732EFF",  "#E9BFDD", "#FFCC9E", "#BCD1EC" ,"#C6C09C","steelblue2","#AABAC2")

#AST不显著，Estimate又太大，先去掉
#forplot1 <- forplot1[!forplot1$pheno%in%c("AST"),]
p1 <- ggplot(forplot1, aes(x = pheno, y = Estimate)) +
  geom_segment( aes(x = pheno, xend = pheno, y = 0, yend = Estimate),color = "grey40")+
  scale_y_reverse(limits = c(0, -0.5), expand = c(0, 0))+ #控制线段的参数，见下
  geom_point(color=  "#7ca6be",size = 12) + 
  geom_text(aes(label=stars), fontface="bold", size=6 ,nudge_y=0.005,nudge_x = 0.15)+
  geom_text(aes(label =signif(Estimate,2)), color = "black", size = 4, nudge_y=0)+
  theme_bw()+theme(panel.grid=element_blank(),legend.position = "none",
                   axis.text = element_text(colour = "black",size=12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(colour = "black",size=13),
                   #panel.grid.major.y = element_blank(),   ##取消y轴的主提示线
                   panel.grid.major.x = element_line(color="grey",linetype = 2),  ##设置z轴的主提示线为灰色的虚线
                   panel.grid.minor.x = element_line(color="grey",linetype = 2))+
  labs(y="Esitimate",title = "Cohort1")+
  coord_flip()

ggsave(p1,file="./bigmeta_phenotype.pdf",width = 6,height = 5)

#========RXX队列表型图============
sur_rxx <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/RenxxInfo_250310.xlsx")
sur_rxx$BMI[is.na(sur_rxx$BMI)]=median(sur_rxx$BMI[!is.na(sur_rxx$BMI)])
sur_rxx$log2AFP <- log2(sur_rxx$AFP)
sur_rxx$BCLC <- ifelse(sur_rxx$BCLC=="0",0,ifelse(sur_rxx$BCLC=="A",1,ifelse(sur_rxx$BCLC=="B",2,3)))
sur_rxx$`Tumor number` <- ifelse(sur_rxx$`Tumor number`==">3",4,sur_rxx$`Tumor number`)
sur_rxx$Edmondson <-  ifelse(sur_rxx$Edmondson=="1",1,ifelse(sur_rxx$Edmondson=="1-2",2,ifelse(sur_rxx$Edmondson=="2",3,
                                                                                               ifelse(sur_rxx$Edmondson=="2-3",4,ifelse(sur_rxx$Edmondson=="3",5,ifelse(sur_rxx$Edmondson=="3-4",6,7))))))


pheno <- c("TBil", "ALT", "AST","log2AFP","TumorSize","NLR","BCLC","TumorNumber","Edmondson")

#====连续型变量====
sur_rxx$BCLC <- as.numeric(sur_rxx$BCLC)
sur_rxx$`Tumor number` <- as.numeric(sur_rxx$`Tumor number` )
sur_rxx$Edmondson <- as.numeric(sur_rxx$Edmondson)
colnames(sur_rxx)[colnames(sur_rxx)=="Tumor size"] <- "TumorSize"
colnames(sur_rxx)[colnames(sur_rxx)=="Tumor number"] <- "TumorNumber"
colnames(sur_rxx)[colnames(sur_rxx)=="Asites"] <- "Ascites"

balance_rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon1.txt")
sur_rxx <- merge(sur_rxx,balance_rxx[,c("sample","balance_value")],by="sample")

sur_rxx <- sur_rxx[!is.na(sur_rxx$AFP)==T,]
res_cor <- list()
res_lm <- list()
res_lm_ad <- list()
for (i in 1:length(pheno)){
  tmp.data <- sur_rxx[,c(pheno[i],"balance_value")]
  tmp.res <-cor.test(tmp.data[,1],tmp.data$balance_value,method = "spearman") 
  res_cor[[i]] <- data.frame(rho=tmp.res[["estimate"]],pvalue=tmp.res[["p.value"]],pheno=pheno[i])
  #tmp.data[[pheno[i]]] <- rmOutlier(tmp.data[[pheno[i]]])
  tmp.data <- na.omit(tmp.data)
  mod1 <- lm(as.formula(paste0(pheno[i],"~balance_value")),data=sur_rxx)
  a <- as.data.frame(summary(mod1)[["coefficients"]])
  a<-  as.data.frame(a[-1,])
  a$pheno <- pheno[i]
  res_lm[[i]] <-a
  tmp.data <-  sur_rxx[,c(pheno[i],"balance_value","Age","Gender","BMI","Smoke","Alcohol")]
  tmp.data$Age <- as.numeric(tmp.data$Age)
  tmp.data$BMI <- as.numeric(tmp.data$BMI)
  tmp.data$Gender <- as.factor(tmp.data$Gender)
  tmp.data$Smoke <- as.factor(tmp.data$Smoke)
  tmp.data$Alcohol <- as.factor(tmp.data$Alcohol)
  #tmp.data[[pheno[i]]] <- rmOutlier(tmp.data[[pheno[i]]])
  tmp.data <- na.omit(tmp.data)
  mod2<-lm(as.formula(paste0(pheno[i],"~balance_value+Age+Gender+BMI+Smoke+Alcohol")),data = tmp.data)
  a <- as.data.frame(summary(mod2)[["coefficients"]])
  a<-  as.data.frame(a[2,])
  a$pheno <- pheno[i]
  res_lm_ad[[i]] <-a
}
res_cor <- do.call("rbind",res_cor)
res_lm <- do.call("rbind",res_lm)
res_lm_ad <- do.call("rbind",res_lm_ad)
res_lm$adj.p <- p.adjust(res_lm$`Pr(>|t|)`,"fdr")
res_lm_ad$adj.p <- p.adjust(res_lm_ad$`Pr(>|t|)`,"fdr")


#========分类型变量==============
sur_rxx$BCLC_group <- ifelse(sur_rxx$BCLC=="0",0,ifelse(sur_rxx$BCLC=="1",0,1))
sur_rxx$Child_Pugh <- ifelse(sur_rxx$Child_Pugh=="A",0,1)
sur_rxx$TumorNumber_group <- ifelse(sur_rxx$TumorNumber=="1",0,ifelse(sur_rxx$TumorNumber=="2",0,1))
sur_rxx$Edmondson_group <- ifelse(sur_rxx$Edmondson=="1",0,ifelse(sur_rxx$Edmondson=="1-2",0,ifelse(sur_rxx$Edmondson=="2",0,1)))
sur_rxx$`Lymphnode Metastasis` <- ifelse(sur_rxx$`Lymphnode Metastasis`=="Yes",1,0)
sur_rxx$`Distant Metastasis` <- ifelse(sur_rxx$`Distant Metastasis`=="Yes",1,0)

pheno <- c("BCLC_group","TumorNumber_group","`Lymphnode Metastasis`","`Distant Metastasis`", "PVTT",                          
           "HVTT" , "Ascites", "Child_Pugh","MVI" ,"Tumor_thrombus" , "Satellite_lesions" ,"Edmondson_group" )

#sur$`Lymphnode Metastasis` <- as.numeric(sur$`Lymphnode Metastasis`)
sur_rxx$HVTT <- factor(sur_rxx$HVTT,levels=c(0,1))
sur_rxx$PVTT <- factor(sur_rxx$PVTT,levels = c(0,1))
sur_rxx$Ascites <- factor(sur_rxx$Ascites,levels = c(0,1))
sur_rxx$Child_Pugh <- factor(sur_rxx$Child_Pugh,levels = c(0,1))
sur_rxx$MVI <- factor(sur_rxx$MVI,levels = c(0,1))
sur_rxx$Tumor_thrombus <- factor(sur_rxx$Tumor_thrombus,levels = c(0,1))
sur_rxx$Satellite_lesions <- factor(sur_rxx$Satellite_lesions,levels = c(0,1))
sur_rxx$BCLC_group <- factor(sur_rxx$BCLC_group,levels = c(0,1))
sur_rxx$TumorNumber_group <- factor(sur_rxx$TumorNumber_group,levels = c(0,1))
sur_rxx$`Lymphnode Metastasis` <- factor(sur_rxx$`Lymphnode Metastasis`,levels = c(0,1))
sur_rxx$`Distant Metastasis` <- factor(sur_rxx$`Distant Metastasis`,levels = c(0,1))
sur_rxx$Edmondson_group <- factor(sur_rxx$Edmondson_group,levels = c(0,1))

res_glm <- list()
res_glm_ad <- list()
for (i in 1:length(pheno)){
  glm.md <- glm(as.formula(paste0(pheno[i],"~balance_value")),data=sur_rxx, family = "binomial")
  a <- summary(glm.md)$coefficients %>% as.data.frame()
  a<-  as.data.frame(a[-1,])
  a$pheno <- pheno[i]
  res_glm[[i]] <-a
  mod2<-glm(as.formula(paste0(pheno[i],"~balance_value+Age+Gender+BMI+Alcohol+Smoke")),data = sur_rxx,family = "binomial")
  a <- summary(mod2)$coefficients %>% as.data.frame()
  a<-  as.data.frame(a[2,])
  a$pheno <- pheno[i]
  res_glm_ad[[i]] <-a
}
res_glm <- do.call("rbind",res_glm)
res_glm_ad <- do.call("rbind",res_glm_ad)
res_glm$adj.p <- p.adjust(res_glm$`Pr(>|z|)`,"fdr")
res_glm_ad$adj.p <- p.adjust(res_glm_ad$`Pr(>|z|)`,"fdr")

#save(sur,file="rxx_sur250107.RData")
#==========画表型棒棒糖图====================
colnames(res_lm_ad)[4] <- "pvalue"
colnames(res_glm_ad)[4] <- "pvalue"
forplot2 <- rbind(res_lm_ad[,c("Estimate","pvalue", "pheno", "adj.p")],res_glm_ad[,c("Estimate","pvalue", "pheno", "adj.p")])
forplot2 <- forplot2[!forplot2$pheno%like%"group",]
forplot2 <- forplot2[!forplot2$pheno%like%"Metastasis",]
forplot2 <- forplot2[!forplot2$Estimate>0,]
forplot2 <- forplot2%>%arrange(Estimate)
forplot2$pheno <- factor(forplot2$pheno,levels = rev(forplot2$pheno))
forplot2$stars <- cut(forplot2$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))

#==pheno保持一致====
pheno <- intersect(forplot1$pheno,forplot2$pheno)
forplot1 <- forplot1[forplot1$pheno%in%pheno,]
#forplot1$pheno <- factor(forplot1$pheno,levels=forplot1$pheno)
forplot2 <- forplot2[forplot2$pheno%in%pheno,]
forplot2 <- forplot2[order(forplot2$Estimate,decreasing = T),]


p1 <- ggplot(forplot1, aes(x = pheno, y = Estimate)) +
  geom_segment( aes(x = pheno, xend = pheno, y = 0, yend = Estimate),color = "grey40")+
  scale_y_reverse(limits = c(0, -0.53), expand = c(0, 0))+ #控制线段的参数，见下
  geom_point(color= "#b96c93",size = 8) + 
  geom_text(aes(label=stars), fontface="bold", size=6 ,nudge_y=0,nudge_x = 0.1,color="white")+
  #geom_text(aes(label =signif(Estimate,2)), color = "black", size = 4, nudge_y=0)+
  theme_bw()+theme(panel.grid=element_blank(),legend.position = "none",
                   axis.text = element_text(colour = "black",size=12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(colour = "black",size=13),
                   #panel.grid.major.y = element_blank(),   ##取消y轴的主提示线
                   panel.grid.major.x = element_line(color="grey",linetype = 2),  ##设置z轴的主提示线为灰色的虚线
                   panel.grid.minor.x = element_line(color="grey",linetype = 2))+
  labs(y="Esitimate",title = "Cohort1")+
  coord_flip()

p2 <- ggplot(forplot2, aes(x = pheno, y = Estimate)) +
  geom_segment( aes(x = pheno, xend = pheno, y = 0, yend = Estimate),color = "grey40")+
  scale_y_reverse(limits = c(0, -0.28), expand = c(0, 0))+ #控制线段的参数，见下
  geom_point(color=  "#b96c93",size = 8) + 
  geom_text(aes(label=stars), fontface="bold", size=6 ,nudge_y=0.005,nudge_x = 0.15,color="white")+
  #geom_text(aes(label =signif(Estimate,2)), color = "black", size = 4, nudge_y=0)+
  theme_bw()+theme(panel.grid=element_blank(),legend.position = "none",
                   axis.text = element_text(colour = "black",size=12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(colour = "black",size=13),
                   #panel.grid.major.y = element_blank(),   ##取消y轴的主提示线
                   panel.grid.major.x = element_line(color="grey",linetype = 2),  ##设置z轴的主提示线为灰色的虚线
                   panel.grid.minor.x = element_line(color="grey",linetype = 2))+
  labs(y="Esitimate",title = "Cohort2")+
  coord_flip()

pdf("./Fig1G260123.pdf",width = 10,height = 5)
p1+p2+plot_layout(heights = c(1,1),widths = c(1,1),ncol = 2,nrow = 1)
dev.off()

save(forplot1,forplot2,file="./Fig1G.RData")
load("./Fig1G.RData")
forplot1$pheno <- factor(forplot1$pheno,levels=rev(c("AST","TumorSize","TBil","NLR","Ascites","Child_Pugh","PVTT","MVI",            
                                                     "BCLC","Satellite_lesions","Edmondson")))
forplot2$pheno <- factor(forplot2$pheno,levels = rev(c("AST","TumorSize","TBil","NLR","Ascites","Child_Pugh","PVTT","MVI",            
                                                       "BCLC","Satellite_lesions","Edmondson")))
library(openxlsx)
write.xlsx(list(cohort1 = forplot1,cohort2 = forplot2),
           file = "/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/Fig1I_sup_data260123.xlsx",
           rowNames = FALSE)


