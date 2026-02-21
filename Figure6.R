rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/FMT/MGS/")
#====library====
library(data.table)
library(pairwiseAdonis)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidydr)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(ggh4x)
library(ape)
library(multcomp)
library(extrafont)
library(patchwork)
library(vegan)
library(ggtern)
library(foreach)
library(ggsci)
library(scales)
source("/groups/ProHuShiX/home/xiashufen/bigmeta/16s/functions.R")

data <- read.table("/groups/microbiome/home/share/Out/Mouse/metagenomic/bigmeta2506/metaphlan/FMT_profiled_metagenome.txt",header = T,check.names = F)
colnames(data) <- gsub(".profiled_metagenome","",colnames(data))
data <- data %>% column_to_rownames(var="clade_name")

sample <- read.table("./input/178quality_sample.txt")

data_f <- data[,colnames(data)%in%sample$sample]

species <- data_f %>% 
  dplyr::filter(str_detect(rownames(data),"s__.*") & !str_detect(rownames(data),"t__.*")) 
species <- species %>% 
  dplyr::mutate(species=str_extract(rownames(species),"s__.*"))
rownames(species) <- NULL
species <- species %>% column_to_rownames(var="species")
rownames(species) <- gsub("s__","",rownames(species)) # 875 species

species_raw <- as.data.frame(t(species))

species <- species[rowSums(species!=0) > (ncol(species)*0.1),] # 0.1: 567 species # 0.2: 360 taxa
species <- apply(species,2,function(x){
  y <- x/sum(x)
  return(y)
})
species <- as.data.frame(t(species))
range(colSums(species))
# 3.724323e-04 4.510877e+01
# 0.001334279 45.238404964

# write.table(species,"/groups/ProHuShiX/home/xiashufen/bigmeta/FMT/MGS/input/metaphlan_0.1trim.txt")

## group info
Group_time <- data.frame(sample=rownames(species))
Group_time$group_t <- sub(".*_","",Group_time$sample)
table(Group_time$group)
# base  abx  fmt  end 
# 40   48   50   40 

Group_time$num <- sub("_.*","",Group_time$sample)
Group_time <- Group_time %>% 
  mutate(num = case_when(
    grepl("A",num) ~ sub("A","A0",num),
    nchar(num) == 1 ~ paste0("A00",num),
    nchar(num) == 2 ~ paste0("A0",num),
    nchar(num) == 3 ~ paste0("A",num),
    TRUE ~ num
  ))

Group_FMT <- read_excel("../FMT_group.xlsx",sheet = "Sheet1")
Group_FMT$group <- gsub("-"," ",Group_FMT$group)

Group_time$group_bs <- Group_FMT$group[match(Group_time$num,Group_FMT$sample)]
table(Group_time$group_t,Group_time$group_bs)

# 平均base在5个组之间的分布
set.seed(666)  
base <- which(Group_time$group_t == "base")
groups <- unique(Group_time$group_bs) 
assigned_groups <- rep(groups, length.out = length(base))
assigned_groups <- sample(assigned_groups)  # 随机打乱顺序
Group_time$group_bs[base] <- assigned_groups

table(Group_time$group_t,Group_time$group_bs)

# writexl::write_xlsx(Group_time,"/groups/ProHuShiX/home/xiashufen/bigmeta/FMT/178group.xlsx")

#====colorset====
col0 <- c("#2482BC","#e73133","#5ba787","#b96c93","#fcd364")
col1 <- c("#7c9680","#7ca6be","#fcd364","#fe9014") # base/abx/fmt/end
col2 <- c("#f99170","#6e7ca5","#b66d32") # donor/pre/post
col3 <- c("#5ba787","#b96c93","#2482BC") # score high/score low/pbs

#====alpha diversity====
shannon <- data.frame(shannon=vegan::diversity(species,index = "shannon"))

Group_time$shannon <- shannon$shannon[match(Group_time$sample,rownames(shannon))]

i="BS-Low CTRL"

for (i in unique(Group_time$group_bs)){
  print(i)
  tmp.shannon <- Group_time[Group_time$group_bs==i,]
  tmp.shannon$group_t <- factor(tmp.shannon$group_t,levels = c("base","abx","fmt","end"))
  comparison <- combn(levels(tmp.shannon$group_t),2,simplify = FALSE)
  print(table(tmp.shannon$group_t))
  p <- ggplot(tmp.shannon,aes(x= group_t, y= shannon,color=group_t))+
    geom_boxplot( )+ 
    geom_jitter(width = 0.2,aes(col = group_t))+ scale_color_manual(values=col1)+
    stat_compare_means(comparisons = comparison)+
    theme_classic()+theme(legend.position = "none",
                          axis.text = element_text(colour = "black",size=12),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(colour = "black",size=13))+
    ylab("Shannon")+ggtitle(i)
  print(p)
  #ggsave(paste0("/groups/ProHuShiX/home/xiashufen/bigmeta/FMT/MGS/output_0.2/",i,"alpha.pdf"),p,width = 6,height = 6)
}


#====beta diveristy====
# prevelance 0.1
species_raw <- species

# prevelance 0.2
# species_raw <- species_raw[,colSums(species_raw!=0) > (nrow(species_raw)*0.2)]
# species_raw <- apply(species_raw,1,function(x){
#   y <- x/sum(x)
#   return(y)
# })
# species_raw <- as.data.frame(t(species_raw))

species_raw <- species_raw[rownames(species_raw)%in%Group_time$sample[Group_time$group_t%in%c("base","fmt")],]

donor <- read_excel("../FMT_group.xlsx",sheet = "donor")
# HCC donor from FMT
FMT <- read.table("/groups/microbiome/home/share/Out/human/metagenomic/HCC/BigMeta2025/250219FMT/metaphlan/FMT_metagenome.txt",header=T,check.names = F)
colnames(FMT) <- gsub("\\.pro.*$","",colnames(FMT))
FMT <- FMT %>% dplyr::filter(str_detect(clade_name,"s__.*") & ! str_detect(clade_name,"t__.*")) %>% 
  dplyr::mutate(species=str_extract(clade_name,"s__.*"))
FMT <- FMT %>% column_to_rownames(var="species")
FMT$clade_name <- NULL
rownames(FMT) <- gsub("s__","",rownames(FMT))
FMT <- FMT[rowSums(FMT!=0) > (ncol(FMT)*0.1),]
FMT <- apply(FMT,2,function(x){
  y <- x/sum(x)
  return(y)
})
FMT <- as.data.frame(t(FMT))

FMT <- FMT[rownames(FMT)%in%donor$sample[donor$group%like%"HCC"],]

# CRTL donor from newMGS
newMGS <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/newMGS_raw_species.txt",check.names = F)
newMGS <- newMGS[rowSums(newMGS!=0) > (ncol(newMGS)*0.1),]
newMGS <- apply(newMGS,2,function(x){
  y <- x/sum(x)
  return(y)
})
newMGS <- as.data.frame(t(newMGS))
newMGS <- newMGS[rownames(newMGS)%in%donor$sample[donor$group%like%"CTRL"],]

length(union(colnames(FMT),union(colnames(newMGS),colnames(species_raw))))# 3035 taxa # 0.1: 995 taxa # 0.2: 633 # 0.3: 316

all.taxa <- bind_rows(FMT,newMGS,species_raw)
all.taxa[is.na(all.taxa)] <- 0

all.taxa_clr <- transform_and_filter_taxa(all.taxa,samples_row = T,method = "clr",missing_filter = 0)

test <- all.taxa[,!colnames(all.taxa)%in%colnames(all.taxa_clr),drop=F]
range(colSums(test))

all <- c("BS High CTRL","BS Low CTRL","BS High HCC","BS Low HCC")

i="BS Low HCC"
# 分开四个组画
for (i in all) {
  
  donor_sample <- donor$sample[donor$group==i]
  pre_sample <- Group_time$sample[Group_time$group_bs==i&Group_time$group_t=="base"]
  post_sample <- Group_time$sample[Group_time$group_bs==i&Group_time$group_t=="fmt"]
  
  tmp.taxa <- all.taxa[rownames(all.taxa)%in%c(donor_sample,pre_sample,post_sample),]
  
  tmp.group <- data.frame(sample=rownames(tmp.taxa))
  tmp.group$batch <- as.factor(ifelse(tmp.group$sample%in%donor_sample,1,2))
  
  data_adj <- matrix(nrow = nrow(tmp.taxa), ncol = ncol(tmp.taxa))
  colnames(data_adj) <- colnames(tmp.taxa)
  rownames(data_adj) <- rownames(tmp.taxa)
  
  # 对每个物种（列）线性回归取残差
  for (n in seq_len(ncol(tmp.taxa))) {
    fit <- lm(tmp.taxa[, n] ~ tmp.group$batch)
    data_adj[, n] <- residuals(fit)
  }
  
  set.seed(666)
  dist_BC <- vegdist(data_adj,method = "euclidian")
  pcoa <- pcoa(dist_BC,correction = "none", rn = NULL) 
  PC1 <- pcoa$vectors[,1]
  PC2 <- pcoa$vectors[,2]
  pcoadata <- data.frame(sample=rownames(pcoa$vectors),PC1=PC1,PC2=PC2)
  pcoadata$group <- Group_time$group_t[match(pcoadata$sample,Group_time$sample)]
  pcoadata$group[is.na(pcoadata$group)] <- "donor"
  pcoadata$group <- gsub("base","pre FMT",pcoadata$group)
  pcoadata$group <- gsub("fmt","post FMT",pcoadata$group)
  
  dist_BC <- as.matrix(dist_BC)
  pcoadata$group <- factor(pcoadata$group,levels = c("donor","pre FMT","post FMT"))
  
  #Adonis test
  res <- pairwise.adonis(dist_BC,factors = pcoadata$group,perm = 999)
  
  
  p <- ggplot(pcoadata, aes(PC1, PC2, color = group)) +
    geom_point(size = 2, alpha = 0.8) +
    stat_ellipse(aes(PC1, PC2), level = 0.95, linetype = 2) +
    #xlim(-0.6, 0.5) + ylim(-0.7, 0.7) +
    labs(
      x = paste0("PCoA1 ( ", floor(pcoa$values$Relative_eig[1] * 100), "% )"),
      y = paste0("PCoA2 ( ", floor(pcoa$values$Relative_eig[2] * 100), "% )"),
      color = "Group"
    ) +
    scale_color_manual(values = col2) +
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
      panel.background = element_rect(fill = NA)) +
    annotate("text", x = -0.15, y = 0.4, 
             label = paste("PERMANOVA:\n", res$pairs[1],"P =",res$p.value[1],"R² =",signif(res$R2[1],2),
                           "\n",res$pairs[2],"P =",res$p.value[2],"R² =",signif(res$R2[2],2)))+
    ggtitle(paste0(i," euclidian"))
  print(p)
  #ggsave(paste0("/groups/ProHuShiX/home/xiashufen/bigmeta/FMT/MGS/output/",i,"noclr.beta.pdf"),p,device = cairo_pdf,width = 8,height = 8)
  
}


#====Fig6B_beta diversity====
#合并四个组并手动加上PBS
donor_sample <- donor$sample[donor$group%in%all]

all1 <- c(all,"PBS")
pre_sample <- Group_time$sample[Group_time$group_bs%in%all1&Group_time$group_t=="base"]
post_sample <- Group_time$sample[Group_time$group_bs%in%all1&Group_time$group_t=="fmt"]

tmp.taxa <- all.taxa_clr[rownames(all.taxa_clr)%in%c(donor_sample,pre_sample,post_sample),]

tmp.group <- data.frame(sample=rownames(tmp.taxa))
tmp.group$batch <- as.factor(ifelse(tmp.group$sample%in%donor_sample,1,2))

data_adj <- matrix(nrow = nrow(tmp.taxa), ncol = ncol(tmp.taxa))
colnames(data_adj) <- colnames(tmp.taxa)
rownames(data_adj) <- rownames(tmp.taxa)

# 对每个物种（列）线性回归取残差
for (n in seq_len(ncol(tmp.taxa))) {
  fit <- lm(tmp.taxa[, n] ~ tmp.group$batch)
  data_adj[, n] <- residuals(fit)
}

set.seed(666)
dist_BC <- vegdist(data_adj,method = "euclidian")
pcoa <- pcoa(dist_BC,correction = "none", rn = NULL) 
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]
pcoadata <- data.frame(sample=rownames(pcoa$vectors),PC1=PC1,PC2=PC2)
pcoadata$group <- Group_time$group_t[match(pcoadata$sample,Group_time$sample)]
pcoadata$group[is.na(pcoadata$group)] <- "donor"
pcoadata$group <- gsub("base","pre FMT",pcoadata$group)
pcoadata$group <- gsub("fmt","post FMT",pcoadata$group)

dist_BC <- as.matrix(dist_BC)
# pcoadata$group <- factor(pcoadata$group,levels = c("donor","pre FMT","post FMT"))

pcoadata$group_d <- NA
pcoadata$group_d[pcoadata$group=="donor"] <- donor$group[match(pcoadata$sample[pcoadata$group=="donor"],donor$sample)]
pcoadata$group_d[pcoadata$group!="donor"] <- Group_time$group_bs[match(pcoadata$sample[pcoadata$group!="donor"],Group_time$sample)]
pcoadata$group_d <- as.factor(pcoadata$group_d)

pcoadata$color_group <- ifelse(pcoadata$group_d=="PBS","PBS",pcoadata$group)
pcoadata$color_group <- factor(pcoadata$color_group,levels = c("donor","pre FMT","post FMT","PBS"))
col_all <- c("donor" = "#1b9e77",
             "pre FMT" = "#d95f02",
             "post FMT" = "#7570b3",
             "PBS" = "#808080")
pcoadata$group <- factor(pcoadata$group,levels = c("donor","pre FMT","post FMT"))
# 点形状：5个实验组
shape_all <- c("BS High CTRL"=16, "BS Low CTRL"=17, 
               "BS High HCC"=21, "BS Low HCC"=24,
               "PBS"=15)  # PBS特殊标记

str(pcoadata)
p<-ggplot(pcoadata, aes(PC1, PC2, color = color_group, shape = group_d)) +
  geom_point(size = 2, alpha = 0.8) +
  
  # 只对4个实验组画椭圆
  stat_ellipse(
    data = subset(pcoadata, color_group %in% c("donor", "pre FMT", "post FMT")),
    aes(group = group, color = color_group),
    level = 0.95,
    linetype = 2,
    size = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = col_all) +
  scale_shape_manual(values = shape_all) +
  labs(
    x = paste0("PCoA1 (", floor(pcoa$values$Relative_eig[1] * 100), "%)"),
    y = paste0("PCoA2 (", floor(pcoa$values$Relative_eig[2] * 100), "%)"),
    color = "Donor type",
    shape = "Experiment Group"
  ) +
  theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = c(.9, .7),  # 修改 ggplot2 3.5.0 兼容性问题
        legend.title = element_text(size = 12, family = "bold", color = "black"),
        legend.text = element_text(size = 9, family = "sans", color = "black"),
        legend.key = element_blank(), 
        axis.text.x = element_text(size = 10, family = "sans", color = "black"),
        axis.text.y = element_text(size = 10, family = "sans", color = "black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


#Adonis test
tmp.dist_BC <- dist_BC[colnames(dist_BC)%in%pcoadata$sample[pcoadata$group_d!="PBS"],rownames(dist_BC)%in%pcoadata$sample[pcoadata$group_d!="PBS"]]
tmp.pcoadata <- pcoadata[pcoadata$group_d!="PBS",]
all(rownames(tmp.dist_BC)==tmp.pcoadata$sample)
set.seed(666)
# res <- pairwise.adonis2(tmp.dist_BC~group,data = tmp.pcoadata,nperm = 999)
# res$`donor_vs_pre FMT`
# res$`donor_vs_post FMT`

res <- pairwise.adonis(tmp.dist_BC,factors = tmp.pcoadata$group,perm=999)

p <- p+annotate("text", x = 0.1, y = 0.3, 
                label = paste("PERMANOVA:\n", res$pairs[1],"P =",res$p.value[1],"R² =",signif(res$R2[1],2),
                              "\n",res$pairs[2],"P =",res$p.value[2],"R² =",signif(res$R2[2],2)))
p
ggsave("./output/Beta_diversity_all_0.1_noclr.pdf",p,device = cairo_pdf,width = 6.5,height = 5)

#====Fig6C_distance=====
ctrl.donor <- newMGS[rownames(newMGS)%in%donor$sample,]
hcc.donor <- FMT[rownames(FMT)%in%donor$sample,]
mouse <- species_raw[rownames(species_raw)%in%Group_time$sample[Group_time$group_t%in%c("base","fmt")],]

length(union(colnames(ctrl.donor),union(colnames(hcc.donor),colnames(mouse))))

data_RDA <- bind_rows(ctrl.donor,hcc.donor,mouse)
data_RDA[is.na(data_RDA)] <- 0
data_RDA <- as.data.frame(t(data_RDA))

data_RDA_clr <- transform_and_filter_taxa(data_RDA,samples_row = F,method = "clr",missing_filter = 0)
dist_EU <- vegdist(data_RDA_clr,method = "euclidean")
distmat <- as.matrix(dist_EU)

# 按照实验分组提取donor和对应的pre FMT/post FMT
i="BS Low CTRL" 
Distance=foreach(i=unique(Group_time$group_bs)[1:4],.combine = rbind) %do% {
  cat(i,"\n")
  d_sample = donor$sample[donor$group==i] # 5 donors
  pre_sample = Group_time$sample[Group_time$group_t=="base"&Group_time$group_bs==i] # 8 base
  post_sample = Group_time$sample[Group_time$group_t=="fmt"&Group_time$group_bs==i] # 10 fmt
  
  # pre FMT提取一个5×8的矩阵
  # post FMT提取一个5×10的矩阵
  pre_dist=distmat[rownames(distmat)%in%d_sample,colnames(distmat)%in%pre_sample]
  post_dist=distmat[rownames(distmat)%in%d_sample,colnames(distmat)%in%post_sample]
  
  pre_distance <- as.data.frame(pre_dist) %>%
    tibble::rownames_to_column("donor") %>%
    pivot_longer(-donor, names_to = "mouse", values_to = "distance")
  pre_distance$donor <- as.character(pre_distance$donor)
  pre_distance$mouse <- as.character(pre_distance$mouse)
  
  post_distance <- as.data.frame(post_dist) %>%
    tibble::rownames_to_column("donor") %>%
    pivot_longer(-donor, names_to = "mouse", values_to = "distance")
  post_distance$donor <- as.character(post_distance$donor)
  post_distance$mouse <- as.character(post_distance$mouse)
  
  tmp_distance = rbind(pre_distance,post_distance)
  
  tmp_distance$group_bs = as.character(i)
  tmp_distance$group_t = ifelse(tmp_distance$mouse%like%"base","preFMT","postFMT")
  
  as.data.frame(tmp_distance)
}


res <- wilcox.test(Distance$distance~Distance$group_t)
pvalue <- res$p.value # 1.25e-56

Distance$group_t <- factor(Distance$group_t,levels = c("preFMT","postFMT"))
ggplot(Distance, aes(x=group_t, y=distance, color=group_t)) + 
  geom_violin(trim=FALSE,aes(fill=group_t,alpha=0.8))+
  geom_boxplot(width=0.15, fill="white")+
  theme_classic()+
  stat_compare_means(comparisons=list(c("preFMT","postFMT")))+
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "none")+
  scale_color_manual(values = col2[2:3])+
  scale_fill_manual(values=col2[2:3])+ylab("Heterogeneity")+ggtitle("Merged groups")
#ggsave("/groups/ProHuShiX/home/xiashufen/bigmeta/FMT/MGS/output/PrePostFMT_distance_merged.pdf",width = 5,height = 6)


i="BS Low CTRL" 
for (i in unique(Distance$group_bs)) {
  tmp.distance = Distance[Distance$group_bs==i,]
  tmp.distance$group_t <- factor(tmp.distance$group_t,levels = c("preFMT","postFMT"))
  print(i)
  print(table(tmp.distance$group_t))
  p <- ggplot(tmp.distance,aes(x= group_t, y= distance,color=group_t))+
    geom_boxplot( )+ 
    geom_jitter(width = 0.2,aes(col = group_t))+ scale_color_manual(values=col2[2:3])+
    stat_compare_means(comparisons = list(c("preFMT","postFMT")))+
    theme_classic()+theme(legend.position = "none",
                          axis.text = element_text(colour = "black",size=12),
                          axis.title.x = element_blank(),
                          axis.title.y = element_text(colour = "black",size=13))+
    ylab("Heterogeneity")+ggtitle(i)
  print(p)
  #ggsave(paste0("/groups/ProHuShiX/home/xiashufen/bigmeta/FMT/MGS/output/",i,"dist.pdf"),p,width = 5,height = 6)
}

#=====Fig6I_pancancer_meta====
load("./for_meta3.RData")
load("../2017Neoplasia/CMD_28923537_2017Neoplasia.RData")

tmp.mela.26 <- balance_df_mela[,c("sample","group","adj_bs")]
tmp.mela.165 <- balance_df_mela165[,c("sample","group3","adj_bs")]
tmp.nsclc.338 <- balance_df_NSCLC[,c("sample","group2","adj_bs")]
tmp.mela.39 <- data.frame(sample=rownames(balance_df4_0.1),group=factor(ifelse(balance_df4_0.1$group%in%c("CR","PR","SD"),"Nonprogressed","Progressed"),levels = c("Nonprogressed","Progressed")),
                          adj_bs=balance_df4_0.1$adj_score)

tmp.mela.26$group <- factor(ifelse(tmp.mela.26$group=="Non progressed","Nonprogressed","Progressed"),levels = c("Nonprogressed","Progressed"))
tmp.mela.165$group3 <- factor(ifelse(tmp.mela.165$group3=="R","Nonprogressed","Progressed"),levels = c("Nonprogressed","Progressed"))
tmp.nsclc.338$group2 <- factor(ifelse(tmp.nsclc.338$group2=="CR/PR/SD","Nonprogressed","Progressed"),levels = c("Nonprogressed","Progressed"))

colnames(tmp.mela.165)[2] <- "group";colnames(tmp.nsclc.338)[2] <- "group"

tmp.mela.26$study <- "McCulloch，2022"
tmp.nsclc.338$study <- "Derosa,2022"
tmp.mela.165$study <- "A.Lee,2022"
tmp.mela.39$study <- "E.Frankel,2017"

all_data <- rbind(tmp.mela.26,tmp.mela.39,tmp.mela.165,tmp.nsclc.338)

all_data$study <- as.factor(all_data$study)

summary_df <- all_data %>%
  group_by(study, group) %>%
  summarise(
    n = n(),
    mean = mean(adj_bs, na.rm = TRUE),
    sd = sd(adj_bs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = c(n, mean, sd),
    names_glue = "{.value}.{group}"
  ) %>%
  dplyr::rename(
    n.e = n.Nonprogressed,
    mean.e = mean.Nonprogressed,
    sd.e = sd.Nonprogressed,
    n.c = n.Progressed,
    mean.c = mean.Progressed,
    sd.c = sd.Progressed
  )

meta_result <- metacont(
  n.e, mean.e, sd.e,
  n.c, mean.c, sd.c,
  data = summary_df,
  studlab = study,
  sm = "SMD",
  method.smd = "Hedges"
)

meta_result$pval.fixed
summary(meta_result)
metainf(meta_result)

z <- meta_result$TE.random / meta_result$seTE.random
p <- 2 * (1 - pnorm(abs(z)))

all.equal(p,meta_result$pval.random)

forest(meta_result)
z <- meta_result$zval.random
p <- meta_result$pval.random
forest(meta_result,comb.random = TRUE,text.random = sprintf("Test for overall effect: Z = %.2f, P = %s",z,p))


pdf("/groups/ProHuShiX/home/xiashufen/bigmeta/BS_pancancer/for_meta/4cohort_metaforest_0.1trim260119.pdf",width = 15,height = 5)
# forest(meta_result)
forest(meta_result,comb.random = TRUE,text.random = sprintf("Test for overall effect: Z = %.2f, P = %s",z,p))
dev.off()




