rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/Fig4/")

#======Fig4A_score_betadiversity==============
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

#====data prepare: R0_data_prepare=====
load("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/input_data/sc_BIGMeta250618.RData")
Cluster <- read.table("../Final_output/Balance.score1.txt",header=T)
sc_hcc <- sc[which(Cluster$group[match(rownames(sc),Cluster$Row.names)]=="HCC"),]
beta_diversity=vegdist(sc_hcc,method = "euclidean")
pcoa_analysis=cmdscale(beta_diversity,k=4,eig=TRUE)
pc<- round(pcoa_analysis$eig/sum(pcoa_analysis$eig)*100,digits=2)#这一步是算解释度的，pc的第一和第二位分别代表X轴和y轴的解释度
pcoa_analysis <- as.data.frame(pcoa_analysis[["points"]])
pcoa_analysis <- merge(pcoa_analysis,Cluster[,c("Row.names","balance_value","Cluster","group")],by.x="row.names",by.y="Row.names",all=F)

cov.data <- Cluster[,c("Row.names","balance_value","Cluster","group")]
cov.data <- cov.data[match(row.names(sc_hcc),cov.data$Row.names),]
cov.data$Cluster <- as.factor(cov.data$Cluster)
cov.data$Enterotype <- as.factor(ifelse(cov.data$Cluster==4,"ETE","ETnonE"))
cov.data$Enterotype <- factor(cov.data$Enterotype,levels=c("ETnonE","ETE"))

otu.adonis=adonis2(sc_hcc ~ cov.data$balance_value, permutations = 999, method = "euclidian")
range(pcoa_analysis$balance_value)
p1 <- ggplot(pcoa_analysis,aes(V1, V2, color = balance_value))+
  geom_point(size = 2)+#ggtitle("France")+
  #stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)
  xlab(paste0("PCoA1(",pc[1],"%)"))+ylab(paste0("PCoA2(",pc[2],"%)"))+
  labs(color="balance_score") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.6,9.5))+
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
  annotate("text",x=-23,y=23,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis[1,5],2),
                                         "\nR2=",signif(otu.adonis[1,3],2)))
p1
# dev.off()
#load("./Fig3A_bigmeta.RData")
save(sc,cov.data,otu.adonis,file = "../Fig4/Fig4S_bigmeta_all.RData")

load("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/input_data/sc_renxx250618.RData")
cli_rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon1.txt")
cli_rxx <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/renxx_BS&shannon_0.1pre.txt")
colnames(cli_rxx)[2] <- "balance_score"
beta_diversity=vegdist(sc_rxx,method = "euclidean")
pcoa_analysis=cmdscale(beta_diversity,k=4,eig=TRUE)
pc<- round(pcoa_analysis$eig/sum(pcoa_analysis$eig)*100,digits=2)#这一步是算解释度的，pc的第一和第二位分别代表X轴和y轴的解释度

pcoa_analysis <- as.data.frame(pcoa_analysis[["points"]])
pcoa_analysis <- merge(pcoa_analysis,cli_rxx[,c("sample","balance_score"),drop=F],by.x="row.names",by.y="sample",all=F)

cov.data <- cli_rxx[match(row.names(sc_rxx),cli_rxx$sample),]

set.seed(111)
otu.adonis=adonis2(sc_rxx ~ cov.data$balance_score, permutations = 999, method = "euclidian")
range(pcoa_analysis$balance_score) # -4.063184  9.420987 # 0.1: -5.208827  9.412680
p2 <- ggplot(pcoa_analysis,aes(V1, V2, color = balance_score))+
  geom_point(size = 2)+#ggtitle("France")+
  #stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)
  xlab(paste0("PCoA1(",pc[1],"%)"))+ylab(paste0("PCoA2(",pc[2],"%)"))+
  labs(color="balance_score") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.6,9.5))+
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
  annotate("text",x=-25,y=20,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis[1,5],2),
                                         "\nR2=",signif(otu.adonis[1,3],2)))

library(patchwork)
pdf("../Fig4/Fig4A_score_betadiversity_HCC+HC250709.pdf",width = 10,height = 5.5)
(p1+ggtitle("Cohort1"))+(p2+ggtitle("Cohort2"))+plot_layout(heights = c(1,1),widths = c(1,1),ncol = 2,nrow = 1)
dev.off()
save(sc,otu.adonis,file="../Fig4/Fig4A_rxx.RData")

#====Fig4S bigmeta_HC+HCC betadiverity in group&bs====
beta_diversity=vegdist(sc,method = "euclidean")
pcoa_analysis=cmdscale(beta_diversity,k=4,eig=TRUE)
pc<- round(pcoa_analysis$eig/sum(pcoa_analysis$eig)*100,digits=2)#这一步是算解释度的，pc的第一和第二位分别代表X轴和y轴的解释度
pcoa_analysis <- as.data.frame(pcoa_analysis[["points"]])
pcoa_analysis <- merge(pcoa_analysis,Cluster[,c("Row.names","balance_value","Cluster","group")],by.x="row.names",by.y="Row.names",all=F)

cov.data <- Cluster[,c("Row.names","balance_value","Cluster","group")]
cov.data <- cov.data[match(row.names(sc),cov.data$Row.names),]
cov.data$Cluster <- as.factor(cov.data$Cluster)
cov.data$Enterotype <- as.factor(ifelse(cov.data$Cluster==4,"ETE","ETnonE"))
cov.data$Enterotype <- factor(cov.data$Enterotype,levels=c("ETnonE","ETE"))

otu.adonis1=adonis2(sc ~ cov.data$group, permutations = 999, method = "euclidian")
otu.adonis2=adonis2(sc ~ cov.data$balance_value, permutations = 999, method = "euclidian")

col1 <- c("#2482BC","#f99170","#e73133")
pcoa_analysis$group <- factor(pcoa_analysis$group,levels = c("Health","HCC"))

p1 <- ggplot(pcoa_analysis,aes(V1, V2, color = group))+
  geom_point(size = 2,alpha = 0.8)+
  stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  #xlim(-0.6,0.6)+ylim(-0.5,0.5)+
  labs(x = paste0("PCoA1(",pc[1],"%)"),
       y = paste0("PCoA2(",pc[2],"%)"),color="group") +
  scale_color_manual(values =c(col1[1], col1[3]))+
  theme_classic()+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position ="top",
        legend.key = element_blank(), 
        # axis.text.x = element_text(size = 10,family ="sans",color="black"),
        #axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))+
  annotate("text",x=-25,y=20,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis1[1,5],2),
                                         "\nR2=",signif(otu.adonis1[1,3],2)))
p1

p2 <- ggplot(pcoa_analysis,aes(V1, V2, color = balance_value))+
  geom_point(size = 2)+#ggtitle("France")+
  #stat_ellipse(aes(V1,V2) ,level=0.95,linetype = 2)+
  # xlim(-0.6,0.6)+ylim(-0.5,0.5)
  xlab(paste0("PCoA1(",pc[1],"%)"))+ylab(paste0("PCoA2(",pc[2],"%)"))+
  labs(color="balance_score") +
  scale_color_gradient2(low="#8B516E",mid="#FEF8A9",high= "#46714C",limit=c(-6.6,9.5))+
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
  annotate("text",x=-25,y=20,label=paste("PERMANOVA:\npvalue=",signif(otu.adonis2[1,5],2),
                                         "\nR2=",signif(otu.adonis2[1,3],2)))

p2

pdf("../Fig4/Fig4S_betadiversity_bigmeta_HCC+HC250709.pdf",width = 10,height = 5.5)
(p1+ggtitle("Group"))+(p2+ggtitle("Balance score"))+plot_layout(heights = c(1,1),widths = c(1,1),ncol = 2,nrow = 1)
dev.off()

#======Fig4B_module_heatmap===============
#WGCNA: R1_WGCNA.R
load("../WGCNA_metabolism//module_4_10_250618.RData")
module_info <- module%>%group_by(moduleColors)
module_info <- module_info[!duplicated(module_info$moduleColors),]%>%arrange(`net[["colors"]]`)
module_info$moduleColors <- paste0("ME",module_info$moduleColors)
module_info$`net[["colors"]]`[1:10] <- paste0("0",module_info$`net[["colors"]]`[1:10] )
module_info$`net[["colors"]]` <- paste0("M",module_info$`net[["colors"]]`)
colnames(module_info) <- c("module_number","module")

library(reshape2)
# Cor0 <- read_csv("../WGCNA_metabolism/moduleTraitCor_HCC.csv")
# Pvalue0 <- read_csv("../WGCNA_metabolism/moduleTraitPvalue_HCC.csv")
# colnames(Cor)[1] <- c("module")
# Cor <- melt(Cor,id.vars=c("module"))
# colnames(Cor) <- c("module","pheno","correlation")
# 
# colnames(Pvalue)[1] <- c("module")
# Pvalue<- melt(Pvalue,id.vars=c("module"))
# colnames(Pvalue) <- c("module","pheno","pvalue")
# 
# data_heatmap <- cbind(Cor,Pvalue$pvalue)
# colnames(data_heatmap)[4] <- c("pvalue")
# 
# forplot <- data_heatmap

#====用lm重新算bigmeta的cor和p====
load("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/cli_p_bigmeta.RData")
ME <- read.csv("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/ME_HCC250618.csv",row.names = 1)
covars <- read_excel("/groups/ProHuShiX/home/share/BIGMeta_clinical_clean/MGSinfo_250310.xlsx") 
covars$sample <- gsub("-",".",covars$sample)
covars <- covars[covars$sample%in%rownames(ME),c("sample","Age","Gender","BMI","Smoke","Drink")]

data <- merge(ME,cli_p,by="row.names")
data <- merge(data,covars,by.x="Row.names",by.y="sample")
data$MVI <- as.factor(data$MVI)

## AFP取下log2，不然不好上色
data$AFP <- log2(data$AFP)

feature <- c("balance_score", "TumorNumber",  "BCLC" , "AFP", "TumorSize")
covars <- covars%>%column_to_rownames(var="sample")
covars$Gender <- as.factor(covars$Gender);covars$Smoke <- as.factor(covars$Smoke);covars$Drink <- as.factor(covars$Drink)
lm_res <- lapply(colnames(cli_p)[c(1,2,3,4,6)], function(feature){
  tmp <- lapply(colnames(ME), function(module){
    glm.md <- lm(as.formula(paste(feature,"~",module ,"+",paste(colnames(covars), collapse = " + "), sep = "" )),
                 data = data)
    tmp <- summary(glm.md)$coefficients %>% as.data.frame()
    tmp$module <- row.names(tmp)
    tmp$pheno <- feature
    tmp <- as.data.frame(tmp[tmp$module==module,])
    return(tmp)
  })
  cat(feature)
  tmp <- do.call("rbind",tmp)
  return(tmp)
})
lm_res <- do.call("rbind",lm_res)

library(lmerTest)
data$MVI2 <- factor(ifelse(data$MVI=="1","Yes","No"))
tmp2 <- lapply(colnames(ME), function(module){
  glm.md <- glm(as.formula(paste("MVI2","~",module ,"+",paste(colnames(covars), collapse = " + "  ), sep = "" )),
                data = data,family = binomial)
  tmp <- summary(glm.md)$coefficients %>% as.data.frame()
  tmp$module <-module
  tmp$pheno <- "MVI"
  tmp <- as.data.frame(tmp[row.names(tmp)==module,])
  return(tmp)
})
lm_res2 <-do.call("rbind",tmp2)
colnames(lm_res) <- c("correlation","Std. Error" ,"t value","pvalue","module", "pheno" )
colnames(lm_res2) <- c("correlation","Std. Error" ,"t value","pvalue", "module","pheno" )
lm_res <-rbind(lm_res,lm_res2)
lm_res$adj.p <- p.adjust(lm_res$pvalue,method = "BH")
forplot <- lm_res

#forplot$module <- module_info$module_number[match(forplot$module,module_info$module)]
forplot$pheno <- factor(forplot$pheno,levels=c("balance_score","TumorSize","TumorNumber","AFP","BCLC","MVI"))
for( i in 1:nrow(forplot)){
  forplot$module[i] <- unique(module_info[module_info$module==forplot$module[i],"module_number"])
}
forplot$module <- as.character(forplot$module)
forplot$adj.p <- p.adjust(forplot$pvalue,method = "BH")
forplot$stars1 <- cut(forplot$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
forplot$stars2 <- cut(forplot$adj.p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("☆☆☆", "☆☆","☆",""))
forplot$stars1 <- as.character(forplot$stars1)
forplot$stars2 <- as.character(forplot$stars2)
forplot$star <- ifelse(forplot$stars2=="",forplot$stars1,forplot$stars2)

library(scales)
#forplot$estimate <- round(forplot$correlation,2)
#correlation本来就很小，可以不rescale了
#forplot$beta <- rescale(forplot$correlation,to=c(-1,1),from=range(forplot$correlation))
forplot$beta <-forplot$correlation
forplot$module <- factor(forplot$module,levels = rev(sort(unique(forplot$module))))
#forplot$adjustp <-signif(forplot$adj.p, 2)
#forplot$pvalue <-signif(forplot$pvalue , 2)
#forplot$pheno <- factor(forplot$pheno,levels = c( "HCC" , "Health", "M.funiformis" ,"E.coli",
#                                                  "balance_score","balance_score_HCC", "TumorNumber" ,"TumorSize","AFP" ,"BCLC","pathology_MVI" ))

fontsize=10
source("../WGCNA_metabolism/theme_nature.r")
c("#4475b4","#eeeeee","#d73027")
#pdf("heatmap.pdf",width = 7,height = 12)
ggplot(data=forplot, aes(x=pheno, y=module)) +
  geom_tile(aes(fill=beta),color="white",linewidth=0.5) +
  coord_equal() + 
  geom_text(aes(label=star), fontface="bold", size=2, nudge_y=0) +
  scale_fill_gradient2(low="#4475b4",mid = "#eeeeee",high= "#d73027")+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        #panel.border=element_rect(fill=NA),
        legend.position = "left", 
        panel.background = element_blank(),
        axis.ticks = element_blank(), ##去掉坐标轴刻度
        axis.title = element_blank(),
        axis.text=element_text(colour="black",size=9))##去掉坐标轴刻度文本
#dev.off()

forplot <- forplot%>%arrange(desc(forplot$module))
plot_wide1 <- forplot[,c("module","pheno","correlation")]
plot_wide1 <- plot_wide1%>%pivot_wider(names_from=pheno, values_from=c(correlation))%>%as.data.frame()
row.names(plot_wide1) <- plot_wide1$module
plot_wide1$module <- NULL
#plot_wide1 <- plot_wide1[match(rev(sort(unique(forplot1$module))),row.names(plot_wide1)),]

plot_star1 <- forplot[,c("module","pheno","star")]
plot_star1 <- plot_star1%>%pivot_wider(names_from=pheno, values_from=c(star))%>%as.data.frame()
row.names(plot_star1) <- plot_star1$module
plot_star1$module <- NULL
#plot_star1 <- plot_star1[match(rev(sort(unique(forplot1$module))),row.names(plot_star1)),]

library(ComplexHeatmap)
library(circlize)
col_fun1 <- colorRamp2(c(min(forplot$correlation), 0, max(forplot$correlation)),c("#4475b4", "#eeeeee", "#d73027"))
#col_fun1 = colorRamp2(c(-0.01, 0, 0.01),c("#4475b4","#eeeeee","#d73027"))
cell_fun1 = function(j, i, x, y, width, height, fill) {
  # 在满足条件时标记星号（比如值大于2）
  if (plot_star1[i, j] !="") {
    grid.text(plot_star1[i, j], x = x, y = y, gp = gpar(fontsize = 9, col = "black",fontface="bold"))
  }
}

plot_wide1 <- plot_wide1[,c("balance_score","TumorSize","TumorNumber","AFP","BCLC","MVI")]
plot_star1 <- plot_star1[,c("balance_score","TumorSize","TumorNumber","AFP","BCLC","MVI")]

H1 <- Heatmap(as.matrix(plot_wide1),
              cluster_columns = F,
              cluster_row = F,
              column_names_rot = 45,
              cell_fun=cell_fun1,
              # clustering_method_rows = "complete",
              show_row_names = T,
              row_names_side = "left",
              row_names_gp =  gpar(fontsize =10),  ##行名大小
              # column_names_side = "none",
              width = unit(5, "cm"),
              show_column_names = T,
              col = col_fun1,
              cluster_rows = F,
              #row_split = bandname,
              #column_split =meta_all$group,
              # row_split =select_for_heatmap$moduleColors,
              #ha = HeatmapAnnotation(df = group, points = anno_points(1:28),
              #bottom_annotation =ha_group ,   ##注释行在顶部
              rect_gp= gpar(col = "white", lty = 1,size=0.5) ##设置边框
              #HeatmapAnnotation(Group=group$group, #设置顶部热图
              #col=groupcol
              #show_legend = T
)
forplot1 <- forplot

#=======Fig4B_rxx_heatmap=============
# load("../WGCNA_metabolism/代谢验证_淑芬/MEs_rxx250618.RData")
load("../WGCNA_metabolism/代谢验证_淑芬/MEs_rxx250915.RData")
# load("../WGCNA_metabolism/代谢验证_淑芬/rxx_cli.RData")
load("../WGCNA_metabolism/代谢验证_淑芬/rxx_cli250915.RData")
data <- merge(MEs,cli,by="row.names")
data <- merge(data,covars,by.x="Row.names",by.y="row.names")

data$MVI <- as.factor(data$MVI)
data$AFP <- log2(data$AFP)

feature <- c("balance_score", "TumorNumber",  "BCLC" , "AFP", "TumorSize","MVI")

lm_res <- lapply(colnames(cli_p)[c(1,2,3,4,6)], function(feature){
  tmp <- lapply(colnames(ME), function(module){
    glm.md <- lm(as.formula(paste(feature,"~",module ,"+",paste(colnames(covars), collapse = " + "), sep = "" )),
                 data = data)
    tmp <- summary(glm.md)$coefficients %>% as.data.frame()
    tmp$module <- row.names(tmp)
    tmp$pheno <- feature
    tmp <- as.data.frame(tmp[tmp$module==module,])
    return(tmp)
  })
  cat(feature)
  tmp <- do.call("rbind",tmp)
  return(tmp)
})
lm_res <- do.call("rbind",lm_res)

library(lmerTest)
data$MVI2 <- factor(ifelse(data$MVI=="1","Yes","No"))
tmp2 <- lapply(colnames(ME), function(module){
  glm.md <- glm(as.formula(paste("MVI2","~",module ,"+",paste(colnames(covars), collapse = " + "  ), sep = "" )),
                data = data,family = binomial)
  tmp <- summary(glm.md)$coefficients %>% as.data.frame()
  tmp$module <-module
  tmp$pheno <- "MVI"
  tmp <- as.data.frame(tmp[row.names(tmp)==module,])
  return(tmp)
})
lm_res2 <-do.call("rbind",tmp2)
colnames(lm_res) <- c("correlation","Std. Error" ,"t value","pvalue","module", "pheno" )
colnames(lm_res2) <- c("correlation","Std. Error" ,"t value","pvalue", "module","pheno" )
lm_res <-rbind(lm_res,lm_res2)
lm_res$adj.p <- p.adjust(lm_res$pvalue,method = "BH")
forplot2 <- lm_res

for( i in 1:nrow(forplot2)){
  forplot2$module[i] <- unique(module_info[module_info$module==forplot2$module[i],"module_number"])
}
forplot2$module <- as.character(forplot2$module)
#forplot2$adj.p <- p.adjust(forplot2$pvalue,method = "fdr")
forplot2$stars1 <- cut(forplot2$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
forplot2$stars2 <- cut(forplot2$adj.p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("☆☆☆", "☆☆","☆",""))
forplot2$stars1 <- as.character(forplot2$stars1)
forplot2$stars2 <- as.character(forplot2$stars2)
forplot2$star <- ifelse(forplot2$stars2=="",forplot2$stars1,forplot2$stars2)

library(scales)
forplot2$beta <-forplot2$correlation
forplot2$module <- factor(forplot2$module,levels = rev(sort(unique(forplot2$module))))
forplot2$pheno <- factor(forplot2$pheno,levels = c("balance_score","BCLC","AFP","TumorNumber","TumorSize","MVI"))

forplot2 <- forplot2%>%arrange(desc(forplot2$module))
plot_wide2 <- forplot2[,c("module","pheno","correlation")]

plot_wide2 <- plot_wide2%>%pivot_wider(names_from=pheno, values_from=c(correlation))%>%as.data.frame()
row.names(plot_wide2) <- plot_wide2$module
plot_wide2$module <- NULL

plot_wide2 <- plot_wide2[,c ("balance_score", "TumorSize", "TumorNumber", "AFP","BCLC","MVI")]

#plot_wide2 <- plot_wide2[match(row.names(plot_wide1),row.names(plot_wide2)),]

plot_star2 <- forplot2[,c("module","pheno","star")]
plot_star2<- plot_star2%>%pivot_wider(names_from=pheno, values_from=c(star))%>%as.data.frame()
row.names(plot_star2) <- plot_star2$module
plot_star2$module <- NULL
plot_star2 <- plot_star2[,c ("balance_score", "TumorSize", "TumorNumber", "AFP","BCLC",         
                             "MVI")]

#plot_star2 <- plot_star2[match(row.names(plot_star1),row.names(plot_star2)),]


library(ComplexHeatmap)
library(circlize)
col_fun2 <- colorRamp2(c(min(forplot2$correlation), 0, max(forplot2$correlation)),c("#4475b4", "#eeeeee", "#d73027"))
#col_fun2 = colorRamp2(c(-0.05, 0, 0.03),c("#4475b4","#eeeeee","#d73027"))
cell_fun2 = function(j, i, x, y, width, height, fill) {
  # 在满足条件时标记星号（比如值大于2）
  if (plot_star2[i, j] !="") {
    grid.text(plot_star2[i, j], x = x, y = y, gp = gpar(fontsize = 9, col = "black",fontface="bold"))
  }
}

H2 <- Heatmap(as.matrix(plot_wide2),
              cluster_columns = F,
              cluster_row = F,
              column_names_rot = 45,
              cell_fun=cell_fun2,
              # clustering_method_rows = "complete",
              show_row_names = F,
              #row_names_side = "left",
              row_names_gp =  gpar(fontsize =10),  ##行名大小
              # column_names_side = "none",
              width = unit(5, "cm"),
              show_column_names = T,
              col = col_fun2,
              cluster_rows = F,
              #row_split = bandname,
              #column_split =meta_all$group,
              # row_split =select_for_heatmap$moduleColors,
              #ha = HeatmapAnnotation(df = group, points = anno_points(1:28),
              #bottom_annotation =ha_group ,   ##注释行在顶部
              rect_gp= gpar(col = "white", lty = 1,size=0.5) ##设置边框
              #HeatmapAnnotation(Group=group$group, #设置顶部热图
              #col=groupcol
              #show_legend = T
)
H2

#pdf("heatmap.pdf",width = 7,height = 12)
# ggplot(data=forplot, aes(x=pheno, y=module)) +
#   geom_tile(aes(fill=beta),color="white",linewidth=0.5) +
#   geom_text(aes(label=star), fontface="bold", size=1.25*fontsize/(14/5), nudge_y=0) +
#   scale_fill_gradient2(low="#0379B5",mid = "#FFFFFF",high= "#E54B8E",limit=c(-0.04,0.04))+
#   theme(axis.text.x = element_text(angle=45, hjust=1),
#         #panel.border=element_rect(fill=NA),
#         legend.position = "left", 
#         panel.background = element_blank(),
#         axis.ticks = element_blank(), ##去掉坐标轴刻度
#         axis.title = element_blank(),
#         axis.text=element_text(colour="black",size=13))

#draw(H1+H2,heatmap_legend_side = "bottom",annotation_legend_side = "right",merge_legend = TRUE)
draw(H1+H2,heatmap_legend_side = "right",merge_legend = TRUE)

pdf("./Fig4B_left_heatmap260118.pdf",width=6,height = 5)
draw(H1+H2,heatmap_legend_side = "right",merge_legend = TRUE)
dev.off()

write.xlsx(forplot1,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/FAHcohort1_WGCNA_lmPvalCor260118.xlsx")
#write.xlsx(forplot2,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/FAHcohort2_WGCNA_lmPvalCor250618.xlsx")
write.xlsx(forplot2,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/FAHcohort2_WGCNA_lmPvalCor260118.xlsx")


#=====两个队列的module cor一致性验证====

#WGCNA Cor值
big_cor <- read.csv("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/moduleTraitCor_HCC.csv")
rxx_cor <- read.csv("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/代谢验证_淑芬/rxx_moduleTraitCor_adj.csv")
i = "MEpurple"
res1=foreach(i=big_cor$X,.combine = rbind) %do% {
  tmp.data <- rbind(bigcor=big_cor[big_cor$X==i,2:7],rxxcor=rxx_cor[rxx_cor$X==i,2:7])
  tmp.data <- as.data.frame(t(tmp.data))
  p <- cor.test(tmp.data$bigcor,tmp.data$rxxcor,method="spearman")[["p.value"]]
  cor <- cor.test(tmp.data$bigcor,tmp.data$rxxcor,method="spearman")[["estimate"]]
  return.string<-data.frame(module=module_info$module_number[match(i,module_info$module)],pvalue=p,cor=cor)
}

#Lm correlation值
table(forplot1$module)
i = "M27"
res2=foreach(i=unique(forplot1$module),.combine = rbind) %do% {
  tmp.data <- merge(forplot1[forplot1$module==i,c("pheno","correlation")],forplot2[forplot2$module==i,c("pheno","correlation")],by="pheno")
  p <- cor.test(tmp.data$correlation.x,tmp.data$correlation.y,method="spearman")[["p.value"]]
  cor <- cor.test(tmp.data$correlation.x,tmp.data$correlation.y,method="spearman")[["estimate"]]
  return.string<-data.frame(module=i,pvalue=p,cor=cor)
}

# 用WGCNA的cor值加一列注释
res1 <- res1%>%arrange(desc(res1$module))
res1$stars <- cut(res1$pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))
res1$stars <- as.character(res1$stars)

plot_wide3 <- data.frame(cor=res1$cor)
rownames(plot_wide3)<-res1$module
plot_wide3 <- as.matrix(plot_wide3)
plot_star3 <- data.frame(star=res1$stars)
rownames(plot_star3) <- res1$module

col_fun3 = colorRamp2(c(-0.95, 0, 0.96),c("#4475b4","#eeeeee","#d73027"))
cell_fun3 = function(j, i, x, y, width, height, fill) {
  # 在满足条件时标记星号（比如值大于2）
  if (plot_star3[i, j] !="") {
    grid.text(plot_star3[i, j], x = x, y = y, gp = gpar(fontsize = 9, col = "black",fontface="bold"))
  }
}


H3 <- Heatmap(plot_wide3,
              cluster_columns = F,
              cluster_row = F,
              column_names_rot = 45,
              cell_fun=cell_fun3,
              # clustering_method_rows = "complete",
              show_row_names = F,
              #row_names_side = "left",
              row_names_gp =  gpar(fontsize =10),  ##行名大小
              # column_names_side = "none",
              width = unit(5, "cm"),
              show_column_names = T,
              col = col_fun3,
              cluster_rows = F,
              #row_split = bandname,
              #column_split =meta_all$group,
              # row_split =select_for_heatmap$moduleColors,
              #ha = HeatmapAnnotation(df = group, points = anno_points(1:28),
              #bottom_annotation =ha_group ,   ##注释行在顶部
              rect_gp= gpar(col = "white", lty = 1,size=0.5) ##设置边框
              #HeatmapAnnotation(Group=group$group, #设置顶部热图
              #col=groupcol
              #show_legend = T
)


#=======整理代谢物信息注释表==============
metab_info <- read_excel("../WGCNA_metabolism/代谢验证_淑芬/已被注释.xlsx", 
                         sheet = "简洁版")
metab <- merge(module,metab_info[,c("MS2 name","SuperClass")],by.x="row.names",by.y="MS2 name",all.x = T)
metab$SuperClass[is.na(metab$SuperClass)==T] <- "Other"

sort(unique(metab$SuperClass))
metab$Category <- ifelse(metab$SuperClass%like%"Lipids and lipid-like molecules","Lipids and lipid-like molecules",
                         ifelse(metab$SuperClass%like%"Nucleosides, nucleotides, and analogues","Nucleosides, nucleotides, and analogues",
                                ifelse(metab$SuperClass%like% "Organic", "Organic compounds",
                                       ifelse(metab$SuperClass%like%"Organoheterocyclic compounds","Organoheterocyclic compounds",
                                              ifelse(metab$SuperClass%like%"Phenylpropanoids and polyketides","Phenylpropanoids and polyketides",metab$SuperClass)))))
sort(unique(metab$Category))


metab$module <- ifelse(metab$`net[["colors"]]`<10,paste0("M0",metab$`net[["colors"]]`),paste0("M",metab$`net[["colors"]]`))
metab <- metab%>%arrange(metab$`net[["colors"]]`)
cl <- c("#F39B7F","#F7DB70","#71D0F5FF" ,"#D2AF81FF","#D5E4A2FF" ,"#197EC0FF","#709AE1FF"  , "#46732EFF",  "#F05C3BFF" , "#E5D7E7", "#ddbea9", "#BCD1EC" ,"#C6C09C","#766a65","#AABAC2","#E6958C","lightgrey","#8FC092")
names(cl) <- sort(unique(metab$Category))

metab$Category <- factor(metab$Category,levels =c(as.character(sort(unique(metab$Category))[1:16]),"Phenylpropanoids and polyketides", "Other") )
pdf("./Fig4B_right_anno_lengend250618.pdf",width=4,height = 5)
anno1 <- ggplot(metab,aes(x =module, fill = Category)) +
  geom_bar(position = "fill")+scale_fill_manual(values = cl)+coord_flip()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "right",
                                   panel.background = element_rect(fill = NA),
                                   axis.line.x = element_line(),
                                   axis.text.y = element_blank())
anno1
dev.off()


a <- rev(c(as.character(sort(unique(metab$Category))[1:16]) ,"Phenylpropanoids and polyketides", "Other"))
metab$Category <- factor(metab$Category,levels =a )
pdf("./Fig4B_right_anno_main250618.pdf",width=4,height = 5)
p<-ggplot(metab,aes(x =module, fill = Category)) +
  geom_bar(position = "fill")+scale_fill_manual(values = cl)+coord_flip()+
  ylab('Pecantage')+xlab("")+theme(legend.position = "right",
                                   panel.background = element_rect(fill = NA),
                                   axis.line.x = element_line(),
                                   axis.text.y = element_blank())
dev.off()

save(forplot1,forplot2,plot_wide1,plot_wide2,plot_star1,plot_star2,cli,forplot1,forplot2,
     metab,file="./Fig4B.RData")
#=====试一下把这三张图拼起来===============
library(cowplot)
library(grid)
library(patchwork)
heatmap_grob <- grid.grabExpr(draw(H1+H2,heatmap_legend_side = "right",merge_legend = TRUE))

final_plot <- plot_grid(
  heatmap_grob,
  anno1,
  ncol = 2,  # 设置列数
  rel_widths = c(1, 1)  # 设置宽度比例
)


wrap_elements(heatmap_grob) + anno1+plot_layout(ncol = 2)

#=======Fig4C_热图=============
load("../WGCNA_metabolism/datKME_HCC250618.RData")
datKME <- datKME_HCC
x=1
hub_gene <- lapply(1:ncol(datKME), function(x){
  module_name <- gsub("KME_","",colnames(datKME)[x]) 
  # 提取对应module的metabolites
  df <- datKME[rownames(datKME)%in%rownames(module[module$moduleColors==module_name,]),]
  #df <-df[df[,x]>0.8,]
  if(nrow(df)>0){
    colnames(df)[x] <- "KME"
    df$hub_gene <- row.names(df)
    df$module <- module_name
    df <- df[,c("hub_gene","module","KME")]}
  return(df)
})
hub_gene_value <- do.call("rbind",hub_gene)

writexl::write_xlsx(hub_gene_value,"/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/KME_28modules.xlsx")

hub_lightgreen <- hub_gene_value[hub_gene_value$module%in%c("lightyellow"),]%>%arrange(-KME)
lightgreen_all <- hub_lightgreen
hub_lightgreen <- hub_lightgreen[1:10,]
hub_lightgreen$hub_gene <- factor(hub_lightgreen$hub_gene,levels = rev(hub_lightgreen$hub_gene))

# 挑出M19的代谢物丰度 吲哚模块
load("../WGCNA_metabolism/input_data/sc_BIGMeta250618.RData")
M18 <-sc[,colnames(sc)%in%hub_lightgreen$hub_gene]
M18_all <-sc[,colnames(sc)%in%lightgreen_all$hub_gene]

##画热图前要手动归一化
sc_M18 <- as.data.frame(t(M18))
exp <- apply(sc_M18, 1,scale)
rownames(exp) <- colnames(sc_M18)
exp <-as.data.frame(t(exp))

load("../WGCNA_metabolism/cli_bigmeta1.RData")
cli$group <- ifelse(cli$HCC==1,"HCC","Health")

meta_all <- cli[match(colnames(exp),row.names(cli)),]

hub_lightgreen <- hub_lightgreen[match(row.names(exp),hub_lightgreen$hub_gene),]%>%arrange(-KME)

#######挑出M07的代谢物丰度 胆汁酸模块
hub_black <- hub_gene_value[hub_gene_value$module=="black",]%>%arrange(-KME)
black_all <- hub_black
hub_black <- hub_black[1:10,]%>%arrange(-KME)
hub_black$hub_gene <- factor(hub_black$hub_gene,levels = rev(hub_black$hub_gene))

M07 <-sc[,colnames(sc)%in%hub_black$hub_gene]
M07_all <- sc[,colnames(sc)%in%black_all$hub_gene]


##画热图前要手动归一化
sc_M07 <- as.data.frame(t(M07))
exp2 <- apply(sc_M07, 1,scale)
rownames(exp2) <- colnames(sc_M07)
exp2 <-as.data.frame( t(exp2))

exp2 <- exp2[,match(colnames(exp),colnames(exp2))]

exp3 <- rbind(exp,exp2)

hub_all <- rbind(hub_lightgreen,hub_black)
hub_all$module_number <- ifelse(hub_all$module=="lightyellow","M19","M07")

meta_all <- cli[match(colnames(exp3),row.names(cli)),]
exp3 <- exp3[match(hub_all$hub_gene,row.names(exp3)),]
hub_all$module_number <- factor(hub_all$module_number,levels=c("M19","M07"))


library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-3.5, 4,8), c("#8B516E", "#FEF5B7","#5D7160"))
col_fun_km=colorRamp2(c(0, 1), c("#eeeeee", "#010985"))
col_fun_sc=colorRamp2(c(-3.5,0, 4.0),c("#4475b4","#eeeeee","#d73027"))

ha_group = HeatmapAnnotation(group=meta_all$group,
                             balance_score=meta_all$balance_score,
                             col = list(
                               group= c("Health" = "#31AEEF", "HCC"="#FF9A9C"),
                               balance_score=col_fun),
                             annotation_name_side = "left",
                             annotation_legend_param = list(
                               group = list(nrow = 1),
                               balance_score = list(direction = "horizontal")))

column_ha = HeatmapAnnotation(balance_score = meta_all$balance_score, col=list(balance_score=col_fun))
row_ha = rowAnnotation(KME = hub_all$KME,col=list(KME=col_fun_km))

exp3<-as.matrix(exp3)

H3 <- Heatmap(exp3,
              cluster_columns = T,
              #cluster_row = F,
              #clustering_method_rows = "complete",
              show_row_names = T,
              cluster_row_slices = FALSE, 
              cluster_column_slices = FALSE,
              #row_order=hub_all$hub_gene,
              row_names_gp =  gpar(fontsize =10, fontface="bold"),  ##行名大小
              # column_names_side = "none",
              width = unit(5, "cm"),
              show_column_names = F,
              # col = col_fun_sc,
              cluster_rows = F,
              border=F,
              # heatmap_legend_param = list(direction = "horizontal"),
              row_split = hub_all$module_number,
              #column_split =meta_all$group,
              #row_split =select_for_heatmap$module,
              #ha = HeatmapAnnotation(df = group, points = anno_points(1:28),
              right_annotation = row_ha,
              bottom_annotation =column_ha, ##注释行在底部
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#d73027","#4475b4")),
                                                               labels = c("M19","M07"), 
                                                               labels_gp = gpar(col = "white", fontsize = 10))),
              row_dend_reorder = F,
              # column_dend_reorder = F
              #rect_gp= gpar(col = "white", lty = 1,size=0.5) ##设置边框
              #HeatmapAnnotation(Group=group$group, #设置顶部热图
              col=col_fun_sc
              #show_legend = T
)

pdf("./Fig4C_module_heatmap250618.pdf",width=7,height=5)
H3
dev.off()
save(exp3,meta_all,hub_all,file="./Fig3.RData")


#======Fig4D_挑代谢物画图===========
load("../WGCNA_metabolism/input_data/sc_BIGMeta250618.RData")
# sc <- sc[row.names(sc)%in%row.names(cli_p),]

Cluster <- read.table("../Final_output/Balance.score1.txt",header=T)
sample <- read.table("../WGCNA_metabolism/input_data/MGS_metabolism_link.txt",header = T)

# for (i in 1:nrow(sample)){
#   row.names(sc) <- gsub(sample$metab_sample[i],sample$mgs_sample[i],row.names(sc))
# }

metab <- colnames(sc)%>%as.data.frame()
metab$num <- paste0("V",1:ncol(sc))
colnames(metab)[1] <- "metabolites"
colnames(sc) <- paste0("V",1:ncol(sc))

tmp.data <- merge(sc,Cluster[,c("Row.names","balance_value")],by.x="row.names",by.y="Row.names",all.x=T)
colnames(tmp.data) <- gsub("balance_value","balance_score",colnames(tmp.data))
tmp_data <- tmp.data

# sur <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/MGSinfo_0514.xlsx")
# tmp_data <- merge(tmp_data,sur[,c("sample", "Age","Gender" , "BMI", "Smoke", "Drink" )],by.x="Row.names",by.y="sample",all.x=T)
# #tmp_data <- tmp_data[tmp_data$Row.names%in%MGS_meta[MGS_meta$group=="HCC","sample"],]
# 
# res_cor <- list()
# res_lm <- list()
# res_lm_ad <- list()
# for (i in 1:nrow(metab)){
#   tmp.data <- tmp_data[,c(metab$num[i],"balance_score")]
#   tmp.res <-cor.test(tmp.data[,1],tmp.data$balance_score,method = "spearman") 
#   res_cor[[i]] <- data.frame(rho=tmp.res[["estimate"]],pvalue=tmp.res[["p.value"]],metab=metab$metabolites[i])
#   # tmp.data[[metab$num[i]]] <- rmOutlier(tmp.data[[metab$num[i]]])
#   #tmp.data <- drop_na(tmp.data, metab$num[i])
#   mod1 <- lm(as.formula(paste0(metab$num[i],"~balance_score")),data=tmp.data)
#   a <- as.data.frame(summary(mod1)[["coefficients"]])
#   a<-  as.data.frame(a[-1,])
#   a$metab <- metab$metabolites[i]
#   res_lm[[i]] <-a
#   tmp.data <-  tmp_data[,c(metab$num[i],"balance_score","Age","Gender","BMI","Smoke","Drink")]
#   tmp.data$Age <- as.numeric(tmp.data$Age)
#   tmp.data$BMI <- as.numeric(tmp.data$BMI)
#   tmp.data$Gender <- as.factor(tmp.data$Gender)
#   tmp.data$Smoke <- as.factor(tmp.data$Smoke)
#   tmp.data$Drink <- as.factor(tmp.data$Drink)
#   # tmp.data[[metab$num[i]]] <- rmOutlier(tmp.data[[metab$num[i]]])
#   #tmp.data <- drop_na(tmp.data, metab$num[i])
#   mod2<-lm(as.formula(paste0(metab$num[i],"~balance_score+Age+Gender+BMI+Smoke+Drink")),data = tmp.data)
#   #mod2<-lm(as.formula(paste0(metab$num[i],"~balance_score+Age+Gender+BMI")),data = tmp.data)
#   a <- as.data.frame(summary(mod2)[["coefficients"]])
#   a<-  as.data.frame(a[2,])
#   a$metab <- metab$metabolites[i]
#   res_lm_ad[[i]] <-a
# }
# res_cor <- do.call("rbind",res_cor)
# res_lm <- do.call("rbind",res_lm)
# res_lm_ad <- do.call("rbind",res_lm_ad)
# 
# res_lm$adj.p <- p.adjust(res_lm$`Pr(>|t|)`,"fdr")
# res_lm_ad$adj.p <- p.adjust(res_lm_ad$`Pr(>|t|)`,"fdr")



load("../WGCNA_metabolism/module_4_10_250618.RData")
module_info <- module%>%group_by(moduleColors)
module_info <- module_info[!duplicated(module_info$moduleColors),]%>%arrange(`net[["colors"]]`)
module_info$moduleColors <- paste0("ME",module_info$moduleColors)
module_info$`net[["colors"]]`[1:10] <- paste0("0",module_info$`net[["colors"]]`[1:10] )
module_info$`net[["colors"]]` <- paste0("M",module_info$`net[["colors"]]`)
colnames(module_info) <- c("module_number","module")


load("../WGCNA_metabolism/datKME_HCC250618.RData")
hub_gene <- lapply(1:ncol(datKME), function(x){
  module_name <- gsub("KME_","",colnames(datKME)[x]) 
  df <- datKME[rownames(datKME)%in%rownames(module[module$moduleColors==module_name,]),]
  #df <-df[df[,x]>0.8,]
  if(nrow(df)>0){
    colnames(df)[x] <- "KME"
    df$hub_gene <- row.names(df)
    df$module <- module_name
    df <- df[,c("hub_gene","module","KME")]}
  return(df)
})

hub_gene_value <- do.call("rbind",hub_gene)
hub_lightcyan <- hub_gene_value[hub_gene_value$module%in%c("lightyellow","black"),]

cor_target <- res_cor[res_cor$metab%in%hub_lightgreen$hub_gene,]
lm_target <- res_lm[res_lm$metab%in%hub_lightgreen$hub_gene,]
lm_ad_target <- res_lm_ad[res_lm_ad$metab%in%hub_lightgreen$hub_gene,]

cor_target <- merge(cor_target,module,by.x="metab",by.y="row.names",all.x=T)
lm_target <- merge(lm_target,module,by.x="metab",by.y="row.names",all.x=T)
lm_ad_target <- merge(lm_ad_target,module,by.x="metab",by.y="row.names",all.x=T)

target <- c("MethylIndole-3-acetate","indolin-2-one","Indoxyl sulfate",
            "Taurocholic acid"  ,"Taurodeoxycholic acid","Taurochenodeoxycholate" )

col_fun_sc=colorRamp2(c(-3.5,0, 4.0),c("#4475b4","#eeeeee","#d73027"))

i=1
plot <- list()
for (i in 1:3){
  num <- metab[metab$metabolites==target[i],"num"]
  plot_df <- tmp_data[,c(num,"balance_score")]
  colnames(plot_df) <- c("value","balance_score")
  plot[[i]]<- ggplot(data = plot_df,aes(x=balance_score,y=value))+
    geom_point(color="#C25759",fill="#C25759",alpha=0.65)+geom_smooth(method = 'lm',color="#B83945",fill="#EDB8B0")+
    # scale_color_manual(values = col3)+
    stat_cor(method = 'spearman', label.y=0.1)+
    theme(panel.background =element_blank(),legend.position ="none",
          axis.line=element_line(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          #axis.text=element_text(face="bold",colour="black",size=20),
          #axis.title = element_text(family = "Cuti",face="bold",colour="black",size=20),
    )+labs(title = target[i],y="Scaled Abundance",x="balance score")
}

# <2.2e-16时获取 spearman 相关系数和 p 值
cor_test_result <- cor.test(plot_df$balance_score, plot_df$value, method = "spearman",exact = FALSE)
R_value <- cor_test_result$estimate  # 相关系数 rho
p_value <- cor_test_result$p.value   # p 值


for (i in 4:6){
  num <- metab[metab$metabolites==target[i],"num"]
  plot_df <- tmp_data[,c(num,"balance_score")]
  colnames(plot_df) <- c("value","balance_score")
  plot[[i]]<- ggplot(data = plot_df,aes(x=balance_score,y=value))+
    geom_point(color="#4475b4",fill="#4475b4",alpha=0.65)+geom_smooth(method = 'lm',color="#115485",fill="#A4D4F2")+
    # scale_color_manual(values = col3)+
    stat_cor(method = 'spearman', label.y=0.1)+
    theme(panel.background =element_blank(),legend.position ="none",
          axis.line=element_line(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          #axis.text=element_text(face="bold",colour="black",size=20),
          #axis.title = element_text(family = "Cuti",face="bold",colour="black",size=20),
    )+labs(title = target[i],y="Scaled Abundance",x="balance score")
}

library(gridExtra)
pdf("./Fig4D250721.pdf",width=9,height=6)
plot[[1]]+plot[[2]]+plot[[3]]+plot[[4]]+plot[[5]]+plot[[6]]+plot_layout(ncol=3,nrow=2)
dev.off()
#ggsave("0.5-0.7.pdf", width = 30, height=10, marrangeGrob(grobs =plist,  ncol = 5,nrow=2))
#marrangeGrob(grobs =plot,  ncol = 3,nrow=2)  ##marrangeGrob无法修改放图顺序
save(tmp_data,metab,file="./Fig4D250618.RData")

#====M19模块中单个代谢物跟BS/ETE做spearman====
# load("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/442MGS_clinical.RData")
cluster <- Cluster[Cluster$Row.names%in%rownames(M18_all),]
M18_all <- M18_all[na.omit(match(Cluster$Row.names,rownames(M18_all))),]
all(cluster$Row.names==rownames(M18_all))
table(cluster$new_group)


res19 <- list()
i=1
for (i in 1:ncol(M18_all)){
  tmp.data <- data.frame(value=M18_all[,i],cluster=cluster$Cluster,BS=cluster$balance_value,new_group=cluster$new_group)
  tmp.data$BS <- as.numeric(tmp.data$BS)
  tmp.data$cluster <- as.factor(ifelse(tmp.data$cluster==4,"ETE","ETnonE"))
  
  # ETE/ETnonE wilcox
  # p.value1 <- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  
  # ETE/ETnonE spearman
  test1 <- cor.test(tmp.data$value,tmp.data$new_group, 
                    method = "spearman", exact = FALSE)
  rho1 = test1$estimate
  p.value2 = test1$p.value
  
  # BS spearman
  test2 <- cor.test(tmp.data$value,tmp.data$BS, 
                    method = "spearman", exact = FALSE)
  rho2 = test2$estimate
  p.value3 = test2$p.value
  
  res19[[i]] <- data.frame(metabolite=colnames(M18_all)[i],rho_ETE=rho1,pvalue_ETE=p.value2,rho_BS=rho2,pvalue_BS=p.value3)
}
res19 <- do.call("rbind",res19)
res19$padj_BS <- p.adjust(res19$pvalue_BS,"BH")
res19$padj_ETE <- p.adjust(res19$pvalue_ETE,"BH")


#====M07模块中单个代谢物跟BS/ETE做spearman====
M07_all <- M07_all[na.omit(match(cluster$Row.names,rownames(M07_all))),]
all(cluster$Row.names==rownames(M07_all))
table(cluster$new_group)

res07 <- list()
i=1
for (i in 1:ncol(M07_all)){
  tmp.data <- data.frame(value=M07_all[,i],cluster=cluster$Cluster,BS=cluster$balance_value,new_group=cluster$new_group)
  tmp.data$BS <- as.numeric(tmp.data$BS)
  tmp.data$cluster <- as.factor(ifelse(tmp.data$cluster==4,"ETE","ETnonE"))
  
  # ETE/ETnonE wilcox
  # p.value1 <- wilcox.test(value~cluster, data=tmp.data)[["p.value"]]
  
  # ETE/ETnonE spearman
  test1 <- cor.test(tmp.data$value,tmp.data$new_group, 
                    method = "spearman", exact = FALSE)
  rho1 = test1$estimate
  p.value2 = test1$p.value
  
  # BS spearman
  test2 <- cor.test(tmp.data$value,tmp.data$BS, 
                    method = "spearman", exact = FALSE)
  rho2 = test2$estimate
  p.value3 = test2$p.value
  
  res07[[i]] <- data.frame(metabolite=colnames(M07_all)[i],rho_ETE=rho1,pvalue_ETE=p.value2,rho_BS=rho2,pvalue_BS=p.value3)
}
res07 <- do.call("rbind",res07)
res07$padj_BS <- p.adjust(res07$pvalue_BS,"BH")
res07$padj_ETE <- p.adjust(res07$pvalue_ETE,"BH")

write_xlsx(list(M19=res19,M07=res07),path="/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/M19M07_metabolites_cor_with_BS&ETE.xlsx")

#======试一下分组比较==============
target_sc <- sc[,colnames(sc)%in%metab[metab$metabolites%in%hub_lightgreen$hub_gene,"num"]]
target_sc <- merge(target_sc,Cluster[,c("Row.names","balance_value")],by.x="row.names",by.y="Row.names",all.x=T)
target_sc$group <- ifelse(target_sc$balance_value>quantile(target_sc$balance_value)[3],"high","low")

res_wilcox <- list()
for (i in 1:nrow(hub_lightgreen)){
  metabolite <- hub_lightgreen$hub_gene[i]
  num <- metab[metab$metabolites==metabolite,"num"]
  mod1 <- lm(as.formula(paste0(num,"~group")),data=target_sc)
  a <- as.data.frame(summary(mod1)[["coefficients"]])
  res_wilcox[[i]] <- data.frame(metab=metabolite,
                                pvalue=wilcox.test(as.formula(paste0(num,"~group")),data=target_sc)[["p.value"]],
                                lm.p=a[2,4]
  )
}
res_wilcox <- do.call("rbind",res_wilcox)

####试一下分组踢离群值
res_wilcox <- list()
for (i in 1:nrow(hub_lightgreen)){
  metabolite <- hub_lightgreen$hub_gene[i]
  num <- metab[metab$metabolites==metabolite,"num"]
  tmp.data <- target_sc[,c(num,"group")]
  a <- tmp.data[tmp.data$group=="high",]
  a[[num]] <- rmOutlier(a[[num]])
  a <- drop_na(a,all_of(num))
  b <- tmp.data[tmp.data$group=="low",]
  b[[num]] <- rmOutlier(b[[num]])
  b <- drop_na(b,all_of(num))
  tmp.data <- rbind(a,b)
  mod1 <- lm(as.formula(paste0(num,"~group")),data=tmp.data)
  a <- as.data.frame(summary(mod1)[["coefficients"]])
  res_wilcox[[i]] <- data.frame(metab=metabolite,
                                pvalue=wilcox.test(as.formula(paste0(num,"~group")),data=tmp.data)[["p.value"]],
                                lm.p=a[2,4]
  )
}
res_wilcox <- do.call("rbind",res_wilcox)

