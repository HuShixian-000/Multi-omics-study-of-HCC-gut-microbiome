rm(list=ls())
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/")

library(WGCNA)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)
set.seed(12345)

#=======================================
#加载数据
#=======================================
load("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/input_data/sc_BIGMeta250618.RData")
MGS <- read.table("./input_data/MGS_metabolism_link.txt",header = T)
MGS$mgs_sample <- gsub("-",".",MGS$mgs_sample)
cli <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt")
#########这部分是合并新的balance_score的
cli <- subset(cli[,c("Row.names","balance_value","group","Cluster"),drop=F])
rownames(cli) <- NULL
cli <- cli %>% column_to_rownames(var="Row.names")
cli$E.coli <- as.numeric(ifelse(cli$Cluster==4,1,0))
cli$HCC <- as.numeric(ifelse(cli$group=="HCC",1,0))
cli$Health <- as.numeric(ifelse(cli$group=="Health",1,0))
cli$group <- NULL
cli$Cluster <- NULL
colnames(cli)[1] <- "balance_score"
# save(cli,file="cli_bigmeta1.RData")
#===========================================
#WGCNA
#===========================================
#基因过滤
gsg=goodSamplesGenes(sc,verbose=3)
gsg$allOK

#基因数目
nGenes = ncol(sc)
#样本数目
nSamples = nrow(sc) 
#提取表达矩阵样本名
metaSamples<-rownames(sc)
#提取表型文件样本名
cliSamples<-rownames(cli)
#表达矩阵样本名顺序匹配表型文件样本名
cliRows<-match(metaSamples,cliSamples)
#获得匹配好样本顺序的表型文件/表型样本名与矩阵样本名对齐
cli<-cli[cliRows,]
# library(openxlsx)
#cli$sample <- rownames(cli)
#write.xlsx(cli,"./bigmeta_cli_forWGCNA.xlsx")


# #绘制树+表型热图
#====样本可视化====
sampleTree2<-hclust(dist(sc),method="average")
traitColors<-numbers2colors(cli,signed=FALSE) # 为所有样本的所有表型值分配颜色
dim(traitColors)
#png("sample-subtype-cluster.png",width = 800,height = 600)
pdf("../Fig4/Fig4S/sample dendrogram.pdf",width = 20,height = 6)
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels=names(cli),
                    main="Sample dendrogram and trait heatmap",
                    cex.colorLabels=1.5,cex.dendroLabels=1,cex.rowText=2)
dev.off() 


#====选择beta====
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
powers = c(seq(1, 20, by = 1))
sft = pickSoftThreshold(sc, powerVector = powers, verbose = 5)
sft$powerEstimate # 3
# sft$fitIndices：包含各个 power 的 fit 


par(mfrow = c(1,2));# 布局参数 1行2列
cex1 = 0.9;
pdf("../Fig4/Fig4S/step2-beta-value.pdf",width = 12,height = 6)
par(mfrow=c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#一步法构建共表达网络
##需要调用WGCNA包自带的cor函数，不然会发生报错奥！
cor<-WGCNA::cor
##在进行共表达网络构建时，power值的选择非常重要，最影响结果的一个参数，需要经过多次尝试，才能找到最适合的。
net = blockwiseModules(sc, power = 4, maxBlockSize = nGenes,
                       TOMType ='unsigned', minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, saveTOMFileBase = "drought",
                       verbose = 3)
table(net$colors)

cor<-stats::cor
#绘制基因聚类树和模块颜色组合
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
pdf("../Fig4/Fig4S/415_step4-genes-modules.pdf",width = 8,height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#计算模块与性状间的相关性及绘制相关性热图
MEs0 = moduleEigengenes(sc, moduleColors)$eigengenes
# Recalculate MEs with color labels

##不同颜色的模块的ME值矩 (样本vs模块)
##模块特征基因（Module Eigengene, ME）
##是一个模块内所有基因表达的第一主成分（PCA第一轴）。
##可以看成是这个模块在每个样本中的“综合表达趋势”。
##是表示一个模块整体表达特征的代表性值。
MEs = orderMEs(MEs0);
moduleTraitCor = cor(MEs, cli, use = "p",method = "spearman"); # 获取cor的矩阵
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples) # 获取p值矩阵 这里Student只适用于pearson
moduleTraitPvalue <- matrix(NA, nrow = ncol(MEs), ncol = ncol(cli),dimnames = list(colnames(MEs), colnames(cli)))
for (i in 1:ncol(MEs)) {
  for (j in 1:ncol(cli)) {
    test_result <- cor.test(MEs[, i], cli[, j], method = "spearman", use = "pairwise.complete.obs")
    moduleTraitPvalue[i, j] <- test_result$p.value
  }
}
#moduleTraitPvalue<- p.adjust(moduleTraitPvalue , method = 'bonferroni')
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#png("step5-Module-trait-relationships.png",width = 800,height = 1200,res = 120)
par(mar = c(5, 6, 1, 1));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(cli),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#dev.off() 
#write.csv(MEs,"module/MEs.csv")
moduleTraitCor <- as.data.frame(moduleTraitCor)
write.csv(moduleTraitCor,"moduleTraitCor_all_spearman250618.csv")
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
write.csv(moduleTraitPvalue,"moduleTraitPvalue_all_spearman250618.csv")
module <- as.data.frame(net[["colors"]])
module <- cbind(module,moduleColors)
save(module,file="module_4_10_250618.RData")
write.csv(MEs,"ME_all250618.csv")

# 模块性状组合图
#par(cex = 0.9)
pdf("../Fig4/Fig4S/415_step7-Eigengene-dendrogram.pdf",width = 8,height = 6)
plotEigengeneNetworks(MEs, "",   ##绘制模块聚类图和热图
                      marDendro =c(0,5,1,5),  ##树类图区间大小谁宿
                      marHeatmap = c(5,6,1,2), cex.lab = 0.8, ##树类图区间大小谁宿
                      xLabelsAngle = 90)  ##轴鲜艿90庿

dev.off()

#批量提取特征代谢物
datKME=signedKME(sc, MEs,outputColumnName="KME_")
hub_gene <- lapply(1:ncol(datKME), function(x){
  module_name <- gsub("KME_","",colnames(datKME)[x]) 
  df <- datKME[rownames(datKME)%in%rownames(module[module$moduleColors==module_name,]),]
  df <-df[df[,x]>0.8,]
  if(nrow(df)>0){
    df$hub_gene <- row.names(df)
    df$module <- module_name
    df <- df[,c("hub_gene","module")]}
  return(df)
})
#names(hub_gene) <- gsub("KME_","",colnames(datKME)) 
hub_gene <- do.call("rbind",hub_gene)
write.table(hub_gene,"hub_gene_4_10_250618.txt",quote = F,row.names = F,sep="\t")

#============================================
#看与HCC表型结果的关系
#============================================
sur <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_input/MGSinfo_0514.xlsx")
cli_p <- sur[sur$sample%in%rownames(cli),] # 361HCC
cli_p <- cli_p %>% select(c("sample","TumorSize","TumorNumber","AFP","BCLC","MVI"))
rownames(cli_p) <- NULL
cli_p <- cli_p %>% column_to_rownames(var="sample")
cli_p$balance_score <- cli$balance_score[match(rownames(cli_p),rownames(cli))]
cli_p$BCLC <- as.numeric(ifelse(cli_p$BCLC=="0",1,ifelse(cli_p$BCLC=="A",2,
                                          ifelse(cli_p$BCLC=="B",3,4))))
cli_p$TumorNumber <- gsub(">3","4",cli_p$TumorNumber);cli_p$TumorNumber <- as.numeric(cli_p$TumorNumber)
cli_p_rownames <- rownames(cli_p)
cli_p$MVI <- as.numeric(cli_p$MVI)
cli_p$TumorSize <- as.numeric(cli_p$TumorSize)
cli_p$balance_score <- as.numeric(cli_p$balance_score)
save(cli_p,file="cli_p_bigmeta.RData")

sc <- sc[row.names(sc)%in%row.names(cli_p),]

#提取表达矩阵样本名
metaSamples<-rownames(sc)
#提取表型文件样本名
cliSamples<-rownames(cli_p)
#表达矩阵样本名顺序匹配表型文件样本名
cliRows<-match(metaSamples,cliSamples)
#获得匹配好样本顺序的表型文件
cli_p<-cli_p[cliRows,]

# Recalculate MEs with color labels
MEs0 <- MEs0[row.names(MEs0)%in%row.names(cli_p),]
MEs0samples <- row.names(MEs0)
MEs0Rows<-match(metaSamples,MEs0samples)
MEs0 <- MEs0[MEs0Rows,]

#MEs0 = moduleEigengenes(sc, moduleColors)$eigengenes
##不同颜色的模块的ME值矩 (样本vs模块)
MEs = orderMEs(MEs0);
MEs=MEs[MEs0Rows,]
moduleTraitCor = cor(MEs, cli_p , use = "p",method = "spearman");
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue <- matrix(NA, nrow = ncol(MEs), ncol = ncol(cli_p),dimnames = list(colnames(MEs), colnames(cli_p)))
for (i in 1:ncol(MEs)) {
  for (j in 1:ncol(cli_p)) {
    test_result <- cor.test(MEs[, i], as.numeric(cli_p[,j]),method = "spearman", exact = FALSE)
    moduleTraitPvalue[i, j] <- test_result$p.value
  }
}
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#png("step5-Module-trait-relationships.png",width = 800,height = 1200,res = 120)
par(mar = c(5.5, 6, 1, 1));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(cli_p),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#dev.off() 

moduleTraitCor <- as.data.frame(moduleTraitCor)
write.csv(moduleTraitCor,"moduleTraitCor_HCC_spearman250618.csv")
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
write.csv(moduleTraitPvalue,"moduleTraitPvalue_HCC_spearman250618.csv")
module <- as.data.frame(net[["colors"]])
module <- cbind(module,moduleColors)
write.csv(MEs,"ME_HCC250618.csv")

datKME_HCC=signedKME(sc, MEs,outputColumnName="KME_")
save(datKME_HCC,file="datKME_HCC250618.RData")

#====ME值和ETE/ETnonE做相关====
cluster <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/Final_output/Balance.score1.txt",header=T)
table(cluster$new_group)
cli_p$Cluster <- as.numeric(cluster$new_group[match(rownames(cli_p),cluster$Row.names)])
table(cli_p$Cluster)
cli_p$Cluster <- as.numeric(ifelse(cli_p$Cluster==0,1,0))

res_list <- lapply(1:ncol(MEs), function(i) {
  test <- cor.test(MEs[, i], as.numeric(as.factor(cli_p$Cluster)), 
                   method = "spearman", exact = FALSE)
  data.frame(
    module = colnames(MEs)[i],
    rho = test$estimate,
    p_value = test$p.value
  )
})

res_df <- do.call(rbind, res_list)
res_df$FDR <- p.adjust(res_df$p_value, method = "BH")
res_df <- res_df[order(res_df$p_value), ]

write_xlsx(res_df,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/FAHcohort1_WGCNA_ME&Cluster_spearman250711.xlsx")
#=======================================
#导出模块代谢物名称做通路富集
#=======================================
library(readxl)
annotation <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/代谢验证_淑芬/已被注释.xlsx",  sheet = "简洁版")
for (i in 1:nrow(annotation)){
  if (is.na(annotation$`MS2 name`[i])==T){
    annotation$`MS2 name`[i] <- annotation$`MS1 name`[i]
  }
}
library(writexl)
module$metabolite <- row.names(module)
colnames(module)[1] <- "module_num"
for(i in unique(module$moduleColors)){
  df <- subset(module,module$moduleColors==i)
  df_mod <- merge(df,annotation,by.x="metabolite",by.y="MS2 name",all.x=T)
  write.xlsx(df_mod,paste0("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/annon3_250618/",i,"_anno.xlsx"))
}

#====合并28个module包含的metabolites====
setwd("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/annon3_250618/")
suffix <- "_anno.xlsx"
file_name <- list.files('/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/annon3_250618/',pattern = suffix)
file_name2 <- map_chr(file_name,function(x){stringr::str_remove_all(x,suffix)})
all_files <- lapply(file_name,function(x){y <- read_excel(x)})
names(all_files) <- file_name2
all_files <- do.call("rbind",all_files)
all_files <- merge(all_files,metab,by.x="metabolite",by.y="Row.names")
write.xlsx(all_files,"/groups/ProHuShiX/home/xiashufen/bigmeta/supplymentary/metabolites_28module250618.xlsx")

# HCC内部与BS正/负相关的module
# BS cor<0 p<0.05 8个
color1 <- intersect(rownames(moduleTraitCor)[which(moduleTraitCor$balance_score<0)],rownames(moduleTraitPvalue)[which(moduleTraitPvalue$balance_score<0.05)])
moduleTraitPvalue$padj <- p.adjust(moduleTraitPvalue$balance_score,method = "fdr")
color11 <- intersect(rownames(moduleTraitCor)[which(moduleTraitCor$balance_score<0)],rownames(moduleTraitPvalue)[which(moduleTraitPvalue$padj<0.05)])
color1 <- gsub("ME","",color1)
color11 <- gsub("ME","",color11)
# BS cor>0 p<0.05 2个
color2 <- intersect(rownames(moduleTraitCor)[which(moduleTraitCor$balance_score>0)],rownames(moduleTraitPvalue)[which(moduleTraitPvalue$balance_score<0.05)])
color22 <- intersect(rownames(moduleTraitCor)[which(moduleTraitCor$balance_score>0)],rownames(moduleTraitPvalue)[which(moduleTraitPvalue$padj<0.05)])
color2 <- gsub("ME","",color2)
color22 <- gsub("ME","",color22)

save(color1,color11,color2,color22,file="pvalue_colors.RData")

df_1 <- subset(module,module$moduleColors%in%color1)
df_1 <- merge(df_1,annotation,by.x="metabolite",by.y="MS2 name",all.x=T)

df_2 <- subset(module,module$moduleColors%in%color2)
df_2 <- merge(df_2,annotation,by.x="metabolite",by.y="MS2 name",all.x=T)

# cor < 0
black <- read_excel("annon2/black_anno.xlsx")
# cor > 0
lightgreen <- read_excel("annon2/lightgreen_anno.xlsx") # power=4

grey60 <- read_excel("annon/grey60_anno.xlsx") # power=3

save(MEs,MEs0,hub_gene,black,lightgreen,moduleTraitCor,moduleTraitPvalue,module,color1,color2,file="annon2/power4_res.RData")

#=======================================
#批量导出模块代谢物丰度，计算相关性
#=======================================
load("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/input_data/sc_BIGMeta250618.RData")
for (i in unique(module$moduleColors)){
  tmp.module <- module[module$moduleColors==i,]
  tmp.color <- unique(tmp.module$moduleColors)
  tmp.num <- unique(tmp.module$module_num)
  tmp.data <- sc[,row.names(tmp.module)]
  write.csv(tmp.data,paste0("/groups/ProHuShiX/home/xiashufen/bigmeta/WGCNA_metabolism/module250618/",tmp.num,"_",i,"_module.csv"))
}

