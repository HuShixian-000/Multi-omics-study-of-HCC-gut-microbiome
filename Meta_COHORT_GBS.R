#=====New York PRJNA541981====
sampleinfo <- read.csv("./2022-NM-melanoma/PRJNA541981_SraRunTable&metadata.csv")
sampleinfo <- sampleinfo%>%filter(Assay.Type=="WGS"&time=="Baseline"&Run!="SRR9033723")
table(sampleinfo$progressed)
# 0  1 
# 15 12 
metadata <- read_excel("./2022-NM-melanoma/metadata.xlsx")
species <- read.table("./2022-NM-melanoma/PRJNA541981_metaphlan_metagenome.txt",header = T,row.names = 1)
colnames(species) <- gsub(".profiled_metagenome","",colnames(species))
species <- species[,colnames(species_s)%in%sampleinfo$Run]

species_s <- species %>% 
  dplyr::filter(str_detect(rownames(species),"s__.*") & !str_detect(rownames(species),"t__.*")) 
species_s <- species_s %>% 
  dplyr::mutate(species=str_extract(rownames(species_s),"s__.*"))
rownames(species_s) <- NULL
species_s <- species_s %>% column_to_rownames(var="species")
rownames(species_s) <- gsub("s__","",rownames(species_s))
species<-species_s
species <- species_s[rowSums(species_s!=0) > (ncol(species_s)*0.1),]
# species <- species[rownames(species)%in%colnames(data_mgs),]
species <- apply(species,2,function(x){
  y <- x/sum(x)
  return(y)
})

species <- as.data.frame(t(species))

taxa_up_index <- na.omit(match(taxa_up$x, colnames(species))) #34个
taxa_down_index <- na.omit(match(taxa_down$x,colnames(species))) #10个 5个

# Compute balance
balance_df2=data.frame(sample=rownames(species),balance_value=NA)
species <- species %>% mutate_all(as.numeric)
species[species==0]=min(species[species>0])
for(sample in rownames(species)){
  balance_df2$balance_value[balance_df2$sample==sample] = log(exp(mean(as.numeric(log(species[sample,taxa_up_index])),na.rm = T))) - log(exp(mean(as.numeric(log(species[sample,taxa_down_index])),na.rm = T)))  
}
balance_df2$balance_value <- as.numeric(balance_df2$balance_value)


balance_df2$group <- as.factor(sampleinfo$progressed[match(balance_df2$sample,sampleinfo$Run)])
balance_df2$gender <- as.factor(sampleinfo$sex[match(balance_df2$sample,sampleinfo$Run)])
balance_df2$age <- as.numeric(sampleinfo$AGE[match(balance_df2$sample,sampleinfo$Run)])
balance_df2$BMI <- as.numeric(sampleinfo$BMI[match(balance_df2$sample,sampleinfo$Run)])

balance_df2$adj_bs <- residuals(lm(balance_value~gender+age+BMI,data=balance_df2))

balance_f<-balance_df2%>%filter(sample!="SRR9033723")
ggboxplot(balance_f, x = "group", y = "balance_value")+
  stat_compare_means(comparisons = list(c("0","1")))
ggboxplot(balance_f, x = "group", y = "adj_bs")+
  stat_compare_means(comparisons = list(c("0","1")))
table(balance_f$group)

balance_f$group <- factor(ifelse(balance_f$group==0,"Non progressed","Progressed"),levels = c("Non progressed","Progressed"))

balance_df_mela <- balance_f


#====PRJNA751792====
species <- read.table("./PRJNA751792_metaphlan_metagenome.txt",row.names = 1,header = T)
metadata <- read_excel("./41591_2021_1655_MOESM4_ESM.xlsx")
link <- read.csv("./SraRunTable.csv")

colnames(species) <- gsub(".profiled_metagenome","",colnames(species))

species_s <- species %>% 
  dplyr::filter(str_detect(rownames(species),"s__.*") & !str_detect(rownames(species),"t__.*")) 
species_s <- species_s %>% 
  dplyr::mutate(species=str_extract(rownames(species_s),"s__.*"))
rownames(species_s) <- NULL
species_s <- species_s %>% column_to_rownames(var="species")
rownames(species_s) <- gsub("s__","",rownames(species_s))

species <- species_s

#species <- species_s[rownames(species_s)%in%colnames(data_taxa_mgs),]
species <- species[rowSums(species!=0) > (ncol(species)*0.1),]
species <- apply(species,2,function(x){
  y <- x/sum(x)
  return(y)
})

species <- as.data.frame(t(species))

shannon <- as.data.frame(vegan::diversity(species, index="shannon"))
colnames(shannon)="shannon"

taxa_up <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/Score.feature.up1.txt")
taxa_down <- read.table("/groups/ProHuShiX/home/xiashufen/bigmeta/BS/Score.feature.down1.txt")
taxa_up_index <- na.omit(match(taxa_up$x, colnames(species))) #34个
taxa_down_index <- na.omit(match(taxa_down$x,colnames(species))) #4个

# Compute balance
balance_df2=data.frame(sample=rownames(species),balance_value=NA)
species <- species %>% mutate_all(as.numeric)
species[species==0]=min(species[species>0])
for(sample in rownames(species)){
  balance_df2$balance_value[balance_df2$sample==sample] = log(exp(mean(as.numeric(log(species[sample,taxa_up_index])),na.rm = T))) - log(exp(mean(as.numeric(log(species[sample,taxa_down_index])),na.rm = T)))  
}
balance_df2$balance_value <- as.numeric(balance_df2$balance_value)

balance_df2$sampleName <- link$Sample.Name[match(balance_df2$sample,link$Run)]
balance_df2$group <- metadata$`Best response`[match(balance_df2$sampleName,metadata$SampleID)]
# CR  PD  PR  SD 
# 7 161  68 102 
balance_df2$group <- factor(balance_df2$group,levels=c("CR","PR","SD","PD"))
balance_df2$group3 <- as.factor(ifelse(balance_df2$group=="CR"|balance_df2$group=="PR","CR/PR","SD/PD"))
ggboxplot(balance_df2, x = "group3", y = "balance_value")+
  stat_compare_means(comparisons = list(c("CR/PR","SD/PD")))

balance_df2$gender <- as.factor(metadata$sex[match(balance_df2$sampleName,metadata$SampleID)])
balance_df2$age <- as.numeric(metadata$age[match(balance_df2$sampleName,metadata$SampleID)])
balance_df2$atb <- as.factor(metadata$atb[match(balance_df2$sampleName,metadata$SampleID)])
balance_df2$shannon <- shannon$shannon[match(balance_df2$sample,rownames(shannon))]

summary(lm(balance_value~group+gender+age+atb,data=balance_df2)) 
summary(lm(balance_value~group2+gender+age+atb,data=balance_df2)) 
summary(lm(balance_value~group3+gender+age+atb,data=balance_df2)) 

balance_df2$adj_bs <- residuals(lm(balance_value~gender+age+atb,data=balance_df2))


balance_df_NSCLC <- balance_df2


#====PRJEB43119====
species <- read.table("/groups/ProHuShiX/home/share/PanCancerData/PRJEB43119/metaphlan/PRJEB43119_metaphlan_metagenome.txt",header=T)
metadata <- read_excel("/groups/ProHuShiX/home/xiashufen/bigmeta/BS_pancancer/PRJEB43119/PRJEB43119-SraRunTable.xlsx")

colnames(species) <- gsub(".profiled_metagenome","",colnames(species))
all(colnames(species[2:166]%in%metadata$Run))

species <- species %>% column_to_rownames(var="clade_name")
species_s <- species %>% 
  dplyr::filter(str_detect(rownames(species),"s__.*") & !str_detect(rownames(species),"t__.*")) 
species_s <- species_s %>% 
  dplyr::mutate(species=str_extract(rownames(species_s),"s__.*"))
rownames(species_s) <- NULL
species_s <- species_s %>% column_to_rownames(var="species")
rownames(species_s) <- gsub("s__","",rownames(species_s))
species<-species_s
species <- species_s[rowSums(species_s!=0) > (ncol(species_s)*0.1),]
# species <- species[rownames(species)%in%colnames(data_mgs),]
species <- apply(species,2,function(x){
  y <- x/sum(x)
  return(y)
})
species <- as.data.frame(t(species))

taxa_up_index <- na.omit(match(taxa_up$x, colnames(species))) #34个
taxa_down_index <- na.omit(match(taxa_down$x,colnames(species))) #14个 / 卡0.1 7个

# Compute balance
balance_df2=data.frame(sample=rownames(species),balance_value=NA)
species <- species %>% mutate_all(as.numeric)
species[species==0]=min(species[species>0])
for(sample in rownames(species)){
  balance_df2$balance_value[balance_df2$sample==sample] = log(exp(mean(as.numeric(log(species[sample,taxa_up_index])),na.rm = T))) - log(exp(mean(as.numeric(log(species[sample,taxa_down_index])),na.rm = T)))  
}
balance_df2$balance_value <- as.numeric(balance_df2$balance_value)

balance_df2$gender <- as.factor(metadata$host_sex[match(balance_df2$sample,metadata$Run)])

metadata <- separate(metadata,host_age_5yr_bin, into = c("age1", "age2"), sep = "-")
metadata$age1 <- as.numeric(metadata$age1);metadata$age2 <- as.numeric(metadata$age2)
metadata$age1[is.na(metadata$age1)] <- median(na.omit(metadata$age1))
metadata$age2[is.na(metadata$age2)] <- median(na.omit(metadata$age2))
metadata$age <- (metadata$age1+metadata$age2)/2
balance_df2$age <- metadata$age[match(balance_df2$sample,metadata$Run)]

metadata$host_body_mass_index[metadata$host_body_mass_index==0] <- median(na.omit(metadata$host_body_mass_index))
balance_df2$BMI <- as.numeric(metadata$host_body_mass_index[match(balance_df2$sample,metadata$Run)])
balance_df2$group <- as.factor(metadata$ORR[match(balance_df2$sample,metadata$Run)])
table(balance_df2$group)
# CR PD PR SD 
# 22 71 42 30

balance_df2$group <- factor(balance_df2$group,levels = c("CR","PR","SD","PD"))

summary(lm(balance_value~group+gender+age+BMI,data=balance_df2)) 
balance_df2$adj_bs <- residuals(lm(balance_value~gender+age+BMI,data=balance_df2))

ggboxplot(balance_df2, x = "group", y = "balance_value")+
  stat_compare_means(comparisons = list(c("CR","PR"),c("CR","SD"),c("CR","PD"),c("PR","SD"),c("PR","PD"),c("SD","PD")))
ggboxplot(balance_df2, x = "group", y = "adj_bs")+
  stat_compare_means(comparisons = list(c("CR","PR"),c("CR","SD"),c("CR","PD"),c("PR","SD"),c("PR","PD"),c("SD","PD")))
balance_df2$group2 <- factor(ifelse(balance_df2$group%in%c("CR","PR"),"R","NR"),levels = c("R","NR"))
table(balance_df2$group2)
ggboxplot(balance_df2, x = "group2", y = "adj_bs")+
  stat_compare_means(comparisons = list(c("R","NR")))
summary(lm(balance_value~group2+gender+age+BMI,data=balance_df2))

balance_df2$group3 <- factor(ifelse(balance_df2$group%in%c("CR","PR","SD"),"R","NR"),levels = c("R","NR"))
table(balance_df2$group3)
ggboxplot(balance_df2, x = "group3", y = "adj_bs")+
  stat_compare_means(comparisons = list(c("R","NR")))
summary(lm(balance_value~group2+gender+age+BMI,data=balance_df2))

balance_df_mela165 <- balance_df2

save(balance_df_mela,balance_df_NSCLC,balance_df_mela165,file="/groups/ProHuShiX/home/xiashufen/bigmeta/BS_pancancer/for_meta/for_meta3.RData")
