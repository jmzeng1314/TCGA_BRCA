rm(list = ls())
### PIK3CA突变与否对免疫特征的影响
load("./Rdata/phe_immune.Rdata")
mutaion_all <- read.delim("./raw_data/mutation_curated_wustl_gene.gz", stringsAsFactors=FALSE)
brca_phe_im$`TCGA Participant Barcode` <- gsub("-",".",brca_phe_im$`TCGA Participant Barcode`)
PIK3CA_all <- as.vector(mutaion_all[mutaion_all$sample == "PIK3CA",-1]) 
names(PIK3CA_all) <- substr(names(PIK3CA_all),1,12)
PIK3CA_mut <- names(PIK3CA_all)[PIK3CA_all==1] #1代表发生非沉默突变

phe <- brca_phe_im[brca_phe_im$`TCGA Participant Barcode` %in% names(PIK3CA_all), 1:64]
phe$PIK3CA <- ifelse(phe$`TCGA Participant Barcode` %in% PIK3CA_mut, "MUT", "WT")
table(phe$PIK3CA )
phe$PIK3CA <- as.numeric(as.factor(phe$PIK3CA))
phe[,5:64] <- apply(phe[,5:64], 2, function(x){as.numeric(x)})
#处理缺失值
n <- apply(phe, 1, function(x){sum(is.na(x))}) #每个患者的缺失值数目
table(n) 
phe <- phe[which(n<15),]
#计算PIK3CA突变与免疫features的相关性是否显著（除去34~37列的生存信息）
p_mut <- apply(phe[,c(5:32,37:64)], 2, function(x){
  t <- cor.test(x, phe$PIK3CA, method = "spearman") 
  return(t$p.value)})
p_mut <- sort(p_mut)
table(p_mut<0.01)
f_mut <- names(p_mut[p_mut<0.05]) #23个features

save(p_mut, PIK3CA_mut, phe, file = "./Rdata/PIK3CA_mut_immune.Rdata")

#相关矩阵图
df <- phe[,f_mut]
rownames(df) <- phe$`TCGA Participant Barcode`
df$PIK3CA <- ifelse(rownames(df) %in% PIK3CA_mut, "MUT", "WT")

library(corrplot)
df2 <- df
df2$PIK3CA <- as.numeric(as.factor(df$PIK3CA))
df2 <- na.omit(df2)
mycor <- cor(df2)
corrplot(mycor,method="color",type="upper",order="hclust",addCoef.col = "black",number.cex=0.52,tl.col = "black",tl.cex = 0.72)

###画图
#boxplot
library(ggpubr)
ggboxplot(df, x="PIK3CA", y="`Homologous Recombination Defects`", 
          color = "PIK3CA", palette = "jco", add = "jitter") +
  stat_compare_means()