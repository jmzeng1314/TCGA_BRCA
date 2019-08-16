## PIK3CA突变与BRCA的免疫特征

下载[pan-cancer的免疫特征数据](https://www.cell.com/immunity/fulltext/S1074-7613(18)30121-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1074761318301213%3Fshowall%3Dtrue#secsectitle0435)TableS1，从中提取出BRCA的免疫信息，与临床信息合并。

```{R}
library(readxl)
pan_immune <- read_excel("./raw_data/Pan Immune Feature Matrix of Immune Characteristics.xlsx")
brca_im <- pan_immune[pan_immune$`TCGA Study`=="BRCA",]
brca_phe <- read.delim("./raw_data/TCGA-BRCA.GDC_phenotype.tsv.gz", stringsAsFactors=FALSE)

###提取需要的表型信息
#表型和免疫共有患者1087个
phe <- brca_phe[brca_phe$sample_type.samples == "Primary Tumor", ]
phe$ID <- substr(phe$submitter_id.samples, 1,12)
phe <- phe[!duplicated(phe$ID),]
rownames(phe) <- phe$ID
phe <- phe[brca_im$`TCGA Participant Barcode`,]

##合并表型与免疫信息
brca_phe_im <- merge(brca_im, phe, by.x = "TCGA Participant Barcode", by.y = "ID")

save(brca_phe_im, file = "./Rdata/phe_immune.Rdata")
```

#### 免疫特征与肿瘤PIK3CA突变的关系

TCGA_BRCA的突变数据[下载地址](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2Fmutation_curated_wustl_gene&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

```{R}
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
```

![Snipaste_2018-12-30_21-48-43.png](https://i.loli.net/2019/01/03/5c2e0ce592158.png)

PIK3CA突变与proliferation相关性最高，与homologous recombination repair相关性也比较高。

![HRD_VS_PIK3CA_mutation.png](https://i.loli.net/2019/01/03/5c2e0d9445002.png)

 



