rm(list = ls())
options(stringsAsFactors = F)
if(F){
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
  
  
  
  
  
}
load(file = "./Rdata/phe_immune.Rdata")
#表型和免疫共有患者1087个
brca_phe_im[1:4,1:8]
library(pheatmap)
grep('Macro',colnames(brca_phe_im))
dat=brca_phe_im[,5:64]
dat=apply(dat,2,as.numeric)
dat=dat[,colMeans(abs(dat),na.rm = T) <1]
pheatmap(dat,show_rownames = F,show_colnames = F) 
dat=scale(dat)
dat[dat>2]=2
dat[dat< -2] = -2 
pheatmap(dat,show_rownames = F,show_colnames = F) 
## 可以看到这些免疫指标分成3类！
im_group=as.data.frame(cutree(hclust(dist(t(dat))),3))
im_group
table(im_group)


###不同免疫指标间确实是相关的，TCGA免疫图谱原文中，作者就是计算了160中免疫特征间的相关系数，画热图，找到了5个模块中的核心特征，进一步进行的免疫分型。
#文章中的图片https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr1.jpg
##补充数据中只给出了这五十多种免疫特征的数据。

#计算BRCA的免疫特征间相关系数，聚类
data <- na.omit(dat)
cor_matrix <- cor(data)
pheatmap(cor_matrix, show_colnames = F) #分为3类，大概有2个比较明显的module
group  <- cutree(hclust(dist(cor_matrix)),3)
table(group)
names(group[group==3]) 
#生存及生存时间在一组
#HRD与would healing, proliferation, Intratumor Heterogeneity等相关性较高，与生存不同组。
#与之前的结果相符：PIK3CA突变不影响生存，与HRD,增殖等相关性高。

