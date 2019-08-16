#####################################
# Some analysis about heterogeneity #
#####################################
rm(list = ls())

load(file = "./Rdata/phe_immune.Rdata")
load(file = "./Rdata/brca_hterogeneity.Rdata")

im_phe <- brca_phe_im
im_phe <- im_phe[,c(4:64)]
im_phe$Sample_ID <- brca_phe_im$submitter_id.samples

info <- merge(heterogeneity, im_phe, by = "Sample_ID")

### the correlation between "MATH" and "Intratumor Heterogeneity"
# the "Intratumor Heterogeneity" is inferred through ABSOLUTE
hets <- info[,c(3,7)]
hets$`Intratumor Heterogeneity`=as.numeric(hets$`Intratumor Heterogeneity`)
hets <- na.omit(hets)
cor(as.numeric(hets$`Intratumor Heterogeneity`), hets$MATH) # cor=0.3, is this coefficient too small?
## 这个时候需要判断，这两个异质性指标，和临床信息的关联，才能判断哪一个更优。
###the relationship between immune/genomic features and heterogeneity in BRCA
#The tumor microenvironment has been reported to contribute to tumour heterogeneity.
dat=info[,c(2,3,5:64)]
dat=apply(dat,2,as.numeric)
data <- na.omit(dat)
# 这些指标的度量衡不一样，可能并不能直接cor它们
cor_matrix <- cor(data)
library(pheatmap)
pheatmap(cor_matrix, show_colnames = F) 
group  <- cutree(hclust(dist(cor_matrix)),5)
table(group)
names(group[group==2])
names(group[group==1])
##As we expected, the "MATH","Intratumor Heterogeneity","SNV Neoantigens","Nonsilent Mutation Rate","Number of Segments", "Aneuploidy Score"... are always in one cluster.
##The "clusters"(number of clones) is more related to survival.
