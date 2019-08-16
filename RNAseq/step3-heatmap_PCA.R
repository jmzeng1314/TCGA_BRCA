## 
### ---------------
###
### Create: Jianming Zeng 
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-08-10  First version
### Update Log: 2019-08-16  second version codes (R version 3.5.1 (2018-07-02))
###
### ---------------


rm(list = ls())
options(stringsAsFactors = F)
load(file = 'tnbc_paired_exprSet.Rdata')
group_list=ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,
                  'tumor','normal')
table(group_list)
dat=exprSet
# 每次都要检测数据
dat[1:4,1:4]
## 下面是画PCA的必须操作，需要看说明书。
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list)
library("FactoMineR")
library("factoextra") 
# The variable group_list (index = ncol(dat)) is removed
# before PCA analysis
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('all_samples_PCA.png')

rm(list = ls())
options(stringsAsFactors = F)
load(file = 'tnbc_paired_exprSet.Rdata')
group_list=ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,
                  'tumor','normal')
table(group_list)
dat=log(edgeR::cpm(exprSet)+1)
# 每次都要检测数据
dat[1:4,1:4]
cg=names(tail(sort(apply(dat,1,sd)),1000))# mad 
mat=dat[cg,]
library(pheatmap)
pheatmap(mat,show_colnames =F,show_rownames = F)
n=t(scale(t(mat)))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
## 可以看到TNBC具有一定的异质性，拿它来区分乳腺癌亚型指导临床治疗还是略显粗糙。
## 首先TNBC本身就可以继续分组，其次那些不是TNBC的病人跟T属于NBC病人的表达量是无法区分的。
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top1000_sd.png')






