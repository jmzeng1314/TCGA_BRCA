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
# https://mp.weixin.qq.com/s/MJLEZPWqzJe4LaKRDtiZQQ
# 挑选出有PIK3CA突变的样本
# 载入下载好的突变信息文件
exprSet[1:4,1:4]
mut <- read.table('../MAF/TCGA-BRCA.mutect2_snv.tsv.gz',
                          header = T,sep = '\t',quote = '')
library(dplyr)
colnames(mut)
mut <- subset(mut,gene=='PIK3CA')
table(mut$effect)
## 取出含有'PIK3CA'信息的行 
PIK3CA_sample <- unique(sort(mut$Sample_ID))
PIK3CA_sample
