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
a=read.table('TCGA-BRCA.survival.tsv.gz',header = T,sep = '\t')
a=read.table('TCGA-BRCA.GDC_phenotype.tsv.gz',header = T,sep = '\t',quote = '')
(tmp=as.data.frame(colnames(a)))
tmp=a[,grepl('her2',colnames(a))]
table(a$breast_carcinoma_estrogen_receptor_status)
table(a$breast_carcinoma_progesterone_receptor_status)
table(a$lab_proc_her2_neu_immunohistochemistry_receptor_status)

# 这个文件在2018年4月Immunity杂志上发表了文章The Immune Landscape of Cancer 附件拿到
b=read.table('TCGA-PAM50-subtype.txt',sep = '\t',header = T)
table(b$TCGA.Subtype)

# 首先简化列名
dat=data.frame(id=a$submitter_id,
               er=a$breast_carcinoma_estrogen_receptor_status,
               pr=a$breast_carcinoma_progesterone_receptor_status,
               her2=a$lab_proc_her2_neu_immunohistochemistry_receptor_status)
dat =dat[with(dat,er %in% c('Negative','Positive')) &
with(dat,pr %in% c('Negative','Positive')) &
with(dat,her2 %in% c('Negative','Positive')),]
table(dat)
# 过滤不确定的IHC信号后， 剩下 747个病人。
# 判断临床亚型
# er或pr为+，且her为+，则IHC_sub为B
# er或pr为+，且her为-，则IHC_sub为A
# er、pr都为-，且her为+，则IHC_sub为HER2过表达
# er、pr、her都为-，则IHC_sub为TNBC
dat$IHC_sub[(dat$er=='Positive' | dat$pr=='Positive') & dat$her2=='Negative'] = 'HR+/HER2–'
dat$IHC_sub[(dat$er=='Positive' | dat$pr=='Positive') & dat$her2=='Positive'] = 'HR+/HER2+'
dat$IHC_sub[(dat$er=='Negative' & dat$pr=='Negative') & dat$her2=='Positive'] = 'HER2'
dat$IHC_sub[(dat$er=='Negative' & dat$pr=='Negative') & dat$her2=='Negative'] = 'TNBC'
table(dat$IHC_sub)
d=merge(dat,b,by.x='id',by.y="TCGA.Participant.Barcode")
colnames(d)
library(ggstatsplot) 
ggpiestats(data = d,
           main = TCGA.Subtype,
           condition = IHC_sub,
           palette = "Set1"
)

