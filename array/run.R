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
# The current release is V7 including 11,688 samples, 53 tissues and 714 donors.
if(F){
  a=read.table('BRCA.medianexp.txt.gz',header = T,sep = '\t',quote = '')
  a[1:4,1:4]
  a=a[-1,]
  rownames(a)=a[,1]
  a=a[,-1]
  genes=rownames(a)
  a=apply(a, 2, as.numeric)
  a=as.data.frame(a)
  rownames(a)= genes
  a=na.omit(a)
  colnames(a)
  group_list=ifelse(as.numeric(substr(colnames(a),14,15)) < 10,'tumor','normal')
  table(group_list)
  breat_array_tcga=a
  save(breat_array_tcga,group_list,file = 'breat_array_tcga.Rdata')
}
load(file = 'breat_array_tcga.Rdata')
# dat=log2(edgeR::cpm(breat_gtex)+1)
dat=breat_array_tcga
dat[1:4,1:4]
if(T){
  ddata=t(dat)
  ddata[1:4,1:4]
  s=colnames(ddata);head(s)
  library(org.Hs.eg.db)
  s2g=toTable(org.Hs.egSYMBOL)
  g=s2g[match(s,s2g$symbol),1];head(g)
  #  probe Gene.symbol Gene.ID
  dannot=data.frame(probe=s,
                    "Gene.Symbol" =s, 
                    "EntrezGene.ID"=g)
  ddata=ddata[,!is.na(dannot$EntrezGene.ID)]
  dannot=dannot[!is.na(dannot$EntrezGene.ID),] 
  head(dannot)
  library(genefu)
  # c("scmgene", "scmod1", "scmod2","pam50", "ssp2006", "ssp2003", "intClust", "AIMS","claudinLow")
  
  s<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                         annot=dannot,do.mapping=TRUE)
  table(s$subtype)
  tmp=as.data.frame(s$subtype)
  subtypes=as.character(s$subtype)
}

library(genefu)
pam50genes=pam50$centroids.map[c(1,3)]
pam50genes[pam50genes$probe=='CDCA1',1]='NUF2'
pam50genes[pam50genes$probe=='KNTC2',1]='NDC80'
pam50genes[pam50genes$probe=='ORC6L',1]='ORC6'
x=dat
x=x[pam50genes$probe[pam50genes$probe %in% rownames(x)] ,]
tmp=data.frame(group=group_list,subtypes=subtypes)
rownames(tmp)=colnames(x)
library(pheatmap)


pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_raw.png') 


x=t(scale(t(x)))
x[x>1.6]=1.6
x[x< -1.6]= -1.6

pheatmap(x,show_rownames = T,show_colnames = F,
         annotation_col = tmp,
         filename = 'ht_by_pam50_scale.png') 






