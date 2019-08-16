rm(list = ls())
options(stringsAsFactors = F)
if(F){ 
  library(data.table)
  ## 下面的文件来自于UCSC的XENA：BRCA.htseq_counts.tsv 
  f='/Users/jmzeng/data/public/TCGA/TCGA-BRCA.htseq_counts.tsv'
  a=fread(f,sep = '\t',header = T)
  a=as.data.frame(a)
  a[1:4,1:4] 
  rownames(a)=a[,1]
  a=a[,-1]
  genes=rownames(a)
  a[1:4,1:4] 
  a=2^a-1
  a[1:4,1:4] 
  save(a,file = 'TCGA-BRCA.htseq_counts.Rdata')
}
load(file = 'TCGA-BRCA.htseq_counts.Rdata')
RNAseq_expr=a
RNAseq_expr=RNAseq_expr[apply(RNAseq_expr,1, function(x) sum(x>1) > 50),]
dim(RNAseq_expr)
RNAseq_expr=log(edgeR::cpm(RNAseq_expr)+1)
head(colnames(RNAseq_expr))

load('../MAF/math.Rdata')
math=as.data.frame(unlist(math))
rownames(math)=substring(rownames(math),1,16)
colnames(math)='math'
RNAseq_expr=RNAseq_expr[,colnames(RNAseq_expr) %in% rownames(math)]
math=math[colnames(RNAseq_expr),,drop=F]
dim(RNAseq_expr)
tmp = apply(RNAseq_expr,1,function(x){
  cor(x,math[,1])
}) 
tmp=as.data.frame(tmp)

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


deg=tmp
ensemb_genes=rownames(deg)
deg$ensemb_genes=rownames(deg)
library(stringr)
class(str_split(ensemb_genes,'[.]',simplify = T))
class(unlist(str_split(ensemb_genes,'[.]')))
deg$ensembl_id=str_split(ensemb_genes,'[.]',simplify = T)[,1]

library(org.Hs.eg.db)
g2s=toTable(org.Hs.egSYMBOL)
g2e=toTable(org.Hs.egENSEMBL)
b=merge(deg,g2e,by='ensembl_id',all.x=T)
d=merge(b,g2s,by='gene_id',all.x=T)

table(d$ensembl_id)[table(d$ensembl_id)>1]

d=d[order(d$ensemb_genes),]
d=d[!duplicated(d$ensemb_genes),]
d=d[match(ensemb_genes,d$ensemb_genes),]


deg=d
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db) 
diffg= deg[abs(deg$tmp)>0.2,1] 
diffg=as.numeric(na.omit(diffg))
kk  <- enrichKEGG(gene         = diffg,
                  organism     = 'hsa', 
                  pvalueCutoff = 0.9,
                  qvalueCutoff =0.9)
head(kk)[,1:6]
dotplot(kk)
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keytype='ENTREZID')
write.csv(kk@result,'kk.csv')

deg=na.omit(deg)
geneList=deg$tmp
names(geneList)=deg$gene_id
geneList=sort(geneList,decreasing = T)
###  GSEA 
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
head(kk_gse)[,1:6]
dotplot(kk_gse)
kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keytype='ENTREZID')
write.csv(kk_gse@result,'kk_gse.csv')
gseaplot(kk_gse, geneSetID = rownames(kk_gse[4,]))















