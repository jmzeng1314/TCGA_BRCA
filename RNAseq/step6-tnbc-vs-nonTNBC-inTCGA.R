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
}
#a[1:4,1:4]
load('tnbc_s.Rdata')
tnbc_p=substring(tnbc_s,1,12)
all_p=substring(colnames(a),1,12)
tnbc_all_expr=a[,all_p %in% tnbc_p]
tnbc_all_expr[1:4,1:4]
nontnbc_all_expr=a[,!all_p %in% tnbc_p]

tnbc_all_expr=tnbc_all_expr[,substring(colnames(tnbc_all_expr),14,15) == '01']
nontnbc_all_expr=nontnbc_all_expr[,substring(colnames(nontnbc_all_expr),14,15) == '01']

exprSet=cbind(tnbc_all_expr,nontnbc_all_expr)
group_list=c(rep('tnbc',ncol(tnbc_all_expr)),rep('nontnbc',ncol(nontnbc_all_expr)))
library(edgeR)
if(T){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  
  dge <- DGEList(counts=exprSet)
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log=TRUE, prior.count=3)
  
  v <- voom(dge,design,plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  
  group_list
  cont.matrix=makeContrasts(contrasts=c('tnbc-nontnbc'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  
  tempOutput = topTable(fit2, coef='tnbc-nontnbc', n=Inf)
  DEG_limma_voom = na.omit(tempOutput)
  head(DEG_limma_voom) 
  
}

deg=DEG_limma_voom
with(deg,plot( logFC,-log10( adj.P.Val)))
diff1=rownames(deg[abs(deg$logFC)>3,])


library(pheatmap)
dat=log2(edgeR::cpm(exprSet+1))
dat[1:4,1:4]
cg=rownames(deg[abs(deg$logFC)>3,]) 
mat=dat[cg,]
library(pheatmap)
pheatmap(mat,show_colnames =F,show_rownames = F)
n=t(scale(t(mat)))
n[n>2]=2
n[n< -2]= -2
n[1:4,1:4]
#pheatmap(n,show_colnames =F,show_rownames = F)
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
pheatmap(n,show_colnames =F,show_rownames = F,
         annotation_col=ac ,filename = 'tnbc-vs-nonTNBC-TCGA-top300-DEG-heatmap.png')

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

write.csv(d,'tnbc-vs-nonTNBC-TCGA-by-limma-voom.csv')

deg=d
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db) 
diffg= deg[abs(deg$logFC)>3,1] 
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
geneList=deg$logFC
names(geneList)=deg$gene_id
geneList=sort(geneList,decreasing = T)
###  GSEA 
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
head(kk_gse)[,1:6]
dotplot(kk_gse)
kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keytype='ENTREZID')
write.csv(kk_gse@result,'kk_gse.csv')
gseaplot(kk_gse, geneSetID = rownames(kk_gse[4,]))








