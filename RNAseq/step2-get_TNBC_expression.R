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
  # 本来是准备在broad的firehose上面下载，发现不方便下游分析就放弃了。
  # a=read.table('BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.data.txt.gz'
  #             ,sep = '\t',header = T,
  #             )
  # a[1:4,1:4]
  library(data.table)
  ## 下面的文件来自于UCSC的XENA：BRCA.htseq_counts.tsv 
  a=fread('TCGA-BRCA.htseq_counts.tsv',sep = '\t',header = T)
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
tmp=apply(tnbc_all_expr,1,function(x){
  sum(x==0) < 10
})
tnbc_all_expr=tnbc_all_expr[tmp,]
tnbc_all_expr=log2(edgeR::cpm(tnbc_all_expr)+1)
# write.csv(tnbc_all_expr,file = 'tnbc_all_expr_log2CPM.csv')

paired_p=names(table(all_p)[table(all_p)==2])
need_p=intersect(tnbc_p,paired_p)
exprSet=a[,all_p %in% need_p]
tmp=apply(exprSet,1,function(x){
  sum(x==0) < 10
})
exprSet=exprSet[tmp,]
save(exprSet,file = 'tnbc_paired_exprSet.Rdata')






