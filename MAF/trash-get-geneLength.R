

# 参考我GitHub项目：https://github.com/jmzeng1314/scRNA_smart_seq2/blob/master/RNA-seq/step7-counts2rpkm.R
# 获取基因长度。
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## 下面是定义基因长度为 非冗余exon长度之和, 即WES的基因长度
if(F){
  exon_txdb=exons(txdb)
  genes_txdb=genes(txdb)
  genes_txdb
  ?GRanges
  # 因为有些基因之间有overlap，所以这个并不是最标准答案。
  o = findOverlaps(exon_txdb,genes_txdb)
  o
  ## exon - 1 : chr1 4807893-4807982
  ## 1        6523
  #  genes_txdb[6523]  # chr1 4807893-4846735 , 18777
  t1=exon_txdb[queryHits(o)]
  t2=genes_txdb[subjectHits(o)]
  t1=as.data.frame(t1)
  t1$geneid=mcols(t2)[,1]
  # 如果觉得速度不够，就参考R语言实现并行计算
  # http://www.bio-info-trainee.com/956.html
  #lapply : 遍历列表向量内的每个元素，并且使用指定函数来对其元素进行处理。返回列表向量。
  #函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组;
  #它的返回值是一个列表，代表分组变量每个水平的观测。
  g_l = lapply(split(t1,t1$geneid),function(x){
    # x=split(t1,t1$geneid)[[1]]
    head(x)
    tmp=apply(x,1,function(y){
      y[2]:y[3]
    })
    length(unique(unlist(tmp)))
    # sum(x[,4])
  })
  head(g_l)
  g_l=data.frame(gene_id=names(g_l),length=as.numeric(g_l))
  head(g_l)
  save(g_l,file = 'gene_length_of_hg38.Rdata')
}
load(file = 'gene_length_of_hg38.Rdata')
sum(g_l$length)




gl=read.table('gene_length.human.txt')
head(gl)
colnames(gl)=c('symbol','length_CCDS')
load(file = 'gene_length_of_hg38.Rdata')
head(g_l)
colnames(g_l)=c('gene_id', 'length_R')
library(org.Hs.eg.db)
columns(org.Hs.eg.db)
e2s=toTable(org.Hs.egSYMBOL)
head(e2s)
g_l=merge(g_l,e2s,by='gene_id')
comp=merge(g_l,gl,by='symbol')
comp[,3]=log(comp[,3])
comp[,4]=log(comp[,4])
plot(comp[,3:4])
head(comp)
library(ggpubr)
p=ggscatter(comp,'length_R','length_CCDS', 
          color = "black", shape = 21, size = 0.3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)
p+geom_abline(intercept = 0, slope = 1, color="red", 
            linetype="dashed", size=1.5)



library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
tx_by_gene = transcriptsBy(txdb, by="gene")
gene_lens = as.data.frame(max(width(tx_by_gene)))
gene_lens$gene_id=rownames(gene_lens)
colnames(gene_lens)=c('length_R','gene_id')

library(org.Hs.eg.db)
columns(org.Hs.eg.db)
e2s=toTable(org.Hs.egSYMBOL)
head(e2s)
gene_lens=merge(gene_lens,e2s,by='gene_id')
comp=merge(gene_lens,gl,by='symbol')
comp[,3]=log(comp[,3])
comp[,4]=log(comp[,4])
plot(comp[,3:4])
head(comp)
library(ggpubr)
p=ggscatter(comp,'length_R','length_CCDS', 
            color = "black", shape = 21, size = 0.3, # Points color, shape and size
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)
p+geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1.5)

