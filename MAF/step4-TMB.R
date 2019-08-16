rm(list = ls()) 
options(stringsAsFactors = F)  
## UCSC xena source 
BRCA.mutect2 = read.table( 'TCGA-BRCA.mutect2_snv.tsv.gz',sep = '\t',header = T)
colnames(BRCA.mutect2)
read.table()
head(BRCA.mutect2)
BRCA.mutect2$pos=paste0(BRCA.mutect2$chrom,':',
                        BRCA.mutect2$start,'-',
                        BRCA.mutect2$end)
TMB=as.numeric(table(BRCA.mutect2$Sample_ID)/38)
dat=data.frame(TMB=log2(TMB+1),BRCA='BRCA')
fivenum(dat$TMB)
library(ggpubr)
ggviolin(dat,x = 'BRCA',  y = 'TMB',ylab = 'log2(TMB+1)',xlab='TCGA')

gl=read.table('gene_length.human.txt')
head(gl) 

TMB=as.data.frame(table(BRCA.mutect2$Sample_ID)/38)
head(TMB)
allgenes=unique(BRCA.mutect2$gene)
tmp=lapply(seq(10,300,by=10), function(size){
  as.numeric(lapply(1:1000, function(x){
    cg=sample(allgenes,size)
    panle_length=sum(gl[gl[,1] %in% cg,2])/1000000
    small.BRCA.mutect2=BRCA.mutect2[BRCA.mutect2$gene %in% cg,]
    small.TMB=as.data.frame(table(small.BRCA.mutect2$Sample_ID)/panle_length)
    comp=merge(small.TMB, TMB,by='Var1')
    cor(comp[,2],comp[,3])
  }))
})
dat=do.call(rbind,tmp)
rownames(dat)=seq(10,300,by=10)
apply(dat, 1, mean )
dat=t(dat)
plot(apply(dat, 2, mean ))
library(reshape2)
df=melt(dat)[,2:3]
colnames(df)=c('number_of_genes','cor')
library(ggpubr)
ggviolin(df,x = 'number_of_genes',  y = 'cor',
         add = "boxplot")

