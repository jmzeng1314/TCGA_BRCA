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
## 根据TCGA样本的命名可以区分正常组织和肿瘤样本的测序结果
group_list=ifelse(as.numeric(substr(colnames(exprSet),14,15)) < 10,
                  'tumor','normal')
exprSet[1:4,1:4]
# strsplit('ENSG00000000003.13','[.]')[[1]][1]
if(T){
  a=data.frame(V1=rownames(exprSet))
  library(stringr)
  class(str_split(a$V1,'[.]',simplify = T))
  class(unlist(str_split(a$V1,'[.]')))
  a$ensembl_id=str_split(a$V1,'[.]',simplify = T)[,1]
  
  library(org.Hs.eg.db)
  g2s=toTable(org.Hs.egSYMBOL)
  g2e=toTable(org.Hs.egENSEMBL)
  b=merge(a,g2e,by='ensembl_id',all.x=T)
  d=merge(b,g2s,by='gene_id',all.x=T)
  
  table(d$ensembl_id)[table(d$ensembl_id)>1]
  
  d=d[order(d$V1),]
  d=d[!duplicated(d$V1),]
  d=d[match(a$V1,d$V1),]
  head(d)
  length(unique(d$symbol))
  length(unique(d$ensembl_id))
  
  dat=exprSet
  dat[1:4,1:4] 
  ids=d
  ids$median=apply(dat,1,median)
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  ids=ids[!duplicated(ids$symbol),]
  dat=dat[ids$V1,]
  ## 没时间搞清楚为什么这里有个基因名字是NA了。
  ids$symbol[is.na(ids$symbol)]='na'
  rownames(dat)=ids$symbol
  dat[1:4,1:4]  
}
exprSet=dat
exprSet[1:4,1:4]


table(group_list)
# https://github.com/jmzeng1314/tcga_example/blob/master/scripts/step02-DEG-3-packages.R
# https://mp.weixin.qq.com/s/hPcqS6M-d2Bun2DDHLvvLg

{
  ## 方法一：DESeq2
  if(T){
    library(DESeq2)
    
    (colData <- data.frame(row.names=colnames(exprSet), 
                           group_list=group_list) )
    dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                  colData = colData,
                                  design = ~ group_list)
    tmp_f=file.path(Rdata_dir,'TCGA-KIRC-miRNA-DESeq2-dds.Rdata')
    if(!file.exists(tmp_f)){
      dds <- DESeq(dds)
      save(dds,file = tmp_f)
    }
    load(file = tmp_f)
    res <- results(dds, 
                   contrast=c("group_list","tumor","normal"))
    resOrdered <- res[order(res$padj),]
    head(resOrdered)
    DEG =as.data.frame(resOrdered)
    DESeq2_DEG = na.omit(DEG)
    
    nrDEG=DESeq2_DEG[,c(2,6)]
    colnames(nrDEG)=c('log2FoldChange','pvalue')  
    draw_h_v(exprSet,nrDEG,'DEseq2',group_list,1)
  }
  
  ### ---------------
  ###
  ### Then run edgeR 
  ###
  ### ---------------
  if(T){
    library(edgeR)
    d <- DGEList(counts=exprSet,group=factor(group_list))
    keep <- rowSums(cpm(d)>1) >= 2
    table(keep)
    d <- d[keep, , keep.lib.sizes=FALSE]
    d$samples$lib.size <- colSums(d$counts)
    d <- calcNormFactors(d)
    d$samples
    dge=d
    design <- model.matrix(~0+factor(group_list))
    rownames(design)<-colnames(dge)
    colnames(design)<-levels(factor(group_list))
    dge=d
    dge <- estimateGLMCommonDisp(dge,design)
    dge <- estimateGLMTrendedDisp(dge, design)
    dge <- estimateGLMTagwiseDisp(dge, design)
    
    fit <- glmFit(dge, design)
    # https://www.biostars.org/p/110861/
    lrt <- glmLRT(fit,  contrast=c(-1,1)) 
    nrDEG=topTags(lrt, n=nrow(dge))
    nrDEG=as.data.frame(nrDEG)
    head(nrDEG)
    edgeR_DEG =nrDEG 
    nrDEG=edgeR_DEG[,c(1,5)]
    colnames(nrDEG)=c('log2FoldChange','pvalue') 
    # draw_h_v(exprSet,nrDEG,'edgeR',group_list,1)
    
  }
  
  
  ### ---------------
  ###
  ### Lastly run voom from limma
  ###
  ### --------------- 
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
    cont.matrix=makeContrasts(contrasts=c('tumor-normal'),levels = design)
    fit2=contrasts.fit(fit,cont.matrix)
    fit2=eBayes(fit2)
    
    tempOutput = topTable(fit2, coef='tumor-normal', n=Inf)
    DEG_limma_voom = na.omit(tempOutput)
    head(DEG_limma_voom)
    nrDEG=DEG_limma_voom[,c(1,4)]
    colnames(nrDEG)=c('log2FoldChange','pvalue') 
    draw_h_v(exprSet,nrDEG,'limma',group_list,1)
    
  }
  
  
  ## 比较一下这三个差异分析的结果
  
  nrDEG1=DEG_limma_voom[,c(1,4)] 
  colnames(nrDEG1)=c('log2FoldChange','pvalue')  
  nrDEG2=edgeR_DEG[,c(1,5)] 
  colnames(nrDEG2)=c('log2FoldChange','pvalue')  
  nrDEG3=DESeq2_DEG[,c(2,6)] 
  colnames(nrDEG3)=c('log2FoldChange','pvalue')  
  
  mi=unique(c(rownames(nrDEG1),rownames(nrDEG1),rownames(nrDEG1)))
  
  lf=data.frame(limma=nrDEG1[mi,1], 
                edgeR=nrDEG2[mi,1], 
                DESeq2=nrDEG3[mi,1]) 
  cor(na.omit(lf))
}



