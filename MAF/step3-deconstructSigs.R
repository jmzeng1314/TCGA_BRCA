require(maftools) 
options(stringsAsFactors = F) 
if(F){
  # 首先需要更改一些镜像配置
  options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
  BiocManager::biocLite('maftools')
  BiocManager::biocLite('deconstructSigs')
  BiocManager::biocLite('BSgenome.Hsapiens.UCSC.hg38')
}
library(maftools)
load(file = 'BRCA_Maf_input.Rdata')
suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
options(stringsAsFactors = F)
mut=laml@data
head(mut)
mut=mut[mut$Variant_Type=='SNP',]
a=mut[,c(16,5,6,12,13)]
colnames(a)=c( "Sample","chr", "pos","ref",  "alt")

a$Sample=as.character(a$Sample)

plot(table(a$Sample),las=2)
sigs.input <- mut.to.sigs.input(mut.ref = a, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
class(sigs.input)
barplot(as.numeric(sigs.input[1,]))
head(t(sigs.input))
w=lapply(unique(a$Sample)[1:10], function(i){
  ## signatures.cosmic signatures.nature2013
  sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')
  print(i)
  return(sample_1$weights)
})
w=do.call(rbind,w)
library(pheatmap)
pheatmap(t(w),cluster_rows = F,cluster_cols = T)
pheatmap(w,cluster_rows = T,cluster_cols = F)

# Determine the signatures contributing to the two example samples
lapply(unique(a$Sample), function(i){
  ## signatures.cosmic signatures.nature2013
  sample_1 = whichSignatures(tumor.ref = sigs.input, 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')
  pdf(paste0(i,'.sig.pdf'))
  plotSignatures(sample_1, sub = i)
  dev.off()
})




