##################################
# Mutational Signatures Analysis #
##################################
rm(list=ls())
require(maftools)
## 我在上层文件夹的MAF里面读取过这个突变文件：
if(F){
  
  #### read maf files
  maf <- "./raw_data/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz"
  laml <- read.maf(maf)
}
## 所以直接加载这个结果即可，避免重复读取文件浪费时间。
load('../MAF/BRCA_Maf_input.Rdata')

suppressPackageStartupMessages(library("deconstructSigs"))
suppressPackageStartupMessages(library("BSgenome"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
options(stringsAsFactors = F)


mut=laml@data
SNP=mut[mut$Variant_Type=='SNP',]
SNP=mut[,c(16,5,6,11,13)]
colnames(SNP)=c("Sample","chr", "pos","ref",  "alt")
SNP$Sample=as.character(SNP$Sample)
# counts the frequency of 96 substitution mutations per sample
sigs.input <- mut.to.sigs.input(mut.ref = SNP, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

# Determines the composition of signatures in each sample
sig <-lapply(unique(SNP$Sample), function(i){
  ## signatures.cosmic signatures.nature2013
  signatures = whichSignatures(tumor.ref = sigs.input, 
                               signatures.ref = signatures.cosmic, 
                               sample.id =  i, 
                               contexts.needed = TRUE,
                               tri.counts.method = 'default')
  return(signatures$weights)
})
sig <- do.call(rbind,sig)
save(sig, file = "./Rdata/brca_cosmic_signatures.Rdata")

library(pheatmap)
png('./pictures/Cosmic_signatures_BRCA_heatmap.png' ,res = 150,width = 1080,height = 643)
pheatmap(t(sig),cluster_rows = F,cluster_cols = T, fontsize = 9,
         show_colnames = F, main = "COSMIC Mutational Signatures in BRCA") 
dev.off()
#signature 1 accounts for a large part in BRCA

