#######################
# Infer Heterogeneity #
#######################
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
laml@data$t_vaf <- laml@data$t_alt_count/laml@data$t_depth ##compute variant allele frequencies(t_vaf)
laml@clinical.data$Tumor_Sample_Barcode <- substr(laml@clinical.data$Tumor_Sample_Barcode,1,16)
laml@data$Tumor_Sample_Barcode <- substr(laml@data$Tumor_Sample_Barcode,1,16)

### infer heterogenity by clustering variant allele frequencies
samples <- as.data.frame(table(laml@data$Tumor_Sample_Barcode))
samples <- as.character(samples[samples$Freq > 10,1]) #samples with more than 10 mutaions
head(samples)
##Ignore variants located on copy-number altered regions
#Copy number alterations results in abnormally high/low variant allele frequencies, 
# which tends to affect clustering. 
## 这里需要写清楚，描述 每个病人的segment文件是如何产生的
seg <- "./raw_data/BRCA_CNV/CNVs_selected_by_maf/masked_segment" #dir of segment files

heterogeneity <- lapply(samples,function(tsb){
  # tsb <- samples[1]
  file=file.path(seg, paste0(tsb,".hg38.seg.txt"))
  if(file.exists(file)){segFile=file}else{segFile=NULL}
  het=inferHeterogeneity(maf=laml, tsb = tsb , segFile=segFile, vafCol="t_vaf")
  mathd=het$clusterMeans 
  mathd=mathd[mathd$cluster !='outlier',] #remove outliers
  mathd=mathd[mathd$cluster !="CN_altered",]  #Ignoring variants in copy number altered regions
  mathd=as.data.frame(table(mathd$Tumor_Sample_Barcode))  #count clusters
  mathd$MATH <- het$clusterData$MATH[1] 
  return(mathd)
})
heterogeneity <- do.call(rbind,heterogeneity)
colnames(heterogeneity) <- c("Sample_ID","clusters","MATH")

save(heterogeneity,file="./Rdata/brca_hterogeneity.Rdata")

load(file="./Rdata/brca_hterogeneity.Rdata")


