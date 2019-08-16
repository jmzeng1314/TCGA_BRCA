rm(list = ls())
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
project='TCGA_BRCA_all'

if(T){
  
  png(paste0('plotmafSummary_',project,'.png'),res = 150,width = 1080,height = 1080)
  plotmafSummary(maf = laml, rmOutlier = TRUE,showBarcodes = T,
                 addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  
  ## ---- fig.align='left',fig.height=5,fig.width=10, fig.align='left'-------
  #We will draw oncoplots for top ten mutated genes.
  png(paste0('oncoplot_top30_',project,'.png') ,res = 150,width = 1080,height = 1080)
  oncoplot(maf = laml, top = 30, fontSize = 12 ,showTumorSampleBarcodes = F )
  dev.off()
  
  
  png(paste0('TMB_',project,'.png') ,res = 150,width = 1080,height = 1080)
  laml.mutload = tcgaCompare(maf = laml, cohortName =project)
  dev.off()
  
  png(paste0('VAF_',project,'.png') ,res = 150,width = 1080,height = 1080)
  plotVaf(maf = laml, top = 20)
  dev.off()
  
}
head(filter_phe)
rownames(filter_phe)=filter_phe$Tumor_Sample_Barcode
png(paste0('oncoplot_top30_','IHC_markers','.png') ,res = 150,width = 1080,height = 1080)
oncoplot(maf = laml, top = 30, fontSize = 12 ,
         annotationDat=filter_phe,clinicalFeatures=c('ER','PR','HER2'),sortByAnnotation=T,
         showTumorSampleBarcodes = F )
dev.off()
filter_phe$tnbc=ifelse(apply(filter_phe[3:5],1, function(x) sum(x=='Negative'))==3,'TNBC','nonTNBC')
head(filter_phe)
png(paste0('oncoplot_top30_','TNBC','.png') ,res = 150,width = 1080,height = 1080)
oncoplot(maf = laml, top = 30, fontSize = 12 ,
         annotationDat=filter_phe,clinicalFeatures='tnbc',sortByAnnotation=T,
         showTumorSampleBarcodes = F )
dev.off()



if(F){
  
  dir.create(paste0('vaf_clust_',project ))
  setwd(paste0('vaf_clust_',project ))
  patients=unique(laml@data$Tumor_Sample_Barcode)
  math=lapply(patients, function(x){
    #x=unique(laml@data$Tumor_Sample_Barcode)[5]
    if(nrow(laml@data[laml@data$Tumor_Sample_Barcode==x,]) > 5 ){
      #png(paste0(x,'_vaf_clust.png'),res=120,width = 1080,height = 1080)
      het = inferHeterogeneity(maf = laml, tsb = x, vafCol = 't_vaf')
      #print(het$clusterMeans) 
      return(het$clusterData$MATH[1])
      #plotClusters(clusters = het)
      #dev.off()
    }

  }) 
  names(math)=patients
  save(math,file = 'math.Rdata')
}
setwd('../')
if(T){
  sort(table(laml@data$Tumor_Sample_Barcode),decreasing = T)
  table(table(laml@data$Tumor_Sample_Barcode)>10)
  length( table(laml@data$Tumor_Sample_Barcode))
  het=inferHeterogeneity(maf = laml, tsb = NULL,top=970) 
  mathd=het$clusterMeans
  mathd=mathd[mathd$cluster !='outlier' ,]
  mathd=as.data.frame(table(mathd$Tumor_Sample_Barcode)) 
  mathd2=unique(as.data.frame(het$clusterData)[c(5,8)])
  mathd=merge(mathd,mathd2,by.x='Var1',by.y='Tumor_Sample_Barcode')
  head(filter_phe)
  mathd$type=filter_phe[match(mathd$Var1,filter_phe$Tumor_Sample_Barcode),'tnbc']
  head(mathd)
  library(ggpubr)
  p <- ggboxplot(mathd, x = "type", y = "MATH",
                 color = "type", 
                 add = "jitter", shape = "type")
  p+stat_compare_means(method = 't.test')
   
  p <- ggboxplot(mathd, x = "type", y = "Freq",
                 color = "type", 
                 add = "jitter", shape = "type")
  p+stat_compare_means(method = 't.test')
}
#We will run mutExclusive on top 10 mutated genes. 
somaticInteractions(maf = laml, top = 25, pvalue = c(0.05, 0.1)) 

