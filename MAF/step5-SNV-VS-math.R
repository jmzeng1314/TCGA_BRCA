rm(list = ls())
require(maftools) 
options(stringsAsFactors = F) 
load(file = 'math.Rdata')
math=as.data.frame(unlist(math))
rownames(math)=substring(rownames(math),1,16)
colnames(math)='math'
load(file = 'BRCA_Maf_input.Rdata')

laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
mut=laml@data[,c("Hugo_Symbol","Chromosome","Start_Position","Tumor_Sample_Barcode","t_vaf")]
mut$pos=paste(mut$Chromosome,mut$Start_Position,sep=':')
mut$Tumor_Sample_Barcode=substring(mut$Tumor_Sample_Barcode,1,16)
mut=as.data.frame(mut)
gs=unique(mut[,c(1,4)])
topg=names(head(sort(table(gs$Hugo_Symbol),decreasing = T),50))
lapply(topg, function(g){
  #g=topg[1]
  mutSamples=unique(mut[mut$Hugo_Symbol == g ,4])
  # mutSamples=unique(mut[mut$Hugo_Symbol %in% c('TP53','CDH1') ,4])
  mutSamples=mutSamples[mutSamples %in% rownames(math)]
  mutorNot=ifelse(rownames(math) %in% mutSamples,'mut','wild')
  this_math=cbind(math,mutorNot)
  library(ggpubr)
  p <- ggboxplot(this_math, x = "mutorNot", y = "math",
                 color = "mutorNot", 
                 add = "jitter", shape = "mutorNot")
  p+stat_compare_means(method = 't.test')
  p
  ggsave(filename = paste0('t.test.',g,'.png'))
})

## only TP53 and CDH1 show significance in all of the BRCA cohort.
 





