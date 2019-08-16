rm(list = ls())
require(maftools) 
options(stringsAsFactors = F) 
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
colnames(laml@data)
## UCSC xena source 
BRCA.mutect2 = read.table( 'TCGA-BRCA.mutect2_snv.tsv.gz',sep = '\t',header = T)
colnames(BRCA.mutect2) 
# 为了避免转换成maf格式代码太麻烦，直接去下载GDC提供的maf文件即可
laml <- read.maf(maf = 'GDC/TCGA.BRCA.mutect.c6a029e5-0ea3-410d-9e67-360bdfee2914.DR-7.0.somatic.maf')
laml 
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
getSampleSummary(laml) 
getGeneSummary(laml) 
getFields(laml)  
colnames(laml@clinical.data)


# 接下来读取UCSC xena 的表型文件。
phe=read.table('../clinical/TCGA-BRCA.GDC_phenotype.tsv.gz',header = T,sep = '\t',quote = "",fill = T)
colnames(phe)
table(phe$metastatic_breast_carcinoma_erbb2_immunohistochemistry_level_result)

table(phe$anatomic_neoplasm_subdivision)
table(phe$breast_carcinoma_estrogen_receptor_status)
table(phe$breast_carcinoma_progesterone_receptor_status)
table(phe$lab_proc_her2_neu_immunohistochemistry_receptor_status)

table(phe$drug_name)
table(phe$margin_status)
table(phe$race.demographic)
table(phe$gender.demographic)
table(phe$state.samples)
table(phe$sample_type_id.samples )
table(phe$vital_status.diagnoses)
table(phe$tumor_stage.diagnoses )
table(phe$metastatic_site_at_diagnosis)

filter_phe=phe[c("submitter_id.samples","anatomic_neoplasm_subdivision",
                 "breast_carcinoma_estrogen_receptor_status", 
                 "breast_carcinoma_progesterone_receptor_status",
                 "lab_proc_her2_neu_immunohistochemistry_receptor_status",
                 "drug_name",
                 "margin_status",
                 "race.demographic",
                 "gender.demographic",
                 "state.samples",
                 "sample_type_id.samples",
                 "vital_status.diagnoses",
                 "tumor_stage.diagnoses",
                 "metastatic_site_at_diagnosis"
                 )]

id=as.character(unique(laml@data$Tumor_Sample_Barcode))
Tumor_Sample_Barcode=data.frame(Tumor_Sample_Barcode=id,
                                submitter_id.samples=substring(id,1,16)
                                )
filter_phe=merge(Tumor_Sample_Barcode,filter_phe,by='submitter_id.samples')
filter_phe =filter_phe[,c(2 ,3:ncol(filter_phe))] 
colnames(filter_phe)[2:6]=c('anatomic','ER','PR','HER2','drug')
filter_phe$tnbc=ifelse(apply(filter_phe[3:5],1, function(x) sum(x=='Negative'))==3,'TNBC','nonTNBC')
head(filter_phe)
write.table(filter_phe,'TCGA_BRCA_filter_phe.tsv',sep = '\t',row.names = F,quote = F)
save(filter_phe,laml,file = 'BRCA_Maf_input.Rdata')
 





