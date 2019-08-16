## 分析影响生存的基因突变

### logrank test方法

载入突变和生存数据：


```{R}
load("./Rdata/brca_phenotype.Rdata")
load("./Rdata/mutation.Rdata")
```
生存数据源自TCGA PanCanAtlas推荐的整理后结果https://gdc.cancer.gov/node/905/  

突变信息源自UCSC Xena整理后的[基因水平非同义体细胞突变结果](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2Fmutation_curated_wustl_gene&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)。1代表基因发生非沉默体细胞突变，包括nonsense, missense, frame-shif indels, splice site mutations, stop codon readthroughs, change of start codon, inframe indels；0代表不存在上述变异。

```{R}
library(survival)
#方案一：按OS
phe <- brca_phe[brca_phe$OS.time != "#N/A",]
phe$OS.time <- as.numeric(phe$OS.time)
phe$OS <- as.numeric(phe$OS)
mut_freq <- apply(mut_info, 1, function(x){sum(x)/ncol(mut_info)}) #每个基因的突变频率
mut_genes <- mut_freq[mut_freq != 0] #15709个基因发生突变
summary(mut_genes) #突变频率的上四分位数为0.005117707

#只保留突变频率≥0.05%的基因
mut <- mut_info[mut_freq >= 0.005 ,phe$bcr_patient_barcode] #4176个基因

mySurv <- with(phe,Surv(OS.time, OS))
log_rank_p <- apply(mut, 1, function(gene){
  phe$group=ifelse(gene == 1,'MUT','WT')
  data.survdiff=survdiff(mySurv~group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

mut_logrank_OS <- sort(log_rank_p)
table(mut_logrank_OS < 0.01) #357个基因

#方案二：按PFI
phe <- brca_phe[brca_phe$PFI.time != "#N/A",]
phe$PFI.time <- as.numeric(phe$PFI.time)
phe$PFI <- as.numeric(phe$PFI)
mut <- mut_info[which(mut_freq >= 0.005) ,phe$bcr_patient_barcode] 

mySurv <- with(phe,Surv(PFI.time, PFI))
log_rank_p <- apply(mut, 1, function(gene){
  phe$group=ifelse(gene == 1,'MUT','WT')
  data.survdiff=survdiff(mySurv~group,data=phe)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})

mut_logrank_PFI <- sort(log_rank_p)
table(mut_logrank_PFI < 0.01) #312个基因

save(mut_logrank_OS, mut_logrank_PFI, mut_genes, file = "./results/survival/log_rank_mut.Rdata")
```

### 用机器学习算法通过基因突变预测生存

```{R}
###lasso回归
library(lars) 
library(glmnet) 
table(colnames(phe) == phe$bcr_patient_barcode) #确保样本顺序一致
x <- t(mut)
y <- phe$OS

model_lasso <- glmnet(x, y, family="binomial", nlambda=50, alpha=1)
print(model_lasso) #在0和1之间，越接近1说明模型的表现越好

cv_fit <- cv.glmnet(x=x, y=y, alpha = 1, nlambda = 1000) #会进行交叉验证，找到最优的lambda值
plot.cv.glmnet(cv_fit) # 两条虚线分别指示了两个特殊的λ值: 可选区间的任意一个
c(cv_fit$lambda.min,cv_fit$lambda.1se) #最小值，1倍离差值

model_lasso <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
rownames(model_lasso$beta)[as.numeric(model_lasso$beta)!=0] #没有可用于预测生存的突变基因

###随机森林
library(randomForest)
table(colnames(phe) == phe$bcr_patient_barcode) #确保样本顺序一致
x <- t(mut)
y <- phe$OS

rf_output <- randomForest(x=x, y=y,importance = TRUE, ntree = 500, proximity=TRUE )
save(rf_output,file = "./results/survival/randomforests_output.Rdata")
rf_importances=importance(rf_output, scale=FALSE) #IncNodePurity代表节点不纯度减少值，值越大越重要
head(rf_importances)
choose_genes <- rownames(tail(rf_importances[order(rf_importances[,2]),],50))
```

### 生存相关基因的交集及突变频率

```{R}
### 高频突变基因与生存相关基因的交集
library(UpSetR)
load("./results/survival/log_rank_mut.Rdata")
OS_mut_genes <- names(mut_logrank_OS[mut_logrank_OS < 0.05]) #574 genes
PFI_mut_genes <- names(mut_logrank_PFI[mut_logrank_PFI < 0.05]) #576 genes
high_mut_100 <- names(head(sort(mut_genes, decreasing = T),100))
rf_genes_100 <- rownames(tail(rf_importances[order(rf_importances[,2]),],100))
genes <- unique(c(OS_mut_genes, PFI_mut_genes, high_mut_100, rf_genes_100)) 

sur_mut_df <- data.frame(genes = genes,
                         OS_genes=ifelse(genes %in% OS_mut_genes, 1, 0),
                         PFI_genes = ifelse(genes %in% PFI_mut_genes, 1, 0),
                         high_mut_100 = ifelse(genes %in% high_mut_100, 1, 0),
                         rf_genes_100 = ifelse(genes %in% rf_genes_100, 1, 0)
)
rownames(sur_mut_df) <- sur_mut_df$genes
sur_mut_df$mut_rate <- mut_freq[rownames(sur_mut_df)]
upset(sur_mut_df[,-1])

intersect(high_mut_100, rf_genes_100) #7个基因，包含PIK3CA和TP53
```

![sur_genes_upset.png](https://i.loli.net/2019/01/03/5c2e0f29374f3.png)

```{R}
###生存相关基因的突变频率boxplot
library(reshape2)
df <- melt(sur_mut_df, measure.vars=c(2:5))
df <- df[df$value==1,]
library(ggpubr)
ggboxplot(df, x = "variable", y = "mut_rate",
          color = "variable", palette = "jco",
          add = "jitter") +
labs(x="",y="mutation rate")  
```

![sur_genes_mut_rate_boxplot.png](https://i.loli.net/2019/01/03/5c2e0eba42f26.png)

logrank生存分析得到的重要基因突变频率都很低。而随机森林得到的重要基因突变率分布较广，包括了一些高突变基因。
