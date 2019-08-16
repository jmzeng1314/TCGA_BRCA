## 癌症免疫图谱

[TOC]

2018年4月Immunity杂志上发表了文章[The Immune Landscape of Cancer](https://www.sciencedirect.com/science/article/pii/S1074761318301213) ，由34个单位共同合作完成。文章对TCGA中33种癌症，超过10,000个肿瘤样本进行了免疫原性分析，将所有肿瘤分成6种免疫亚型，进一步分析：

* 不同亚型间巨噬细胞或淋巴细胞特征、Th1:Th2细胞比例、肿瘤异质性程度、非整倍性、新抗原负荷程度、细胞增殖、免疫调节基因的表达、预后等指标的差异。
* 与免疫相关的驱动突变
* 参与肿瘤免疫的细胞内和细胞间调控网络（调控网络包含转录、microRNA、拷贝数、表观遗传信息）

### 分析流程

![](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-fx1.jpg)

- step1. 下载TCGA平台33种癌症类型的6种分子平台数据（mRNA、microRNA、外显子测序数据；甲基化、拷贝数、蛋白质芯片数据）

- step2. 为全部非血液类癌症进行**免疫表达特征**（160种）评分，将TCGA全部肿瘤聚类成6种免疫亚型

- step3. 根据DNA甲基化、mRNA和H&E图像分析确定白细胞类型及比例（用于评估免疫浸润细胞的组成）

- step4. 生存分析评估免疫亚型与预后的关系

- step5. 分别以TCGA肿瘤类型、TCGA定义的分子亚型、6种免疫亚型分组，比较新抗原、病毒RNA表达、T细胞受体(TCR)和B细胞受体(BCR)、免疫调节因子(IM)表达与调控等特征（评估影响免疫原性和免疫浸润的因素间的关系）

> 免疫原性（immunogenicity)：能引起机体产生免疫应答的物质
>
> 肿瘤浸润免疫细胞包括：T cells, B cells, natural killer cells, macrophages, neutrophils, dendritic cells, mast cells, eosinophils, basophils等

- step6. 鉴定与免疫浸润相关的体细胞变异（包括pathway, 拷贝数变异，driver mutation）

- step7. 评估性别与祖源是否对特定肿瘤类型的免疫反应有影响

- step8. 鉴定控制肿瘤免疫应答的细胞内调控网络、参与形成肿瘤免疫微环境(TME)的细胞外通讯网络

### 癌症的免疫亚型

计算160种免疫特征间的相关系数，聚类分析，得到5个核心模块免疫表达特征（巨噬细胞 , 免疫浸润淋巴细胞, TGF-β response , IFN-γ response, wound healing ) → 据此将TCGA中的非血液肿瘤聚成6种免疫亚型(C1-C6) → 6种免疫亚型中这5种免疫特征的分布情况。（图1A）

6种免疫亚型中肿瘤样本表现出的其它关键免疫特征（图1B,C）

每种肿瘤类型中免疫亚型的组成情况，图1D。

![Figure1](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr1.jpg)

Figure 1. Immune Subtypes in Cancer 

### 肿瘤免疫浸润的组成

不同免疫亚型的主要免疫浸润细胞分布，图2A

白细胞比例(leukocyte fraction,LF)在不同免疫亚型间（图1C）和不同肿瘤类型间（图2B）存在差异。

肿瘤基质比例与白细胞比例的相关性在肿瘤类型中不同，图2C显示2种代表性的肿瘤。（该相关性在免疫亚型的差别见补充图2B）

> [基质细胞](https://en.wikipedia.org/wiki/Stromal_cell)：器官中的结缔组织细胞，为器官中实质细胞提供支持和营养。
>
> 已知基质细胞和肿瘤细胞间的相互作用在肿瘤生长和进展中起主要作用。基质细胞可以提供肿瘤细胞可以生长的细胞外基质，调节细胞生长因子，还可能通过一氧化氮的产生限制T细胞增殖，从而阻碍免疫能力。

通过TCGA的H&E染色图片分析估计出具有浸润淋巴细胞的肿瘤区域的空间比例，按免疫亚型查看，发现C2亚型的淋巴细胞浸润区比例最高，图2D。

![](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr2.jpg)

Figure 2. Composition of the Tumor Immune Infiltrate 

### 肿瘤免疫反应与预后

6种免疫亚型OS生存分析（图3A）：C3预后最佳；C1,C2虽然免疫成分多，但预后并不好；C4,C6免疫特征混杂，预后最差。（PFI结果见补充图3A）

计算5种免疫表达特征分数与OS的一致性指数，红色代表高风险。按免疫亚型和肿瘤类型分组（图3B），淋巴细胞特征与C1,C2的预后改善相关；5种特征中任何一个增加都会导致C3预后更差。
计算辅助T细胞（Th1,Th2,Th17）与OS的一致性指数，按免疫亚型分组（图3C），Th17细胞的增加与多数亚型的预后改善相关。

> 辅助T细胞（T helper cells, Th）,是一种T细胞（白细胞的一种），它的表面有抗原受体，可以识别抗原呈递细胞的MHC-II 类分子呈献的抗原片段。 辅助T细胞主要可区分为Th1 Th2 Th17 及 Thαβ等四种。

使用弹性网络CoxPH建立免疫特征分数高低的生存模型，分数高低与生存显著相关，图3D是验证集的KM-polt，图3E显示模型的预测准确度。图3F说明添加肿瘤stage和tpye信息后模型预测更准确。

> 疑问：图F的y轴Delta Log-Likelihood是什么意思？为什么tissue+stage组的值在stage组和tissue组的中间呢？

补充图3E,F讨论了其它免疫特征与OS的关系，某些肿瘤类型中免疫亚型与OS的关系

![](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr3.jpg)

Figure 3. Immune Response and Prognostics 

### 免疫反应与体细胞变异

6中免疫亚型的免疫浸润（白细胞比例LF）与DNA损伤（包括CNV负荷、非整倍性、杂合性丢失LOH、同源重组缺陷HRD、肿瘤内异质性ITH）的相关系数热图，图4A。发现LF与C6和C2的相关性最强，从整体来看与非整倍性、LOH、HRD及突变负荷成正相关。

图B显示平均LF值在1-22号染色体上的观察值与预估值间的差异，左右图分别为拷贝数扩增和缺失，黑色部分代表差异显著。揭示了LF与拷贝数变异的相关性，如：chr1p扩增，LF高于预估值，缺失则LF低于预估值，据此进一步将该位置上的重要基因的CNV与免疫浸润相联系。

> 基因的平均LF预估值：每种肿瘤的平均IF的平均值由每种疾病类型中存在的“扩增/缺失”样品的数量加权。

补充图4探讨了肿瘤内异质性ITH与预后、LF、免疫亚型间的关系。

将299个癌症驱动基因的突变与免疫亚型相关联，发现33个显著相关的基因，图C。

图D为LF-驱动基因变异火山图，x轴为基因变异与LF的多因素相关性（考虑了肿瘤类型和突变数量），y轴为显著性，橙色即为与LF显著相关的驱动基因变异。其中23个驱动基因变异与LF增加相关，12个与LF降低相关。

补充图4C分析了体细胞变异（突变+CNVs）对8种癌信号通路的影响。

![](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr4.jpg)

Figure 4. Immune Response and Genome State 

### 人口统计学因素及遗传变异与免疫反应

考察了不同肿瘤类型中免疫细胞与性别和祖源的相关性，图4E。发现在某些肿瘤中，女性的PD-L1表达高于男性，非洲血统患者PD-L1表达较低。

### 免疫原性调查

文章中的免疫原性主要讨论的是新抗原负荷，即SNV和Indel突变产生的预测能与pMHC结合并诱导免疫反应的肽。

> [MHC](https://en.wikipedia.org/wiki/Major_histocompatibility_complex) proteins(pMHC)：指一些列细胞表面蛋白，主要作用是结合抗原并将其呈递至细胞表面，使其被T-cell识别，即帮助免疫系统识别外来物质。

突变与pMHC：图5A为与突变数量相关的pMHCs数量分布，标注了源于超过40个突变的pMHCs。99.8的pMHC由单突变产生，0.2%是由不同突变产生的相同pMHC。

常见pMHC：图5B为y轴表达共有pMHC的肿瘤数，标注了产生最常见pMHC的已知癌基因(PIK3CA, KRAS, BRAF等)

补充图5A,B讨论了pMHC数与预后的关系：在多数肿瘤中pMHC数与生存无关，但在某些免疫亚型中，pMHC数与PFI预后相关。

补充图5C讨论了病毒与肿瘤中免疫特征的关系，病毒含量会影响免疫细胞的组成。

![](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr5.jpg)

Figure 5. The Tumor-Immune Interface 

### 癌症中的适应性免疫受体库

> 适应性免疫[Adaptive immune](https://en.wikipedia.org/wiki/Adaptive_immune_system)：又称获得性免疫，对特定病原体具有高度特异性。获得性免疫系统的细胞是T和B淋巴细胞。B细胞主要参与体液免疫，而T细胞则参与细胞免疫。B细胞和T细胞都携带能够识别特定靶标的受体分子。 
>
> T细胞负责识别“非自身”靶标，如病原体，但需要一种被称为主要组织相容性复合体（MHC）呈递抗原之后才能实现。B细胞抗原特异性受体则是位于B细胞表面的抗体，负责识别病原体，不需要经过抗原呈递。 
>
> 免疫球蛋白[immunoglobulin](https://en.wikipedia.org/wiki/Antibody)，又称抗体，由适应性免疫系统的B细胞分泌。

根据RNA-seq数据评估TCR(T cell receptor)α和β、免疫球蛋白重链和轻链库。

发现TCR多样性值因免疫亚型而不同(图5C)，C6和C2的多样性最高。（补充图5D展示了不同肿瘤类型的TCR多样性）

评估CDR3的α和β链，以确定具有相同TCR的患者频率。图5C为CDR3α-CDR3β pairs，颜色代表p值，圆的大小代表sample数目，至少在2种肿瘤中出现的α-β pairs有2812个。

检验SNV特异性pMHC-CDR3 pairs，图5D展示了206个pMHC-CDR3α pairs和196个pMHC-CDR3β pairs。（说明这些患者中抗原受体多，抗原少，多数T细胞抗原反应是由公共抗原介导）

> 抗原受体通常由位于重链和轻链上的两个可变结构域组成，每个可变结构域上有三个非连续排列的CDR（CDR1，CDR2和CDR3）
>
> [Complementarity-determining regions, CDR](https://en.wikipedia.org/wiki/Complementarity-determining_region): 互补决定区，是免疫球蛋白（抗体）和T细胞受体中可变链的一部分，分别由B细胞和T细胞产生，与其特异性抗原结合。作为分子中变化最大的部分，CDR对淋巴细胞产生的抗原特异性的多样性至关重要。 

在某些肿瘤类型中，TCR多样性与PFI改善相关，见补充图5F。

免疫球蛋白重链和轻链的多样性模式与TCR类似。

### 免疫调节剂的调控

免疫调节剂(Immunomodulators, IM)基因的情况可以反映癌细胞对肿瘤微环境(TME)的调节。

IM基因的表达在免疫亚型中不同（图6A），图6B为亚型间差异最大的基因EDNRB和CXCL10的表达值。

许多IM的甲基化程度与表达值负相关，即表观遗传沉默（图6A），如图6C的C3亚型中的CD40基因。

294个miRNA可能是IM基因的调节因子，补充图6C。

拷贝数改变影响多个IM，并在免疫亚型中存在差异。C1,C2亚型中IM基因频繁扩增和缺失，C3和C5中IM基因变异较少（图6A）。KIR2DL3等基因在C5中明显缺失，CD40等基因常在C1中发生扩增（图6D）。

> 图6A热图中的扩增/缺失频率：obs是指特定亚型中IM基因扩增/缺失频率，exp是指全部样本中扩增/缺失频率的差异，热图反映的是obs与exp的差异

![](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr6.jpg)

Figure 6. Regulation of Immunomodulators 

### 免疫反应的网络调控

免疫应答由肿瘤细胞、免疫细胞、基质细胞中的细胞内分子网络状态和细胞外通讯网络（包括细胞间直接互作、通过可溶性蛋白或细胞因子通讯）构成。

根据[FANTOM5](http://fantom.gsc.riken.jp/5/suppl/Ramilowski_et_al_2015/)中记录的信息获得配体 - 受体，细胞 - 受体和细胞 - 配体对的网络，图A,B,C分别是提取的IFM-γ, TGF-β, T细胞的细胞外通讯子网络网络。展示了免疫亚型间的一些关联，例如C2,C3中，CD4 T细胞、CD8 T细胞和NK细胞都与IFM-γ及CCL5的表达相关；C2,C3,C6中，多种细胞都与TGF-β的表达相关……

构建泛癌症转录调控网络，关联基因组事件 - 转录调控因子 - 下游靶基因 - 免疫浸润-预后 。使用Master Regulator(MR)和SYGNAL两种互补的方法，产生了2个转录调控网络：

* 图D是Master Regulator(MR)-Pan-Immune网络，26个MRs(橙色填充)代表基因表达与LF相关的hub，连接15个上游driver events(橙色环)。
* SYGNAL-PamImmune网络，包含了171个富集且与LF相关的的IMs，主要展示下游调控网络（突变 - 转录因子 - miRNA间的联系）。（图7E展示部分SYGNAL-PamImmune网络）

> MRs网络图的构建过程：根据下游靶基因的表达变化，推断影响蛋白质活性的转录因子，为MRs → 根据蛋白互作寻找与MRs相关的发生变异的蛋白cMRs(与LF的相关系数高的优先) → 寻找与cMRs相关的mutation和CNVs事件。 
>
> SYGNAL网络构建过程：用高表达的基因作为input，得到具有转录因子和miRNA共调节的基因(bicluster) → 过滤筛选得到171个与免疫浸润和免疫调节均显著相关的biclusters→ 添加biclusters的转录因子和miRNA调控信息 → 整合这些调控关系，连接成网络
>
> （网络图构建步骤比较复杂，大致应该是这个意思）

综合2个调控网络，MR及SYGNAL中有7个共享的TFs，调节27个IM基因，推断AKAP9, HRAS, KRAS和PREX2基因突变影响IMs。

图7E展示免疫亚型-IM基因表达调控-肿瘤类型间的关系，可以鉴定不同免疫亚型内的IM基因调控因子。例如，C1和C2中的IM基因共同受BCL5B, ETV7, IRF1, IRF2, IRF4, PRDM1, 和SPIB调控；C3中IM受KLF15 and miR-141-3p调控……

某些肿瘤中，bicluster的表达增加与预后相关，在KIRC, LGG, LUSC和READ中预后差，在SKCM中预后好。

![](https://ars.els-cdn.com/content/image/1-s2.0-S1074761318301213-gr7.jpg)

Figure 7. Predicted Networks Modulating the Immune Response to Tumors 

### 结果总结

文章利用TCGA数据构建了癌症的免疫基因组学图谱：

* 根据5种免疫表达特征将TCGA的肿瘤分为6种免疫亚型
* 通过多种方法评估肿瘤样本中的免疫组成，包括：根据基因表达和甲基化数据估计免疫细胞组成、预测来自突变和HLA-typing的新抗原-MHC pairs、根据RNA测序数据评估BCR和TCR库。
* 比较癌症和免疫亚型的免疫组成
* 鉴定与肿瘤微环境相关的体细胞变异
* 构建影响肿瘤微环境的分子调控网络及细胞间通讯网络
* 联系免疫特征与生存，OS和PFI在肿瘤间或一些肿瘤内的免疫亚型中存在差异

数据和结果的补充表格下载地址：https://gdc.cancer.gov/about-data/publications/panimmune





> 文章很清晰，做了pan-cancer的免疫亚型分类，然后把能用到的所有基因组指标和免疫指标都进行了分组比较。
>
> 文章的数据量很大，用到的方法很多，方法介绍也写得很详细，需要时可以再详细了解；补充表格的数据有很多信息可以直接拿来用
>
> 免疫相关的知识我几乎不懂，对文章的理解也不一定准确，仅供参考


