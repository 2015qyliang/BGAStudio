# BGAStudio
**Bacterial genomes analysis on Windows** (基于 windows 平台的细菌基因组分析流程)
## **Just One Step** & **Batch Processing**

在科研工作中, 对于非生物信息信息分析专业的人来讲, 经常接触的是 windows 操作系统, 然而多数的细菌基因组分析工具在该系统下并不能够"完美"运行; 因此整合了可在 windows 系统平台下运行的细菌基因组预测软件 - [prodigal v2.6.3](https://github.com/hyattpd/Prodigal/releases) 和蛋白比对软件 - [diamond v0.9.28](https://github.com/bbuchfink/diamond/releases), 同时梳理了 [uniprot](https://www.ebi.ac.uk/uniprot/download-center) & [COG](https://www.ncbi.nlm.nih.gov/COG/) & [orthoMcl_OG5](https://orthomcl.org/common/downloads/release-5/) & [dbCAN](http://bcb.unl.edu/dbCAN2/) 数据库, 从而使得能够在该平台下对细菌基因组进行简单的功能注释. 

系统发育分析: **(1)** 核酸序列水平 ([UBCG](https://www.ncbi.nlm.nih.gov/pubmed/29492869) & [GTDB](https://doi.org/10.1038/nbt.4229)) , **(2)** 氨基酸序列水平 (*rpoB* & Ribosomal proteins) *氨基酸序列水平的系统发育分析的[参考文献](https://doi.org/10.1038/s41564-019-0588-1)* ; 后续系统发育树的构建过程可参考 https://github.com/2015qyliang/PhyTA 中的相关内容

在该流程中, 使用者仅仅提供已完成基因组测序细菌基因组核酸序列即可, 将其存放在 01_GenomeFna 文件夹中, 序列文件的命名需要使用下划线将物种名进行连接, 同时文件的后缀为 .fasta, 文件名示例: Bdellovibrio_bacteriovorus_109J.fasta. 

序列数据准备完成之后, 仅需双击 Step_1_Prodigal.bat 即可等待程序启动开始进行细菌基因预测, 待其分析任务完成之后, 双击 Step_2_Diamond.bat 即可等待程序启动开始进行基因功能注释.

**! 需要在 windows 系统中配置 R 语言的编译环境** 以及 **需要配置 python 语言的编译环境 (pip install biopython 安装需要的分析模块)**

---

### 文献推荐

1 [Cultivation and functional characterization of 79 planctomycetes uncovers their unique biology -- 2019](https://doi.org/10.1038/s41564-019-0588-1)

2 [A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life -- 2018](https://www.nature.com/articles/nbt.4229)

3 [Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life -- 2017](https://www.nature.com/articles/s41564-017-0012-7)

4 [UBCG: Up-to-date bacterial core gene set and pipeline for phylogenomic tree reconstruction -- 2018](https://www.ncbi.nlm.nih.gov/pubmed/29492869)

5 [dbCAN2: a meta server for automated carbohydrate-active enzyme annotation -- 2018](https://academic.oup.com/nar/article/46/W1/W95/4996582)

---

2019-11-17

├─01_GenomeFna

├─02_ProteinsFaa

├─03_GenesFna

├─04_DiamondTmpDir

├─05_ResultUniprot

├─06_ResultOrthoMclOG5

├─07_ResultCOG

├─08_ResultCAZy

├─09_ResultGTDB

├─PhyloAnalysisResult

├─Database

├─Lib

└─Scripts

---
