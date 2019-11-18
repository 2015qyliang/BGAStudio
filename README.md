# BGAStudio
Bacterial genomes analysis on Windows (**基于 windows 平台的细菌基因组分析流程**)

## **Just One Step**

在科研工作中, 对于非生物信息信息分析专业的人来讲, 经常接触的是 windows 操作系统, 然而多数的细菌基因组分析工具在该系统下并不能够"完美"运行; 因此整合了可在 windows 系统平台下运行的细菌基因组预测软件 - [prodigal v2.6.3](https://github.com/hyattpd/Prodigal/releases) 和蛋白比对软件 - [diamond v0.9.28](https://github.com/bbuchfink/diamond/releases), 同时梳理了 [uniprot](https://www.ebi.ac.uk/uniprot/download-center) & [COG](https://www.ncbi.nlm.nih.gov/COG/) & [orthoMcl_OG5](https://orthomcl.org/common/downloads/release-5/) 数据库, 从而使得能够在该平台下对细菌基因组进行简单的功能注释. 

在该流程中, 使用者仅仅提供已完成基因组测序细菌基因组核酸序列即可, 将其存放在 01_GenomeFna 文件夹中, 序列文件的命名需要使用下划线将物种名进行连接, 同时文件的后缀为 .fasta, 文件名示例: Bdellovibrio_bacteriovorus_109J.fasta. 

序列数据准备完成之后, 仅需双击 Step_1_Prodigal.bat 即可等待程序启动开始进行细菌基因预测, 待其分析任务完成之后, 双击 Step_2_Diamond.bat 即可等待程序启动开始进行基因功能注释.

---

2019-11-17

├─01_GenomeFna

├─02_ProteinsFaa

├─03_GenesFna

├─04_DiamondTmpDir

├─05_ResultUniprot

├─06_ResultOrthoMclOG5

├─07_ResultCOG

├─Database

├─Lib

└─Scripts

---
