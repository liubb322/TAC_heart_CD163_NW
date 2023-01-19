#差异的富集，TAC组中ko导致的基因差异。
rm(list = ls())
options(stringsAsFactors = F)

library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(ggplot2)
library(tidyverse)

# 读取差异分析结果
load(file = "Data/Step03-WTvsKT.Rdata")
ls()

# 提取所有差异表达的基因名
DEG <- DEG_edgeR_symbol[DEG_edgeR_symbol$regulated!="normal",2]
head(DEG)

## ===GO数据库, 输出所有结果，后续可根据pvalue挑选结果
ego_CC <- enrichGO(gene=DEG, OrgDb= 'org.Mm.eg.db', keyType='SYMBOL', ont="CC", pvalueCutoff= 1,qvalueCutoff= 1)
ego_MF <- enrichGO(gene=DEG, OrgDb= 'org.Mm.eg.db', keyType='SYMBOL', ont="MF", pvalueCutoff= 1,qvalueCutoff= 1)
ego_BP <- enrichGO(gene=DEG, OrgDb= 'org.Mm.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff= 1,qvalueCutoff= 1)

p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("Biological process")
p_CC <- barplot(ego_CC,showCategory = 5) + ggtitle("Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("Molecular function")
plotc <- p_CC
plotc
ggsave('6.enrichGO.png', plotc, width = 10,height = 16)

ego_BP <- data.frame(ego_BP)
ego_CC <- data.frame(ego_CC)
ego_MF <- data.frame(ego_MF)
write.csv(ego_BP,'6.enrichGO_BP.csv')
write.csv(ego_CC,'6.enrichGO_CC.csv')
write.csv(ego_MF,'6.enrichGO_MF.csv')



## === KEGG
genelist <- bitr(gene=DEG, fromType="SYMBOL", toType="ENTREZID", OrgDb='org.Mm.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'mmu', pvalueCutoff = 1, qvalueCutoff = 1)#human是hsa
p1 <- barplot(ekegg, showCategory=10)
p2 <- dotplot(ekegg, showCategory=10)
plotc = p1/p2
plotc
ggsave('6.enrichKEGG.png', plot = plotc, width = 8, height = 10)

ekegg <- data.frame(ekegg)
write.csv(ekegg,'6.enrichKEGG.csv')



## === 其他数据库通路
geneset <- read.gmt("data/MsigDB/v7.4/h.all.v7.4.symbols.gmt")
table(geneset$term)
geneset$term <- gsub(pattern = "HALLMARK_","", geneset$term)
geneset$term <- str_to_title(geneset$term)

my_path <- enricher(gene=DEG, pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE=geneset)
p1 <- barplot(my_path, showCategory=15,color = "pvalue")
p1
ggsave("result/6.enrich_HALLMARK.png", plot = p1, width = 8, height = 7)
  
my_path <- data.frame(my_path)
write.csv(my_path,"result/6.enrich_HALLMARK.csv")
  

