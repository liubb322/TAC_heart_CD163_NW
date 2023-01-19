
rm(list = ls())
options(stringsAsFactors = F)
#tac后变化的基因被cd163调控

##比较WT小鼠TAC之后的基因表达的变化的情况 down tac < sham up tac > sham
#Tac 上调的基因
TAC_sham_up <- read.csv("Data/WTvsWS_Sig_up.csv")
#去重
TAC_sham_up = TAC_sham_up[!is.na(TAC_sham_up$SYMBOL),]
TAC_sham_up = TAC_sham_up[!duplicated(TAC_sham_up$SYMBOL),]
#Tac 下调的基因
TAC_sham_down <- read.csv("Data/WTvsWS_Sig_down.csv")
#去重
TAC_sham_down = TAC_sham_down[!is.na(TAC_sham_down$SYMBOL),]
TAC_sham_down = TAC_sham_down[!duplicated(TAC_sham_down$SYMBOL),]

##TAC后 KO小鼠相对于WT变化的基因 up KO>WT in TAC; down KO< WT in TAC
#KO 上调的的基因
KO_WT_up<- read.csv("Data/WTvsKT_Sig_up.csv")
#去重
KO_WT_up = KO_WT_up[!is.na(KO_WT_up$SYMBOL),]
KO_WT_up = KO_WT_up[!duplicated(KO_WT_up$SYMBOL),]
#KO 下调的的基因
KO_WT_down  <- read.csv("Data/WTvsKT_Sig_down.csv")
#去重
KO_WT_down = KO_WT_down[!is.na(KO_WT_down$SYMBOL),]
KO_WT_down = KO_WT_down[!duplicated(KO_WT_down$SYMBOL),]


###交集加并集
d_d = TAC_sham_down[TAC_sham_down$SYMBOL%in%KO_WT_down$SYMBOL,]
u_u = TAC_sham_up  [TAC_sham_up$SYMBOL%in%KO_WT_up$SYMBOL,]  

DEGm=c(d_d$SYMBOL,u_u$SYMBOL)


##====== 检查是否上下调设置错了
# 挑选一个差异表达基因
head(d_d)

# 读取基因表达矩阵信息并查看分组信息和表达矩阵数据
lname <- load(file = "Data/01.Rdata")
lname
exp <- c(t(express_cpm[match("ENSMUSG00000000441",rownames(express_cpm)),]))
test <- data.frame(value=exp, group=group$sample_title)
library(ggplot2)

ggplot(data=test,aes(x=group,y=value,fill=group)) + geom_boxplot()

##同源转化
library(homologene)
#这里以小鼠的三个基因为例
#更多基因方法是一样的
#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人
y=homologene(DEGm, inTax = 10090, outTax = 9606)
y = y[!is.na(y$`9606`),]
y = y[!duplicated(y$`9606`),]

#y2=DEG[DEG$SYMBOL%in%y$`10090`,]
#y3=y[y$`10090`%in%y2$SYMBOL,]
#y4 <- merge(DEG,y,by.x = "SYMBOL",by.y = "10090")#合并同源转化的基因和表达矩阵
DEG <- y$`9606`


## ===GO数据库, 输出所有结果，后续可根据pvalue挑选结果

library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)

ego_CC <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="CC", pvalueCutoff= 1,qvalueCutoff= 1)
ego_MF <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="MF", pvalueCutoff= 1,qvalueCutoff= 1)
ego_BP <- enrichGO(gene=DEG, OrgDb= 'org.Hs.eg.db', keyType='SYMBOL', ont="BP", pvalueCutoff= 1,qvalueCutoff= 1)

p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("Biological process")
p_CC <- barplot(ego_CC,showCategory = 5) + ggtitle("Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("Molecular function")
plotc <- p_BP
plotc
ggsave('Result/6.enrichGO.png', plotc, width = 10,height = 16)

ego_BP <- data.frame(ego_BP)
ego_CC <- data.frame(ego_CC)
ego_MF <- data.frame(ego_MF)
write.csv(ego_BP,'Result/6.enrichGO_BP.csv')
write.csv(ego_CC,'Result/6.enrichGO_CC.csv')
write.csv(ego_MF,'Result/6.enrichGO_MF.csv')


#取通路的基因名字
a=read.csv(file = "Result/6.enrichGO_BP.csv")
b=a[8,]$geneID
library(stringr)
c=str_split(b,"/")#字符串拆分
d=c[[1]]#基因名字的向量

##同源转化
#这里以小鼠的三个基因为例
#更多基因方法是一样的
#使用homologene函进行转换
#@genelist是要转换的基因列表
#@inTax是输入的基因列表所属的物种号，10090是小鼠
#@outTax是要转换成的物种号，9606是人
dy =y[y$`9606`%in%d,]

cytom <- dy$`10090`

###画热图

# 加载原始表达矩阵
lname <- load(file = "Data/01.Rdata")
lname

express_cpm1 <- rownames_to_column(as.data.frame(express_cpm) ,var = "ID")

# 读取差异分析结果
lname <- load(file = "data/Step03-WTvsKT.Rdata")
lname
#取要画的基因的差异分析的结果
cc=DEG_edgeR_symbol[DEG_edgeR_symbol$SYMBOL%in%cytom,]

#构建count的文件
data <- merge(cc,express_cpm1,by.x = "ENSEMBL",by.y = "ID")
data <- na.omit(data)
data <- data[!duplicated(data$SYMBOL),]



# 绘制热图
dat <- select(data,starts_with(c("K","W")))
rownames(dat) <- data$SYMBOL
dat[1:4,1:12]
anno <- data.frame(group=group$sample_title,row.names = group$run_accession)
xx=group$run_accession
class(xx)

# 按照指定顺序绘制热图
WT <- express_cpm[,match(rownames(anno)[which(anno$group=="WT")],
                         colnames(express_cpm))]

KT <- express_cpm[,match(rownames(anno)[which(anno$group=="KT")],
                         colnames(express_cpm))]
WS <- express_cpm[,match(rownames(anno)[which(anno$group=="WS")],
                         colnames(express_cpm))]
KS <- express_cpm[,match(rownames(anno)[which(anno$group=="KS")],
                         colnames(express_cpm))]

data_new <- cbind(WS,WT, KT,KS)
dat1 <- data_new[match(cc$ENSEMBL,rownames(data_new)),]
dat1
rownames(dat1)=cc$SYMBOL
#求均值
WTm = apply(WT,1,mean)
KTm = apply(KT,1,mean)
WSm = apply(WS,1,mean)
KSm = apply(KS,1,mean)
data_newm <- cbind(WSm,WTm, KTm,KSm)
dat1m <- data_newm[match(cc$ENSEMBL,rownames(data_new)),]
dat1m
rownames(dat1m)=cc$SYMBOL
dat1m
#画图
library(pheatmap)
pheatmap(dat1m, scale = "row",show_colnames =T,show_rownames = T, 
         cluster_cols = F,
         main = "edgeR's DEG")

