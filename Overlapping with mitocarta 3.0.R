## overlapping between DEG between WT and KO after TAC and mitocarta 3.0 in mice
rm(list = ls())
options(stringsAsFactors = F)


#tac后变化的基因被cd163调控
TAC <- read.csv("Data/WTvsKT_Sig.csv")
mito  <- read.csv("Mitocart/mito 3.0.csv")

# 去掉没有配对上symbol的行
TAC <- TAC[!is.na(TAC$SYMBOL),]

# 去掉重复行
TAC <- TAC[!duplicated(TAC$SYMBOL),]

###交集加并集
OV=TAC[TAC$SYMBOL%in%mito$Symbol,]

TACu=TAC[!TAC$SYMBOL%in%mito$Symbol,]

mitou=mito[!mito$Symbol%in%TAC$SYMBOL,]

dim(OV)

write.csv(OV,"Data/OV mito.csv", row.names = F)
write.table(OV$SYMBOL,
            file="Data/OV mito.txt",
            row.names = F,
            col.names = F,
            quote = F)
