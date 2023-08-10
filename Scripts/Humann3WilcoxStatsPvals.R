#humann3 unstratified ttest and wilcoxon results(gvhd vs no gvhd)
rm(list = ls())
library(stringr)
unloadNamespace("scales")
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/unstratified.txt",
                sep = "\t", header = T, row.names = 1, check.names = F)
myT<-myT[ , -which(colSums(myT) == 0)]
metaData<-metaData[colnames(myT), , drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData)]
rownames(myT)<-sapply(str_split(rownames(myT), ":", n = 2 ), `[`, 1)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
NormT<-myT
for (i in 1:ncol(NormT)) {
  NormT[,i]<-NormT[,i]/n[i]
}
NormT<-log10(NormT*(sumx/ncol(NormT))+1)
NormT<-data.frame(NormT, check.names = F)
NormT<-na.omit(NormT)
#Filter
lowAbundance<-which(rowMeans(NormT)<1)
NormT<-NormT[-lowAbundance, ]

NormT<-NormT[, rownames(metaData), drop = F]
#Stats and pvals
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(NormT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(NormT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(NormT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/Humann3GvNUnstratifiedResults.txt", sep = "\t", row.names = T)

#humann3 unstratified ttest and wilcoxon results(sensitive vs stratified)
rm(list = ls())
unloadNamespace("scales")
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/unstratified.txt",
                sep = "\t", header = T, row.names = 1, check.names = F)
myT<-myT[ , -which(colSums(myT) == 0)]
metaData<-metaData[colnames(myT), , drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData)]
rownames(myT)<-sapply(str_split(rownames(myT), ":", n = 2 ), `[`, 1)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
NormT<-myT
for (i in 1:ncol(NormT)) {
  NormT[,i]<-NormT[,i]/n[i]
}
NormT<-log10(NormT*(sumx/ncol(NormT))+1)
NormT<-data.frame(NormT, check.names = F)
NormT<-na.omit(NormT)
#Filter
lowAbundance<-which(rowMeans(NormT)<1)
NormT<-NormT[-lowAbundance, ]

NormT<-NormT[, rownames(metaData), drop = F]
#Stats and pvals
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(NormT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(NormT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(NormT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/Humann3SvRUnstratifiedResults.txt", sep = "\t", row.names = T)