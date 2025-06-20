#humann3 unstratified ttest and wilcoxon results(gvhd vs no gvhd)
rm(list = ls())
library(stringr)
unloadNamespace("scales")
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/unstratified.txt",
                sep = "\t", header = T, row.names = 1, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
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
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(NormT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(NormT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(NormT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/Humann3GvNUnstratifiedResults.txt", sep = "\t", row.names = T)

#humann3 unstratified ttest and wilcoxon results(sensitive vs stratified)
rm(list = ls())
unloadNamespace("scales")
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/unstratified.txt",
                sep = "\t", header = T, row.names = 1, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
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
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(NormT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(NormT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(NormT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/Humann3SvRUnstratifiedResults.txt", sep = "\t", row.names = T)