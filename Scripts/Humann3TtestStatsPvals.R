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
NormT<-log10((myT/n)*(sumx/nrow(myT))+1)
NormT<-na.omit(NormT)
#Filter
lowAbundance<-which(rowMeans(NormT)<1)
NormT<-NormT[-lowAbundance, ]

#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
for (i in 1:dim(NormT)[1]){
  t_stats[i]<-t.test(unlist(NormT[i,])~metaData$dx)$statistic
  t_test_p[i]<-t.test(unlist(NormT[i,])~metaData$dx)$p.value
}
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
myT<-myT[ , -which(colSums(myT) == 0)]
metaData<-metaData[colnames(myT), , drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData)]
rownames(myT)<-sapply(str_split(rownames(myT), ":", n = 2 ), `[`, 1)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
NormT<-log10((myT/n)*(sumx/nrow(myT))+1)
NormT<-na.omit(NormT)
#Filter
lowAbundance<-which(rowMeans(NormT)<1)
NormT<-NormT[-lowAbundance, ]

#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
for (i in 1:dim(NormT)[1]){
  t_stats[i]<-t.test(unlist(NormT[i,])~metaData$dx)$statistic
  t_test_p[i]<-t.test(unlist(NormT[i,])~metaData$dx)$p.value
}
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(NormT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/Humann3SvRUnstratifiedResults.txt", sep = "\t", row.names = T)