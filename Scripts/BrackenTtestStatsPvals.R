#Bracken phylum ttest results(gvhd vs no gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(myT)<2)
myT<-myT[-lowAbundance, ]

metaData<-metaData[colnames(myT),, drop = F]
#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(myT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/BrackenGvNPhylumResults.txt", sep = "\t", row.names = T)

#Bracken genus ttest results(gvhd vs no gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(myT)<2)
myT<-myT[-lowAbundance, ]

metaData<-metaData[colnames(myT),, drop = F]
#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(myT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/BrackenGvNGenusResults.txt", sep = "\t", row.names = T)

#Bracken species ttest results(gvhd vs no gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(myT)<2)
myT<-myT[-lowAbundance, ]

metaData<-metaData[colnames(myT),, drop = F]
#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(myT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/BrackenGvNSpeciesResults.txt", sep = "\t", row.names = T)

#Bracken phylum ttest results(sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(myT)<2)
myT<-myT[-lowAbundance, ]

metaData<-metaData[colnames(myT),, drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData), drop = F]
#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(myT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/BrackenSvRPhylumResults.txt", sep = "\t", row.names = T)

#Bracken genus ttest results(sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(myT)<2)
myT<-myT[-lowAbundance, ]

metaData<-metaData[colnames(myT),, drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData), drop = F]
#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(myT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/BrackenSvRGenusResults.txt", sep = "\t", row.names = T)

#Bracken species ttest results(sensitive vs refractory)
rm(list = ls())
library("coin")
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)

#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(myT)<2)
myT<-myT[-lowAbundance, ]

metaData<-metaData[colnames(myT),, drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData), drop = F]
#Stats and pvals
t_stats<-vector()
t_test_p<-vector()
t_stats<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$stat})
t_test_p<-apply(myT, 1, function(x){t.test(unlist(x)~metaData$dx)$p.value})
T_stat_table<-cbind(t_stats,t_test_p)
row.names(T_stat_table)<-rownames(myT)
T_stat_table<-data.frame(T_stat_table)
T_stat_table$padj<-p.adjust(T_stat_table$t_test_p, method = "fdr")

write.table(T_stat_table, "ttest/BrackenSvRSpeciesResults.txt", sep = "\t", row.names = T)

