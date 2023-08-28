#Bracken phylum ttest results(gvhd vs no gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
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
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(myT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/BrackenGvNPhylumResults.txt", sep = "\t", row.names = T)

#Bracken genus ttest results(gvhd vs no gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
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
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(myT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/BrackenGvNGenusResults.txt", sep = "\t", row.names = T)

#Bracken species ttest results(gvhd vs no gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
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
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(myT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/BrackenGvNSpeciesResults.txt", sep = "\t", row.names = T)

#Bracken phylum ttest results(sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
metaData<-metaData[colnames(myT),, drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData), drop = F]
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
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(myT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/BrackenSvRPhylumResults.txt", sep = "\t", row.names = T)

#Bracken genus ttest results(sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
metaData<-metaData[colnames(myT),, drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData), drop = F]
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
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(myT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/BrackenSvRGenusResults.txt", sep = "\t", row.names = T)

#Bracken species ttest results(sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
metaData<-metaData[colnames(myT),, drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData), drop = F]
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
wilcox_stats<-vector()
wilcox_p<-vector()
wilcox_stats<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$stat})
wilcox_p<-apply(myT, 1, function(x){wilcox.test(unlist(x)~metaData$dx)$p.value})
wilcox_stat_table<-cbind(wilcox_stats,wilcox_p)
row.names(wilcox_stat_table)<-rownames(myT)
wilcox_stat_table<-data.frame(wilcox_stat_table)
wilcox_stat_table$padj<-p.adjust(wilcox_stat_table$wilcox_p, method = "fdr")

write.table(wilcox_stat_table, "wilcox/BrackenSvRSpeciesResults.txt", sep = "\t", row.names = T)

