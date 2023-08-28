#aldex 2 analysis
rm(list = ls())
library(ALDEx2)
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.csv("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[, intersect(rownames(metaData), colnames(myT)), drop = F]
metaData<-metaData[intersect(rownames(metaData), colnames(myT)), , drop = F]

x <- aldex.clr(myT, metaData$dx, mc.samples=128)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)

write.table(x.tt[, c(1,2), drop = F], "ttest/BrackenGvNPhylumALDEx2Results.txt", sep = "\t", row.names = T)
write.table(x.tt[, c(3,4), drop = F], "wilcox/BrackenGvNPhylumALDEx2Results.txt", sep = "\t", row.names = T)

metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.csv("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[, intersect(rownames(metaData), colnames(myT)), drop = F]
metaData<-metaData[intersect(rownames(metaData), colnames(myT)), , drop = F]

x <- aldex.clr(myT, metaData$dx, mc.samples=128)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)

write.table(x.tt[, c(1,2), drop = F], "ttest/BrackenGvNGenusALDEx2Results.txt", sep = "\t", row.names = T)
write.table(x.tt[, c(3,4), drop = F], "wilcox/BrackenGvNGenusALDEx2Results.txt", sep = "\t", row.names = T)

metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.csv("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[, intersect(rownames(metaData), colnames(myT)), drop = F]
metaData<-metaData[intersect(rownames(metaData), colnames(myT)), , drop = F]

x <- aldex.clr(myT, metaData$dx, mc.samples=128)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)

write.table(x.tt[, c(1,2), drop = F], "ttest/BrackenGvNSpeciesALDEx2Results.txt", sep = "\t", row.names = T)
write.table(x.tt[, c(3,4), drop = F], "wilcox/BrackenGvNSpeciesALDEx2Results.txt", sep = "\t", row.names = T)

#SvR
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.csv("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[, intersect(rownames(metaData), colnames(myT)), drop = F]
metaData<-metaData[intersect(rownames(metaData), colnames(myT)), , drop = F]

x <- aldex.clr(myT, metaData$dx, mc.samples=128)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)

write.table(x.tt[, c(1,2), drop = F], "ttest/BrackenSvRPhylumALDEx2Results.txt", sep = "\t", row.names = T)
write.table(x.tt[, c(3,4), drop = F], "wilcox/BrackenSvRPhylumALDEx2Results.txt", sep = "\t", row.names = T)

metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.csv("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[, intersect(rownames(metaData), colnames(myT)), drop = F]
metaData<-metaData[intersect(rownames(metaData), colnames(myT)), , drop = F]

x <- aldex.clr(myT, metaData$dx, mc.samples=128)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)

write.table(x.tt[, c(1,2), drop = F], "ttest/BrackenSvRGenusALDEx2Results.txt", sep = "\t", row.names = T)
write.table(x.tt[, c(3,4), drop = F], "wilcox/BrackenSvRGenusALDEx2Results.txt", sep = "\t", row.names = T)

metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.csv("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[, intersect(rownames(metaData), colnames(myT)), drop = F]
metaData<-metaData[intersect(rownames(metaData), colnames(myT)), , drop = F]

x <- aldex.clr(myT, metaData$dx, mc.samples=128)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)

write.table(x.tt[, c(1,2), drop = F], "ttest/BrackenSvRSpeciesALDEx2Results.txt", sep = "\t", row.names = T)
write.table(x.tt[, c(3,4), drop = F], "wilcox/BrackenSvRSpeciesALDEx2Results.txt", sep = "\t", row.names = T)
