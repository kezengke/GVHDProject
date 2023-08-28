#Filtering calculation (only prints, no saving for the results)

#phylum
rm(list = ls())
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
normT<-myT
for (i in 1:ncol(normT)) {
  normT[,i]<-normT[,i]/n[i]
}
normT<-log10(normT*(sumx/ncol(normT))+1)
normT<-data.frame(normT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(normT)<2)

#percentage of remained taxa from total
sum(myT[-lowAbundance,])/sum(myT)
#percentage of filtered taxa from total
sum(myT[lowAbundance,])/sum(myT)

#relative abundance of the top 1 taxa that got filtered
mean(as.numeric(myT[names(sort(rowSums(myT[lowAbundance, ]), decreasing = T)[1]), ])/sum(myT))

#genus
rm(list = ls())
myT<-read.table("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
normT<-myT
for (i in 1:ncol(normT)) {
  normT[,i]<-normT[,i]/n[i]
}
normT<-log10(normT*(sumx/ncol(normT))+1)
normT<-data.frame(normT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(normT)<2)

#percentage of remained taxa from total
sum(myT[-lowAbundance,])/sum(myT)
#percentage of filtered taxa from total
sum(myT[lowAbundance,])/sum(myT)

#relative abundance of the top 1 taxa that got filtered
mean(as.numeric(myT[names(sort(rowSums(myT[lowAbundance, ]), decreasing = T)[1]), ])/sum(myT))

#species
rm(list = ls())
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
normT<-myT
for (i in 1:ncol(normT)) {
  normT[,i]<-normT[,i]/n[i]
}
normT<-log10(normT*(sumx/ncol(normT))+1)
normT<-data.frame(normT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(normT)<2)

#percentage of remained taxa from total
sum(myT[-lowAbundance,])/sum(myT)
#percentage of filtered taxa from total
sum(myT[lowAbundance,])/sum(myT)

#relative abundance of the top 1 taxa that got filtered
mean(as.numeric(myT[names(sort(rowSums(myT[lowAbundance, ]), decreasing = T)[1]), ])/sum(myT))

#pathway
rm(list = ls())
myT<-read.delim("CountsTables/unstratified.txt", 
                row.name = 1, sep = "\t", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[ , -which(colSums(myT) == 0)]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
normT<-myT
for (i in 1:ncol(normT)) {
  normT[,i]<-normT[,i]/n[i]
}
normT<-log10(normT*(sumx/ncol(normT))+1)
normT<-data.frame(normT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(normT)<1)

#percentage of remained taxa from total
sum(myT[-lowAbundance,])/sum(myT)
#percentage of filtered taxa from total
sum(myT[lowAbundance,])/sum(myT)

#relative abundance of the top 1 taxa that got filtered
mean(as.numeric(myT[names(sort(rowSums(myT[lowAbundance, ]), decreasing = T)[1]), ])/sum(myT))
