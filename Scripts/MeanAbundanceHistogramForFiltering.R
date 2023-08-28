pdf("Plots/MeanAbundance(Log).pdf", width=21, height=14)
par(mfrow=c(2,3))
#Bracken
#phylum
rm(list = ls())
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
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "Bracken(phylum)", 
     xlab = "taxa mean abundance", ylab = "log10(Frequency)")
#genus
rm(list = ls())
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
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "Bracken(genus)", 
     xlab = "taxa mean abundance", ylab = "log10(Frequency)")
#species
rm(list = ls())
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
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "Bracken(species)", 
     xlab = "taxa mean abundance", ylab = "log10(Frequency)")

#humann
rm(list = ls())
myT<-read.delim("CountsTables/unstratified.txt",
                sep = "\t", header = T, row.names = 1, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
myT<-myT[ , -which(colSums(myT) == 0)]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
        myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
myT<-na.omit(myT)
hist.data<-hist(rowMeans(myT), breaks = 50, plot=F)
hist.data$counts<-log10(hist.data$counts+1)
plot(hist.data, cex.lab=1.5,
     main = "Humann3(unstratified)", 
     xlab = "pathway mean abundance", ylab = "log10(Frequency)")

dev.off()

