#Normalized counts box plots with sorted pvals (gvhd vs no-gvhd)
rm(list = ls())

metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-metaData[colnames(myT),, drop = F]
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

wilcox_res<-read.table("wilcox/BrackenGvNPhylumResults.txt", sep = "\t", header = T, row.names = 1)
wilcox_res<-wilcox_res[rownames(myT), , drop = F]

wilcox_res<-wilcox_res[order(wilcox_res$wilcox_p), , drop = F]
wilcox_res<-na.omit(wilcox_res)

sigWT<-myT[rownames(wilcox_res), , drop = F]
sigWT<-na.omit(sigWT)

pdf("Plots/BrackenPhylumGvNSortedTaxaBoxPlots(Wilcox).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("GVHD", "No GVHD")
for (i in 1:dim(sigWT)[1]){
  boxplot(unlist(sigWT[i,])~metaData$dx, outline = F,
          ylab=rownames(sigWT)[i], xlab="Treatment", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("royalblue3","gold1"))
  stripchart(unlist(sigWT[i,])~metaData$dx, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("Wilcoxon P-value:",signif(wilcox_res[rownames(sigWT[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Wilcoxon P-value:",signif(wilcox_res[rownames(sigWT[i,]),3],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()


#Normalized counts box plots with sorted pvals (sensitive vs refractory)
rm(list = ls())

metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-metaData[colnames(myT),, drop = F]
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

wilcox_res<-read.table("wilcox/BrackenSvRPhylumResults.txt", sep = "\t", header = T, row.names = 1)
wilcox_res<-wilcox_res[rownames(myT), , drop = F]

wilcox_res<-wilcox_res[order(wilcox_res$wilcox_p), , drop = F]
wilcox_res<-na.omit(wilcox_res)

sigWT<-myT[rownames(wilcox_res), , drop = F]
sigWT<-na.omit(sigWT)

pdf("Plots/BrackenPhylumSvRSortedTaxaBoxPlots(Wilcox).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("Steroid Refractory", "Steroid Sensitive")
for (i in 1:dim(sigWT)[1]){
  boxplot(unlist(sigWT[i,])~metaData$dx, outline = F,
          ylab=rownames(sigWT)[i], xlab="Treatment", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("royalblue3","gold1"))
  stripchart(unlist(sigWT[i,])~metaData$dx, method="jitter", 
             vertical=T, pch=20, add=T)
  
  mtext(paste("Wilcoxon P-value:",signif(wilcox_res[rownames(sigWT[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. Wilcoxon P-value:",signif(wilcox_res[rownames(sigWT[i,]),3],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()

