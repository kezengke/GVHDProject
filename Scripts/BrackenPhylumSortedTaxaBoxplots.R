#Normalized counts box plots with sorted pvals (gvhd vs no-gvhd)
rm(list = ls())

metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
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

ttest_res<-read.table("ttest/BrackenGvNPhylumResults.txt", sep = "\t", header = T, row.names = 1)
ttest_res<-ttest_res[rownames(myT), , drop = F]

ttest_res<-ttest_res[order(ttest_res$t_test_p), , drop = F]
ttest_res<-na.omit(ttest_res)

sigTT<-myT[rownames(ttest_res), , drop = F]
sigTT<-na.omit(sigTT)

pdf("Plots/BrackenPhylumGvNSortedTaxaBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("GVHD", "No GVHD")
for (i in 1:dim(sigTT)[1]){
  boxplot(unlist(sigTT[i,])~metaData$dx, outline = F,
          ylab=rownames(sigTT)[i], xlab="Treatment", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("royalblue3","gold1"))
  stripchart(unlist(sigTT[i,])~metaData$dx, method="jitter", 
             vertical=T, pch=19, add=T)
  
  mtext(paste("T-Test P-value:",signif(ttest_res[rownames(sigTT[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. T-Test P-value:",signif(ttest_res[rownames(sigTT[i,]),3],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()


#Normalized counts box plots with sorted pvals (sensitive vs refractory)
rm(list = ls())

metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
myT<-myT[, !(names(myT) %in% c("3921_07-18-17", "6673_12-07-16"))]
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

ttest_res<-read.table("ttest/BrackenSvRPhylumResults.txt", sep = "\t", header = T, row.names = 1)
ttest_res<-ttest_res[rownames(myT), , drop = F]

ttest_res<-ttest_res[order(ttest_res$t_test_p), , drop = F]
ttest_res<-na.omit(ttest_res)

sigTT<-myT[rownames(ttest_res), , drop = F]
sigTT<-na.omit(sigTT)

pdf("Plots/BrackenPhylumSvRSortedTaxaBoxPlots(Ttest).pdf", width=14, height=7)
par(mfrow=c(2, 4))
par(mar=c(5,6,4,1)+.1)
label = c("Steroid Refractory", "Steroid Sensitive")
for (i in 1:dim(sigTT)[1]){
  boxplot(unlist(sigTT[i,])~metaData$dx, outline = F,
          ylab=rownames(sigTT)[i], xlab="Treatment", 
          names = label,
          cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
          col=c("royalblue3","gold1"))
  stripchart(unlist(sigTT[i,])~metaData$dx, method="jitter", 
             vertical=T, pch=20, add=T)
  
  mtext(paste("T-Test P-value:",signif(ttest_res[rownames(sigTT[i,]),2],digits = 3)),3, 0.9, cex = 0.6, adj = 0)
  mtext(paste("Adj. T-Test P-value:",signif(ttest_res[rownames(sigTT[i,]),3],digits = 3)),3, 0.9, cex = 0.6, adj = 1)
}
dev.off()

