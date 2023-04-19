#MDS plot for humann3 unstratified table (gvhd vs no gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/unstratified.txt",
                sep = "\t", header = T, row.names = 1, check.names = F)
myT<-myT[ , -which(colSums(myT) == 0)] #get rid of sample with all 0s.
metaData<-metaData[colnames(myT), , drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData)]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
NormT<-myT
for (i in 1:ncol(NormT)) {
  NormT[,i]<-NormT[,i]/n[i]
}
NormT<-log10(NormT*(sumx/ncol(NormT))+1)
NormT<-data.frame(NormT, check.names = F)

library("vegan")
pdf("Plots/Humann3GvNUnstratifiedMDS.pdf",onefile = T)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circlecol<-c("red","blue")
#status
MDS<-capscale(t(NormT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$dx)]
pval<-adonis(t(NormT)~metaData$dx, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("GVHD/No GVHD \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="no_gvhd",label=T,font=2,cex=1) 
legend("topright",c("gvhd","no gvhd"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()

#MDS plot for humann3 unstratified table (sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/unstratified.txt",
                sep = "\t", header = T, row.names = 1, check.names = F)
myT<-myT[ , -which(colSums(myT) == 0)] #get rid of sample with all 0s.
metaData<-metaData[colnames(myT), , drop = F]
metaData<-na.omit(metaData)
myT<-myT[, rownames(metaData)]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
NormT<-myT
for (i in 1:ncol(NormT)) {
  NormT[,i]<-NormT[,i]/n[i]
}
NormT<-log10(NormT*(sumx/ncol(NormT))+1)
NormT<-data.frame(NormT, check.names = F)

library("vegan")
pdf("Plots/Humann3SvRUnstratifiedMDS.pdf",onefile = T)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circlecol<-c("red","blue")
#status
MDS<-capscale(t(NormT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$dx)]
pval<-adonis(t(NormT)~metaData$dx, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("Steroid Sensitive/Steroid Refractory \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="steroid_sensitive_gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="steroid_refractory_gvhd",label=T,font=2,cex=1) 
legend("topright",c("sensitive","refractory"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()
