#MDS plot at genus level (gvhd vs no-gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-metaData[colnames(myT),]
metaData<-data.frame(metaData)
rownames(metaData)<-colnames(myT)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)

library("vegan")
pdf("Plots/GvNGenusMDS.pdf",onefile = T)
par(mar=c(5,6,4,1)+.1)
par(mfrow=c(1,1))
circlecol<-c("red","blue")

MDS<-capscale(t(myT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$metaData)]
pval<-adonis(t(myT)~metaData$metaData, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("GVHD/No GVHD \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$metaData, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$metaData, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="no_gvhd",label=T,font=2,cex=1) 
legend("topright",c("gvhd","no gvhd"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()

#MDS plot at species level (gvhd vs no-gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-metaData[colnames(myT),]
metaData<-data.frame(metaData)
rownames(metaData)<-colnames(myT)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)

library("vegan")
pdf("Plots/GvNSpeciesMDS.pdf",onefile = T)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circlecol<-c("red","blue")

MDS<-capscale(t(myT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$metaData)]
pval<-adonis(t(myT)~metaData$metaData, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("GVHD/No GVHD \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$metaData, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$metaData, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="no_gvhd",label=T,font=2,cex=1) 
legend("topright",c("gvhd","no gvhd"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()

#MDS plot at phylum level (gvhd vs no-gvhd)
rm(list = ls())
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-metaData[colnames(myT),]
metaData<-data.frame(metaData)
rownames(metaData)<-colnames(myT)
metaData<-na.omit(metaData)
myT<-myT[,rownames(metaData)]
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
for (i in 1:ncol(myT)) {
  myT[,i]<-myT[,i]/n[i]
}
myT<-log10(myT*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)

library("vegan")
pdf("Plots/GvNPhylumMDS.pdf",onefile = T)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circlecol<-c("red","blue")
#status
MDS<-capscale(t(myT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$metaData)]
pval<-adonis(t(myT)~metaData$metaData, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("GVHD/No GVHD \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$metaData, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$metaData, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="no_gvhd",label=T,font=2,cex=1) 
legend("topright",c("gvhd","no gvhd"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()

#MDS plot at genus level (sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
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

library("vegan")
pdf("Plots/SvRGenusMDS.pdf",onefile = T)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circlecol<-c("red","blue")

MDS<-capscale(t(myT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$dx)]
pval<-adonis(t(myT)~metaData$dx, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("Steroid Sensitive/Steroid Refractory \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="steroid_sensitive_gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="steroid_refractory_gvhd",label=T,font=2,cex=1) 
legend("topright",c("sensitive","refractory"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()

#MDS plot at genus level (sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
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

library("vegan")
pdf("Plots/SvRSpeciesMDS.pdf",onefile = T)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circlecol<-c("red","blue")

MDS<-capscale(t(myT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$dx)]
pval<-adonis(t(myT)~metaData$dx, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("Steroid Sensitive/Steroid Refractory \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="steroid_sensitive_gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="steroid_refractory_gvhd",label=T,font=2,cex=1) 
legend("topright",c("sensitive","refractory"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()

#MDS plot at phylum level (sensitive vs refractory)
rm(list = ls())
metaData<-read.table("metaSvR.txt", sep = "\t", header = T, row.names = 1)
myT<-read.table("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
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

library("vegan")
pdf("Plots/SvRPhylumMDS.pdf",onefile = T)
par(mfrow=c(1,1))
par(mar=c(5,6,4,1)+.1)
circlecol<-c("red","blue")

MDS<-capscale(t(myT)~1,distance = "bray")
percentVariance<-MDS$CA$eig/sum(eigenvals(MDS))*100
cols<-c("red","blue")
col1<-circlecol[factor(metaData$dx)]
pval<-adonis(t(myT)~metaData$dx, method="bray")$aov.tab$`Pr(>F)`[1]
statusPlot<-ordiplot(MDS,choices=c(1,2),type="none",cex.lab=2.5, cex.axis = 2, cex.main = 2, 
                     xlab=paste("MDS1  ", format(percentVariance[1], digits = 4), "%", sep = ""),
                     ylab=paste("MDS2  ", format(percentVariance[2], digits = 4), "%", sep = ""),
                     main=paste("Steroid Sensitive/Steroid Refractory \nP-value:",pval))
points(statusPlot,"sites",col=adjustcolor(col1, alpha.f = 0.5),pch=16,cex=2.5)
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="red",show.groups="steroid_sensitive_gvhd",label=T,font=2,cex=1) 
ordiellipse(statusPlot, metaData$dx, kind="se", conf=0.95, lwd=4, draw = "lines", col="blue",show.groups="steroid_refractory_gvhd",label=T,font=2,cex=1) 
legend("topright",c("sensitive","refractory"),col=circlecol[1:2],cex=1.8,pch=16,bty = "n")

dev.off()
