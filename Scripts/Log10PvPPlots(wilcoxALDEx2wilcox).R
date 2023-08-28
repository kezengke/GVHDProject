#log10 pvp for different statistics methods (wilcoxon vs. ALDEx2 wilcoxon)
rm(list = ls())
library(scales)

pvpPlotFun <- function(result1, result2, testType1, testType2, condition, taxaLevel) {
  result1$direction<-result1[, 2]/abs(result1[, 2])
  result1$logP<-log10(result1[, 2])
  result1$p_directed<-result1$direction*result1$logP
  
  result2$direction<-result1$direction
  result2$logP<-log10(result2[, 1])
  result2$p_directed<-result2$direction*result2$logP
  
  plot(result1$p_directed[result1[, 3]<0.05 & result2[, 3]<0.05], 
       result2$p_directed[result1[, 3]<0.05 & result2[, 3]<0.05], 
       xlim=range(na.omit(result1$p_directed)),
       ylim=range(na.omit(result2$p_directed)),
       main=paste(testType1, "vs.", testType2, taxaLevel, condition), 
       xlab=paste0(testType1, "(log10)"), 
       ylab=paste0(testType2, "(log10)"), 
       col=alpha("turquoise1", 0.3), pch=19,
       cex.lab = 1.4, cex.main = 1.5, cex.axis = 1.4)
  points(result1$p_directed[result1[, 3]>=0.05 & result2[, 2]>=0.05], 
         result2$p_directed[result1[, 3]>=0.05 & result2[, 2]>=0.05],
         col=alpha("green", 0.3), pch=19)
  points(result1$p_directed[result1[, 3]>=0.05 & result2[, 2]<0.05], 
         result2$p_directed[result1[, 3]>=0.05 & result2[, 2]<0.05],
         col=alpha("red", 0.3), pch=19)
  points(result1$p_directed[result1[, 3]<0.05 & result2[, 2]>=0.05], 
         result2$p_directed[result1[, 3]<0.05 & result2[, 2]>=0.05],
         col=alpha("purple", 0.3), pch=19)
  abline(0,1, col = "gray")
  lgnames<-c(paste0("Both(", length(result1$p_directed[result1[, 3]<0.05 & result2[, 2]<0.05]), ")"), 
             paste0("Neither(", length(result1$p_directed[result1[, 3]>=0.05 & result2[, 2]>=0.05]), ")"), 
             paste0(testType2, "(", length(result1$p_directed[result1[, 3]>=0.05 & result2[, 2]<0.05]), ")"), 
             paste0(testType1, "(", length(result1$p_directed[result1[, 3]<0.05 & result2[, 2]>=0.05]), ")"))
  legend("bottomright", lgnames,
         col = c("turquoise1", "green", "red", "purple"),
         inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", pch = 19, cex = 0.8)
}

pdf("Plots/Log10PvalPlots(TtestAldex2Wilcoxon).pdf", width=15, height=10)
par(mfrow=c(2,3))

wilcoxRES<-read.table("wilcox/BrackenGvNPhylumResults.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-read.table("wilcox/BrackenGvNPhylumALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-aldex2Wilcox[rownames(wilcoxRES), , drop = F]
pvpPlotFun(wilcoxRES, aldex2Wilcox, "Wilcoxon", "ALDEx2Wilcox", "GvN", "Phylum")

wilcoxRES<-read.table("wilcox/BrackenGvNGenusResults.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-read.table("wilcox/BrackenGvNGenusALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-aldex2Wilcox[rownames(wilcoxRES), , drop = F]
pvpPlotFun(wilcoxRES, aldex2Wilcox, "Wilcoxon", "ALDEx2Wilcox", "GvN", "Genus")

wilcoxRES<-read.table("wilcox/BrackenGvNSpeciesResults.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-read.table("wilcox/BrackenGvNSpeciesALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-aldex2Wilcox[rownames(wilcoxRES), , drop = F]
pvpPlotFun(wilcoxRES, aldex2Wilcox, "Wilcoxon", "ALDEx2Wilcox", "GvN", "Species")

wilcoxRES<-read.table("wilcox/BrackenSvRPhylumResults.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-read.table("wilcox/BrackenSvRPhylumALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-aldex2Wilcox[rownames(wilcoxRES), , drop = F]
pvpPlotFun(wilcoxRES, aldex2Wilcox, "Wilcoxon", "ALDEx2Wilcox", "SvR", "Phylum")

wilcoxRES<-read.table("wilcox/BrackenSvRGenusResults.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-read.table("wilcox/BrackenSvRGenusALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-aldex2Wilcox[rownames(wilcoxRES), , drop = F]
pvpPlotFun(wilcoxRES, aldex2Wilcox, "Wilcoxon", "ALDEx2Wilcox", "SvR", "Genus")

wilcoxRES<-read.table("wilcox/BrackenSvRSpeciesResults.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-read.table("wilcox/BrackenSvRSpeciesALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
aldex2Wilcox<-aldex2Wilcox[rownames(wilcoxRES), , drop = F]
pvpPlotFun(wilcoxRES, aldex2Wilcox, "Wilcoxon", "ALDEx2Wilcox", "SvR", "Species")

dev.off()