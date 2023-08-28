#log10 pvp for different statistics methods (ttest vs. wilcoxon)
rm(list = ls())
library(scales)

pvpPlotFun <- function(result1, result2, testType1, testType2, condition, taxaLevel) {
  result1$direction<-result1[, 1]/abs(result1[, 1])
  result1$logP<-log10(result1[, 2])
  result1$p_directed<-result1$direction*result1$logP
  
  result2$direction<-result1$direction
  result2$logP<-log10(result2[, 2])
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
  points(result1$p_directed[result1[, 3]>=0.05 & result2[, 3]>=0.05], 
         result2$p_directed[result1[, 3]>=0.05 & result2[, 3]>=0.05],
         col=alpha("green", 0.3), pch=19)
  points(result1$p_directed[result1[, 3]>=0.05 & result2[, 3]<0.05], 
         result2$p_directed[result1[, 3]>=0.05 & result2[, 3]<0.05],
         col=alpha("red", 0.3), pch=19)
  points(result1$p_directed[result1[, 3]<0.05 & result2[, 3]>=0.05], 
         result2$p_directed[result1[, 3]<0.05 & result2[, 3]>=0.05],
         col=alpha("purple", 0.3), pch=19)
  abline(0,1, col = "gray")
  lgnames<-c(paste0("Both(", length(result1$p_directed[result1[, 3]<0.05 & result2[, 3]<0.05]), ")"), 
             paste0("Neither(", length(result1$p_directed[result1[, 3]>=0.05 & result2[, 3]>=0.05]), ")"), 
             paste0(testType2, "(", length(result1$p_directed[result1[, 3]>=0.05 & result2[, 3]<0.05]), ")"), 
             paste0(testType1, "(", length(result1$p_directed[result1[, 3]<0.05 & result2[, 3]>=0.05]), ")"))
  legend("bottomright", lgnames,
         col = c("turquoise1", "green", "red", "purple"),
         inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n", pch = 19, cex = 0.8)
}

pdf("Plots/Log10PvalPlots(TtestWilcox).pdf", width=15, height=10)
par(mfrow=c(2,3))

ttestRES<-read.table("ttest/BrackenGvNPhylumResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-read.table("wilcox/BrackenGvNPhylumResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-wilcoxRES[rownames(ttestRES), , drop = F]
pvpPlotFun(ttestRES, wilcoxRES, "T-test", "Wilcoxon", "GvN", "Phylum")

ttestRES<-read.table("ttest/BrackenGvNGenusResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-read.table("wilcox/BrackenGvNGenusResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-wilcoxRES[rownames(ttestRES), , drop = F]
pvpPlotFun(ttestRES, wilcoxRES, "T-test", "Wilcoxon", "GvN", "Genus")

ttestRES<-read.table("ttest/BrackenGvNSpeciesResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-read.table("wilcox/BrackenGvNSpeciesResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-wilcoxRES[rownames(ttestRES), , drop = F]
pvpPlotFun(ttestRES, wilcoxRES, "T-test", "Wilcoxon", "GvN", "Species")

ttestRES<-read.table("ttest/BrackenSvRPhylumResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-read.table("wilcox/BrackenSvRPhylumResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-wilcoxRES[rownames(ttestRES), , drop = F]
pvpPlotFun(ttestRES, wilcoxRES, "T-test", "Wilcoxon", "SvR", "Phylum")

ttestRES<-read.table("ttest/BrackenSvRGenusResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-read.table("wilcox/BrackenSvRGenusResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-wilcoxRES[rownames(ttestRES), , drop = F]
pvpPlotFun(ttestRES, wilcoxRES, "T-test", "Wilcoxon", "SvR", "Genus")

ttestRES<-read.table("ttest/BrackenSvRSpeciesResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-read.table("wilcox/BrackenSvRSpeciesResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxRES<-wilcoxRES[rownames(ttestRES), , drop = F]
pvpPlotFun(ttestRES, wilcoxRES, "T-test", "Wilcoxon", "SvR", "Species")

dev.off()