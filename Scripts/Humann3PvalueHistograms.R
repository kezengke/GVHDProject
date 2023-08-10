#humann3 results p-value histograms
rm(list = ls())
ttestGvN<-read.delim("ttest/Humann3GvNUnstratifiedResults.txt", sep = "\t", header = T, row.names = 1)
ttestSvR<-read.delim("ttest/Humann3SvRUnstratifiedResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxGvN<-read.delim("wilcox/Humann3GvNUnstratifiedResults.txt", sep = "\t", header = T, row.names = 1)
wilcoxSvR<-read.delim("wilcox/Humann3SvRUnstratifiedResults.txt", sep = "\t", header = T, row.names = 1)

pdf("Plots/Humann3PvalHistograms.pdf", width=12, height=12)
par(mfrow=c(2, 2))
par(mar=c(5, 6, 4, 1)+.1)

hist(ttestGvN$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Humann3 GvN Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(ttestSvR$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Humann3 SvR Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(wilcoxGvN$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Humann3 GvN Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(wilcoxSvR$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Humann3 SvR Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()