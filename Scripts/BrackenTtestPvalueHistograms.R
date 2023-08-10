#bracken results p-value histograms
rm(list = ls())
brackenGvNPhylumT<-read.delim("ttest/BrackenGvNPhylumResults.txt", sep = "\t", header = T, row.names = 1)
brackenGvNGenusT<-read.delim("ttest/BrackenGvNGenusResults.txt", sep = "\t", header = T, row.names = 1)
brackenGvNSpeciesT<-read.delim("ttest/BrackenGvNSpeciesResults.txt", sep = "\t", header = T, row.names = 1)

brackenSvRPhylumT<-read.delim("ttest/BrackenSvRPhylumResults.txt", sep = "\t", header = T, row.names = 1)
brackenSvRGenusT<-read.delim("ttest/BrackenSvRGenusResults.txt", sep = "\t", header = T, row.names = 1)
brackenSvRSpeciesT<-read.delim("ttest/BrackenSvRSpeciesResults.txt", sep = "\t", header = T, row.names = 1)

pdf("Plots/BrackenTtestPvalHistograms.pdf", width=18, height=12)
par(mfrow=c(2, 3))
par(mar=c(5, 6, 4, 1)+.1)
#bracken
hist(brackenGvNPhylumT$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Phylum GvN Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNGenusT$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Genus GvN Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNSpeciesT$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Species GvN Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRPhylumT$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Phylum SvR Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRGenusT$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Genus SvR Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRSpeciesT$t_test_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Species SvR Ttest P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
