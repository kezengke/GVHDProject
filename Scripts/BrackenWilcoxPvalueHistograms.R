#bracken results p-value histograms
rm(list = ls())
brackenGvNPhylumT<-read.delim("wilcox/BrackenGvNPhylumResults.txt", sep = "\t", header = T, row.names = 1)
brackenGvNGenusT<-read.delim("wilcox/BrackenGvNGenusResults.txt", sep = "\t", header = T, row.names = 1)
brackenGvNSpeciesT<-read.delim("wilcox/BrackenGvNSpeciesResults.txt", sep = "\t", header = T, row.names = 1)

brackenSvRPhylumT<-read.delim("wilcox/BrackenSvRPhylumResults.txt", sep = "\t", header = T, row.names = 1)
brackenSvRGenusT<-read.delim("wilcox/BrackenSvRGenusResults.txt", sep = "\t", header = T, row.names = 1)
brackenSvRSpeciesT<-read.delim("wilcox/BrackenSvRSpeciesResults.txt", sep = "\t", header = T, row.names = 1)

pdf("Plots/BrackenWilcoxPvalHistograms.pdf", width=18, height=12)
par(mfrow=c(2, 3))
par(mar=c(5, 6, 4, 1)+.1)
#bracken
hist(brackenGvNPhylumT$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Phylum GvN Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNGenusT$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Genus GvN Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNSpeciesT$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Species GvN Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRPhylumT$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Phylum SvR Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRGenusT$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Genus SvR Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRSpeciesT$wilcox_p, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Species SvR Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
