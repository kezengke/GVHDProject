#bracken results p-value histograms
pdf("Plots/BrackenALDEx2PvalHistograms.pdf", width=18, height=12)
par(mfrow=c(2, 3))
par(mar=c(5, 6, 4, 1)+.1)

rm(list = ls())
brackenGvNPhylumT<-read.delim("ttest/BrackenGvNPhylumALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenGvNGenusT<-read.delim("ttest/BrackenGvNGenusALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenGvNSpeciesT<-read.delim("ttest/BrackenGvNSpeciesALDEx2Results.txt", sep = "\t", header = T, row.names = 1)

brackenSvRPhylumT<-read.delim("ttest/BrackenSvRPhylumALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenSvRGenusT<-read.delim("ttest/BrackenSvRGenusALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenSvRSpeciesT<-read.delim("ttest/BrackenSvRSpeciesALDEx2Results.txt", sep = "\t", header = T, row.names = 1)

pdf("Plots/BrackenALDEx2PvalHistograms.pdf", width=18, height=12)
par(mfrow=c(2, 3))
par(mar=c(5, 6, 4, 1)+.1)
#bracken
hist(brackenGvNPhylumT$we.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Phylum GvN ALDEx2 Welch-T P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNGenusT$we.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Genus GvN ALDEx2 Welch-T P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNSpeciesT$we.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Species GvN ALDEx2 Welch-T P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRPhylumT$we.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Phylum SvR ALDEx2 Welch-T P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRGenusT$we.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Genus SvR ALDEx2 Welch-T P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRSpeciesT$we.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "tan2",
     main = "Bracken-Species SvR ALDEx2 Welch-T P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

rm(list = ls())
brackenGvNPhylumT<-read.delim("wilcox/BrackenGvNPhylumALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenGvNGenusT<-read.delim("wilcox/BrackenGvNGenusALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenGvNSpeciesT<-read.delim("wilcox/BrackenGvNSpeciesALDEx2Results.txt", sep = "\t", header = T, row.names = 1)

brackenSvRPhylumT<-read.delim("wilcox/BrackenSvRPhylumALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenSvRGenusT<-read.delim("wilcox/BrackenSvRGenusALDEx2Results.txt", sep = "\t", header = T, row.names = 1)
brackenSvRSpeciesT<-read.delim("wilcox/BrackenSvRSpeciesALDEx2Results.txt", sep = "\t", header = T, row.names = 1)

#bracken
hist(brackenGvNPhylumT$wi.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Phylum GvN ALDEx2 Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNGenusT$wi.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Genus GvN ALDEx2 P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenGvNSpeciesT$wi.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Species GvN ALDEx2 Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRPhylumT$wi.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Phylum SvR ALDEx2 Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRGenusT$wi.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Genus SvR ALDEx2 Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

hist(brackenSvRSpeciesT$wi.ep, breaks=seq(0, 1, 0.05), xlab = "p-value", col = "cornflowerblue",
     main = "Bracken-Species SvR ALDEx2 Wilcox P-vals", 
     xlim=c(0,1), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.4)
abline(v=0.05, col=gray(.5), lty=2)
mtext(at=0.05, side=3, text=0.05, col=gray(.5))

dev.off()
