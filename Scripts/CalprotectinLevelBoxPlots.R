library("coin")
pdf("Plots/CalprotectinLevelBoxplots.pdf", width=10, height=5)
par(mar=c(5,6,4,1)+.1)
par(mfrow=c(1,2))
rm(list = ls())
calprotectin<-read.table("CountsTables/gvhd_calprotectin_results.tsv", sep = "\t", header = T)
calprotectin<-na.omit(calprotectin)
rownames(calprotectin)<-calprotectin$sample_code
calprotectin<-calprotectin[, 3, drop = F]
metaData<-read.delim("metaGvN.txt", sep = "\t", row.names = 1, header = T)
t_test_p<-t.test(unlist(calprotectin)~factor(metaData$dx))$p.value
wilcox_p<-pvalue(wilcox_test(unlist(calprotectin)~factor(metaData$dx)))
label = c("GVHD", "No GVHD")
boxplot(unlist(calprotectin)~factor(metaData$dx), outline = F,
        cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Calprotectin levels(ug/mg)", xlab = "conditions", 
        names = label,
        main = paste0("P-value:", signif(wilcox_p, 4)))
stripchart(unlist(calprotectin)~factor(metaData$dx), method="jitter", 
           vertical=T, pch=19, add=T)

rm(list = ls())
calprotectin<-read.table("CountsTables/gvhd_calprotectin_results.tsv", sep = "\t", header = T)
calprotectin<-na.omit(calprotectin)
rownames(calprotectin)<-calprotectin$sample_code
calprotectin<-calprotectin[, 3, drop = F]
metaData<-read.delim("metaSvR.txt", sep = "\t", row.names = 1, header = T)
calprotectin<-calprotectin[rownames(metaData), , drop = F]
t_test_p<-t.test(unlist(calprotectin)~factor(metaData$dx))$p.value
wilcox_p<-pvalue(wilcox_test(unlist(calprotectin)~factor(metaData$dx)))
label = c("Steroid Refractory", "Steroid Sensitive")
boxplot(unlist(calprotectin)~factor(metaData$dx), outline = F,
        cex.lab = 1.5, cex.axis = 1, cex.main = 1.8, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Calprotectin levels(ug/mg)", xlab = "conditions", 
        names = label,
        main = paste0("P-value:", signif(wilcox_p, 4)))
stripchart(unlist(calprotectin)~factor(metaData$dx), method="jitter", 
           vertical=T, pch=19, add=T)

dev.off()
