#phylum
rm(list = ls())
calprotectin<-read.table("CountsTables/gvhd_calprotectin_results.tsv", sep = "\t", header = T)
calprotectin<-na.omit(calprotectin)
rownames(calprotectin)<-calprotectin$sample_code
calprotectin<-calprotectin[, 3, drop = F]
myT<-read.delim("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-read.delim("metaSvR.txt", sep = "\t", row.names = 1, header = T)

metaData<-metaData[intersect(colnames(myT), rownames(metaData)), , drop = F]
calprotectin<-calprotectin[intersect(colnames(myT), rownames(calprotectin)), , drop = F]
myT<-myT[, metaData$dx == "steroid_refractory_gvhd", drop = F]
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
calprotectin<-calprotectin[intersect(colnames(myT), rownames(calprotectin)), , drop = F]

cor_pvals<-vector()
lm_pvals<-vector()
for (i in 1:nrow(myT)) {
  cor_pvals[i]<-cor.test(unlist(myT[i, ]), calprotectin[,1], method = "kendall")$p.value
  lm_pvals[i]<-anova(lm(unlist(myT[i, ])~calprotectin[,1]))$"Pr(>F)"[1]
}
adj_corP<-p.adjust(cor_pvals, method = "BH")
adj_lmP<-p.adjust(lm_pvals, method = "BH")

pdf("Plots/BrackenPhylumRefractoryCalprotectinCorLmANOVAPvalScatterPlots(SortedByKendall).pdf", width=24, height=12)
par(mfrow=c(2,4))
par(mar=c(5,6,4,1)+.1)
plot_order<-order(cor_pvals)
for (i in 1:nrow(myT)) {
  mainText<-paste0("Kendall pval=", signif(cor_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_corP[plot_order[i]], 3),
                   '\n', "ANOVA pval=", signif(lm_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_lmP[plot_order[i]], 3))
  plot(unlist(myT[plot_order[i], ]), calprotectin[,1], main = mainText,
       cex.lab = 2.5, cex.axis = 2, cex.main = 2, 
       xlab = rownames(myT)[plot_order[i]], ylab = "Calprotectin Level",
       col = "palevioletred", pch = 19)
}
dev.off()

pdf("Plots/BrackenPhylumRefractoryCalprotectinCorLmANOVAPvalScatterPlots(SortedByANOVA).pdf", width=24, height=12)
par(mfrow=c(2,4))
par(mar=c(5,6,4,1)+.1)
plot_order<-order(lm_pvals)
for (i in 1:nrow(myT)) {
  mainText<-paste0("Kendall pval=", signif(cor_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_corP[plot_order[i]], 3),
                   '\n', "ANOVA pval=", signif(lm_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_lmP[plot_order[i]], 3))
  plot(unlist(myT[plot_order[i], ]), calprotectin[,1], main = mainText,
       cex.lab = 2.5, cex.axis = 2, cex.main = 2, 
       xlab = rownames(myT)[plot_order[i]], ylab = "Calprotectin Level",
       col = "palevioletred", pch = 19)
}
dev.off()

#genus
rm(list = ls())
calprotectin<-read.table("CountsTables/gvhd_calprotectin_results.tsv", sep = "\t", header = T)
calprotectin<-na.omit(calprotectin)
rownames(calprotectin)<-calprotectin$sample_code
calprotectin<-calprotectin[, 3, drop = F]
myT<-read.delim("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-read.delim("metaSvR.txt", sep = "\t", row.names = 1, header = T)

metaData<-metaData[intersect(colnames(myT), rownames(metaData)), , drop = F]
calprotectin<-calprotectin[intersect(colnames(myT), rownames(calprotectin)), , drop = F]
myT<-myT[, metaData$dx == "steroid_refractory_gvhd", drop = F]
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
calprotectin<-calprotectin[intersect(colnames(myT), rownames(calprotectin)), , drop = F]

cor_pvals<-vector()
lm_pvals<-vector()
for (i in 1:nrow(myT)) {
  cor_pvals[i]<-cor.test(unlist(myT[i, ]), calprotectin[,1], method = "kendall")$p.value
  lm_pvals[i]<-anova(lm(unlist(myT[i, ])~calprotectin[,1]))$"Pr(>F)"[1]
}
adj_corP<-p.adjust(cor_pvals, method = "BH")
adj_lmP<-p.adjust(lm_pvals, method = "BH")

pdf("Plots/BrackenGenusRefractoryCalprotectinCorLmANOVAPvalScatterPlots(SortedByKendall).pdf", width=24, height=12)
par(mfrow=c(2,4))
par(mar=c(5,6,4,1)+.1)
plot_order<-order(cor_pvals)
for (i in 1:nrow(myT)) {
  mainText<-paste0("Kendall pval=", signif(cor_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_corP[plot_order[i]], 3),
                   '\n', "ANOVA pval=", signif(lm_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_lmP[plot_order[i]], 3))
  plot(unlist(myT[plot_order[i], ]), calprotectin[,1], main = mainText,
       cex.lab = 2.5, cex.axis = 2, cex.main = 2, 
       xlab = rownames(myT)[plot_order[i]], ylab = "Calprotectin Level",
       col = "palevioletred", pch = 19)
}
dev.off()

pdf("Plots/BrackenGenusRefractoryCalprotectinCorLmANOVAPvalScatterPlots(SortedByANOVA).pdf", width=24, height=12)
par(mfrow=c(2,4))
par(mar=c(5,6,4,1)+.1)
plot_order<-order(lm_pvals)
for (i in 1:nrow(myT)) {
  mainText<-paste0("Kendall pval=", signif(cor_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_corP[plot_order[i]], 3),
                   '\n', "ANOVA pval=", signif(lm_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_lmP[plot_order[i]], 3))
  plot(unlist(myT[plot_order[i], ]), calprotectin[,1], main = mainText,
       cex.lab = 2.5, cex.axis = 2, cex.main = 2, 
       xlab = rownames(myT)[plot_order[i]], ylab = "Calprotectin Level",
       col = "palevioletred", pch = 19)
}
dev.off()

#species
rm(list = ls())
calprotectin<-read.table("CountsTables/gvhd_calprotectin_results.tsv", sep = "\t", header = T)
calprotectin<-na.omit(calprotectin)
rownames(calprotectin)<-calprotectin$sample_code
calprotectin<-calprotectin[, 3, drop = F]
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-read.delim("metaSvR.txt", sep = "\t", row.names = 1, header = T)

metaData<-metaData[intersect(colnames(myT), rownames(metaData)), , drop = F]
calprotectin<-calprotectin[intersect(colnames(myT), rownames(calprotectin)), , drop = F]
myT<-myT[, metaData$dx == "steroid_refractory_gvhd", drop = F]
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
calprotectin<-calprotectin[intersect(colnames(myT), rownames(calprotectin)), , drop = F]

cor_pvals<-vector()
lm_pvals<-vector()
for (i in 1:nrow(myT)) {
  cor_pvals[i]<-cor.test(unlist(myT[i, ]), calprotectin[,1], method = "kendall")$p.value
  lm_pvals[i]<-anova(lm(unlist(myT[i, ])~calprotectin[,1]))$"Pr(>F)"[1]
}
adj_corP<-p.adjust(cor_pvals, method = "BH")
adj_lmP<-p.adjust(lm_pvals, method = "BH")

pdf("Plots/BrackenSpeciesRefractoryCalprotectinCorLmANOVAPvalScatterPlots(SortedByKendall).pdf", width=24, height=12)
par(mfrow=c(2,4))
par(mar=c(5,6,4,1)+.1)
plot_order<-order(cor_pvals)
for (i in 1:nrow(myT)) {
  mainText<-paste0("Kendall pval=", signif(cor_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_corP[plot_order[i]], 3),
                   '\n', "ANOVA pval=", signif(lm_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_lmP[plot_order[i]], 3))
  plot(unlist(myT[plot_order[i], ]), calprotectin[,1], main = mainText,
       cex.lab = 2.5, cex.axis = 2, cex.main = 2, 
       xlab = rownames(myT)[plot_order[i]], ylab = "Calprotectin Level",
       col = "palevioletred", pch = 19)
}
dev.off()

pdf("Plots/BrackenSpeciesRefractoryCalprotectinCorLmANOVAPvalScatterPlots(SortedByANOVA).pdf", width=24, height=12)
par(mfrow=c(2,4))
par(mar=c(5,6,4,1)+.1)
plot_order<-order(lm_pvals)
for (i in 1:nrow(myT)) {
  mainText<-paste0("Kendall pval=", signif(cor_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_corP[plot_order[i]], 3),
                   '\n', "ANOVA pval=", signif(lm_pvals[plot_order[i]], 3), ", adj.p=", signif(adj_lmP[plot_order[i]], 3))
  plot(unlist(myT[plot_order[i], ]), calprotectin[,1], main = mainText,
       cex.lab = 2.5, cex.axis = 2, cex.main = 2, 
       xlab = rownames(myT)[plot_order[i]], ylab = "Calprotectin Level",
       col = "palevioletred", pch = 19)
}
dev.off()
