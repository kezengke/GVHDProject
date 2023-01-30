pdf("Plots/BrackenShannonDiversity.pdf", width=8, height=12)
par(mar=c(5,6,4,1)+.1)
par(mfrow=c(3,2))
library("vegan")
library("coin")

#phylum
rm(list = ls())
meta<-read.delim("metaGvN.txt", sep = "\t", header=T, row.names = 1)
myT<-read.delim("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
meta<-meta[intersect(rownames(meta), colnames(myT)), , drop = F]
myT<-myT[, intersect(rownames(meta), colnames(myT)), drop = F]
group1<-rownames(meta)[meta$dx == "no_gvhd"]
group2<-rownames(meta)[meta$dx == "gvhd"]
T1<-myT[, group1, drop = F]
T2<-myT[, group2, drop = F]
H1<-diversity(T1, MARGIN = 2)
H2<-diversity(T2, MARGIN = 2)
H3<-diversity(myT, MARGIN = 2)
wilcoxP<-pvalue(wilcox_test(H3~factor(meta$dx)))
x<-list("No GVHD"=H1, "GVHD"=H2)
boxplot(x, outline = F,
        cex.lab = 2.5, cex.axis = 1.4, cex.main = 2, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        main = paste0("(Phylum)P-value=",signif(wilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

rm(list = ls())
meta<-read.delim("metaSvR.txt", sep = "\t", header=T, row.names = 1)
myT<-read.delim("CountsTables/bracken_phylum_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
meta<-meta[intersect(rownames(meta), colnames(myT)), , drop = F]
myT<-myT[, intersect(rownames(meta), colnames(myT)), drop = F]
group1<-rownames(meta)[meta$dx == "steroid_refractory_gvhd"]
group2<-rownames(meta)[meta$dx == "steroid_sensitive_gvhd"]
T1<-myT[, group1, drop = F]
T2<-myT[, group2, drop = F]
H1<-diversity(T1, MARGIN = 2)
H2<-diversity(T2, MARGIN = 2)
H3<-diversity(myT, MARGIN = 2)
wilcoxP<-pvalue(wilcox_test(H3~factor(meta$dx)))
x<-list("Steroid Refractory"=H1, "Steroid Sensitive"=H2)
boxplot(x, outline = F,
        cex.lab = 2.5, cex.axis = 1.4, cex.main = 2, 
        col=c("gold1", "coral3"),
        ylab = "Shannon Diversity",
        main = paste0("(Phylum)P-value=",signif(wilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)


#genus
rm(list = ls())
meta<-read.delim("metaGvN.txt", sep = "\t", header=T, row.names = 1)
myT<-read.delim("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
meta<-meta[intersect(rownames(meta), colnames(myT)), , drop = F]
myT<-myT[, intersect(rownames(meta), colnames(myT)), drop = F]
group1<-rownames(meta)[meta$dx == "no_gvhd"]
group2<-rownames(meta)[meta$dx == "gvhd"]
T1<-myT[, group1, drop = F]
T2<-myT[, group2, drop = F]
H1<-diversity(T1, MARGIN = 2)
H2<-diversity(T2, MARGIN = 2)
H3<-diversity(myT, MARGIN = 2)
wilcoxP<-pvalue(wilcox_test(H3~factor(meta$dx)))
x<-list("No GVHD"=H1, "GVHD"=H2)
boxplot(x, outline = F,
        cex.lab = 2.5, cex.axis = 1.4, cex.main = 2, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        main = paste0("(Genus)P-value=",signif(wilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

rm(list = ls())
meta<-read.delim("metaSvR.txt", sep = "\t", header=T, row.names = 1)
myT<-read.delim("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
meta<-meta[intersect(rownames(meta), colnames(myT)), , drop = F]
myT<-myT[, intersect(rownames(meta), colnames(myT)), drop = F]
group1<-rownames(meta)[meta$dx == "steroid_refractory_gvhd"]
group2<-rownames(meta)[meta$dx == "steroid_sensitive_gvhd"]
T1<-myT[, group1, drop = F]
T2<-myT[, group2, drop = F]
H1<-diversity(T1, MARGIN = 2)
H2<-diversity(T2, MARGIN = 2)
H3<-diversity(myT, MARGIN = 2)
wilcoxP<-pvalue(wilcox_test(H3~factor(meta$dx)))
x<-list("Steroid Refractory"=H1, "Steroid Sensitive"=H2)
boxplot(x, outline = F,
        cex.lab = 2.5, cex.axis = 1.4, cex.main = 2, 
        col=c("gold1", "coral3"),
        ylab = "Shannon Diversity",
        main = paste0("(Genus)P-value=",signif(wilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

#species
rm(list = ls())
meta<-read.delim("metaGvN.txt", sep = "\t", header=T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
meta<-meta[intersect(rownames(meta), colnames(myT)), , drop = F]
myT<-myT[, intersect(rownames(meta), colnames(myT)), drop = F]
group1<-rownames(meta)[meta$dx == "no_gvhd"]
group2<-rownames(meta)[meta$dx == "gvhd"]
T1<-myT[, group1, drop = F]
T2<-myT[, group2, drop = F]
H1<-diversity(T1, MARGIN = 2)
H2<-diversity(T2, MARGIN = 2)
H3<-diversity(myT, MARGIN = 2)
wilcoxP<-pvalue(wilcox_test(H3~factor(meta$dx)))
x<-list("No GVHD"=H1, "GVHD"=H2)
boxplot(x, outline = F,
        cex.lab = 2.5, cex.axis = 1.4, cex.main = 2, 
        col=c("olivedrab3", "cornflowerblue"),
        ylab = "Shannon Diversity",
        main = paste0("(Species)P-value=",signif(wilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

rm(list = ls())
meta<-read.delim("metaSvR.txt", sep = "\t", header=T, row.names = 1)
myT<-read.delim("CountsTables/bracken_species_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
meta<-meta[intersect(rownames(meta), colnames(myT)), , drop = F]
myT<-myT[, intersect(rownames(meta), colnames(myT)), drop = F]
group1<-rownames(meta)[meta$dx == "steroid_refractory_gvhd"]
group2<-rownames(meta)[meta$dx == "steroid_sensitive_gvhd"]
T1<-myT[, group1, drop = F]
T2<-myT[, group2, drop = F]
H1<-diversity(T1, MARGIN = 2)
H2<-diversity(T2, MARGIN = 2)
H3<-diversity(myT, MARGIN = 2)
wilcoxP<-pvalue(wilcox_test(H3~factor(meta$dx)))
x<-list("Steroid Refractory"=H1, "Steroid Sensitive"=H2)
boxplot(x, outline = F,
        cex.lab = 2.5, cex.axis = 1.4, cex.main = 2, 
        col=c("gold1", "coral3"),
        ylab = "Shannon Diversity",
        main = paste0("(Species)P-value=",signif(wilcoxP, 4)))
stripchart(x, method="jitter", 
           vertical=T, pch=19, add=T)

dev.off()
