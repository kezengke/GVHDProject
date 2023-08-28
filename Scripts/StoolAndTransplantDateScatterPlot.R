rm(list= ls())
library(stringr)
meta1<-read.delim("metaGvN_index.txt", header = T, row.names = 1)
meta2<-read.delim("metaSvR_index.txt", header = T, row.names = 1)

meta1[rownames(meta2), 1]<-meta2$dx
rownames(meta1)<-sapply(str_split(rownames(meta1), "_", n = 2), `[`, 1)
rownames(meta1)<-as.numeric(rownames(meta1))

myT<-read.csv("TransplantStoolDate.csv", header = T)
myT<-na.omit(myT)

myT$stool_date<-as.Date(myT$stool_date, format="%Y-%m-%d")
myT$bmt_date<-as.Date(myT$bmt_date, format="%Y-%m-%d")

myT$days<-as.numeric(myT$stool_date - myT$bmt_date)

myT$condition<-meta1[as.character(myT$code), 1]
myT$index<-meta1[as.character(myT$code), 2]
myT<-myT[myT$days >= 0, , drop = F]
myT<-myT[myT$days < 400, , drop = F]

myT<-myT[order(myT$condition), , drop = F]

pdf("Plots/StoolTransplantDatePlot.pdf", width = 12, height = 15)
par(mar=c(6,7,4,1)+.1)
plot(1, 1, xlim = c(0,220), 
     ylim = c(1, nrow(myT)), 
     type="n", xaxt="n", xlab="", ylab="", main = "Stool Collection Date relative to Transplant Date",
     yaxt="n", las=1)
color_mapping<-c("no_gvhd" = "cornflowerblue", 
                 "steroid_refractory_gvhd" = "coral3",
                 "steroid_sensitive_gvhd" = "tan2")

points(myT$days, 1:nrow(myT), col=color_mapping[myT$condition], pch=19)
axis(2, at=1:nrow(myT), labels=myT$index, las=1)
axis(1, at=seq(0, 220, by=5), las = 3)
mtext("Days", side=1, line=3)
mtext("Sample Type", side = 2, line = 5)
legend("topright", legend = names(color_mapping), fill = color_mapping, border = "white", bty="n")


dev.off()
