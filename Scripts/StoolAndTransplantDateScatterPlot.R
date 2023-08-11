rm(list= ls())

myT<-read.csv("DatesOfStoolAndTP.csv", header = T)
myT<-na.omit(myT)

myT$Stool.Date<-as.Date(myT$Stool.Date, format="%Y-%m-%d")
myT$Transplant.date<-as.Date(myT$Transplant.date, format="%Y-%m-%d")

myT$days<-as.numeric(myT$Stool.Date - myT$Transplant.date)

pdf("Plots/StoolTransplantDatePlot.pdf", width = 15, height = 15)
par(mar=c(6,6,4,1)+.1)
plot(1, 1, xlim = range(myT$days), 
     ylim = c(1, nrow(myT)), 
     type="n", xaxt="n", xlab="", ylab="Sample Code", main = "Stool Collection Date relative to Transplant Date",
     yaxt="n", las=1)
points(myT$days, 1:nrow(myT), col="cornflowerblue", pch=19)
axis(2, at=1:nrow(myT), labels=myT$Code, las=1)
axis(1, at=seq(min(myT$days), max(myT$days), by=50), las = 3)
mtext("Days", side=1, line=3)

dev.off()
