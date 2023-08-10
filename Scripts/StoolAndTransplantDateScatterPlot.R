rm(list= ls())

myT<-read.csv("DatesOfStoolAndTP.csv", header = T)
myT<-na.omit(myT)

myT$Stool.Date<-as.Date(myT$Stool.Date, format="%Y-%m-%d")
myT$Transplant.date<-as.Date(myT$Transplant.date, format="%Y-%m-%d")

pdf("Plots/StoolTransplantDatePlot.pdf", width = 13, height = 13)
par(mar=c(6,6,4,1)+.1)
plot(1, 1, xlim = range(c(myT$Stool.Date, myT$Transplant.date)), 
     ylim = c(1, nrow(myT)), 
     type="n", xaxt="n", xlab="", ylab="Sample Code", main = "Stool Collection and Transplant Date",
     yaxt="n", las=1)
points(myT$Stool.Date, 1:nrow(myT), col="blue", pch=19)
points(myT$Transplant.date, 1:nrow(myT), col="red", pch=19)

axis(2, at=1:nrow(myT), labels=myT$Code, las=1)
axis.Date(1, at=seq(min(c(myT$Stool.Date, myT$Transplant.date)), max(c(myT$Stool.Date, myT$Transplant.date)), by="months"), 
          format="%m-%Y", las = 3)
mtext("Date", side=1, line=5)
legend("bottomleft", legend=c("stool date", "transplant date"), fill=c("blue", "red"), bty="n")


dev.off()
