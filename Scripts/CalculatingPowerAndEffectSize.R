rm(list = ls())
#estimate effect size
metaData<-read.table("metaphlanGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("metaphlan/countTables/GVHD_species.txt", 
                row.name = 1, sep = "\t", header = T, check.names = F)
colnames(myT)<-gsub("X", "", colnames(myT))
metaData<-metaData[colnames(myT), , drop = F]
zeroP<-rowSums(myT==0)/ncol(myT)
myT<-myT[-which(zeroP>0.75),]

test<-t.test(unlist(myT["Enterococcus_faecium", ])~metaData$dx)
D<-(test$estimate[1]-test$estimate[2])/sd(myT["Enterococcus_faecium",])
pvals<-test$p.value

D<-vector()
pvals<-vector()
for (i in 1:nrow(myT)){
  test<-t.test(unlist(myT[i,])~metaData$dx)
  D[i]<-(test$estimate[1]-test$estimate[2])/sd(myT[i,])
  pvals[i]<-test$p.value
}
mean(abs(D))
length(which(abs(D)>0.32))/length(D)

t_func <- function(simNum, N, d) {
  x1 <- rnorm(N, 0, 1)
  x2 <- rnorm(N, d, 1)
  t <- t.test(x1, x2, var.equal=TRUE)  # run t-test on generated data
  stat <- t$statistic
  p <- t$p.value
  return(c(t=stat, p=p, sig=(p < .05)))
  # return a named vector with the results we want to keep
}
library(paramtest)
power_ttest <- run_test(t_func, n.iter=500, output='data.frame', N=36, d=0.312)  # simulate data
fd=p.adjust(results(power_ttest)$p,method="fdr")
length(which(fd<0.05))/500

power_ttest <- run_test(t_func, n.iter=500, output='data.frame', N=100, d=0.312)  # simulate data
fd=p.adjust(results(power_ttest)$p,method="fdr")
length(which(fd<0.05))/500

power_ttest <- run_test(t_func, n.iter=500, output='data.frame', N=180, d=0.312)  # simulate data
fd=p.adjust(results(power_ttest)$p,method="fdr")
length(which(fd<0.05))/500
