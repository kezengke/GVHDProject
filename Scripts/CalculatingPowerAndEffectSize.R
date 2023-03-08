rm(list = ls())
#estimate effect size
metaData<-read.table("metaGvN.txt", sep = "\t", header = T, row.names = 1)
myT<-read.delim("CountsTables/bracken_genus_reads.csv", 
                row.name = 1, sep = ",", header = T, check.names = F)
metaData<-metaData[colnames(myT),]
metaData<-data.frame(metaData)
rownames(metaData)<-colnames(myT)
#Normalize
n<-colSums(myT)
sumx<-sum(myT)
myT<-log10((myT/n)*(sumx/ncol(myT))+1)
myT<-data.frame(myT, check.names = F)
#Filter
lowAbundance<-which(rowMeans(myT)<2)
myT<-myT[-lowAbundance, ]

test<-t.test(unlist(myT["Saccharomyces", ])~metaData$metaData)
#effect size D 
D<-(test$estimate[1]-test$estimate[2])/sd(myT["Saccharomyces",])
D
pvals<-test$p.value

D<-vector()
pvals<-vector()
for (i in 1:nrow(myT)){
  test<-t.test(unlist(myT[i,])~metaData$metaData)
  D[i]<-(test$estimate[1]-test$estimate[2])/sd(myT[i,])
  pvals[i]<-test$p.value
}
mean(abs(D))
length(which(abs(D)>0.2))/length(D)

t_func <- function(simNum, N, d) {
  x1 <- rnorm(N, 0, 1)
  x2 <- rnorm(N, d, 1)
  t <- t.test(x1, x2, var.equal=TRUE)  # run t-test on generated data
  stat <- t$statistic
  p <- t$p.value
  return(c(t=stat, p=p, sig=(p < .05)))
  # return a named vector with the results we want to keep
}

#try out different N with previously calculated D to find the optimal sample size for desired power (for one condition)
#output is the power
library(paramtest)
power_ttest <- run_test(t_func, n.iter=500, output='data.frame', N=300, d=0.24)  # simulate data
fd=p.adjust(results(power_ttest)$p,method="fdr")
length(which(fd<0.05))/500

