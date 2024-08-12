rm(list = ls())
d <- read.table("D:/中英文论文/data/M1(1).txt", header=TRUE)
M1 <- ts(d[["M1SL"]], start=c(1959,1), frequency=12)
length(M1)  #736
plot(M1,main="M1",col="red",xlab="",ylab="")
minor.tick(nx = 10,  tick.ratio = 0.5)


M1_d<- diff(M1)
T=length(M1_d)  #735
plot(M1_d,main="M1_d",col="red",xlab="",ylab="")
minor.tick(nx = 10,  tick.ratio = 0.5)




M1_lrt<-diff(log(M1))
T=length(M1_lrt)  #735
plot(M1_lrt,main="M1_lrt",col="red",xlab="",ylab="")
minor.tick(nx = 10,  tick.ratio = 0.5)



d1 <- read.table("D:/中英文论文/data/PPI(1).txt", header=TRUE)
PPI <- ts(d1[["PPIACO"]], start=c(1959,1), frequency=12)
length(PPI)  #736
plot(PPI,main="PPI",col="red",xlab="",ylab="")
minor.tick(nx = 10,  tick.ratio = 0.5)


PPI_d<- diff(PPI)
T=length(PPI_d)  #735
plot(PPI_d,main="PPI_d",col="red",xlab="",ylab="",ylim=c(-12,8))
minor.tick(nx = 10,  tick.ratio = 0.5)




PPI_lrt<-diff(log(PPI))
T=length(PPI_lrt)  #735
plot(PPI_lrt,main="PPI_lrt",col="red",xlab="",ylab="")
minor.tick(nx = 10,  tick.ratio = 0.5)


opar <- par(mfrow=c(2,1))


d.bp <- tibble(
  M1_lrt,PPI_lrt
  ) 
knitr::kable(d.bp)

with(d.bp, {plot(M1_lrt,main="M1_lrt",col="red",xlab="",ylab="");
  minor.tick(nx = 10,  tick.ratio = 0.5);
  plot(PPI_lrt,main="PPI_lrt",col="red",xlab="",ylab="");
  minor.tick(nx = 10,  tick.ratio = 0.5)})




