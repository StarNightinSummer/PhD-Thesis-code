library(rstudioapi)
library(Rsolnp)
library(numDeriv)
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("final functions.R")
source("p_match.R")


system("R CMD SHLIB garch_mu.c")
dyn.load("garch_mu.dll")
#dyn.unload("garch_mu.dll")

system("R CMD SHLIB garch_ar.c")
dyn.load("garch_ar.dll")
#dyn.unload("garch_ar.dll")
#diffh = h-ht
#diffe = heps[[2]]-epsit



N=2
include.AR = T
type="diagonal"
startval = NULL



####################################   input data    ############################################

HSI<-read.csv("data.csv", header = T)[,2]
SSE<-read.csv("data.csv", header = T)[,3]

plot(HSI,type="l" )
plot(SSE,type="l")

z<-cbind(HSI,SSE)
MTSplot(z)



#compute log return
HSI_lrt<-diff(log(HSI))*100
SSE_lrt<-diff(log(SSE))*100

plot(HSI_lrt,type="l" )
plot(SSE_lrt,type="l")



z_lrt<-cbind(HSI_lrt,SSE_lrt)  # create a vector series
MTSplot(z_lrt)
dim(z_lrt)
# 3333 2
MTSplot(z_lrt)





#######################################    estimation         ########################################
HSI<-HSI[313:3312]
SSE<-SSE[313:3312]
z<-cbind(HSI,SSE)
MTSplot(z)


HSI_lrt<-HSI_lrt[313:3312]
SSE_lrt<-SSE_lrt[313:3312]



z_lrt<-z_lrt[313:3312,]
MTSplot(z_lrt)

nobs=dim(z_lrt)[1]
weight = weights(data=z_lrt)
start.time = Sys.time()
parEst1 = CCCfit(data=z_lrt, w=weight, include.AR = T, type = type, solver="gosolnp", tr=1)
print(as.vector(parEst1$pars))
end.time = Sys.time()
print(difftime(end.time, start.time, units="mins"))


para.mat = p.match(parEst1$pars, include.AR, type, N)
print(para.mat$mu)
print(para.mat$phi)
print(para.mat$W)
print(para.mat$A)
print(para.mat$B)
print(para.mat$Gamma)


start.time = Sys.time()
AV = AVSWMC(data=z_lrt, w=weight, type=type, par=parEst1$pars, include.AR=T, lag = 6)
parEst2 = Localfit(data=z_lrt, par=parEst1$pars, par1=parEst1$pars, type=type, include.AR=T)
print(as.vector(parEst2$fit))
para.mat = p.match(parEst2$fit, include.AR, type, N)
print(para.mat$mu)
print(para.mat$phi)
print(para.mat$W)
print(para.mat$A)
print(para.mat$B)
print(para.mat$Gamma)




AV1   =  LocalVAMC(data=z_lrt, par=parEst2$fit, type=type, include.AR=include.AR, lag = 6)

end.time = Sys.time()
print(difftime(end.time, start.time, units="mins"))


print(AV$vari)
print(AV$QMM)
print(AV1$vari)
print(AV1$QMM)


pv=1-pchisq(AV1$QMM,6)
pv







C<- rep(1,13)
tao<- diag(C)
r<- rep(0,13)

W<- c(T*t(tao%*%parEst2$fit-r)%*%solve(tao%*%AV1$vari1%*%t(tao))%*%(tao%*%parEst2$fit-r))

pv=1-pchisq(W,13)
pv   ### 1.431285e-09





tao1<-matrix(numeric(52),nrow=4,ncol=13)
tao1[1,1]<-1
tao1[2,2]<-1
tao1[3,4]<-1
tao1[4,6]<-1
tao1

r1<- rep(0,4)

W1<- c(T*t(tao1%*%parEst2$fit-r1)%*%solve(tao1%*%AV1$vari1%*%t(tao1))%*%(tao1%*%parEst2$fit-r1))

pv=1-pchisq(W1,4)
pv  ##    


