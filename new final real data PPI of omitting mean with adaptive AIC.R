#####################   first difference of PPI     #########################
install.packages("tsibble")
library(tsibble)



rm(list = ls())
d <- read.table("D:/中英文论文/data/PPI(1).txt", header=TRUE)
PPI <- ts(d[["PPIACO"]], start=c(1959,1), frequency=12)
length(PPI)  #736
plot(PPI,main="PPI",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)


PPI_d<- diff(PPI)
T=length(PPI_d)  #735
plot(PPI_d,main="PPI_d",col="red",xlab="",ylab="",ylim=c(-12,8))
# minor.tick(nx = 10,  tick.ratio = 0.5)




PPI_lrt<-diff(log(PPI))
PPI_lrt<-PPI_lrt-mean(PPI_lrt)
T=length(PPI_lrt)  #735
plot(PPI_lrt,main="PPI_lrt",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)




forecast::Acf(PPI_lrt,xlab="",ylab="",col="red",main="ACF of PPI_lrt")
forecast::Pacf(PPI_lrt,xlab="",ylab="",col="red",main="PACF of PPI_lrt")
forecast::Acf(PPI_lrt^2,xlab="",ylab="",col="red",main="ACF of PPI_lrt^2")
forecast::Pacf(PPI_lrt^2,xlab="",ylab="",col="red",main="PACF of PPI_lrt^2")



p=10
q=10
table.aic<-matrix(numeric(p*q), nrow=p,ncol=q)
for (i in 1:p) {
  for (j in 1:q) {
    table.aic[i,j]<- arima(PPI_lrt,order = c(i-1,0,j-1))$aic 
    j=j+1
  }
  i=i+1
}
table.aic
which.min(table.aic)   # 90 




ar(PPI_lrt,method = "ols")



############  ARMA(14,0) ############################

#Test mean being zero.
t.test(PPI_lrt)   


m1<-arima(PPI_lrt, order = c(14,0,0),include.mean = FALSE)
m1

coef1<-m1$coef
coef1   ##  0.3195  0.1484  0.0537  -0.0233  -0.0560  0.0302  0.0289  0.0045  -0.0427  0.0306  0.1431  0.0069  -0.1228  0.1053  aic = -5015.05

mean(PPI_lrt)
# -6.925618e-20




########################  calculate residuals  ###################################

epsilonhat=numeric(T)


epsilonhat[1]<- PPI_lrt[1] 

epsilonhat[2]<- PPI_lrt[2] - coef1[1]*PPI_lrt[1]

epsilonhat[3]<- PPI_lrt[3] - coef1[1]*PPI_lrt[2] - coef1[2]*PPI_lrt[1]

epsilonhat[4]<- PPI_lrt[4] - coef1[1]*PPI_lrt[3] - coef1[2]*PPI_lrt[2] - coef1[3]*PPI_lrt[1] 

epsilonhat[5]<- PPI_lrt[5] - coef1[1]*PPI_lrt[4] - coef1[2]*PPI_lrt[3] - coef1[3]*PPI_lrt[2] - coef1[4]*PPI_lrt[1]

epsilonhat[6]<- PPI_lrt[6] - coef1[1]*PPI_lrt[5] - coef1[2]*PPI_lrt[4] - coef1[3]*PPI_lrt[3] - coef1[4]*PPI_lrt[2] - coef1[5]*PPI_lrt[1]

epsilonhat[7]<- PPI_lrt[7] - coef1[1]*PPI_lrt[6] - coef1[2]*PPI_lrt[5] - coef1[3]*PPI_lrt[4] - coef1[4]*PPI_lrt[3] - coef1[5]*PPI_lrt[2] - coef1[6]*PPI_lrt[1]

epsilonhat[8]<- PPI_lrt[8] - coef1[1]*PPI_lrt[7] - coef1[2]*PPI_lrt[6] - coef1[3]*PPI_lrt[5] - coef1[4]*PPI_lrt[4] - coef1[5]*PPI_lrt[3] - coef1[6]*PPI_lrt[2] - coef1[7]*PPI_lrt[1]

epsilonhat[9]<- PPI_lrt[9] - coef1[1]*PPI_lrt[8] - coef1[2]*PPI_lrt[7] - coef1[3]*PPI_lrt[6] - coef1[4]*PPI_lrt[5] - coef1[5]*PPI_lrt[4] - coef1[6]*PPI_lrt[3] - coef1[7]*PPI_lrt[2] - coef1[8]*PPI_lrt[1]

epsilonhat[10]<- PPI_lrt[10] - coef1[1]*PPI_lrt[9] - coef1[2]*PPI_lrt[8] - coef1[3]*PPI_lrt[7] - coef1[4]*PPI_lrt[6] - coef1[5]*PPI_lrt[5] - coef1[6]*PPI_lrt[4] - coef1[7]*PPI_lrt[3] - coef1[8]*PPI_lrt[2] - coef1[9]*PPI_lrt[1]

epsilonhat[11]<- PPI_lrt[11] - coef1[1]*PPI_lrt[10] - coef1[2]*PPI_lrt[9] - coef1[3]*PPI_lrt[8] - coef1[4]*PPI_lrt[7] - coef1[5]*PPI_lrt[6] - coef1[6]*PPI_lrt[5] - coef1[7]*PPI_lrt[4] - coef1[8]*PPI_lrt[3] - coef1[9]*PPI_lrt[2] - coef1[10]*PPI_lrt[1]

epsilonhat[12]<- PPI_lrt[12] - coef1[1]*PPI_lrt[11] - coef1[2]*PPI_lrt[10] - coef1[3]*PPI_lrt[9] - coef1[4]*PPI_lrt[8] - coef1[5]*PPI_lrt[7] - coef1[6]*PPI_lrt[6] - coef1[7]*PPI_lrt[5] - coef1[8]*PPI_lrt[4] - coef1[9]*PPI_lrt[3] - coef1[10]*PPI_lrt[2] - coef1[11]*PPI_lrt[1]

epsilonhat[13]<- PPI_lrt[13] - coef1[1]*PPI_lrt[12] - coef1[2]*PPI_lrt[11] - coef1[3]*PPI_lrt[10] - coef1[4]*PPI_lrt[9] - coef1[5]*PPI_lrt[8] - coef1[6]*PPI_lrt[7] - coef1[7]*PPI_lrt[6] - coef1[8]*PPI_lrt[5] - coef1[9]*PPI_lrt[4] - coef1[10]*PPI_lrt[3] - coef1[11]*PPI_lrt[2] - coef1[12]*PPI_lrt[1]

epsilonhat[14]<- PPI_lrt[14] - coef1[1]*PPI_lrt[13] - coef1[2]*PPI_lrt[12] - coef1[3]*PPI_lrt[11] - coef1[4]*PPI_lrt[10] - coef1[5]*PPI_lrt[9] - coef1[6]*PPI_lrt[8] - coef1[7]*PPI_lrt[7] - coef1[8]*PPI_lrt[6] - coef1[9]*PPI_lrt[5] - coef1[10]*PPI_lrt[4] - coef1[11]*PPI_lrt[3] - coef1[12]*PPI_lrt[2] - coef1[13]*PPI_lrt[1]


for (i in 15:T) {
  epsilonhat[i]<- PPI_lrt[i] - coef1[1]*PPI_lrt[i-1] - coef1[2]*PPI_lrt[i-2] - coef1[3]*PPI_lrt[i-3] - coef1[4]*PPI_lrt[i-4] - coef1[5]*PPI_lrt[i-5] - coef1[6]*PPI_lrt[i-6] - coef1[7]*PPI_lrt[i-7] - coef1[8]*PPI_lrt[i-8] - coef1[9]*PPI_lrt[i-9] - coef1[10]*PPI_lrt[i-10] - coef1[11]*PPI_lrt[i-11] - coef1[12]*PPI_lrt[i-12] - coef1[13]*PPI_lrt[i-13] - coef1[14]*PPI_lrt[i-14]  
  
}


##  adaptive AIC and BIC


options(np.messages = FALSE)
X<-seq(1,T,length=T)
X<-X/T
Y<-epsilonhat^2
bw0 <- np::npregbw(formula = Y ~ X)
kre0 <- np::npreg(bw0)
# h[j] <- bw0$bw
plot(X,Y,pch=20)
lines(kre0$eval$X,kre0$mean, lwd=4, col=2)

LALS = mean(log(kre0$mean)+(Y/kre0$mean))
AIC  = LALS + 2*(14+0)/T
BIC  = LALS + (14+0)*log(T)/T
AIC
#  -9.590053

BIC
# -9.502437







########################    model checking M=16  (S(2)) ##########################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)
rho10<- sum(epsilonhat[1:(T-10)]*epsilonhat[11:T])/sum(epsilonhat^2)
rho11<- sum(epsilonhat[1:(T-11)]*epsilonhat[12:T])/sum(epsilonhat^2)
rho12<- sum(epsilonhat[1:(T-12)]*epsilonhat[13:T])/sum(epsilonhat^2)
rho13<- sum(epsilonhat[1:(T-13)]*epsilonhat[14:T])/sum(epsilonhat^2)
rho14<- sum(epsilonhat[1:(T-14)]*epsilonhat[15:T])/sum(epsilonhat^2)
rho15<- sum(epsilonhat[1:(T-15)]*epsilonhat[16:T])/sum(epsilonhat^2)
rho16<- sum(epsilonhat[1:(T-16)]*epsilonhat[17:T])/sum(epsilonhat^2)



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12,rho13,rho14,rho15,rho16),nrow = 1,ncol = 16)

C<-rep(1,16)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   0.4357927


pv=1-pchisq(P,2)
pv  #  0.8042088



############   M=20 (S(6)) #############################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)
rho10<- sum(epsilonhat[1:(T-10)]*epsilonhat[11:T])/sum(epsilonhat^2)
rho11<- sum(epsilonhat[1:(T-11)]*epsilonhat[12:T])/sum(epsilonhat^2)
rho12<- sum(epsilonhat[1:(T-12)]*epsilonhat[13:T])/sum(epsilonhat^2)
rho13<- sum(epsilonhat[1:(T-13)]*epsilonhat[14:T])/sum(epsilonhat^2)
rho14<- sum(epsilonhat[1:(T-14)]*epsilonhat[15:T])/sum(epsilonhat^2)
rho15<- sum(epsilonhat[1:(T-15)]*epsilonhat[16:T])/sum(epsilonhat^2)
rho16<- sum(epsilonhat[1:(T-16)]*epsilonhat[17:T])/sum(epsilonhat^2)
rho17<- sum(epsilonhat[1:(T-17)]*epsilonhat[18:T])/sum(epsilonhat^2)
rho18<- sum(epsilonhat[1:(T-18)]*epsilonhat[19:T])/sum(epsilonhat^2)
rho19<- sum(epsilonhat[1:(T-19)]*epsilonhat[20:T])/sum(epsilonhat^2)
rho20<- sum(epsilonhat[1:(T-20)]*epsilonhat[21:T])/sum(epsilonhat^2)



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12,rho13,rho14,rho15,rho16,rho17,rho18,rho19,rho20),nrow = 1,ncol = 20)

C<-rep(1,20)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   1.141566


pv=1-pchisq(P,6)
pv   ###  0.9796725



############   M=24  (S(10)) #############################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)
rho10<- sum(epsilonhat[1:(T-10)]*epsilonhat[11:T])/sum(epsilonhat^2)
rho11<- sum(epsilonhat[1:(T-11)]*epsilonhat[12:T])/sum(epsilonhat^2)
rho12<- sum(epsilonhat[1:(T-12)]*epsilonhat[13:T])/sum(epsilonhat^2)
rho13<- sum(epsilonhat[1:(T-13)]*epsilonhat[14:T])/sum(epsilonhat^2)
rho14<- sum(epsilonhat[1:(T-14)]*epsilonhat[15:T])/sum(epsilonhat^2)
rho15<- sum(epsilonhat[1:(T-15)]*epsilonhat[16:T])/sum(epsilonhat^2)
rho16<- sum(epsilonhat[1:(T-16)]*epsilonhat[17:T])/sum(epsilonhat^2)
rho17<- sum(epsilonhat[1:(T-17)]*epsilonhat[18:T])/sum(epsilonhat^2)
rho18<- sum(epsilonhat[1:(T-18)]*epsilonhat[19:T])/sum(epsilonhat^2)
rho19<- sum(epsilonhat[1:(T-19)]*epsilonhat[20:T])/sum(epsilonhat^2)
rho20<- sum(epsilonhat[1:(T-20)]*epsilonhat[21:T])/sum(epsilonhat^2)
rho21<- sum(epsilonhat[1:(T-21)]*epsilonhat[22:T])/sum(epsilonhat^2)
rho22<- sum(epsilonhat[1:(T-22)]*epsilonhat[23:T])/sum(epsilonhat^2)
rho23<- sum(epsilonhat[1:(T-23)]*epsilonhat[24:T])/sum(epsilonhat^2)
rho24<- sum(epsilonhat[1:(T-24)]*epsilonhat[25:T])/sum(epsilonhat^2)


rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12,rho13,rho14,rho15,rho16,rho17,rho18,rho19,rho20,rho21,rho22,rho23,rho24),nrow = 1,ncol = 24)

C<-rep(1,24)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #  2.362917


pv=1-pchisq(P,10)
pv   ###   0.9927264




###########################     Wald test  ########################################

##########  derivatives  ######################

derivativeofphi<- matrix(numeric(14*T),nrow = 14, ncol = T)


derivativeofphi1<- numeric(T)

for (i in 2:T) {
  derivativeofphi1[i]<- - PPI_lrt[i-1] 
}


derivativeofphi2<- numeric(T)

for (i in 3:T) {
  derivativeofphi2[i]<- - PPI_lrt[i-2]
}


derivativeofphi3<- numeric(T)

for (i in 4:T) {
  derivativeofphi3[i]<- - PPI_lrt[i-3]  
}



derivativeofphi4<- numeric(T)

for (i in 5:T) {
  derivativeofphi4[i]<- - PPI_lrt[i-4]  
}


derivativeofphi5<- numeric(T)

for (i in 6:T) {
  derivativeofphi5[i]<- - PPI_lrt[i-5]  
}



derivativeofphi6<- numeric(T)

for (i in 7:T) {
  derivativeofphi6[i]<- - PPI_lrt[i-6]  
}



derivativeofphi7<- numeric(T)

for (i in 8:T) {
  derivativeofphi7[i]<- - PPI_lrt[i-7]  
}


derivativeofphi8<- numeric(T)

for (i in 9:T) {
  derivativeofphi8[i]<- - PPI_lrt[i-8]  
}


derivativeofphi9<- numeric(T)

for (i in 10:T) {
  derivativeofphi9[i]<- - PPI_lrt[i-9]  
}


derivativeofphi10<- numeric(T)

for (i in 11:T) {
  derivativeofphi10[i]<- - PPI_lrt[i-10]  
}


derivativeofphi11<- numeric(T)

for (i in 12:T) {
  derivativeofphi11[i]<- - PPI_lrt[i-11]  
}


derivativeofphi12<- numeric(T)

for (i in 13:T) {
  derivativeofphi12[i]<- - PPI_lrt[i-12]  
}


derivativeofphi13<- numeric(T)

for (i in 14:T) {
  derivativeofphi13[i]<- - PPI_lrt[i-13]
}


derivativeofphi14<- numeric(T)

for (i in 15:T) {
  derivativeofphi14[i]<- - PPI_lrt[i-14]
}





for (i in 1:T) {
  derivativeofphi[,i]<- c(derivativeofphi1[i],derivativeofphi2[i],derivativeofphi3[i],derivativeofphi4[i],derivativeofphi5[i],derivativeofphi6[i],derivativeofphi7[i],derivativeofphi8[i],derivativeofphi9[i],derivativeofphi10[i],derivativeofphi11[i],derivativeofphi12[i],derivativeofphi13[i],derivativeofphi14[i])
  
}

####################################################

omega<- matrix(numeric(196*T),nrow=196,ncol=T)
sigma<- matrix(numeric(196*T),nrow=196,ncol=T)


for (j in 1:T) {
  
  omega[,j]<- as.vector(4*epsilonhat[j]^2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  sigma[,j]<- as.vector(2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  
}

Omega<-matrix(apply(omega,1,mean),nrow=14,ncol=14)
Sigma<-matrix(apply(sigma,1,mean),nrow=14,ncol=14)

Gamma<-solve(Sigma)%*%Omega%*%solve(Sigma)
var<- diag(Gamma)/T
round(sqrt(var),4)
##   0.0814 0.0758 0.0654 0.0586 0.0650 0.0487 0.0470 0.0483 0.0516 0.0512 0.0649 0.0636 0.0471 0.0497



##########  test all  are significantly not equal to 0 ################################
C<- rep(1,14)
tao<- diag(C)
r<- rep(0,14)

W<- c(T*t(tao%*%coef1-r)%*%solve(tao%*%Gamma%*%t(tao))%*%(tao%*%coef1-r))

pv=1-pchisq(W,14)
pv   ###  8.68121e-10


#########  test theta2 and theta4  equal to 0 ##################

tao1<-matrix(numeric(126),nrow=9,ncol=14)
tao1[1,3]<-1
tao1[2,4]<-1
tao1[3,5]<-1
tao1[4,6]<-1
tao1[5,7]<-1
tao1[6,8]<-1
tao1[7,9]<-1
tao1[8,10]<-1
tao1[9,12]<-1
tao1

r1<- rep(0,9)

W1<- c(T*t(tao1%*%coef1-r1)%*%solve(tao1%*%Gamma%*%t(tao1))%*%(tao1%*%coef1-r1))

pv=1-pchisq(W1,9)
pv  ##    0.9498942



############  fix parameters ############################

#Test mean being zero.
t.test(PPI_lrt)   


m1<-arima(PPI_lrt, order = c(14,0,0),include.mean = FALSE, fixed = c(NA,NA,0,0,0,0,0,0,0,0,NA,0,NA,NA),transform.pars=FALSE)
m1

coef1<-m1$coef
coef1   ##  0.3236  0.1519    0    0    0    0    0    0    0     0  0.1465     0  -0.1208  0.1159  aic = -5026.81


mean(PPI_lrt)
# -5.04827e-19




########################  calculate residuals  ###################################

epsilonhat=numeric(T)


epsilonhat[1]<- PPI_lrt[1] 

epsilonhat[2]<- PPI_lrt[2] - coef1[1]*PPI_lrt[1]

epsilonhat[3]<- PPI_lrt[3] - coef1[1]*PPI_lrt[2] - coef1[2]*PPI_lrt[1]

epsilonhat[4]<- PPI_lrt[4] - coef1[1]*PPI_lrt[3] - coef1[2]*PPI_lrt[2] 

epsilonhat[5]<- PPI_lrt[5] - coef1[1]*PPI_lrt[4] - coef1[2]*PPI_lrt[3] 

epsilonhat[6]<- PPI_lrt[6] - coef1[1]*PPI_lrt[5] - coef1[2]*PPI_lrt[4] 

epsilonhat[7]<- PPI_lrt[7] - coef1[1]*PPI_lrt[6] - coef1[2]*PPI_lrt[5] 

epsilonhat[8]<- PPI_lrt[8] - coef1[1]*PPI_lrt[7] - coef1[2]*PPI_lrt[6] 

epsilonhat[9]<- PPI_lrt[9] - coef1[1]*PPI_lrt[8] - coef1[2]*PPI_lrt[7] 

epsilonhat[10]<- PPI_lrt[10] - coef1[1]*PPI_lrt[9] - coef1[2]*PPI_lrt[8]

epsilonhat[11]<- PPI_lrt[11] - coef1[1]*PPI_lrt[10] - coef1[2]*PPI_lrt[9]

epsilonhat[12]<- PPI_lrt[12] - coef1[1]*PPI_lrt[11] - coef1[2]*PPI_lrt[10] - coef1[11]*PPI_lrt[1]

epsilonhat[13]<- PPI_lrt[13] - coef1[1]*PPI_lrt[12] - coef1[2]*PPI_lrt[11] - coef1[11]*PPI_lrt[2]

epsilonhat[14]<- PPI_lrt[14] - coef1[1]*PPI_lrt[13] - coef1[2]*PPI_lrt[12] - coef1[11]*PPI_lrt[3] - coef1[13]*PPI_lrt[1]


for (i in 15:T) {
  epsilonhat[i]<- PPI_lrt[i] - coef1[1]*PPI_lrt[i-1] - coef1[2]*PPI_lrt[i-2] - coef1[11]*PPI_lrt[i-11] - coef1[13]*PPI_lrt[i-13] - coef1[14]*PPI_lrt[i-14]  
  
}



##  adaptive AIC and BIC


options(np.messages = FALSE)
X<-seq(1,T,length=T)
X<-X/T
Y<-epsilonhat^2
bw0 <- np::npregbw(formula = Y ~ X)
kre0 <- np::npreg(bw0)
# h[j] <- bw0$bw
plot(X,Y,pch=20)
lines(kre0$eval$X,kre0$mean, lwd=4, col=2)

LALS = mean(log(kre0$mean)+(Y/kre0$mean))
AIC  = LALS + 2*(5+0)/T
BIC  = LALS + (5+0)*log(T)/T
AIC
#-9.607896

BIC
# -9.576604





########################    model checking M=7  (S(2)) ##########################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7),nrow = 1,ncol = 7)

C<-rep(1,7)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   1.136627


pv=1-pchisq(P,2)
pv  #  0.5664801



############   M=11 (S(6)) #############################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)
rho10<- sum(epsilonhat[1:(T-10)]*epsilonhat[11:T])/sum(epsilonhat^2)
rho11<- sum(epsilonhat[1:(T-11)]*epsilonhat[12:T])/sum(epsilonhat^2)



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11),nrow = 1,ncol = 11)

C<-rep(1,11)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   1.385819


pv=1-pchisq(P,6)
pv   ###  0.9667154



############   M=15  (S(10)) #############################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)
rho10<- sum(epsilonhat[1:(T-10)]*epsilonhat[11:T])/sum(epsilonhat^2)
rho11<- sum(epsilonhat[1:(T-11)]*epsilonhat[12:T])/sum(epsilonhat^2)
rho12<- sum(epsilonhat[1:(T-12)]*epsilonhat[13:T])/sum(epsilonhat^2)
rho13<- sum(epsilonhat[1:(T-13)]*epsilonhat[14:T])/sum(epsilonhat^2)
rho14<- sum(epsilonhat[1:(T-14)]*epsilonhat[15:T])/sum(epsilonhat^2)
rho15<- sum(epsilonhat[1:(T-15)]*epsilonhat[16:T])/sum(epsilonhat^2)


rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12,rho13,rho14,rho15),nrow = 1,ncol = 15)

C<-rep(1,15)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #  1.520445


pv=1-pchisq(P,10)
pv   ###   0.9988702




###########################     Wald test  ########################################

##########  derivatives  ######################

derivativeofphi<- matrix(numeric(5*T),nrow = 5, ncol = T)


derivativeofphi1<- numeric(T)

for (i in 2:T) {
  derivativeofphi1[i]<- - PPI_lrt[i-1] 
}


derivativeofphi2<- numeric(T)

for (i in 3:T) {
  derivativeofphi2[i]<- - PPI_lrt[i-2]
}



derivativeofphi11<- numeric(T)

for (i in 12:T) {
  derivativeofphi11[i]<- - PPI_lrt[i-11]  
}



derivativeofphi13<- numeric(T)

for (i in 14:T) {
  derivativeofphi13[i]<- - PPI_lrt[i-13]
}


derivativeofphi14<- numeric(T)

for (i in 15:T) {
  derivativeofphi14[i]<- - PPI_lrt[i-14]
}





for (i in 1:T) {
  derivativeofphi[,i]<- c(derivativeofphi1[i],derivativeofphi2[i],derivativeofphi11[i],derivativeofphi13[i],derivativeofphi14[i])
  
}


####################################################

omega<- matrix(numeric(25*T),nrow=25,ncol=T)
sigma<- matrix(numeric(25*T),nrow=25,ncol=T)


for (j in 1:T) {
  
  omega[,j]<- as.vector(4*epsilonhat[j]^2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  sigma[,j]<- as.vector(2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  
}

Omega<-matrix(apply(omega,1,mean),nrow=5,ncol=5)
Sigma<-matrix(apply(sigma,1,mean),nrow=5,ncol=5)

Gamma<-solve(Sigma)%*%Omega%*%solve(Sigma)
var<- diag(Gamma)/T
round(sqrt(var),4)
##   0.0807 0.0687 0.0531 0.0426 0.0469



##########  test all  are significantly not equal to 0 ################################
C<- rep(1,5)
tao<- diag(C)
r<- rep(0,5)

W<- c(T*t(tao%*%coef1[c(1,2,11,13,14)]-r)%*%solve(tao%*%Gamma%*%t(tao))%*%(tao%*%coef1[c(1,2,11,13,14)]-r))

pv=1-pchisq(W,5)
pv   ### 9.057755e-12








#######################   log return of PPI  #####################################

rm(list = ls())
d <- read.table("D:/中英文论文/data/PPI(1).txt", header=TRUE)
PPI <- ts(d[["PPIACO"]], start=c(1959,1), frequency=12)
length(PPI)  #736
plot(PPI,main="PPI",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)


PPI_d<- diff(PPI)
T=length(PPI_d)  #735
plot(PPI_d,main="PPI_d",col="red",xlab="",ylab="",ylim=c(-12,8))
# minor.tick(nx = 10,  tick.ratio = 0.5)




PPI_lrt<-diff(log(PPI))
PPI_lrt<-PPI_lrt-mean(PPI_lrt)
T=length(PPI_lrt)  #735
plot(PPI_lrt,main="PPI_lrt",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)




forecast::Acf(PPI_lrt,xlab="",ylab="",col="red",main="ACF of PPI_lrt")
forecast::Pacf(PPI_lrt,xlab="",ylab="",col="red",main="PACF of PPI_lrt")
forecast::Acf(PPI_lrt^2,xlab="",ylab="",col="red",main="ACF of PPI_lrt^2")
forecast::Pacf(PPI_lrt^2,xlab="",ylab="",col="red",main="PACF of PPI_lrt^2")



p=10
q=10
table.aic<-matrix(numeric(p*q), nrow=p,ncol=q)
for (i in 1:p) {
  for (j in 1:q) {
    table.aic[i,j]<- arima(PPI_lrt,order = c(i-1,0,j-1))$aic 
    j=j+1
  }
  i=i+1
}
table.aic
which.min(table.aic)   # 90 



#Test mean being zero.
t.test(PPI_lrt)   



#######################   ARMA(3,2)  model     #############################

m1<-arima(PPI_lrt, order = c(3,0,2),include.mean = FALSE)
m1

coef1<-m1$coef
coef1   ##  0.6807  -0.9485  0.3353  -0.3356  0.9425    aic = -5009.59


mean(PPI_lrt)
# -6.925618e-20



########################  calculate residuals  ###################################

epsilonhat=numeric(T)


epsilonhat[1]<- PPI_lrt[1] 

epsilonhat[2]<- PPI_lrt[2] - coef1[1]*PPI_lrt[1] - coef1[4]*epsilonhat[1] 

epsilonhat[3]<- PPI_lrt[3] - coef1[1]*PPI_lrt[2] - coef1[2]*PPI_lrt[1] - coef1[4]*epsilonhat[2] - coef1[5]*epsilonhat[1] 



for (i in 4:T) {
  epsilonhat[i]<- PPI_lrt[i] - coef1[1]*PPI_lrt[i-1] - coef1[2]*PPI_lrt[i-2] - coef1[3]*PPI_lrt[i-3] - coef1[4]*epsilonhat[i-1] - coef1[5]*epsilonhat[i-2]
  
}



##  adaptive AIC and BIC


options(np.messages = FALSE)
X<-seq(1,T,length=T)
X<-X/T
Y<-epsilonhat^2
bw0 <- np::npregbw(formula = Y ~ X)
kre0 <- np::npreg(bw0)
# h[j] <- bw0$bw
plot(X,Y,pch=20)
lines(kre0$eval$X,kre0$mean, lwd=4, col=2)

LALS = mean(log(kre0$mean)+(Y/kre0$mean))
AIC  = LALS + 2*(3+2)/T
BIC  = LALS + (3+2)*log(T)/T
AIC
#   -9.545492

BIC
# -9.5142




########################    model checking M=7  (S(2)) ##########################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)


rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7),nrow = 1,ncol = 7 )

C<-rep(1,7)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #  1.409575


pv=1-pchisq(P,2)
pv  #  0.4942137



############   M=11 (S(6)) #############################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)
rho10<- sum(epsilonhat[1:(T-10)]*epsilonhat[11:T])/sum(epsilonhat^2)
rho11<- sum(epsilonhat[1:(T-11)]*epsilonhat[12:T])/sum(epsilonhat^2)



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11),nrow = 1,ncol = 11)

C<-rep(1,11)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   3.506677


pv=1-pchisq(P,6)
pv   ###  0.7430812



############   M=15  (S(10)) #############################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)
rho10<- sum(epsilonhat[1:(T-10)]*epsilonhat[11:T])/sum(epsilonhat^2)
rho11<- sum(epsilonhat[1:(T-11)]*epsilonhat[12:T])/sum(epsilonhat^2)
rho12<- sum(epsilonhat[1:(T-12)]*epsilonhat[13:T])/sum(epsilonhat^2)
rho13<- sum(epsilonhat[1:(T-13)]*epsilonhat[14:T])/sum(epsilonhat^2)
rho14<- sum(epsilonhat[1:(T-14)]*epsilonhat[15:T])/sum(epsilonhat^2)
rho15<- sum(epsilonhat[1:(T-15)]*epsilonhat[16:T])/sum(epsilonhat^2)



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12,rho13,rho14,rho15),nrow = 1,ncol = 15)

C<-rep(1,15)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   7.254234


pv=1-pchisq(P,10)
pv   ###  0.7012455




###########################     Wald test  ########################################

##########  derivatives  ######################

derivativeofphi<- matrix(numeric(5*T),nrow = 5, ncol = T)



derivativeofphi1<- numeric(T)

derivativeofphi1[2]<- - PPI_lrt[1] -  coef1[4]*derivativeofphi1[1]

for (i in 3:T) {
  derivativeofphi1[i]<- - PPI_lrt[i-1] - coef1[4]*derivativeofphi1[i-1] - coef1[5]*derivativeofphi1[i-2]
}



derivativeofphi2<- numeric(T)

derivativeofphi2[2]<- -  coef1[4]*derivativeofphi2[1]

for (i in 3:T) {
  derivativeofphi2[i]<- - PPI_lrt[i-2] - coef1[4]*derivativeofphi2[i-1] - coef1[5]*derivativeofphi2[i-2]
}


derivativeofphi3<- numeric(T)

derivativeofphi3[2]<- -  coef1[4]*derivativeofphi3[1]

derivativeofphi3[3]<- -  coef1[4]*derivativeofphi3[2] - coef1[5]*derivativeofphi3[1]

for (i in 4:T) {
  derivativeofphi3[i]<- - PPI_lrt[i-3] - coef1[4]*derivativeofphi3[i-1] - coef1[5]*derivativeofphi3[i-2]
}



derivativeoftheta1<- numeric(T)

derivativeoftheta1[2]<- - epsilonhat[1] - coef1[4]*derivativeoftheta1[1]


for (i in 3:T) {
  derivativeoftheta1[i]<- - epsilonhat[i-1] - coef1[4]*derivativeoftheta1[i-1] - coef1[5]*derivativeoftheta1[i-2] 
}



derivativeoftheta2<- numeric(T)

derivativeoftheta2[2]<-  - coef1[4]*derivativeoftheta2[1]

for (i in 3:T) {
  derivativeoftheta2[i]<- - epsilonhat[i-2] - coef1[4]*derivativeoftheta2[i-1] - coef1[5]*derivativeoftheta2[i-2]
}




for (i in 1:T) {
  derivativeofphi[,i]<- c(derivativeofphi1[i],derivativeofphi2[i],derivativeofphi3[i],derivativeoftheta1[i],derivativeoftheta2[i])
  
}

####################################################

omega<- matrix(numeric(25*T),nrow=25,ncol=T)
sigma<- matrix(numeric(25*T),nrow=25,ncol=T)


for (j in 1:T) {
  
  omega[,j]<- as.vector(4*epsilonhat[j]^2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  sigma[,j]<- as.vector(2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  
}

Omega<-matrix(apply(omega,1,mean),nrow=5,ncol=5)
Sigma<-matrix(apply(sigma,1,mean),nrow=5,ncol=5)

Gamma<-solve(Sigma)%*%Omega%*%solve(Sigma)
var<- diag(Gamma)/T
round(sqrt(var),4)
##   0.0834 0.0612 0.0730 0.0311 0.0266



##########  test all  are significantly not equal to 0 ################################
C<- rep(1,5)
tao<- diag(C)
r<- rep(0,5)

W<- c(T*t(tao%*%coef1-r)%*%solve(tao%*%Gamma%*%t(tao))%*%(tao%*%coef1-r))

pv=1-pchisq(W,5)
pv   ### 0


