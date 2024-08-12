#####################   first difference of M1    #########################
install.packages("tsibble")
library(tsibble)


rm(list = ls())
d <- read.table("D:/中英文论文/data/M1(1).txt", header=TRUE)
M1 <- ts(d[["M1SL"]], start=c(1959,1), frequency=12)
length(M1)  #736
plot(M1,main="M1",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)


M1_d<- diff(M1)
T=length(M1_d)  #735
plot(M1_d,main="M1_d",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)




M1_lrt<-diff(log(M1))
M1_lrt<-M1_lrt-mean(M1_lrt)
T=length(M1_lrt)  #735
plot(M1_lrt,main="M1_lrt",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)


############  fit ARMA(3,2) model  ############################

#Test mean being zero.
t.test(M1_lrt)   


m1<-arima(M1_lrt, order = c(3,0,2),include.mean = FALSE)
m1

coef1<-m1$coef
coef1   ##  0.2188  0.5559  0.1082  -0.0232  -0.5138   aic = -5071.39


mean(M1_lrt)
# -5.04827e-19



########################  calculate residuals  ###################################

epsilonhat=numeric(T)


epsilonhat[1]<- M1_lrt[1] 

epsilonhat[2]<- M1_lrt[2] - coef1[1]*M1_lrt[1] - coef1[4]*epsilonhat[1] 

epsilonhat[3]<- M1_lrt[3] - coef1[1]*M1_lrt[2] - coef1[2]*M1_lrt[1] - coef1[4]*epsilonhat[2] - coef1[5]*epsilonhat[1] 



for (i in 4:T) {
  epsilonhat[i]<- M1_lrt[i] - coef1[1]*M1_lrt[i-1] - coef1[2]*M1_lrt[i-2] - coef1[3]*M1_lrt[i-3] - coef1[4]*epsilonhat[i-1] - coef1[5]*epsilonhat[i-2] 
  
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
#  -10.94523

BIC
# -10.91394



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
P  #  0.5689472


pv=1-pchisq(P,2)
pv  #  0.7524102



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
P  #   1.303229


pv=1-pchisq(P,6)
pv   ###   0.9714794



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
P  # 2.821233


pv=1-pchisq(P,10)
pv   ###  0.9853235




###########################     Wald test  ########################################

##########  derivatives  ######################

derivativeofphi<- matrix(numeric(5*T),nrow = 5, ncol = T)



derivativeofphi1<- numeric(T)

derivativeofphi1[2]<- - M1_lrt[1] -  coef1[4]*derivativeofphi1[1]

for (i in 3:T) {
  derivativeofphi1[i]<- - M1_lrt[i-1] - coef1[4]*derivativeofphi1[i-1] - coef1[5]*derivativeofphi1[i-2]
}



derivativeofphi2<- numeric(T)

derivativeofphi2[2]<- -  coef1[4]*derivativeofphi2[1]

for (i in 3:T) {
  derivativeofphi2[i]<- - M1_lrt[i-2] - coef1[4]*derivativeofphi2[i-1] - coef1[5]*derivativeofphi2[i-2]
}



derivativeofphi3<- numeric(T)

derivativeofphi3[2]<- -  coef1[4]*derivativeofphi3[1]

derivativeofphi3[3]<- -  coef1[4]*derivativeofphi3[2] - coef1[5]*derivativeofphi3[1]

for (i in 4:T) {
  derivativeofphi3[i]<- - M1_lrt[i-3] - coef1[4]*derivativeofphi3[i-1] - coef1[5]*derivativeofphi3[i-2]
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
##   0.2213 0.1736 0.1225 0.1807 0.1498



##########  test all  are significantly not equal to 0 ################################
C<- rep(1,5)
tao<- diag(C)
r<- rep(0,5)

W<- c(T*t(tao%*%coef1-r)%*%solve(tao%*%Gamma%*%t(tao))%*%(tao%*%coef1-r))

pv=1-pchisq(W,5)
pv   ### 0



##################  phi1, theta1 #####################

tao1<-matrix(numeric(10),nrow=2,ncol=5)
tao1[1,1]<-1
tao1[2,4]<-1
tao1

r1<- rep(0,2)

W1<- c(T*t(tao1%*%coef1-r1)%*%solve(tao1%*%Gamma%*%t(tao1))%*%(tao1%*%coef1-r1))

pv=1-pchisq(W1,2)
pv  ##    0.4958995






#################   fix parameter  ###############################


#Test mean being zero.
t.test(M1_lrt)   


m1<-arima(M1_lrt, order = c(3,0,2),include.mean = FALSE, fixed = c(0,NA,NA,0,NA), transform.pars=FALSE)
m1

coef1<-m1$coef
coef1   ##   0  0.6436  0.2189    0  -0.5667   aic = -5052.87


mean(M1_lrt)
# -5.04827e-19



########################  calculate residuals  ###################################

epsilonhat=numeric(T)


epsilonhat[1]<- M1_lrt[1] 

epsilonhat[2]<- M1_lrt[2]

epsilonhat[3]<- M1_lrt[3] - coef1[2]*M1_lrt[1] - coef1[5]*epsilonhat[1] 



for (i in 4:T) {
  epsilonhat[i]<- M1_lrt[i] - coef1[2]*M1_lrt[i-2] - coef1[3]*M1_lrt[i-3] - coef1[5]*epsilonhat[i-2] 
  
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
AIC  = LALS + 2*(2+1)/T
BIC  = LALS + (2+1)*log(T)/T
AIC
#  -10.90134

BIC
# -10.88256



########################    model checking M=5  (S(2)) ##########################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)


rho<- matrix(c(rho1,rho2,rho3,rho4,rho5),nrow = 1,ncol = 5)

C<-rep(1,5)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #  1.224488


pv=1-pchisq(P,2)
pv  #  0.5421328



############   M=9 (S(6)) #############################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)
rho9<- sum(epsilonhat[1:(T-9)]*epsilonhat[10:T])/sum(epsilonhat^2)



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9),nrow = 1,ncol = 9)

C<-rep(1,9)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   2.104668


pv=1-pchisq(P,6)
pv   ###   0.9098249



############   M=13  (S(10)) #############################


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




rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12,rho13),nrow = 1,ncol = 13)

C<-rep(1,13)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  # 2.97499


pv=1-pchisq(P,10)
pv   ###  0.9820065




###########################     Wald test  ########################################

##########  derivatives  ######################

derivativeofphi<- matrix(numeric(3*T),nrow = 3, ncol = T)



derivativeofphi2<- numeric(T)


for (i in 3:T) {
  derivativeofphi2[i]<- - M1_lrt[i-2] - coef1[5]*derivativeofphi2[i-2]
}



derivativeofphi3<- numeric(T)

derivativeofphi3[3]<- - coef1[5]*derivativeofphi3[1]

for (i in 4:T) {
  derivativeofphi3[i]<- - M1_lrt[i-3] - coef1[5]*derivativeofphi3[i-2]
}




derivativeoftheta2<- numeric(T)

for (i in 3:T) {
  derivativeoftheta2[i]<- - epsilonhat[i-2] - coef1[5]*derivativeoftheta2[i-2]
}



for (i in 1:T) {
  derivativeofphi[,i]<- c(derivativeofphi2[i],derivativeofphi3[i],derivativeoftheta2[i])
  
}

####################################################

omega<- matrix(numeric(9*T),nrow=9,ncol=T)
sigma<- matrix(numeric(9*T),nrow=9,ncol=T)


for (j in 1:T) {
  
  omega[,j]<- as.vector(4*epsilonhat[j]^2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  sigma[,j]<- as.vector(2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  
}

Omega<-matrix(apply(omega,1,mean),nrow=3,ncol=3)
Sigma<-matrix(apply(sigma,1,mean),nrow=3,ncol=3)

Gamma<-solve(Sigma)%*%Omega%*%solve(Sigma)
var<- diag(Gamma)/T
round(sqrt(var),4)
# 0.1101 0.0726 0.1283



##########  test all  are significantly not equal to 0 ################################
C<- rep(1,3)
tao<- diag(C)
r<- rep(0,3)

W<- c(T*t(tao%*%coef1[c(2,3,5)]-r)%*%solve(tao%*%Gamma%*%t(tao))%*%(tao%*%coef1[c(2,3,5)]-r))

pv=1-pchisq(W,3)
pv   ### 0







#######################   log return of M1  #####################################

rm(list = ls())
d <- read.table("D:/中英文论文/data/M1(1).txt", header=TRUE)
M1 <- ts(d[["M1SL"]], start=c(1959,1), frequency=12)
length(M1)  #736
plot(M1,main="M1",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)


M1_d<- diff(M1)
T=length(M1_d)  #735
plot(M1_d,main="M1_d",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)




M1_lrt<-diff(log(M1))
M1_lrt<-M1_lrt-mean(M1_lrt)
T=length(M1_lrt)  #735
plot(M1_lrt,main="M1_lrt",col="red",xlab="",ylab="")
# minor.tick(nx = 10,  tick.ratio = 0.5)



forecast::Acf(M1_lrt,xlab="",ylab="",col="red",main="ACF of M1_lrt")
forecast::Pacf(M1_lrt,xlab="",ylab="",col="red",main="PACF of M1_lrt")
forecast::Acf(M1_lrt^2,xlab="",ylab="",col="red",main="ACF of M1_lrt^2")
forecast::Pacf(M1_lrt^2,xlab="",ylab="",col="red",main="PACF of M1_lrt^2")



p=9
q=9
table.aic<-matrix(numeric(p*q), nrow=p,ncol=q)
for (i in 1:p) {
  for (j in 1:q) {
    table.aic[i,j]<- arima(M1_lrt,order = c(i-1,0,j-1))$aic 
    j=j+1
  }
  i=i+1
}
table.aic
which.min(table.aic)   # 


ar(M1_lrt,method = "ols")



############  ARMA(14,0) ############################

#Test mean being zero.
t.test(M1_lrt)   


m1<-arima(M1_lrt, order = c(14,0,0),include.mean = FALSE)
m1

coef1<-m1$coef
coef1   ##  0.2216  0.0781  0.1666  -0.0211  0.1251  0.1321  0.0070  -0.0381  0.1510  -0.0239  -0.0309  -0.0475  -0.1299  0.1548    aic = -5101.6

mean(M1_lrt)
# -5.04827e-19




########################  calculate residuals  ###################################

epsilonhat=numeric(T)


epsilonhat[1]<- M1_lrt[1] 

epsilonhat[2]<- M1_lrt[2] - coef1[1]*M1_lrt[1]

epsilonhat[3]<- M1_lrt[3] - coef1[1]*M1_lrt[2] - coef1[2]*M1_lrt[1]

epsilonhat[4]<- M1_lrt[4] - coef1[1]*M1_lrt[3] - coef1[2]*M1_lrt[2] - coef1[3]*M1_lrt[1] 

epsilonhat[5]<- M1_lrt[5] - coef1[1]*M1_lrt[4] - coef1[2]*M1_lrt[3] - coef1[3]*M1_lrt[2] - coef1[4]*M1_lrt[1]

epsilonhat[6]<- M1_lrt[6] - coef1[1]*M1_lrt[5] - coef1[2]*M1_lrt[4] - coef1[3]*M1_lrt[3] - coef1[4]*M1_lrt[2] - coef1[5]*M1_lrt[1]

epsilonhat[7]<- M1_lrt[7] - coef1[1]*M1_lrt[6] - coef1[2]*M1_lrt[5] - coef1[3]*M1_lrt[4] - coef1[4]*M1_lrt[3] - coef1[5]*M1_lrt[2] - coef1[6]*M1_lrt[1]

epsilonhat[8]<- M1_lrt[8] - coef1[1]*M1_lrt[7] - coef1[2]*M1_lrt[6] - coef1[3]*M1_lrt[5] - coef1[4]*M1_lrt[4] - coef1[5]*M1_lrt[3] - coef1[6]*M1_lrt[2] - coef1[7]*M1_lrt[1]

epsilonhat[9]<- M1_lrt[9] - coef1[1]*M1_lrt[8] - coef1[2]*M1_lrt[7] - coef1[3]*M1_lrt[6] - coef1[4]*M1_lrt[5] - coef1[5]*M1_lrt[4] - coef1[6]*M1_lrt[3] - coef1[7]*M1_lrt[2] - coef1[8]*M1_lrt[1]

epsilonhat[10]<- M1_lrt[10] - coef1[1]*M1_lrt[9] - coef1[2]*M1_lrt[8] - coef1[3]*M1_lrt[7] - coef1[4]*M1_lrt[6] - coef1[5]*M1_lrt[5] - coef1[6]*M1_lrt[4] - coef1[7]*M1_lrt[3] - coef1[8]*M1_lrt[2] - coef1[9]*M1_lrt[1]

epsilonhat[11]<- M1_lrt[11] - coef1[1]*M1_lrt[10] - coef1[2]*M1_lrt[9] - coef1[3]*M1_lrt[8] - coef1[4]*M1_lrt[7] - coef1[5]*M1_lrt[6] - coef1[6]*M1_lrt[5] - coef1[7]*M1_lrt[4] - coef1[8]*M1_lrt[3] - coef1[9]*M1_lrt[2] - coef1[10]*M1_lrt[1]

epsilonhat[12]<- M1_lrt[12] - coef1[1]*M1_lrt[11] - coef1[2]*M1_lrt[10] - coef1[3]*M1_lrt[9] - coef1[4]*M1_lrt[8] - coef1[5]*M1_lrt[7] - coef1[6]*M1_lrt[6] - coef1[7]*M1_lrt[5] - coef1[8]*M1_lrt[4] - coef1[9]*M1_lrt[3] - coef1[10]*M1_lrt[2] - coef1[11]*M1_lrt[1]

epsilonhat[13]<- M1_lrt[13] - coef1[1]*M1_lrt[12] - coef1[2]*M1_lrt[11] - coef1[3]*M1_lrt[10] - coef1[4]*M1_lrt[9] - coef1[5]*M1_lrt[8] - coef1[6]*M1_lrt[7] - coef1[7]*M1_lrt[6] - coef1[8]*M1_lrt[5] - coef1[9]*M1_lrt[4] - coef1[10]*M1_lrt[3] - coef1[11]*M1_lrt[2] - coef1[12]*M1_lrt[1]

epsilonhat[14]<- M1_lrt[14] - coef1[1]*M1_lrt[13] - coef1[2]*M1_lrt[12] - coef1[3]*M1_lrt[11] - coef1[4]*M1_lrt[10] - coef1[5]*M1_lrt[9] - coef1[6]*M1_lrt[8] - coef1[7]*M1_lrt[7] - coef1[8]*M1_lrt[6] - coef1[9]*M1_lrt[5] - coef1[10]*M1_lrt[4] - coef1[11]*M1_lrt[3] - coef1[12]*M1_lrt[2] - coef1[13]*M1_lrt[1]


for (i in 15:T) {
  epsilonhat[i]<- M1_lrt[i] - coef1[1]*M1_lrt[i-1] - coef1[2]*M1_lrt[i-2] - coef1[3]*M1_lrt[i-3] - coef1[4]*M1_lrt[i-4] - coef1[5]*M1_lrt[i-5] - coef1[6]*M1_lrt[i-6] - coef1[7]*M1_lrt[i-7] - coef1[8]*M1_lrt[i-8] - coef1[9]*M1_lrt[i-9] - coef1[10]*M1_lrt[i-10] - coef1[11]*M1_lrt[i-11] - coef1[12]*M1_lrt[i-12] - coef1[13]*M1_lrt[i-13] - coef1[14]*M1_lrt[i-14]  
  
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
#  -11.068

BIC
# -10.98039



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
P  #   0.2588194


pv=1-pchisq(P,2)
pv  #  0.8786139



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
P  #   0.6307024


pv=1-pchisq(P,6)
pv   ###  0.9958663



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
P  #  1.823816


pv=1-pchisq(P,10)
pv   ###   0.9975208




###########################     Wald test  ########################################

##########  derivatives  ######################

derivativeofphi<- matrix(numeric(14*T),nrow = 14, ncol = T)


derivativeofphi1<- numeric(T)

for (i in 2:T) {
  derivativeofphi1[i]<- - M1_lrt[i-1] 
}


derivativeofphi2<- numeric(T)

for (i in 3:T) {
  derivativeofphi2[i]<- - M1_lrt[i-2]
}


derivativeofphi3<- numeric(T)

for (i in 4:T) {
  derivativeofphi3[i]<- - M1_lrt[i-3]  
}



derivativeofphi4<- numeric(T)

for (i in 5:T) {
  derivativeofphi4[i]<- - M1_lrt[i-4]  
}


derivativeofphi5<- numeric(T)

for (i in 6:T) {
  derivativeofphi5[i]<- - M1_lrt[i-5]  
}



derivativeofphi6<- numeric(T)

for (i in 7:T) {
  derivativeofphi6[i]<- - M1_lrt[i-6]  
}



derivativeofphi7<- numeric(T)

for (i in 8:T) {
  derivativeofphi7[i]<- - M1_lrt[i-7]  
}


derivativeofphi8<- numeric(T)

for (i in 9:T) {
  derivativeofphi8[i]<- - M1_lrt[i-8]  
}


derivativeofphi9<- numeric(T)

for (i in 10:T) {
  derivativeofphi9[i]<- - M1_lrt[i-9]  
}


derivativeofphi10<- numeric(T)

for (i in 11:T) {
  derivativeofphi10[i]<- - M1_lrt[i-10]  
}


derivativeofphi11<- numeric(T)

for (i in 12:T) {
  derivativeofphi11[i]<- - M1_lrt[i-11]  
}


derivativeofphi12<- numeric(T)

for (i in 13:T) {
  derivativeofphi12[i]<- - M1_lrt[i-12]  
}


derivativeofphi13<- numeric(T)

for (i in 14:T) {
  derivativeofphi13[i]<- - M1_lrt[i-13]
}


derivativeofphi14<- numeric(T)

for (i in 15:T) {
  derivativeofphi14[i]<- - M1_lrt[i-14]
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
##   0.1709 0.0610 0.0681 0.0632 0.0520 0.0464 0.0425 0.0434 0.0570 0.0564 0.0458 0.0511 0.0699 0.0424



##########  test all  are significantly not equal to 0 ################################
C<- rep(1,14)
tao<- diag(C)
r<- rep(0,14)

W<- c(T*t(tao%*%coef1-r)%*%solve(tao%*%Gamma%*%t(tao))%*%(tao%*%coef1-r))

pv=1-pchisq(W,14)
pv   ### 0


#########  test theta2 and theta4  equal to 0 ##################

tao1<-matrix(numeric(112),nrow=8,ncol=14)
tao1[1,1]<-1
tao1[2,2]<-1
tao1[3,4]<-1
tao1[4,7]<-1
tao1[5,8]<-1
tao1[6,10]<-1
tao1[7,11]<-1
tao1[8,12]<-1
tao1

r1<- rep(0,8)

W1<- c(T*t(tao1%*%coef1-r1)%*%solve(tao1%*%Gamma%*%t(tao1))%*%(tao1%*%coef1-r1))

pv=1-pchisq(W1,8)
pv  ##    0.6892303



############  fix parameters ############################

#Test mean being zero.
t.test(M1_lrt)   


m1<-arima(M1_lrt, order = c(14,0,0),include.mean = FALSE, fixed = c(0,0,NA,0,NA,NA,0,0,NA,0,0,0,NA,NA),transform.pars=FALSE)
m1

coef1<-m1$coef
coef1   ##  0    0  0.1844    0  0.1407  0.1535    0    0  0.1513     0     0     0  -0.1267  0.1235   aic = -5081.23


mean(M1_lrt)
# -5.04827e-19




########################  calculate residuals  ###################################

epsilonhat=numeric(T)


epsilonhat[1]<- M1_lrt[1] 

epsilonhat[2]<- M1_lrt[2] 

epsilonhat[3]<- M1_lrt[3] 

epsilonhat[4]<- M1_lrt[4] - coef1[3]*M1_lrt[1] 

epsilonhat[5]<- M1_lrt[5] - coef1[3]*M1_lrt[2]

epsilonhat[6]<- M1_lrt[6] - coef1[3]*M1_lrt[3] - coef1[5]*M1_lrt[1]

epsilonhat[7]<- M1_lrt[7] - coef1[3]*M1_lrt[4] - coef1[5]*M1_lrt[2] - coef1[6]*M1_lrt[1]

epsilonhat[8]<- M1_lrt[8] - coef1[3]*M1_lrt[5] - coef1[5]*M1_lrt[3] - coef1[6]*M1_lrt[2]

epsilonhat[9]<- M1_lrt[9] - coef1[3]*M1_lrt[6] - coef1[5]*M1_lrt[4] - coef1[6]*M1_lrt[3]

epsilonhat[10]<- M1_lrt[10] - coef1[3]*M1_lrt[7] - coef1[5]*M1_lrt[5] - coef1[6]*M1_lrt[4] - coef1[9]*M1_lrt[1]

epsilonhat[11]<- M1_lrt[11] - coef1[3]*M1_lrt[8] - coef1[5]*M1_lrt[6] - coef1[6]*M1_lrt[5] - coef1[9]*M1_lrt[2]

epsilonhat[12]<- M1_lrt[12] - coef1[3]*M1_lrt[9] - coef1[5]*M1_lrt[7] - coef1[6]*M1_lrt[6] - coef1[9]*M1_lrt[3]

epsilonhat[13]<- M1_lrt[13] - coef1[3]*M1_lrt[10] - coef1[5]*M1_lrt[8] - coef1[6]*M1_lrt[7] - coef1[9]*M1_lrt[4]

epsilonhat[14]<- M1_lrt[14] - coef1[3]*M1_lrt[11] - coef1[5]*M1_lrt[9] - coef1[6]*M1_lrt[8] - coef1[9]*M1_lrt[5] - coef1[13]*M1_lrt[1]



for (i in 15:T) {
  epsilonhat[i]<- M1_lrt[i] - coef1[3]*M1_lrt[i-3] - coef1[5]*M1_lrt[i-5] - coef1[6]*M1_lrt[i-6] - coef1[9]*M1_lrt[i-9] - coef1[13]*M1_lrt[i-13] - coef1[14]*M1_lrt[i-14]
  
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
AIC  = LALS + 2*(6+0)/T
BIC  = LALS + (6+0)*log(T)/T
AIC
#  -10.89221

BIC
# -10.85466





########################    model checking M=8  (S(2)) ##########################


rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
rho7<- sum(epsilonhat[1:(T-7)]*epsilonhat[8:T])/sum(epsilonhat^2)
rho8<- sum(epsilonhat[1:(T-8)]*epsilonhat[9:T])/sum(epsilonhat^2)


rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8),nrow = 1,ncol = 8 )

C<-rep(1,8)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   1.832499


pv=1-pchisq(P,2)
pv  #   0.4000166



############   M=12 (S(6)) #############################


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



rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10,rho11,rho12),nrow = 1,ncol = 12)

C<-rep(1,12)
I<-diag(C)


integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
q<- integral1*I

P<- T*rho %*% solve(q) %*% t(rho)
P
P<- as.numeric(P)
P  #   2.131315


pv=1-pchisq(P,6)
pv   ###  0.9072339



############   M=16  (S(10)) #############################


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
P  #  2.215656


pv=1-pchisq(P,10)
pv   ###   0.9944041




###########################     Wald test  ########################################

##########  derivatives  ######################

derivativeofphi<- matrix(numeric(6*T),nrow = 6, ncol = T)


derivativeofphi3<- numeric(T)

for (i in 4:T) {
  derivativeofphi3[i]<- - M1_lrt[i-3]  
}


derivativeofphi5<- numeric(T)

for (i in 6:T) {
  derivativeofphi5[i]<- - M1_lrt[i-5]  
}



derivativeofphi6<- numeric(T)

for (i in 7:T) {
  derivativeofphi6[i]<- - M1_lrt[i-6]  
}



derivativeofphi9<- numeric(T)

for (i in 10:T) {
  derivativeofphi9[i]<- - M1_lrt[i-9]  
}



derivativeofphi13<- numeric(T)

for (i in 14:T) {
  derivativeofphi13[i]<- - M1_lrt[i-13]  
}


derivativeofphi14<- numeric(T)

for (i in 15:T) {
  derivativeofphi14[i]<- - M1_lrt[i-14]  
}



for (i in 1:T) {
  derivativeofphi[,i]<- c(derivativeofphi3[i],derivativeofphi5[i],derivativeofphi6[i],derivativeofphi9[i],derivativeofphi13[i],derivativeofphi14[i])
  
}

####################################################

omega<- matrix(numeric(36*T),nrow=36,ncol=T)
sigma<- matrix(numeric(36*T),nrow=36,ncol=T)


for (j in 1:T) {
  
  omega[,j]<- as.vector(4*epsilonhat[j]^2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  sigma[,j]<- as.vector(2*derivativeofphi[,j,drop=FALSE]%*%t(derivativeofphi[,j,drop=FALSE]))
  
}

Omega<-matrix(apply(omega,1,mean),nrow=6,ncol=6)
Sigma<-matrix(apply(sigma,1,mean),nrow=6,ncol=6)

Gamma<-solve(Sigma)%*%Omega%*%solve(Sigma)
var<- diag(Gamma)/T
round(sqrt(var),4)
##   0.0672 0.0476 0.0506 0.0579 0.0676 0.0379



##########  test all  are significantly not equal to 0 ################################
C<- rep(1,6)
tao<- diag(C)
r<- rep(0,6)

W<- c(T*t(tao%*%coef1[c(3,5,6,9,13,14)]-r)%*%solve(tao%*%Gamma%*%t(tao))%*%(tao%*%coef1[c(3,5,6,9,13,14)]-r))

pv=1-pchisq(W,6)
pv   ### 0

