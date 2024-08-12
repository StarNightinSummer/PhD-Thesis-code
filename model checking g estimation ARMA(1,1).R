rm(list=ls())##回看数据不用运行，产生新数据需要运行来清理运行环境


simulationTimes=10000


################# parameter setting ##################################
T=600  #sample size
phi01=0.3
phi02= 0.5
theta0=0.7 # 0.5 0.7
tao=0.5 # 0.5 0.9
#m=2
delta=5 # 1 5
sigma0=1 #1
sigma1=sigma0*delta
r0=0
epsilon0=0


C<-rep(1,6)
I<-diag(C)

##### calculation of (the intergral of g^4)/(the intergral of g^2)^2

integral<- (tao+(1-tao)*delta^4)/((tao+(1-tao)*delta^2)^2)##abrupt



#C1<- 1/(1-phi0^2)   #C1 of AR(1) model# 
C1<- matrix(c(1/(1-phi01^2),1/(1+phi01*theta0),1/(1+phi01*theta0),1/(1-theta0^2)),nrow=2,ncol=2)#C1 of ARMA(1,1) model#


#C2<-matrix(c(-1, -phi0,-phi0^2,-phi0^3,-phi0^4,-phi0^5),nrow = 1, ncol = 6)#C2 of AR(1) model#
C2<- matrix(c(-1,-1,-phi01,theta0,-phi01^2,-theta0^2,-phi01^3,theta0^3,-phi01^4,-theta0^4,-phi01^5,theta0^5),nrow=2,ncol=6)#C2 of ARMA(1,1)model#

q<- integral*I - integral*t(C2)%*%solve(C1)%*%C2
q1<- solve(q)
q2<- integral*I
q5<- solve(q2)


Q=numeric(simulationTimes)

P=numeric(simulationTimes)


#EstimateOfcoef<-numeric(simulationTimes)#estimation for coef of AR(1) model#
EstimateOfcoef<-matrix(numeric(2*simulationTimes),nrow = 2,ncol = simulationTimes)#estimation for coef#


#epsilonhat=numeric(T)


######### data generation function #############
data=function(n,phi1,phi2,theta){
  
  noise=rnorm(n)
  #noise=rlaplace(n, location=0, scale=sqrt(2)/2)
  r=numeric(n)
  b=numeric(n)
  g=numeric(n)
  epsilon=numeric(n) 
  
  
  for (i in 1:n) { 
    if((i/n)>=tao) b[i]<-1
    #b[i]<-(i/n)^m
    g[i]<-(sigma0)^2+((sigma1)^2-(sigma0)^2)*b[i]
    epsilon[i]<-sqrt(g[i])*noise[i]
  }
  
  
  
  r[1]<-phi1*r0+epsilon[1]+theta*epsilon0
  r[2]<-phi1*r[1]+phi2*r0+epsilon[2]+theta*epsilon[1]
  
  for (i in 3:n) {
    r[i]<-phi1*r[i-1]+phi2*r[i-2]+epsilon[i]+theta*epsilon[i-1]
  }
  
  
  return(r)
  
}


for (j in 1:simulationTimes)
{
  
  
  sample=data(T, phi01, phi02, theta0)
  #print(sample)
  #print(length(sample))
  #plot(sample,main="sample",ylab="",xlab="",type="l")
  
  
  #m0<-arima(sample,order=c(1,0,0), include.mean = FALSE)# fitted AR(1) model#
  EstimateOfcoef0<-arima(sample,order=c(1,0,1), include.mean = FALSE)$coef
  
  
  epsilonhat=numeric(T)
  
  
  epsilonhat[1]<- sample[1]
  # epsilonhat[2:T]<- sample[2:T] - EstimateOfcoef0*sample[1:(T-1)] 
  
  for (i in 2:T) {
    epsilonhat[i]<- sample[i] - EstimateOfcoef0[1]*sample[i-1] - EstimateOfcoef0[2]*epsilonhat[i-1]
  }
  #lines(epsilonhat, lwd=1,lty=2, col='blue')
  # plot(epsilonhat,main="epsilonhat",ylab="",xlab="",type="l")
  
  
  ########## checking whether "m0$residuals" and "epsilonhat" are identical 
  
  # m0$residuals %in% epsilonhat
  # table(m0$residuals %in% epsilonhat)
  # identical(m0$residuals,epsilonhat)
  
  
  rho1<- sum(epsilonhat[1:(T-1)]*epsilonhat[2:T])/sum(epsilonhat^2)
  rho2<- sum(epsilonhat[1:(T-2)]*epsilonhat[3:T])/sum(epsilonhat^2)
  rho3<- sum(epsilonhat[1:(T-3)]*epsilonhat[4:T])/sum(epsilonhat^2)
  rho4<- sum(epsilonhat[1:(T-4)]*epsilonhat[5:T])/sum(epsilonhat^2)
  rho5<- sum(epsilonhat[1:(T-5)]*epsilonhat[6:T])/sum(epsilonhat^2)
  rho6<- sum(epsilonhat[1:(T-6)]*epsilonhat[7:T])/sum(epsilonhat^2)
  rho<- matrix(c(rho1,rho2,rho3,rho4,rho5,rho6),nrow = 1,ncol = 6)
  
  #################  Ljung-Box test rho ##################
  
  #Rho<- matrix(c(rho1/sqrt(T-1),rho2/sqrt(T-2),rho3/sqrt(T-3),rho4/sqrt(T-4),rho5/sqrt(T-5),rho6/sqrt(T-6)),nrow = 1,ncol = 6)
  
  
  
  ############################################################ the estimation of g #################################### 
  
  
  integral1<- (mean(epsilonhat[2:T]^2*epsilonhat[1:(T-1)]^2))/((mean(epsilonhat^2))^2)
  q3<- integral1*I - integral1*t(C2)%*%solve(C1)%*%C2
  q4<- integral1*I
  
  
  Q[j]<- T*rho %*% solve(q2) %*% t(rho)  
  
  P[j]<- T*rho %*% solve(q4) %*% t(rho)
  
  
  
  EstimateOfcoef[,j]<- as.vector(EstimateOfcoef0)
  
 
  
  
  
}
EstimateOfcoef1<- EstimateOfcoef[1,]
EstimateOfcoef2<- EstimateOfcoef[2,]

meancoef1<- mean(EstimateOfcoef1,na.rm = TRUE)
meancoef2<- mean(EstimateOfcoef2,na.rm = TRUE)


sum(Q>9.49)/simulationTimes # chi distribution of df=4#
sum(P>9.49)/simulationTimes # chi distribution of df=4#

forecast::Acf(sample,xlab="",ylab="",col="red",main="ACF of sample")
forecast::Pacf(sample,xlab="",ylab="",col="red",main="PACF of sample")
forecast::Acf(sample^2,xlab="",ylab="",col="red",main="ACF of sample^2")
forecast::Pacf(sample^2,xlab="",ylab="",col="red",main="PACF of sample^2")



