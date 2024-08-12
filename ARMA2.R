######################### in R format #############################
rm(list=ls())##回看数据不用运行，产生新数据需要运行来清理运行环境  
# install.packages("LaplacesDemon")
# library(LaplacesDemon)


simulationTimes=10000
T=400 # sample size

EstimateOfcoef<-matrix(numeric(2*simulationTimes),nrow = 2,ncol = simulationTimes)#estimation for coef#


Omega_11<-numeric(simulationTimes)
Omega_12<-numeric(simulationTimes)
Omega_22<-numeric(simulationTimes)
Sigma_11<-numeric(simulationTimes)
Sigma_12<-numeric(simulationTimes)
Sigma_22<-numeric(simulationTimes)

Gamma<-matrix(numeric(4*simulationTimes),nrow = 4, ncol = simulationTimes)


############################### generation function #####################################
data=function(n,phi,theta){
  
  noise=rnorm(n)
  #noise=rlaplace(n, location=0, scale=1)
  r=numeric(n)
  b=numeric(n)
  g=numeric(n)
  epsilon=numeric(n)
  
  
  for (i in 1:n) {
    
    if((i/n)>=tao) b[i]<-1 else b[i]<-0 #else语句可以去掉
    #b[i]<-(i/n)^m
    g[i]<-(sigma0)^2+((sigma1)^2-(sigma0)^2)*b[i]
    epsilon[i]<-sqrt(g[i])*noise[i]
    
  }
  
  
  r[1]<-phi*r0+epsilon[1]+theta*epsilon0
  
  for (i in 2:n) {
    r[i]<-phi*r[i-1]+epsilon[i]+theta*epsilon[i-1]
  }

  #length(r)
  return(r)
  
}

############################################# Estimation function #########################

for (j in 1:simulationTimes)
{
  
  
  phi0=0.3 # 0.5 0.9
  theta0=0.7
  tao=0.9 # 0.5 0.9
  #m=2 # 2
  delta=5 # 1 5
  sigma0=1
  sigma1=sigma0*delta
  r0=0
  epsilon0=0
  sample=data(T, phi0, theta0)
  #print(sample)
  #print(length(sample))
  #plot(sample,main=" ",ylab="",xlab="",type="l")
  
  #m0<-arima(sample,order=c(1,0,1), include.mean = FALSE)
  EstimateOfcoef0=arima(sample,order=c(1,0,1), include.mean = FALSE)$coef
  
  
  epsilonhat=numeric(T)
  
  derivativeofphi=numeric(T)
  derivativeoftheta=numeric(T)
  
  omega_11<-numeric(T)
  omega_12<-numeric(T)
  omega_22<-numeric(T)
  sigma_11<-numeric(T)
  sigma_12<-numeric(T)
  sigma_22<-numeric(T)
  
  
  epsilonhat[1]<- sample[1]
  for (i in 2:T) {
    epsilonhat[i]<- sample[i] - EstimateOfcoef0[1]*sample[i-1] - EstimateOfcoef0[2]*epsilonhat[i-1]
  }

  
  derivativeofphi[1]=0
  derivativeoftheta[1]=0
  
  for (i in 2:T) {
    derivativeofphi[i]= - sample[i-1] - EstimateOfcoef0[2]*derivativeofphi[i-1]
    derivativeoftheta[i]= - epsilonhat[i-1] - EstimateOfcoef0[2]*derivativeoftheta[i-1]
  }

  
  omega_11<- 4*epsilonhat^2*derivativeofphi^2
  omega_12<- 4*epsilonhat^2*derivativeofphi*derivativeoftheta
  omega_22<- 4*epsilonhat^2*derivativeoftheta^2

  
  Omega_11[j]<- mean(omega_11)
  Omega_12[j]<- mean(omega_12)
  Omega_22[j]<- mean(omega_22)
  
  
  sigma_11<- 2*derivativeofphi^2
  sigma_12<- 2*derivativeofphi*derivativeoftheta
  sigma_22<- 2*derivativeoftheta^2
  
  
  Sigma_11[j]<- mean(sigma_11)
  Sigma_12[j]<- mean(sigma_12)
  Sigma_22[j]<- mean(sigma_22)
  
  # Omega[j]<-matrix(c(Omega_11[j],Omega_12[j],Omega_12[j], Omega_22[j]),nrow=2,ncol=2)
  # Sigma[j]<-matrix(c(Sigma_11[j],Sigma_12[j],Sigma_12[j], Sigma_22[j]),nrow=2,ncol=2)
  
  Gamma[,j]<-as.vector(solve(matrix(c(Sigma_11[j],Sigma_12[j],Sigma_12[j], Sigma_22[j]),nrow=2,ncol=2)) %*% (matrix(c(Omega_11[j],Omega_12[j],Omega_12[j], Omega_22[j]),nrow=2,ncol=2)) %*% solve(matrix(c(Sigma_11[j],Sigma_12[j],Sigma_12[j], Sigma_22[j]),nrow=2,ncol=2)))

  EstimateOfcoef[,j]<- as.vector(EstimateOfcoef0)

  
}


#print(EstimateOfcoef)
#print(EstimateOfvar)


#meancoef<-round(apply(EstimateOfcoef,1,mean),digits=6)


EstimateOfcoef1<- EstimateOfcoef[1,]
EstimateOfcoef2<- EstimateOfcoef[2,]

meancoef1<- mean(EstimateOfcoef1,na.rm = TRUE)
meancoef2<- mean(EstimateOfcoef2,na.rm = TRUE) 



#meanvar110<-var(EstimateOfcoef1,na.rm=TRUE)*T
#meanvar220<-var(EstimateOfcoef2,na.rm=TRUE)*T
#meanvar120<-cov(EstimateOfcoef1, EstimateOfcoef2,use = "pairwise.complete.obs")*T
meanvar<-var(t(EstimateOfcoef))*T


#################################### calculate the limitation of Omega #########################


### calculation of (the integral of g^4)/(the integral of g^2)^2

integral<- (tao+(1-tao)*delta^4)/((tao+(1-tao)*delta^2)^2)##abrupt
#integral<- (1+2*(delta^2-1)/(m+1)+(delta^2-1)^2/(2*m+1))/((1+(delta^2-1)/(m+1))^2)#trending

Gamma0<-matrix(c(1/(1-phi0^2),1/(1+phi0*theta0),1/(1+phi0*theta0),1/(1-theta0^2)),nrow=2,ncol=2)
Gamma00<-integral*solve(Gamma0)

# use the estimators to calculate the asymptotic variance
#Gamma000<-matrix(c(1/(1-meancoef1^2),1/(1+meancoef1*meancoef2),1/(1+meancoef1*meancoef2),1/(1-meancoef2^2)),nrow=2,ncol=2)
#Gamma0000<-integral*solve(Gamma000)


################################################################################################
### use iteration to obtain the estimation of Omega, Sigma and Gamma

Gamma1<-apply(Gamma,1,mean)


plot(sample,main=" ",ylab="",xlab="",type="l")

as.vector(Gamma00)
as.vector(meanvar)
Gamma1


Gamma00[1,1]
Gamma00[2,2]

#Gamma0000[1,1]
#Gamma0000[2,2]


meanvar[1,1]
meanvar[2,2]

Gamma1[1]
Gamma1[4]
