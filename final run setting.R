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




###############     second stationary for VAR(1)-GARCH(1,1) model   ############
N=2   #  dimension


mu = rep(0.5,N)
phi = matrix(c(0.5,0.05,0.05,0.5),N,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.225,N)
B = matrix(0,N,N)
diag(B) = rep(0.5,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model2 = list(
  mu = mu,
  phi = phi,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)





N=3   #  dimension


mu = rep(0.5,N)
phi = matrix(c(0.5,0.05,0.05,0.5),N,N)
W = rep(0.1,N)
A = matrix(0.15,N,N)
diag(A) = rep(0.2,N)
B = matrix(0,N,N)
diag(B) = rep(0.4,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model2 = list(
  mu = mu,
  phi = phi,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)

################################################################################

###############    the IGARCH case for VAR(1)-GARCH(1,1) model   ###############

N=2   #  dimension

mu = rep(0.5,N)
phi = matrix(0.05,N,N)
diag(phi) = rep(0.5,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.6,N)
B = matrix(0,N,N)
diag(B) = rep(0.2,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model2 = list(
  mu = mu,
  phi = phi,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)




N=3   #  dimension

mu = rep(0.5,N)
phi = matrix(0.05,N,N)
diag(phi) = rep(0.5,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.4,N)
B = matrix(0,N,N)
diag(B) = rep(0.2,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model2 = list(
  mu = mu,
  phi = phi,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)


################################################################################


###############  second stationary for diagonal VAR(1)-GARCH(1,1) model   ######

N=2
# N=3

mu = rep(0.5,N)
phi = matrix(c(0.5,0.05,0.05,0.5),N,N)
W = rep(0.1,N)
A = matrix(0,N,N)
diag(A) = rep(0.4,N)
B = matrix(0,N,N)
diag(B) = rep(0.5,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)



model2 = list(
  mu = mu,
  phi = phi,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)


################################################################################



###############  the IGARCH case for diagonal VAR(1)-GARCH(1,1) model   ########

N=2
# N=3


mu = rep(0.5,N)
phi = matrix(0.05,N,N)
diag(phi) = rep(0.5,N)
W = rep(0.1,N)
A = matrix(0,N,N)
diag(A) = rep(0.7,N)
B = matrix(0,N,N)
diag(B) = rep(0.3,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)




model2 = list(
  mu = mu,
  phi = phi,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)


################################################################################


#######   model checking for second stationary VARMA(1,1)-GARCH(1,1) model  ####


N=2


mu = rep(0.5,N)
phi = matrix(c(0.2,0.05,0.05,0.2),N,N)
theta = matrix(c(0.7,0.05,0.05,0.7),N,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.225,N)
B = matrix(0,N,N)
diag(B) = rep(0.5,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model3 = list(
  mu = mu,
  phi = phi,
  theta = theta,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)




N=3


mu = rep(0.5,N)
phi = matrix(c(0.2,0.05,0.05,0.05,0.2,0.05,0.05,0.05,0.2),N,N)
theta = matrix(c(0.7,0.05,0.05,0.05,0.7,0.05,0.05,0.05,0.7),N,N)
W = rep(0.1,N)
A = matrix(0.15,N,N)
diag(A) = rep(0.2,N)
B = matrix(0,N,N)
diag(B) = rep(0.4,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model3 = list(
  mu = mu,
  phi = phi,
  theta = theta,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)



################################################################################



#######   model checking for the IGARCH case VARMA(1,1)-GARCH(1,1) model  ######


N=2


mu = rep(0.5,N)
phi = matrix(c(0.2,0.05,0.05,0.2),N,N)
theta = matrix(c(0.7,0.05,0.05,0.7),N,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.6,N)
B = matrix(0,N,N)
diag(B) = rep(0.2,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model3 = list(
  mu = mu,
  phi = phi,
  theta = theta,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)




N=3


mu = rep(0.5,N)
phi = matrix(c(0.2,0.05,0.05,0.05,0.2,0.05,0.05,0.05,0.2),N,N)
theta = matrix(c(0.7,0.05,0.05,0.05,0.7,0.05,0.05,0.05,0.7),N,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.4,N)
B = matrix(0,N,N)
diag(B) = rep(0.2,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model3 = list(
  mu = mu,
  phi = phi,
  theta = theta,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)


################################################################################


#######   model checking for second stationary VAR(2)-GARCH(1,1) model  ########


N=2


mu = rep(0.5,N)
phi1 = matrix(c(0.2,0.05,0.05,0.2),N,N)
phi2 = matrix(c(0.5,0.05,0.05,0.5),N,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.225,N)
B = matrix(0,N,N)
diag(B) = rep(0.5,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model4 = list(
  mu = mu,
  phi1 = phi1,
  phi2 = phi2,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)




N=3


mu = rep(0.5,N)
phi1 = matrix(c(0.2,0.05,0.05,0.05,0.2,0.05,0.05,0.05,0.2),N,N)
phi2 = matrix(c(0.5,0.05,0.05,0.05,0.5,0.05,0.05,0.05,0.5),N,N)
W = rep(0.1,N)
A = matrix(0.15,N,N)
diag(A) = rep(0.2,N)
B = matrix(0,N,N)
diag(B) = rep(0.4,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model4 = list(
  mu = mu,
  phi1 = phi1,
  phi2 = phi2,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)



################################################################################



########   model checking for the IGARCH case VAR(2)-GARCH(1,1) model  #########


N=2


mu = rep(0.5,N)
phi1 = matrix(c(0.2,0.05,0.05,0.2),N,N)
phi2 = matrix(c(0.5,0.05,0.05,0.5),N,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.6,N)
B = matrix(0,N,N)
diag(B) = rep(0.2,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model4 = list(
  mu = mu,
  phi1 = phi1,
  phi2 = phi2,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)




N=3


mu = rep(0.5,N)
phi1 = matrix(c(0.2,0.05,0.05,0.05,0.2,0.05,0.05,0.05,0.2),N,N)
phi2 = matrix(c(0.5,0.05,0.05,0.05,0.5,0.05,0.05,0.05,0.5),N,N)
W = rep(0.1,N)
A = matrix(0.2,N,N)
diag(A) = rep(0.4,N)
B = matrix(0,N,N)
diag(B) = rep(0.2,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)


model4 = list(
  mu = mu,
  phi1 = phi1,
  phi2 = phi2,
  W = W,
  A = A,
  B = B,
  Gamma = Gamma
)


################################################################################




#######   model checking for second stationary VAR(1)-GARCH(2,1) model  ########


N=2

mu = rep(0.5,N)
phi = matrix(0.05,N,N)
diag(phi) = rep(0.5,N)
W = rep(0.1,N)
A1 = matrix(0.05,N,N)
diag(A1) = rep(0.2,N)
A2 = matrix(0.05,N,N)
diag(A2) = rep(0.2,N)
B = matrix(0.05,N,N)
diag(B) = rep(0.4,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)



model4 = list(
  mu = mu,
  phi = phi,
  W = W,
  A1 = A1,
  A2 = A2,
  B = B,
  Gamma = Gamma
)





N=3

mu = rep(0.5,N)
phi = matrix(0.05,N,N)
diag(phi) = rep(0.5,N)
W = rep(0.1,N)
A1 = matrix(0.05,N,N)
diag(A1) = rep(0.1,N)
A2 = matrix(0.05,N,N)
diag(A2) = rep(0.1,N)
B = matrix(0.05,N,N)
diag(B) = rep(0.4,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)



model4 = list(
  mu = mu,
  phi = phi,
  W = W,
  A1 = A1,
  A2 = A2,
  B = B,
  Gamma = Gamma
)


################################################################################





#######   model checking for the IGARCH case VAR(1)-GARCH(2,1) model  ##########


N=2

mu = rep(0.5,N)
phi = matrix(0.05,N,N)
diag(phi) = rep(0.5,N)
W = rep(0.1,N)
A1 = matrix(0.05,N,N)
diag(A1) = rep(0.25,N)
A2 = matrix(0.05,N,N)
diag(A2) = rep(0.2,N)
B = matrix(0.05,N,N)
diag(B) = rep(0.4,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)



model4 = list(
  mu = mu,
  phi = phi,
  W = W,
  A1 = A1,
  A2 = A2,
  B = B,
  Gamma = Gamma
)





N=3

mu = rep(0.5,N)
phi = matrix(0.05,N,N)
diag(phi) = rep(0.5,N)
W = rep(0.1,N)
A1 = matrix(0.05,N,N)
diag(A1) = rep(0.1,N)
A2 = matrix(0.05,N,N)
diag(A2) = rep(0.2,N)
B = matrix(0.05,N,N)
diag(B) = rep(0.4,N)
Gamma = matrix(0.3,N,N)
diag(Gamma)=rep(1,N)



model4 = list(
  mu = mu,
  phi = phi,
  W = W,
  A1 = A1,
  A2 = A2,
  B = B,
  Gamma = Gamma
)


################################################################################



################################################################################


nobs = 1000               # sample size
skip = 1000               # skip size 
include.AR = T            # include AR part or not 
type="half-diagonal"      # type of the GARCH part: "diagonal", "half-diagonal", or "extended"
startval = NULL

set.seed(6)

replicates = 1000         # replication time







##########  type=="diagonal" only means GARCH part is diagonal   ###############


if(include.AR){
  if(type=="diagonal"){
    pars=matrix(0,replicates,(N*(N+4)+N*(N-1)/2))
    pars1=matrix(0,replicates,(N*(N+4)+N*(N-1)/2))
    AVs=matrix(0,replicates,(3*(N*(N+4)+N*(N-1)/2)))
    AVs1=matrix(0,replicates,(3*(N*(N+4)+N*(N-1)/2)))
    nums=rep(0,replicates)
  }else if(type=="half-diagonal"){
    pars=matrix(0,replicates,(N*(2*N+3)+N*(N-1)/2))
    pars1=matrix(0,replicates,(N*(2*N+3)+N*(N-1)/2))
    AVs=matrix(0,replicates,(3*(N*(2*N+3)+N*(N-1)/2)))
    AVs1=matrix(0,replicates,(3*(N*(2*N+3)+N*(N-1)/2)))
    nums=rep(0,replicates)
  }else{
    pars=matrix(0,replicates,(N*(3*N+2)+N*(N-1)/2))
    pars1=matrix(0,replicates,(N*(3*N+2)+N*(N-1)/2))
    AVs=matrix(0,replicates,(3*(N*(3*N+2)+N*(N-1)/2)))
    AVs1=matrix(0,replicates,(3*(N*(3*N+2)+N*(N-1)/2)))
    nums=rep(0,replicates)
  }
}else{
  if(type=="diagonal"){
    pars=matrix(0,replicates,(4*N+N*(N-1)/2))
    pars1=matrix(0,replicates,(4*N+N*(N-1)/2))
    AVs=matrix(0,replicates,(3*(4*N+N*(N-1)/2)))
    AVs1=matrix(0,replicates,(3*(4*N+N*(N-1)/2)))
    nums=rep(0,replicates)
  }else if(type=="half-diagonal"){
    pars=matrix(0,replicates,(N*(N+3)+N*(N-1)/2))
    pars1=matrix(0,replicates,(N*(N+3)+N*(N-1)/2))
    AVs=matrix(0,replicates,(3*(N*(N+3)+N*(N-1)/2)))
    AVs1=matrix(0,replicates,(3*(N*(N+3)+N*(N-1)/2)))
    nums=rep(0,replicates)
  }else{
    pars=matrix(0,replicates,(N*(2*N+2)+N*(N-1)/2))
    pars1=matrix(0,replicates,(N*(2*N+2)+N*(N-1)/2))
    AVs=matrix(0,replicates,(3*(N*(2*N+2)+N*(N-1)/2)))
    AVs1=matrix(0,replicates,(3*(N*(2*N+2)+N*(N-1)/2)))
    nums=rep(0,replicates)
  }
}



######################  the wald statistics  ###################################
# W11 - the ARMA part of self-weighted QMLE
# W22 - the ARMA part of local QMLE 
# W33 - the GARCH part of self-weighted QMLE
# W44 - the GARCH part of local QMLE 


W11=rep(0,replicates)
W22=rep(0,replicates)
W33=rep(0,replicates)
W44=rep(0,replicates)


######################  the model checking  ###################################
# M33 - the self-weighted QMLE
# M44 - the local QMLE 


M33=rep(0,replicates)
M44=rep(0,replicates)



#############  test whether the sum of the diagonal elements of (A+B) equals to 1 or not  ##############

if(include.AR){
  if(type=="diagonal"){
    tao=matrix(0,1,(N*(N+4)+N*(N-1)/2))
    tao[1,(N*(N+2)+1)]=1
    tao[1,(N*(N+3)+1)]=1
  }else if(type=="half-diagonal"){
    tao=matrix(0,1,(N*(2*N+3)+N*(N-1)/2))
    tao[1,(N*(N+2)+1)]=1
    tao[1,(2*N*(N+1)+1)]=1
  }else{
    tao=matrix(0,1,(N*(3*N+2)+N*(N-1)/2))
    tao[1,(N*(N+2)+1)]=1
    tao[1,(2*N*(N+1)+1)]=1
  }
}else{
  if(type=="diagonal"){
    tao=matrix(0,1,(4*N+N*(N-1)/2))
    tao[1,(2*N+1)]=1
    tao[1,(3*N+1)]=1
  }else if(type=="half-diagonal"){
    tao=matrix(0,1,(N*(N+3)+N*(N-1)/2))
    tao[1,(2*N+1)]=1
    tao[1,(N*(N+2)+1)]=1
  }else{
    tao=matrix(0,1,(2*N*(N+1)+N*(N-1)/2))
    tao[1,(2*N+1)]=1
    tao[1,(N*(N+2)+1)]=1
  }
}






count=0
i=1
while (i <= replicates)  {
  options(warn = 1)
  
  cat("====================================================================\n")
  cat("Round ", i, "\n")
  
  
  
  ###############   parameter   estimation    ##################################
  
  data <- CCCsim(model = model2, nobs = nobs, skip = skip, include.AR = include.AR, dist = "Gaussian")
  start.time = Sys.time()
  weight = weights(data=data)
  end.time = Sys.time()
  print(difftime(end.time, start.time, units="mins"))
  start.time = Sys.time()
  parEst1 = CCCfitSW4(data=data, w=weight, type=type, include.AR = include.AR)
  # print(as.vector(parEst1$pars))
  # print(parEst1$convergence)
  end.time = Sys.time()
  print(difftime(end.time, start.time, units="mins"))
  pars[i,] = parEst1$pars
  nums[i] = parEst1$convergence
  
  
  
  ######################   variance     ########################################
  
  
  start.time = Sys.time()
  if(parEst1$convergence==0){
    AV = AVSWMC(data=data, w=weight, type=type, par=parEst1$pars, include.AR=include.AR)
    parEst2 = Localfit(data=data, par=parEst1$pars, par1=parEst1$pars, type=type, include.AR=include.AR)
    # print(as.vector(parEst2$fit))
 
   
    j=1
    while (parEst2$error==3&j<=4) {
      initial  = parEst2$fit
      parEst2  = Localfit(data=data, par=initial, par1=parEst1$pars, type=type, include.AR=include.AR)
      print(as.vector(parEst2$fit))
      j=j+1
    }
    
    
    
    AV1   =  LocalVAMC(data=data, par=parEst2$fit, type=type, include.AR=include.AR, lag = 6)
    
    
    if(AV$error==0&(parEst2$error==0|parEst2$error==3)&AV1$error==0){
      # if(AV$error==0&parEst1$convergence==0){
      
      AVs[i,]   = as.vector(AV$vari)
      pars1[i,] = parEst2$fit
      AVs1[i,]  = as.vector(AV1$vari)
      
      W11[i] = nobs*parEst1$pars[5]*(AV$vari1[5,5])^{-1}*parEst1$pars[5]
      W22[i] = nobs*parEst2$fit[5]*(AV1$vari1[5,5])^{-1}*parEst2$fit[5]
      W33[i] = nobs*t(tao%*%parEst1$pars-1)*(tao%*%AV$vari1%*%t(tao))^{-1}*(tao%*%parEst1$pars-1)
      W44[i] = nobs*t(tao%*%parEst2$fit-1)*(tao%*%AV1$vari1%*%t(tao))^{-1}*(tao%*%parEst2$fit-1)

      M33[i] = AV$QMM
      M44[i] = AV1$QMM
      
      
      i=i+1
    }else{
      i=i
      count=count+1
    }
  }else{
    i=i
    count=count+1
  }
  end.time = Sys.time()
  print(difftime(end.time, start.time, units="mins"))
  

  
  
  
  cat("====================================================================\n")
  
}

par = colMeans(pars)
AV11 = colMeans(AVs)
par1 = colMeans(pars1)
AV22 = colMeans(AVs1)

if(include.AR){
  if(type=="diagonal"){
    
    ########################   self-weighted QMLE  #############################
    
    para.mat = p.match(par, include.AR, type, N)
    
    H.se   = AV11[1:(N*(N+4)+N*(N-1)/2)]
    out.se = AV11[(N*(N+4)+N*(N-1)/2+1):(2*(N*(N+4)+N*(N-1)/2))]
    rob.se = AV11[(2*(N*(N+4)+N*(N-1)/2)+1):(3*(N*(N+4)+N*(N-1)/2))]
    
    #########################   local QMLE  ####################################
    
    para.mat1 = p.match(par1, include.AR, type, N)
    
    H.se1   = AV22[1:(N*(N+4)+N*(N-1)/2)]
    out.se1 = AV22[(N*(N+4)+N*(N-1)/2+1):(2*(N*(N+4)+N*(N-1)/2))]
    rob.se1 = AV22[(2*(N*(N+4)+N*(N-1)/2)+1):(3*(N*(N+4)+N*(N-1)/2))]
    
  }else if(type=="half-diagonal"){
    
    ########################   self-weighted QMLE  #############################
    
    para.mat = p.match(par, include.AR, type, N)
    
    H.se   = AV11[1:(N*(2*N+3)+N*(N-1)/2)]
    out.se = AV11[(N*(2*N+3)+N*(N-1)/2+1):(2*(N*(2*N+3)+N*(N-1)/2))]
    rob.se = AV11[(2*(N*(2*N+3)+N*(N-1)/2)+1):(3*(N*(2*N+3)+N*(N-1)/2))]
    
    
    #################   local QMLE  ###############################
    para.mat1 = p.match(par1, include.AR, type, N)
    
    H.se1   = AV22[1:(N*(2*N+3)+N*(N-1)/2)]
    out.se1 = AV22[(N*(2*N+3)+N*(N-1)/2+1):(2*(N*(2*N+3)+N*(N-1)/2))]
    rob.se1 = AV22[(2*(N*(2*N+3)+N*(N-1)/2)+1):(3*(N*(2*N+3)+N*(N-1)/2))]
    
  }else{
    
    ########################   self-weighted QMLE  #############################
    
    para.mat = p.match(par, include.AR, type, N)
    
    H.se   = AV11[1:(N*(3*N+2)+N*(N-1)/2)]
    out.se = AV11[(N*(3*N+2)+N*(N-1)/2+1):(2*(N*(3*N+2)+N*(N-1)/2))]
    rob.se = AV11[(2*(N*(3*N+2)+N*(N-1)/2)+1):(3*(N*(3*N+2)+N*(N-1)/2))]
    
    
    #################   local QMLE  ###############################
    para.mat1 = p.match(par1, include.AR, type, N)
    
    H.se1   = AV22[1:(N*(3*N+2)+N*(N-1)/2)]
    out.se1 = AV22[(N*(3*N+2)+N*(N-1)/2+1):(2*(N*(3*N+2)+N*(N-1)/2))]
    rob.se1 = AV22[(2*(N*(3*N+2)+N*(N-1)/2)+1):(3*(N*(3*N+2)+N*(N-1)/2))]
    
  }
  
  ########################   self-weighted QMLE  #############################
  
  
  print(para.mat$mu)
  print(para.mat$phi)
  print(para.mat$W)
  print(para.mat$A)
  print(para.mat$B)
  print(para.mat$Gamma)
  
  
  ####################   local QMLE  #########################################
  
  print(para.mat1$mu)
  print(para.mat1$phi)
  print(para.mat1$W)
  print(para.mat1$A)
  print(para.mat1$B)
  print(para.mat1$Gamma)
  
  
  
}else{
  if(type=="diagonal"){
    
    ########################   self-weighted QMLE  #############################
    
    para.mat = p.match(par, include.AR, type, N)
    
    H.se   = AV11[1:(4*N+N*(N-1)/2)]
    out.se = AV11[(4*N+N*(N-1)/2+1):(2*(4*N+N*(N-1)/2))]
    rob.se = AV11[(2*(4*N+N*(N-1)/2)+1):(3*(4*N+N*(N-1)/2))]
    
    
    
    #######################   local  QMLE   ####################################
    
    para.mat1 = p.match(par1, include.AR, type, N)
    
    H.se1   = AV22[1:(4*N+N*(N-1)/2)]
    out.se1 = AV22[(4*N+N*(N-1)/2+1):(2*(4*N+N*(N-1)/2))]
    rob.se1 = AV22[(2*(4*N+N*(N-1)/2)+1):(3*(4*N+N*(N-1)/2))]
    
    
  }else if(type=="half-diagonal"){
    
    ########################   self-weighted QMLE  #############################
    
    para.mat = p.match(par, include.AR, type, N)
    
    H.se   = AV11[1:(N*(N+3)+N*(N-1)/2)]
    out.se = AV11[(N*(N+3)+N*(N-1)/2+1):(2*(N*(N+3)+N*(N-1)/2))]
    rob.se = AV11[(2*(N*(N+3)+N*(N-1)/2)+1):(3*(N*(N+3)+N*(N-1)/2))]
    
    
    ####################   local QMLE  #########################################
    
    
    para.mat1 = p.match(par1, include.AR, type, N)
    
    H.se1   = AV22[1:(N*(N+3)+N*(N-1)/2)]
    out.se1 = AV22[(N*(N+3)+N*(N-1)/2+1):(2*(N*(N+3)+N*(N-1)/2))]
    rob.se1 = AV22[(2*(N*(N+3)+N*(N-1)/2)+1):(3*(N*(N+3)+N*(N-1)/2))]
    
    
  }else{
    
    ########################   self-weighted QMLE  #############################
    
    para.mat = p.match(par, include.AR, type, N)
    
    H.se   = AV11[1:(N*(2*N+2)+N*(N-1)/2)]
    out.se = AV11[(N*(2*N+2)+N*(N-1)/2+1):(2*(N*(2*N+2)+N*(N-1)/2))]
    rob.se = AV11[(2*(N*(2*N+2)+N*(N-1)/2)+1):(3*(N*(2*N+2)+N*(N-1)/2))]
    
    
    ####################   local QMLE  #########################################
    
    
    para.mat1 = p.match(par1, include.AR, type, N)
    
    H.se1   = AV22[1:(N*(2*N+2)+N*(N-1)/2)]
    out.se1 = AV22[(N*(2*N+2)+N*(N-1)/2+1):(2*(N*(2*N+2)+N*(N-1)/2))]
    rob.se1 = AV22[(2*(N*(2*N+2)+N*(N-1)/2)+1):(3*(N*(2*N+2)+N*(N-1)/2))]
    
    
  }
  
  ########################   self-weighted QMLE  #############################
  
  
  print(para.mat$mu)
  print(para.mat$W)
  print(para.mat$A)
  print(para.mat$B)
  print(para.mat$Gamma)
  
  
  ####################   local QMLE  #########################################
  
  
  print(para.mat1$mu)
  print(para.mat1$W)
  print(para.mat1$A)
  print(para.mat1$B)
  print(para.mat1$Gamma)
}




SV = sqrt(diag(var(pars)))
# print(SV)
cbind(H.se,out.se,rob.se,SV)


SV1 = sqrt(diag(var(pars1)))
# print(SV1)
cbind(H.se1,out.se1,rob.se1,SV1)



print(para.mat$mu-rep(0.5,N))
print(para.mat$phi-matrix(c(0.5,0.05,0.05,0.5),nrow = 2,ncol=2))
print(para.mat$W-rep(0.1,N))
print(para.mat$A-matrix(c(0.225,0.2,0.2,0.225),nrow = 2,ncol=2))
print(para.mat$B-matrix(c(0.5,0,0,0.5),nrow = 2,ncol=2))
print(para.mat$Gamma-matrix(c(1,0.3,0.3,1),nrow = 2,ncol=2))
print(para.mat1$mu-rep(0.5,N))
print(para.mat1$phi-matrix(c(0.5,0.05,0.05,0.5),nrow = 2,ncol=2))
print(para.mat1$W-rep(0.1,N))
print(para.mat1$A-matrix(c(0.225,0.2,0.2,0.225),nrow = 2,ncol=2))
print(para.mat1$B-matrix(c(0.5,0,0,0.5),nrow = 2,ncol=2))
print(para.mat1$Gamma-matrix(c(1,0.3,0.3,1),nrow = 2,ncol=2))


sum(W11>3.84)*100/replicates
sum(W22>3.84)*100/replicates
sum(W33>3.84)*100/replicates
sum(W44>3.84)*100/replicates


sum(M33>12.592)*100/replicates
sum(M44>12.592)*100/replicates
print(count)




