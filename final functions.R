############################    data generating   ##########################################################################################

"CCCsim" <- function (model = list(),nobs = NULL,skip = NULL, dist = "Gaussian", include.AR = T)
{
  # Description:
  # Generates a series of VARMA-GARCH model specified by Ling & McAleer (2003)
  # using Gaussian innovations.
  # AR(1)-GARCH(1,1)
  
  # Parameters:
  # model -- a list with the model parameters as entries
  #     mu - a vector as the mean value of ARMA specification
  #     phi - a matrix of AR coefficients [phi1,phi2,...]
  #     W - a vector as the drift of GARCH specification
  #     A - a matrix of ARCH coefficients [A1,A2,...]
  #     B - a matrix of GARCH coefficients [B1,B2,...]
  #     Gamma - a covariance matrix of innovations
  # nobs -- number of observations generated
  # skip -- number of warm-up observations
  # dist -- Multivariate "Gaussian" or Student's "t" innovations
  # include.AR -- include AR or not
  
  ########## Passin the arguments ##########
  if (include.AR) {
    phi = model$phi
  }
  mu = model$mu
  W  = model$W
  A  = model$A
  B  = model$B
  Gamma = model$Gamma
  
  ########## Check parameters ##########
  
  # Model orders:
  N = dim(Gamma)[1] # dimension of vectors
  
  # Check if Gamma is a correlation matrix
  if (!identical(diag(Gamma),rep(1,N)))
    stop("Diagonal elements of Gamma are not all ones!")
  
  # Check if Gamma is a correlation matrix
  eigenGamma <- eigen(Gamma)$values
  if (max(abs(Gamma[lower.tri(Gamma)]))>1.0||min(eigenGamma)<0||!is.double(eigenGamma))
    stop("Gamma is not a real correlation matrix!")
  
  # # Check if the model is second-order stationary
  # if (max(abs(eigen(A+B)$values))>=1)
  #   stop("GARCH: Stationary condition is violated!")
  
  # Check if the model is strictly stationary
  if (max(abs(eigen(A+B)$values))>1)
    stop("GARCH: Stationary condition is violated!")
  
  
  # check if the AR coefficient is stationary
  if (include.AR){
    if (max(abs(eigen(phi)$values))>=1){
      stop("AR: Stationary condition is violated!")
    }
  }
  
  # # Check if the second-order moment exists
  # At = rbind(cbind(A,B),cbind(A,B))
  # if (max(abs(eigen(At)$values))>=1)
  #   stop("Moment condition is violated!")
  
  ########## Check parameters ##########
  
  # Generate observations
  nT = nobs+skip
  
  # Initial values
  r.ini = t(mu)   # using unconditional mean of r_t
  h.ini = matrix(0,1,N)
  # h.ini = t(solve(diag(N)-A-B)%*%W)   # using unconditional mean of H_t
  if (dist == "Gaussian"){
    etat = MASS::mvrnorm(nT+1,rep(0,N),Gamma)
  }else if (dist == "t"){
    df   = 5
    etat = t(t(MASS::mvrnorm(nT+1, rep(0,N), Gamma/df*(df-2)) * sqrt(df / rchisq(nT+1, df))))
  }else if (dist == "laplace"){
    etat = LaplacesDemon::rmvl(nT+1,rep(0,N),Gamma)
  }
  
  eta.ini = etat[1,]
  epsi.ini = sqrt(h.ini)*eta.ini
  
  # Simulate observations
  rt = rbind(r.ini,matrix(0,nT,N))
  ht = rbind(h.ini,matrix(0,nT,N))
  epsit = rbind(epsi.ini,matrix(0,nT,N))
  
  for (it in 2:(nT+1)){
    ht[it,] = W+A%*%(epsit[it-1,]**2)+B%*%ht[it-1,]
    epsit[it,] = sqrt(ht[it,])*etat[it,]
    if (include.AR){
      rt[it,] = mu+phi%*%(rt[it-1,]-mu)+epsit[it,]
    }else{
      rt[it,] = mu+epsit[it,]
    }
  }
  
  # skip the first "skip" points
  rts = rt[(2+skip):(nT+1),]
  CCCsim <- rts
}






"CCCsim1" <- function (model = list(),nobs = NULL,skip = NULL, dist = "Gaussian", include.AR = T)
{
  # Description:
  # Generates a series of VARMA-GARCH model specified by Ling & McAleer (2003)
  # using Gaussian innovations.
  # ARMA(1,1)-GARCH(1,1)
  
  # Parameters:
  # model -- a list with the model parameters as entries
  #     mu - a vector as the mean value of ARMA specification
  #     phi - a matrix of AR coefficients [phi1,phi2,...]
  #     theta - a matrix of MA coefficients [theta1,theta2,...]
  #     W - a vector as the drift of GARCH specification
  #     A - a matrix of ARCH coefficients [A1,A2,...]
  #     B - a matrix of GARCH coefficients [B1,B2,...]
  #     Gamma - a covariance matrix of innovations
  # nobs -- number of observations generated
  # skip -- number of warm-up observations
  # dist -- Multivariate "Gaussian" or Student's "t" innovations
  # include.AR -- include ARMA or not
  
  ########## Passin the arguments ##########
  if (include.AR) {
    phi = model$phi
    theta = model$theta
    
  }
  mu = model$mu
  W = model$W
  A = model$A
  B = model$B
  Gamma = model$Gamma
  
  ########## Check parameters ##########
  
  # Model orders:
  N = dim(Gamma)[1] # dimension of vectors
  
  # Check if Gamma is a correlation matrix
  if (!identical(diag(Gamma),rep(1,N)))
    stop("Diagonal elements of Gamma are not all ones!")
  
  # Check if Gamma is a correlation matrix
  eigenGamma <- eigen(Gamma)$values
  if (max(abs(Gamma[lower.tri(Gamma)]))>1.0||min(eigenGamma)<0||!is.double(eigenGamma))
    stop("Gamma is not a real correlation matrix!")
  
  # # Check if the model is second-order stationary
  # if (max(abs(eigen(A+B)$values))>=1)
  #   stop("GARCH: Stationary condition is violated!")
  
  # Check if the model is strictly stationary
  if (max(abs(eigen(A+B)$values))>1)
    stop("GARCH: Stationary condition is violated!")
  
  
  # check if the AR coefficient is stationary
  if (include.AR){
    if (max(abs(eigen(phi)$values))>=1){
      stop("AR: Stationary condition is violated!")
    }
  }
  
  
  # check if the MA coefficient is invertible
  if (include.AR){
    if (max(abs(eigen(theta)$values))>=1){
      stop("MA: Invertible condition is violated!")
    }
  }
  
  
  # # Check if the second-order moment exists
  # At = rbind(cbind(A,B),cbind(A,B))
  # if (max(abs(eigen(At)$values))>=1)
  #   stop("Moment condition is violated!")
  
  ########## Check parameters ##########
  
  # Generate observations
  nT = nobs+skip
  
  # Initial values
  r.ini = t(mu)   # using unconditional mean of r_t
  h.ini = matrix(0,1,N)
  # h.ini = t(solve(diag(N)-A-B)%*%W)   # using unconditional mean of H_t
  if (dist == "Gaussian"){
    etat = MASS::mvrnorm(nT+1,rep(0,N),Gamma)
  }else if (dist == "t"){
    df = 5
    etat = t(t(MASS::mvrnorm(nT+1, rep(0,N), Gamma/df*(df-2)) * sqrt(df / rchisq(nT+1, df))))
  }else if (dist == "laplace"){
    etat = LaplacesDemon::rmvl(nT+1,rep(0,N),Gamma)
  }
  
  eta.ini = etat[1,]
  epsi.ini = sqrt(h.ini)*eta.ini
  
  # Simulate observations
  rt = rbind(r.ini,matrix(0,nT,N))
  ht = rbind(h.ini,matrix(0,nT,N))
  epsit = rbind(epsi.ini,matrix(0,nT,N))
  
  for (it in 2:(nT+1)){
    ht[it,] = W+A%*%(epsit[it-1,]**2)+B%*%ht[it-1,]
    epsit[it,] = sqrt(ht[it,])*etat[it,]
    if (include.AR){
      rt[it,] = mu+phi%*%(rt[it-1,]-mu)+epsit[it,]+theta%*%(epsit[it-1,])
    }else{
      rt[it,] = mu+epsit[it,]
    }
  }
  
  # skip the first "skip" points
  rts = rt[(2+skip):(nT+1),]
  CCCsim <- rts
}




"CCCsim2" <- function (model = list(),nobs = NULL,skip = NULL, dist = "Gaussian", include.AR = T)
{
  # Description:
  # Generates a series of VARMA-GARCH model specified by Ling & McAleer (2003)
  # using Gaussian innovations.
  # AR(1)-GARCH(2,1)
  
  # Parameters:
  # model -- a list with the model parameters as entries
  #     mu - a vector as the mean value of ARMA specification
  #     phi - a matrix of AR coefficients [phi1,phi2,...]
  #     W - a vector as the drift of GARCH specification
  #     A - a matrix of ARCH coefficients [A1,A2,...]
  #     B - a matrix of GARCH coefficients [B1,B2,...]
  #     Gamma - a covariance matrix of innovations
  # nobs -- number of observations generated
  # skip -- number of warm-up observations
  # dist -- Multivariate "Gaussian" or Student's "t" innovations
  # include.AR -- include AR or not
  
  ########## Passin the arguments ##########
  if (include.AR) {
    phi = model$phi
  }
  mu = model$mu
  W = model$W
  A1 = model$A1
  A2 = model$A2
  B = model$B
  Gamma = model$Gamma
  
  ########## Check parameters ##########
  
  # Model orders:
  N = dim(Gamma)[1] # dimension of vectors
  
  # Check if Gamma is a correlation matrix
  if (!identical(diag(Gamma),rep(1,N)))
    stop("Diagonal elements of Gamma are not all ones!")
  
  # Check if Gamma is a correlation matrix
  eigenGamma <- eigen(Gamma)$values
  if (max(abs(Gamma[lower.tri(Gamma)]))>1.0||min(eigenGamma)<0||!is.double(eigenGamma))
    stop("Gamma is not a real correlation matrix!")
  
  # # Check if the model is second-order stationary
  # if (max(abs(eigen(A1+A2+B)$values))>=1)
  #   stop("GARCH: Stationary condition is violated!")
  
  # Check if the model is strictly stationary
  MatAAB=rbind(cbind((A1+B),A2),cbind(diag(N),matrix(0,N,N)))
  if (max(abs(eigen(MatAAB)$values))>1)
    stop("GARCH: Stationary condition is violated!")
  
  
  # check if the AR coefficient is stationary
  if (include.AR){
    if (max(abs(eigen(phi)$values))>=1){
      stop("AR: Stationary condition is violated!")
    }
  }
  
  
  ########## Check parameters ##########
  
  # Generate observations
  nT = nobs+skip
  
  # Initial values
  r.inii = t(mu)   # using unconditional mean of r_t
  r.ini  = t(mu)   # using unconditional mean of r_t
  h.inii = matrix(0,1,N)
  h.ini  = matrix(0,1,N)
  # h.inii = t(solve(diag(N)-A1-A2-B)%*%W)   # using unconditional mean of H_t
  # h.ini = t(solve(diag(N)-A1-A2-B)%*%W)   # using unconditional mean of H_t
  if (dist == "Gaussian"){
    etat = MASS::mvrnorm(nT+2,rep(0,N),Gamma)
  }else if (dist == "t"){
    df = 5
    etat = t(t(MASS::mvrnorm(nT+2, rep(0,N), Gamma/df*(df-2)) * sqrt(df / rchisq(nT+2, df))))
  }else if (dist == "laplace"){
    etat = LaplacesDemon::rmvl(nT+2,rep(0,N),Gamma)
  }
  
  eta.inii = etat[1,]
  epsi.inii = sqrt(h.inii)*eta.inii
  eta.ini = etat[2,]
  epsi.ini = sqrt(h.ini)*eta.ini
  
  # Simulate observations
  rt = rbind(r.inii,r.ini,matrix(0,nT,N))
  ht = rbind(h.inii,h.ini,matrix(0,nT,N))
  epsit = rbind(epsi.inii,epsi.ini,matrix(0,nT,N))
  
  for (it in 3:(nT+2)){
    ht[it,] = W+A1%*%(epsit[it-1,]**2)+A2%*%(epsit[it-2,]**2)+B%*%ht[it-1,]
    epsit[it,] = sqrt(ht[it,])*etat[it,]
    if (include.AR){
      rt[it,] = mu+phi%*%(rt[it-1,]-mu)+epsit[it,]
    }else{
      rt[it,] = mu+epsit[it,]
    }
  }
  
  # skip the first "skip" points
  rts = rt[(3+skip):(nT+2),]
  CCCsim <- rts
}






"CCCsim3" <- function (model = list(),nobs = NULL,skip = NULL, dist = "Gaussian", include.AR = T)
{
  # Description:
  # Generates a series of VARMA-GARCH model specified by Ling & McAleer (2003)
  # using Gaussian innovations.
  # AR(2)-GARCH(1,1)
  
  # Parameters:
  # model -- a list with the model parameters as entries
  #     mu - a vector as the mean value of ARMA specification
  #     phi1,phi2 - a matrix of AR coefficients [phi1,phi2,...]
  #     W - a vector as the drift of GARCH specification
  #     A - a matrix of ARCH coefficients [A1,A2,...]
  #     B - a matrix of GARCH coefficients [B1,B2,...]
  #     Gamma - a covariance matrix of innovations
  # nobs -- number of observations generated
  # skip -- number of warm-up observations
  # dist -- Multivariate "Gaussian" or Student's "t" innovations
  # include.AR -- include ARMA or not
  
  ########## Passin the arguments ##########
  if (include.AR) {
    phi1 = model$phi1
    phi2 = model$phi2
    
  }
  mu = model$mu
  W = model$W
  A = model$A
  B = model$B
  Gamma = model$Gamma
  
  ########## Check parameters ##########
  
  # Model orders:
  N = dim(Gamma)[1] # dimension of vectors
  
  # Check if Gamma is a correlation matrix
  if (!identical(diag(Gamma),rep(1,N)))
    stop("Diagonal elements of Gamma are not all ones!")
  
  # Check if Gamma is a correlation matrix
  eigenGamma <- eigen(Gamma)$values
  if (max(abs(Gamma[lower.tri(Gamma)]))>1.0||min(eigenGamma)<0||!is.double(eigenGamma))
    stop("Gamma is not a real correlation matrix!")
  
  # # Check if the model is second-order stationary
  # if (max(abs(eigen(A+B)$values))>=1)
  #   stop("GARCH: Stationary condition is violated!")
  
  # Check if the model is strictly stationary
  if (max(abs(eigen(A+B)$values))>1)
    stop("GARCH: Stationary condition is violated!")
  
  
  # check if the AR coefficient is stationary
  if (include.AR){
    if (max(abs(eigen(phi1+phi2)$values))>=1){
      stop("AR: Stationary condition is violated!")
    }
  }
  
  
  
  # # Check if the second-order moment exists
  # At = rbind(cbind(A,B),cbind(A,B))
  # if (max(abs(eigen(At)$values))>=1)
  #   stop("Moment condition is violated!")
  
  ########## Check parameters ##########
  
  # Generate observations
  nT = nobs+skip
  
  # Initial values
  r.inii = t(mu)   # using unconditional mean of r_t
  r.ini  = t(mu)   # using unconditional mean of r_t
  h.inii = matrix(0,1,N)
  h.ini  = matrix(0,1,N)
  # h.inii = t(solve(diag(N)-A-B)%*%W)   # using unconditional mean of H_t
  # h.ini = t(solve(diag(N)-A-B)%*%W)   # using unconditional mean of H_t
  if (dist == "Gaussian"){
    etat = MASS::mvrnorm(nT+2,rep(0,N),Gamma)
  }else if (dist == "t"){
    df = 5
    etat = t(t(MASS::mvrnorm(nT+2, rep(0,N), Gamma/df*(df-2)) * sqrt(df / rchisq(nT+2, df))))
  }else if (dist == "laplace"){
    etat = LaplacesDemon::rmvl(nT+2,rep(0,N),Gamma)
  }
  
  eta.inii = etat[1,]
  epsi.inii = sqrt(h.inii)*eta.inii
  eta.ini = etat[2,]
  epsi.ini = sqrt(h.ini)*eta.ini
  
  # Simulate observations
  rt = rbind(r.inii,r.ini,matrix(0,nT,N))
  ht = rbind(h.inii,h.ini,matrix(0,nT,N))
  epsit = rbind(epsi.inii,epsi.ini,matrix(0,nT,N))
  
  for (it in 3:(nT+2)){
    ht[it,] = W+A%*%(epsit[it-1,]**2)+B%*%ht[it-1,]
    epsit[it,] = sqrt(ht[it,])*etat[it,]
    if (include.AR){
      rt[it,] = mu+phi1%*%(rt[it-1,]-mu)+phi2%*%(rt[it-2,]-mu)+epsit[it,]
    }else{
      rt[it,] = mu+epsit[it,]
    }
  }
  
  # skip the first "skip" points
  rts = rt[(3+skip):(nT+2),]
  CCCsim <- rts
}



############################################################################################################################################



############################   weight    ###################################################################################################


"weights" <- function(data)
{
  nT    = dim(data)[1]
  N     = dim(data)[2]
  data  = matrix(as.numeric(data),nrow=nT,ncol=N)
  data1 = rep(0,nT)
  for (t in 1:nT) {
    data1[t] = norm(data[t,,drop=FALSE],"F")
  }
  C = quantile(data1,0.9)
  m = 1:(nT-1)
  m = m^(-9)
  # m[1]=2
  de = as.numeric(data1 > C)
  w  = rep(1,nT)
  S  = rep(0,(nT-1))
  for (i in 2:nT) {
    S[i-1] = sum(rev(data1[1:(i-1)]*de[1:(i-1)])*m[1:(i-1)])
    w[i]   = (max(1,S[i-1]/C))^(-4)
  }
  
  # w = rep(1,nT)
  return(w)
}



############################################################################################################################################


# making a symmetric matrix positive definite. From r-help 2003.12.27
# except normalization. 

make.pd <- function(x, tol=1e-6) {
  eig <- eigen(x, symmetric=TRUE)
  rtol <- tol * eig$values[1]
  if(min(eig$values) < rtol) {
    vals <- eig$values
    vals[vals < rtol] <- rtol
    srev <- eig$vectors %*% (vals * t(eig$vectors))
    dimnames(srev) <- dimnames(x)
    srev <- diag(1/sqrt(diag(srev)))%*%srev%*%diag(1/sqrt(diag(srev)))
    return(srev) 
  } else {return(x)}
}





# constructing parameter vectors and matrices
# input: parameters, include.AR, type and dimension
# "diagonal": A and B are diagonal
# "half-diagonal": A is full and B is diagonal
# "extended": A and B are full



p.match <- function(par, include.AR, type, ndim){
  if(include.AR){
    if(type=="diagonal"){
      mu  = par[1:ndim]
      phi = matrix(par[(ndim+1):(ndim*(ndim+1))],ndim,ndim)
      W   = par[(ndim*(ndim+1)+1):(ndim*(ndim+2))]
      A   = diag(par[(ndim*(ndim+2)+1):(ndim*(ndim+3))])
      B   = diag(par[(ndim*(ndim+3)+1):(ndim*(ndim+4))])
      Gammas = par[(ndim*(ndim+4)+1):(ndim*(ndim+4)+ndim*(ndim-1)/2)]
    }else if(type=="half-diagonal"){
      mu  = par[1:ndim]
      phi = matrix(par[(ndim+1):(ndim*(ndim+1))],ndim,ndim)
      W   = par[(ndim*(ndim+1)+1):(ndim*(ndim+2))]
      A   = matrix(par[(ndim*(ndim+2)+1):(2*ndim*(ndim+1))],ndim,ndim)
      B   = diag(par[(2*ndim*(ndim+1)+1):(ndim*(2*ndim+3))])
      Gammas = par[(ndim*(2*ndim+3)+1):(ndim*(2*ndim+3)+ndim*(ndim-1)/2)]
    }else{
      mu  = par[1:ndim]
      phi = matrix(par[(ndim+1):(ndim*(ndim+1))],ndim,ndim)
      W   = par[(ndim*(ndim+1)+1):(ndim*(ndim+2))]
      A   = matrix(par[(ndim*(ndim+2)+1):(2*ndim*(ndim+1))],ndim,ndim)
      B   = matrix(par[(2*ndim*(ndim+1)+1):(ndim*(3*ndim+2))],ndim,ndim)
      Gammas = par[(ndim*(3*ndim+2)+1):(ndim*(3*ndim+2)+ndim*(ndim-1)/2)]
    }
    Gamma                      = matrix(0, ndim, ndim)
    Gamma[lower.tri(Gamma)]    = Gammas
    Gamma                      = t(Gamma)
    Gamma[lower.tri(Gamma)]    = Gammas
    diag(Gamma)                = rep(1,ndim)
    
    list(mu=mu, phi=phi, W=W, A=A, B=B, Gamma=Gamma)
  }else{
    if(type=="diagonal"){
      mu  = par[1:ndim]
      W   = par[(ndim+1):(ndim*2)]
      A   = diag(par[(ndim*2+1):(3*ndim)])
      B   = diag(par[(3*ndim+1):(4*ndim)])
      Gammas = par[(4*ndim+1):(4*ndim+ndim*(ndim-1)/2)]
    }else if(type=="half-diagonal"){
      mu  = par[1:ndim]
      W   = par[(ndim+1):(ndim*2)]
      A   = matrix(par[(ndim*2+1):(ndim*(ndim+2))],ndim,ndim)
      B   = diag(par[(ndim*(ndim+2)+1):(ndim*(ndim+3))])
      Gammas = par[(ndim*(ndim+3)+1):(ndim*(ndim+3)+ndim*(ndim-1)/2)]
    }else{
      mu  = par[1:ndim]
      W   = par[(ndim+1):(ndim*2)]
      A   = matrix(par[(ndim*2+1):(ndim*(ndim+2))],ndim,ndim)
      B   = matrix(par[(ndim*(ndim+2)+1):(ndim*(2*ndim+2))],ndim,ndim)
      Gammas = par[(ndim*(2*ndim+2)+1):(ndim*(2*ndim+2)+ndim*(ndim-1)/2)]
    }
    Gamma                   = matrix(0, ndim, ndim)
    Gamma[lower.tri(Gamma)] = Gammas
    Gamma                   = t(Gamma)
    Gamma[lower.tri(Gamma)] = Gammas
    diag(Gamma)             = rep(1,ndim)
    
    list(mu=mu, W=W, A=A, B=B, Gamma=Gamma)
  }
}





############################    estimation    ##############################################################################################
############################    the true values are the initial values  ####################################################################

############################    This is for the simulation   ###############################################################################



"CCCfitSW4" <- function (data, w = weight, include.AR = T, type="extended", startval = NULL, solver = "Solnp")
{
  ################################################################################
  ##                                                                            ##
  ## Description:                                                               ##
  ##    QMLE for VARMA-GARCH model specified by Ling & McAleer (2003)           ##
  ##    Specified as AR(1)-GARCH(1,1) or mu-GARCH(1,1)                          ##
  ##                                                                            ##
  ## Arguments:                                                                 ##
  ##    data -- a nT-by-N matrix of data                                        ##
  ##    include.AR -- include AR or not                                         ##
  ##    solver -- optimization solvers for QMLE, including R function "Solnp"   ##
  ##    from R package "Rsolnp"                                                 ##
  ##    type -- only GARCH part is "diagonal", "half-diagonal", or "extended"   ##
  ##                                                                            ##
  ################################################################################
  
  
  # ---------------------------------------------------------------------------
  # Passin the arguments
  
  nT = dim(data)[1]
  N  = dim(data)[2]
  
  
  # ---------------------------------------------------------------------------
  # Starting Values
  
  # # Conditional mean coefficients (self-weighted LSE)
  # if (include.AR){
  #   ist     = 2
  #   ne      = nT-1
  #   xmtx    = matrix(1,ne,1)
  #   ymtx    = as.matrix(data[ist:nT,])
  #   xmtx    = cbind(xmtx,data[(ist-1):(nT-1),])
  #   xmtx    = as.matrix(xmtx)
  #   xmtxt   = t(xmtx)
  #   xtx     = matrix(0,(N+1),(N+1))
  #   xty     = matrix(0,(N+1),N)
  #   for (i in 1:(nT-1)) {
  #     xtx = xtx + w[i+1]*xmtxt[,i,drop=FALSE]%*%xmtx[i,,drop=FALSE]
  #     xty = xty + w[i+1]*(xmtxt[,i,drop=FALSE]%*%ymtx[i,,drop=FALSE])
  #   }
  #   xtxinv  = solve(xtx)
  #   ARcoff  = t(xtxinv%*%xty)
  #   mu.ini  = ARcoff[,1]
  #   phi.ini = ARcoff[,2:(N+1)]
  #   mu.ini  = solve(diag(N)-phi.ini)%*%mu.ini
  # }else{
  #   mu.ini  = colMeans(data)
  # }
  
  
  
  #  typical values as initial values
  
  if(include.AR){
    mu.ini        = rep(0.5,N)
    phi.ini       = matrix(0.05,N,N)
    diag(phi.ini) = rep(0.5,N)
  }else{
    mu.ini        = rep(0.5,N)
  }
  
  # (LSE)
  # if (include.AR){
  #   ist     = 2
  #   ne      = nT-1
  #   xmtx    = matrix(1,ne,1)
  #   ymtx    = as.matrix(data[ist:nT,])
  #   xmtx    = cbind(xmtx,data[(ist-1):(nT-1),])
  #   xmtx    = as.matrix(xmtx)
  #   xtx     = crossprod(xmtx,xmtx)
  #   xty     = crossprod(xmtx,ymtx)
  #   xtxinv  = solve(xtx)
  #   ARcoff  = t(xtxinv%*%xty)
  #   mu.ini  = ARcoff[,1]
  #   phi.ini = ARcoff[,2:(N+1)]
  #   mu.ini  = solve(diag(N)-phi.ini)%*%mu.ini
  # }else{
  #   mu.ini  = colMeans(data)
  # }
  
  
  #   Conditional variance coefficients
  if(type=="diagonal"){
    W.ini = rep(0.1,N)
    A.ini = rep(0.7,N)
    B.ini = rep(0.3,N)
    
    # #   Correlation coefficients
    # r.ini    = matrix(0,1,N)
    # epsi.ini = matrix(0,1,N)
    # rt       = rbind(r.ini,data)
    # epsit    = rbind(epsi.ini,matrix(0,nT,N)) 
    # 
    # for (it in 1:nT){
    #   if (include.AR) {
    #     epsit[it+1,] = (rt[it+1,]-mu.ini)-phi.ini%*%(rt[it,]-mu.ini)
    #   }else{
    #     epsit[it+1,] = rt[it+1,]-mu.ini
    #   }
    # }
  }else if(type=="half-diagonal"){
    W.ini       = rep(0.1,N)
    A.ini       = matrix(0.2,N,N)
    diag(A.ini) = rep(0.6,N)
    B.ini       = rep(0.2,N)
    
    # #   Correlation coefficients
    # r.ini    = matrix(0,1,N)
    # epsi.ini = matrix(0,1,N)
    # rt       = rbind(r.ini,data)
    # epsit    = rbind(epsi.ini,matrix(0,nT,N)) 
    # 
    # for (it in 1:nT){
    #   if (include.AR) {
    #     epsit[it+1,] = (rt[it+1,]-mu.ini)-phi.ini%*%(rt[it,]-mu.ini)
    #   }else{
    #     epsit[it+1,] = rt[it+1,]-mu.ini
    #   }
    # }
    
  }else{
    W.ini       = rep(0.1,N)
    A.ini       = matrix(0.05,N,N)
    diag(A.ini) = rep(0.2,N)
    B.ini       = matrix(0.05,N,N)
    diag(B.ini) = rep(0.7,N) # Typical values as starting values
    
    # #   Correlation coefficients
    # r.ini    = matrix(0,1,N)
    # epsi.ini = matrix(0,1,N)
    # rt       = rbind(r.ini,data)
    # epsit    = rbind(epsi.ini,matrix(0,nT,N)) 
    # 
    # for (it in 1:nT){
    #   if (include.AR) {
    #     epsit[it+1,] = (rt[it+1,]-mu.ini)-phi.ini%*%(rt[it,]-mu.ini)
    #   }else{
    #     epsit[it+1,] = rt[it+1,]-mu.ini
    #   }
    # }
  }
  # epsit      = epsit[2:(nT+1),]
  # Gamma.ini  = cor(epsit)
  # Gammas.ini = Gamma.ini[lower.tri(Gamma.ini)] # Correlation matrix of epsi as starting values
  
  Gammas.ini = rep(0.3,N*(N-1)/2) # typical values as starting values
  
  #   Combining coefficients as starting value for optimization
  if (include.AR){
    if (length(startval) == 0){
      x0 = c(as.vector(cbind(mu.ini,phi.ini,W.ini,A.ini,B.ini)),Gammas.ini)
    }else{
      x0 = c(as.vector(cbind(mu.ini,phi.ini,startval)),Gammas.ini)
    }
  }else{
    if (length(startval) == 0){
      x0 = c(as.vector(cbind(mu.ini,W.ini,A.ini,B.ini)),Gammas.ini)
    }else{
      x0 = c(as.vector(cbind(mu.ini,startval)),Gammas.ini)
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # Negative conditional log likelihood function
  
  
  LLT <- function (par){
    if (include.AR) {
      para.mat = p.match(par, include.AR, type, N)
      heps     = .Call("garch_ar", data, para.mat$mu, para.mat$phi, para.mat$W, para.mat$A, para.mat$B)
    }else{
      para.mat = p.match(par, include.AR, type, N)
      heps     = .Call("garch_mu", data, para.mat$mu, para.mat$W, para.mat$A, para.mat$B)
    }
    
    # check if R is positive definite. If not, making R positive definite
    para.mat$Gamma <- make.pd(para.mat$Gamma)
    
    
    h      = heps[[1]]
    eta    = heps[[2]]/sqrt(h)
    logh   = log(h)
    dlogh  = rowSums(logh)
    lndetR = log(det(para.mat$Gamma))
    invR   = solve(para.mat$Gamma)
    llt    = 0.5*N*log(2*pi)*sum(w)+0.5*sum(w*dlogh)+0.5*lndetR*sum(w)+0.5*sum(w*((eta%*%invR)*eta))
    # llt    = 0.5*N*nT*log(2*pi)+0.5*sum(log(h))+0.5*nT*lndetR+0.5*sum((eta%*%invR)*eta) (log-likelihood function)
    llt    = llt/nT
    return(llt)
  }
  
  # ---------------------------------------------------------------------------
  # Nonlinear constraints to shape the parameter space
  # A small value is added to avoid the "solnp" takes the values on the boundary
  # 10^-8 is chosen since the feasibility tolerance of "solnp" is 10^-8 by default
  
  LLTineq <- function (par){
    para.mat = p.match(par, include.AR, type, N)
    c1 = max(abs(eigen(para.mat$B)$values))-(1-1e-8)
    c2 = -max(abs(eigen(para.mat$Gamma)$values))+1e-8 # Spectral radius of Gamma has a positive lower bound
    c3 = -min(eigen(para.mat$Gamma)$values)+1e-8 # Gamma is positive definite
    if (include.AR){
      c4 = max(abs(eigen(para.mat$phi)$values))-(1-1e-8) # AR stationarity
      return(c(c1,c2,c3,c4))
    }else{
      return(c(c1,c2,c3))
    }
  }
  
  # ---------------------------------------------------------------------------
  # Optimization
  
  if (include.AR) {
    if(type=="diagonal"){
      lb     = c(rep(-Inf, N*(N+1)), rep(1e-8, N), rep(1e-8, 2*N), rep(-1, (N-1)*N/2))     # Positive lower and upper bounds for W
      ub     = c(rep(Inf, N*(N+1)), rep(0.3, N),  rep(Inf, 2*N), rep(1, (N-1)*N/2))
    }else if(type=="half-diagonal"){
      lb     = c(rep(-Inf, N*(N+1)), rep(1e-8, N), rep(1e-8, (N*N+N)), rep(-1, (N-1)*N/2)) # Positive lower and upper bounds for W
      ub     = c(rep(Inf, N*(N+1)), rep(0.3, N),  rep(Inf, (N*N+N)), rep(1, (N-1)*N/2))
    }else{
      lb     = c(rep(-Inf, N*(N+1)), rep(1e-8, N), rep(1e-8, 2*N*N), rep(-1, (N-1)*N/2))   # Positive lower and upper bounds for W
      ub     = c(rep(Inf, N*(N+1)), rep(0.3, N),  rep(Inf, 2*N*N), rep(1, (N-1)*N/2))
    }
    lbineq = c(-Inf, -Inf, -Inf, -Inf)
    ubineq = c(0, 0, 0, 0)
  }else{
    if(type=="diagonal"){
      lb     = c(rep(-Inf, N), rep(1e-8, N), rep(1e-8, 2*N), rep(-1, (N-1)*N/2))     # Positive lower and upper bounds for W
      ub     = c(rep(Inf, N), rep(0.3, N),  rep(Inf, 2*N), rep(1, (N-1)*N/2))
    }else if(type=="half-diagonal"){
      lb     = c(rep(-Inf, N), rep(1e-8, N), rep(1e-8, (N*N+N)), rep(-1, (N-1)*N/2)) # Positive lower and upper bounds for W
      ub     = c(rep(Inf, N), rep(0.3, N),  rep(Inf, (N*N+N)), rep(1, (N-1)*N/2))
    }else{
      lb     = c(rep(-Inf, N), rep(1e-8, N), rep(1e-8, 2*N*N), rep(-1, (N-1)*N/2))   # Positive lower and upper bounds for W
      ub     = c(rep(Inf, N), rep(0.3, N),  rep(Inf, 2*N*N), rep(1, (N-1)*N/2))
    }
    lbineq = c(-Inf, -Inf, -Inf)
    ubineq = c(0, 0, 0)
  }
  
  opts = list("trace" = 0,
              # "delta"=1e-3,
              "rho"=100)
  
  est  = solnp(x0, LLT, ineqfun = LLTineq, ineqLB = lbineq, ineqUB = ubineq, LB = lb, UB = ub, control = opts)
  return(est)
}




#############################   This is for the application   ##############################################################################



"CCCfit" <- function (data, w = weight, include.AR = T, type = extended, startval = NULL, solver = "solnp", tr = 1)
{
  ################################################################################
  ##                                                                            ##
  ## Description:                                                               ##
  ##    QMLE for VARMA-GARCH model specified by Ling & McAleer (2003)           ##
  ##    Specified as AR(1)-GARCH(1,1) or mu-GARCH(1,1)                          ##
  ##                                                                            ##
  ## Arguments:                                                                 ##
  ##    data -- a nT-by-N matrix of data                                        ##
  ##    include.AR -- include AR or not                                         ##
  ##    startval -- a N-by-N(2N+1) matrix as starting values for [W, A, B]      ##
  ##    solver -- optimization solvers for QMLE, including R function "solnp"   ##
  ##    and "gosolnp", both from R package "Rsolnp"                             ##
  ##    type -- only GARCH part is "diagonal", "half-diagonal", or "extended"   ##
  ##                                                                            ##
  ################################################################################
  
  # ---------------------------------------------------------------------------
  # Passin the arguments
  
  nT = dim(data)[1]
  N  = dim(data)[2]
  
  # ---------------------------------------------------------------------------
  # Starting Values
  
  # Conditional mean coefficients (self-weighted LSE)
  if (include.AR){
    ist     = 2
    ne      = nT-1
    xmtx    = matrix(1,ne,1)
    ymtx    = as.matrix(data[ist:nT,])
    xmtx    = cbind(xmtx,data[(ist-1):ne,])
    xmtx    = as.matrix(xmtx)
    xmtxt   = t(xmtx)
    xtx     = matrix(0,(N+1),(N+1))
    xty     = matrix(0,(N+1),N)
    for (i in 1:(nT-1)) {
      xtx = xtx + w[i+1]*xmtxt[,i,drop=FALSE]%*%xmtx[i,,drop=FALSE]
      xty = xty + w[i+1]*(xmtxt[,i,drop=FALSE]%*%ymtx[i,,drop=FALSE])
    }
    xtxinv  = solve(xtx)
    ARcoff  = t(xtxinv%*%xty)
    mu.ini  = ARcoff[,1]
    phi.ini = ARcoff[,2:(N+1)]
    mu.ini  = solve(diag(N)-phi.ini)%*%mu.ini
  }else{
    mu.ini  = colMeans(data)
  }
  
  
  
  #   Conditional variance coefficients
  W.ini       = rep(0.1,N)
  A.ini       = rep(0.1,N)   # if A is diagonal matrix
  # A.ini       = rep(0.1,N+1)   #  m=2, and A is semi-diagonal
  B.ini       = rep(0.6,N)   # if B is diagonal matrix
  # diag(B.ini) = rep(0.6,N) # Typical values as starting values
  
  #   Correlation coefficients
  r.ini    = matrix(0,1,N)
  epsi.ini = matrix(0,1,N)
  rt       = rbind(r.ini,data)
  epsit    = rbind(epsi.ini,matrix(0,nT,N))
  
  for (it in 1:nT){
    if (include.AR) {
      epsit[it+1,] = (rt[it+1,]-mu.ini)-phi.ini%*%(rt[it,]-mu.ini)
    }else{
      epsit[it+1,] = rt[it+1,]-mu.ini
    }
  }
  epsit      = epsit[2:(nT+1),]
  Gamma.ini  = cor(epsit)
  Gammas.ini = Gamma.ini[lower.tri(Gamma.ini)] # Correlation matrix of epsi as starting values
  
  
  #   Combining coefficients as starting value for optimization
  if (include.AR){
    if (length(startval) == 0){
      x0 = c(as.vector(cbind(mu.ini,phi.ini,W.ini,A.ini,B.ini)),Gammas.ini)
    }else{
      x0 = c(as.vector(cbind(mu.ini,phi.ini,startval)),Gammas.ini)
    }
  }else{
    if (length(startval) == 0){
      x0 = c(as.vector(cbind(mu.ini,W.ini,A.ini,B.ini)),Gammas.ini)
    }else{
      x0 = c(as.vector(cbind(mu.ini,startval)),Gammas.ini)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Negative conditional log likelihood function
  
  LLT <- function (par){
    if (include.AR) {
      para.mat = p.match(par, include.AR, type, N)
      heps = .Call("garch_ar", data, para.mat$mu, para.mat$phi, para.mat$W, para.mat$A, para.mat$B)
    }else{
      para.mat = p.match(par, include.AR, type, N)
      heps = .Call("garch_mu", data, para.mat$mu, para.mat$W, para.mat$A, para.mat$B)
    }
    
    # check if R is positive definite. If not, making R positive definite
    para.mat$Gamma <- make.pd(para.mat$Gamma)
    
    
    h      = heps[[1]]
    eta    = heps[[2]]/sqrt(h)
    logh   = log(h)
    dlogh  = rowSums(logh)
    lndetR = log(det(para.mat$Gamma))
    invR   = solve(para.mat$Gamma)
    llt    = 0.5*N*log(2*pi)*sum(w)+0.5*sum(w*dlogh)+0.5*lndetR*sum(w)+0.5*sum(w*((eta%*%invR)*eta))
    # llt    = 0.5*N*nT*log(2*pi)+0.5*sum(log(h))+0.5*nT*lndetR+0.5*sum((eta%*%invR)*eta) (log-likelihood function)
    llt    = llt/nT
    return(llt)
  }
  
  # ---------------------------------------------------------------------------
  # Nonlinear constraints to shape the parameter space
  
  LLTineq <- function (par){
    para.mat = p.match(par, include.AR, type, N)
    c1 = max(abs(eigen(para.mat$B)$values))-(1-1e-8)
    c2 = -max(abs(eigen(para.mat$Gamma)$values))+1e-8 # Spectral radius of Gamma has a positive lower bound
    c3 = -min(eigen(para.mat$Gamma)$values)+1e-8 # Gamma is positive definite
    if (include.AR){
      c4 = max(abs(eigen(para.mat$phi)$values))-(1-1e-8) # AR stationarity
      return(c(c1,c2,c3,c4))
    }else{
      return(c(c1,c2,c3))
    }
  }
  
  # ---------------------------------------------------------------------------
  # Optimization
  
  
  
  
  meanW       = rep(1e-4, N)
  sdW         = rep(1e-2, N)
  meanA       = rep(0.1,N)
  sdA         = rep(0.15,N)
  meanB       = rep(0.4,N)
  # diag(meanB) = rep(0.4,N)  # if B is a diagonal matrix
  sdB         = rep(0.1,N)
  # diag(sdB)   = rep(0.1,N)  # if B is a diagonal matrix
  
  
  
  if (include.AR) {
    lb     = c(rep(-Inf, N*(N+1)), rep(1e-8, N), rep(1e-8, 2*N), rep(-1, (N-1)*N/2))     # Positive lower and upper bounds for W
    ub     = c(rep(Inf, N*(N+1)), rep(0.3, N),  rep(Inf, 2*N), rep(1, (N-1)*N/2))
    mean = c(rep(0, N*(N+1)), as.vector(cbind(meanW,meanA,meanB)), rep(1, (N-1)*N/2))
    sd = c(rep(0, N*(N+1)), as.vector(cbind(sdW,sdA,sdB)), rep(1, (N-1)*N/2))
    lbineq = c(-Inf, -Inf, -Inf, -Inf)
    ubineq = c(0, 0, 0, 0)
    fixed  = c(1:(N*(N+1)), (N*(N+4)+1):(N*(N+4)+N*(N-1)/2))
  }else{
    lb     = c(rep(-Inf, N), rep(1e-8, N), rep(1e-8, 2*N), rep(-1, (N-1)*N/2))     # Positive lower and upper bounds for W
    ub     = c(rep(Inf, N), rep(0.3, N),  rep(Inf, 2*N), rep(1, (N-1)*N/2))
    mean = c(rep(0, N), as.vector(cbind(meanW,meanA,meanB)), rep(1, (N-1)*N/2))
    sd = c(rep(0, N), as.vector(cbind(sdW,sdA,sdB)), rep(1, (N-1)*N/2))
    lbineq = c(-Inf, -Inf, -Inf)
    ubineq = c(0, 0, 0)
    fixed  = c(1:N, (4*N+1):(4*N+N*(N-1)/2))
  }
  
  distr.opt = apply(rbind(mean, sd), 2, as.list)
  
  if (solver == "solnp"){
    opts = list("trace" = 0)
    est  = solnp(x0, LLT, ineqfun = LLTineq, ineqLB = lbineq, ineqUB = ubineq, LB = lb, UB = ub, control = opts)
  } else if (solver == "gosolnp"){
    est  = gosolnp(x0, fixed, fun = LLT, ineqfun = LLTineq, ineqLB = lbineq,
                   ineqUB = ubineq, LB = lb, UB = ub, distr = rep(2, length(lb)), distr.opt = distr.opt,
                   n.restarts = 10, n.sim = 10000, rseed = 1996, control = list("trace" = tr))
  }
  return(est)
}



############################################################################################################################################


##############   Asymptotic variance of self-weighted QMLE   ############################################################################### 


"AVSW"<- function(data, par, type, w = weight, include.AR = F)
{
  ###########################################################################################
  ##                                                                                       ##
  ##  Description:                                                                         ##
  ##  Asymptotic variance and model checking for VARMA-GARCH model specified               ##
  ##  by LING & McAleer (2003), specified as AR(1)-GARCH(1,1) or mu-GARCH(1,1)             ##                            
  ##  par -- QMLE or self-weighted QMLE                                                    ##
  ##  data -- a nT-by-N matrix of data                                                     ##
  ##  type -- only GARCH part is "diagonal" or "extended"                                  ##
  ##  adj -- number of coefficient parameters in the ARMA part in the fitted model.        ##
  ##        (without counting those in the mean and covariance matrix)                     ##
  ##  type -- only GARCH part is "diagonal", "half-diagonal", or "extended"                ##
  ##                                                                                       ##
  ##                                                                                       ##
  ###########################################################################################
  
  # -----------------------------------------------------------------------------------------
  # passin the arguments
  
  nT     = dim(data)[1]    # number of observations
  N      = dim(data)[2]    # number of dimensions
  n      = length(par)     # length of parameters 
  
  
  
  if (include.AR) {
    para.mat = p.match(par, include.AR, type, N)
    heps     = .Call("garch_ar", data, para.mat$mu, para.mat$phi, para.mat$W, para.mat$A, para.mat$B)
  }else{
    para.mat = p.match(par, include.AR, type, N)
    heps     = .Call("garch_mu", data, para.mat$mu, para.mat$W, para.mat$A, para.mat$B)
  }
  
  
  para.mat$Gamma <- make.pd(para.mat$Gamma)
  
  # constructing volatilities
  h       = heps[[1]]             # estimated volatility
  epsit   = heps[[2]]
  invh    = 1/h
  sq.h    = sqrt(h)
  sq.invh = 1/sq.h
  eta     = heps[[2]]/sq.h     # eta or standardized residuals
  invR    = solve(para.mat$Gamma)
  Cn      = invR*para.mat$Gamma + diag(N)
  
  
  # -----------------------------------------------------------------------
  # the score function and information function( i.e. the first derivatives and second derivatives)
  
  
  # partial derivatives of Gamma w.r.t. correlation coefficients
  
  DTGammaDdelta = matrix(0,(N*(N-1)/2),N*N)
  
  for (i in 1:(N*(N-1)/2)){
    delta = matrix(0,(N*(N-1)/2),1)
    delta[i] = 1
    Gamma1 = matrix(0, N, N)
    Gamma1[lower.tri(Gamma1)] = delta
    Gamma1 = t(Gamma1)
    Gamma1[lower.tri(Gamma1)] = delta
    DTGammaDdelta[i,]=as.vector(Gamma1)
  }
  
  
  P = (diag(N)%x%invR)%*%t(DTGammaDdelta)
  C1=matrix(0,N,(N*N))
  for (i in 1:N) {
    C1[i,((i-1)*N+i)]=1
  }
  
  
  
  r.ini=matrix(0,1,N)
  h.ini=matrix(0,1,N)
  epsi.ini=matrix(0,1,N)
  rt=rbind(r.ini,data)
  hthat=rbind(h.ini,h)
  epsithat=rbind(epsi.ini,epsit)
  
  
  
  if(include.AR){
    
    DTepsitDphi  = array(0,dim=c((N*(N+1)),N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,(N*(N+1)),(nT+1)))
    
    if(type=="diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+4)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+4)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+4)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2),nT))
      G      = matrix(0,(N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11  = matrix(0, (N*(N+4)), (N*(N+4)))
      H12  = matrix(0, (N*(N+4)), (N*(N-1)/2))
      H22  = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
      
    }else if(type=="half-diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(2*N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(2*N+3)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(2*N+3)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(2*N+3)), (N*(2*N+3)))
      H12 = matrix(0, (N*(2*N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(3*N+2)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(3*N+2)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(3*N+2)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(3*N+2)), (N*(3*N+2)))
      H12 = matrix(0, (N*(3*N+2)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }else{
    
    DTepsitDphi  = array(0,dim=c(N,N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,N,(nT+1)))
    
    if(type=="diagonal"){
      
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((4*N),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(4*N),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      
      grad   = matrix(0,(4*N+N*(N-1)/2),nT)
      grads  = array(0, dim = c((4*N+N*(N-1)/2),(4*N+N*(N-1)/2),nT))
      G      = matrix(0, (4*N+N*(N-1)/2),(4*N+N*(N-1)/2))
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (4*N), (4*N))
      H12 = matrix(0, (4*N), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else if(type=="half-diagonal"){
      
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+3)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+3)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(N+3)), (N*(N+3)))
      H12 = matrix(0, (N*(N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((2*N*(N+1)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(2*N*(N+1)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(2*N*(N+1)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2),nT))
      G      = matrix(0, (2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (2*N*(N+1)), (2*N*(N+1)))
      H12 = matrix(0, (2*N*(N+1)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }
  
  
  
  DTepsitDphi  = DTepsitDphi[,,(2:(nT+1))]
  DhtDTphi     = DhtDTphi[,,(2:(nT+1))]
  DhtDTlambda  = DhtDTlambda[,,(2:(nT+1))]
  DTepsitDpl   = DTepsitDpl[,,(2:(nT+1))]
  DhtDTpl      = DhtDTpl[,,(2:(nT+1))]
  
  
  invDRD  =  array(0,dim=c(N,N,nT))
  for (i in 1:nT){
    D              = diag(sq.h[i,])
    invD           = diag(sq.invh[i,])
    invH           = diag(invh[i,])
    z              = eta[i,]%o%eta[i,]
    invDRD[,,i]    = invD%*%invR%*%invD
    zeta           = matrix(1,N,1)-diag(eta[i,])%*%invR%*%eta[i,]
    DlDphi[,i]     = DTepsitDphi[,,i]%*%invDRD[,,i]%*%epsit[i,]+0.5*t(DhtDTphi[,,i])%*%invH%*%zeta
    DlDlambda[,i]  = 0.5*t(DhtDTlambda[,,i])%*%invH%*%zeta
    DlDdelta[,i]   = 0.5*DTGammaDdelta%*%as.vector(invR - invR%*%z%*%invR)
    
    
    
    H11           = H11 + w[i]*(DTepsitDpl[,,i]%*%invDRD[,,i]%*%t(DTepsitDpl[,,i]) + 0.25*t(DhtDTpl[,,i])%*%invH%*%Cn%*%invH%*%DhtDTpl[,,i])
    H12           = H12 + w[i]*(0.5*t(DhtDTpl[,,i])%*%invH%*%C1%*%P)
    H22           = H22 + w[i]*0.5*t(P)%*%P
    
    
    grad[,i]      = c(DlDphi[,i],DlDlambda[,i],DlDdelta[,i])
    grads[,,i]    = w[i]*w[i]*(grad[,i]%*%t(grad[,i]))
    G             = G + grads[,,i] 
    
    
  }
  
  
  H   = rbind(cbind(H11,H12),cbind(t(H12),H22))         # analytical Hessian
  
  
  
  
  # ---------------------------------------------------------------------------
  # calculation of asymptotic variance 
  
  
  gradi      = G/nT                                    # quasi-information matrix calculated by the procudt of gradient
  invHe      = try(solve(H/nT),silent = TRUE)
  if(class(invHe)[1]=="try-error")
  {
    message("Caught an error 1")
    return(list(error=1, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0, n, n)))
  }else{
    Avar    = (invHe%*%gradi%*%invHe)/nT
    rob.se  = sqrt(diag(Avar))                                       # robust standard errors
    out.se  = try(sqrt(diag(solve(gradi)/nT)),silent = TRUE)         # standard errors based on the outer product of the gradient
    if(class(out.se)[1]=="try-error")
    {
      message("Caught an error 1")
      return(list(error=1, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0, n, n)))
    }else{
      H.se    = sqrt(diag(invHe/nT))               #  standard errors based on the inverted Hessian
      Avar1   = invHe%*%gradi%*%invHe              #  for wald test
      return(list(error=0,vari = cbind(H.se, out.se, rob.se), vari1 = Avar1))
    }
  }
}



##############   Asymptotic variance of self-weighted QMLE and model checking  #############################################################


"AVSWMC"<- function(data, par, type, w = weight, include.AR = F, lag = 6)
{
  ###########################################################################################
  ##                                                                                       ##
  ##  Description:                                                                         ##
  ##  Asymptotic variance and model checking for VARMA-GARCH model specified               ##
  ##  by LING & McAleer (2003), specified as AR(1)-GARCH(1,1) or mu-GARCH(1,1)             ##                            
  ##  par -- QMLE or self-weighted QMLE                                                    ##
  ##  data -- a nT-by-N matrix of data                                                     ##
  ##  type -- only GARCH part is "diagonal" or "extended"                                  ##
  ##  adj -- number of coefficient parameters in the ARMA part in the fitted model.        ##
  ##        (without counting those in the mean and covariance matrix)                     ##
  ##  type -- only GARCH part is "diagonal", "half-diagonal", or "extended"                ##
  ##                                                                                       ##
  ##                                                                                       ##
  ###########################################################################################
  
  # -----------------------------------------------------------------------------------------
  # passin the arguments
  
  nT     = dim(data)[1]    # number of observations
  N      = dim(data)[2]    # number of dimensions
  n      = length(par)     # length of parameters 
  n1     = n - (N*(N-1)/2)   # length of parameters except correlation coefficients
  
  
  if (include.AR) {
    para.mat = p.match(par, include.AR, type, N)
    heps = .Call("garch_ar", data, para.mat$mu, para.mat$phi, para.mat$W, para.mat$A, para.mat$B)
  }else{
    para.mat = p.match(par, include.AR, type, N)
    heps = .Call("garch_mu", data, para.mat$mu, para.mat$W, para.mat$A, para.mat$B)
  }
  
  
  para.mat$Gamma <- make.pd(para.mat$Gamma)
  
  # constructing volatilities
  h       = heps[[1]]             # estimated volatility
  epsit   = heps[[2]]
  invh    = 1/h
  sq.h    = sqrt(h)
  sq.invh = 1/sq.h
  eta     = heps[[2]]/sq.h     # eta or standardized residuals
  invR    = solve(para.mat$Gamma)
  Cn      = invR*para.mat$Gamma + diag(N)
  
  
  # -----------------------------------------------------------------------
  # the score function and information function( i.e. the first derivatives and second derivatives)
  
  
  # partial derivatives of Gamma w.r.t. correlation coefficients
  
  DTGammaDdelta = matrix(0,(N*(N-1)/2),N*N)
  
  for (i in 1:(N*(N-1)/2)){
    delta = matrix(0,(N*(N-1)/2),1)
    delta[i] = 1
    Gamma1 = matrix(0, N, N)
    Gamma1[lower.tri(Gamma1)] = delta
    Gamma1 = t(Gamma1)
    Gamma1[lower.tri(Gamma1)] = delta
    DTGammaDdelta[i,]=as.vector(Gamma1)
  }
  
  
  P = (diag(N)%x%invR)%*%t(DTGammaDdelta)
  C1=matrix(0,N,(N*N))
  for (i in 1:N) {
    C1[i,((i-1)*N+i)]=1
  }
  
  
  
  r.ini=matrix(0,1,N)
  h.ini=matrix(0,1,N)
  epsi.ini=matrix(0,1,N)
  rt=rbind(r.ini,data)
  hthat=rbind(h.ini,h)
  epsithat=rbind(epsi.ini,epsit)
  
  
  
  if(include.AR){
    
    DTepsitDphi  = array(0,dim=c((N*(N+1)),N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,(N*(N+1)),(nT+1)))
    
    if(type=="diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+4)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+4)),(nT+1)))
      
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(N+4)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta    =  array(0,dim=c((N*N),((N*(N+4))+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+4)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2),nT))
      G      = matrix(0,(N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11  = matrix(0, (N*(N+4)), (N*(N+4)))
      H12  = matrix(0, (N*(N+4)), (N*(N-1)/2))
      H22  = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
      
    }else if(type=="half-diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(2*N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(2*N+3)),(nT+1)))
      
      
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(2*N+3)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),((N*(2*N+3))+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(2*N+3)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(2*N+3)), (N*(2*N+3)))
      H12 = matrix(0, (N*(2*N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
      
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(3*N+2)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(3*N+2)),(nT+1)))
      
      
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(3*N+2)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),((N*(3*N+2))+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(3*N+2)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(3*N+2)), (N*(3*N+2)))
      H12 = matrix(0, (N*(3*N+2)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }else{
    
    DTepsitDphi  = array(0,dim=c(N,N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,N,(nT+1)))
    
    if(type=="diagonal"){
      
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((4*N),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(4*N),(nT+1)))
      
      
      DvDTpl      =  array(0,dim=c((N*N),(4*N),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),((4*N)+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      
      grad   = matrix(0,(4*N+N*(N-1)/2),nT)
      grads  = array(0, dim = c((4*N+N*(N-1)/2),(4*N+N*(N-1)/2),nT))
      G      = matrix(0, (4*N+N*(N-1)/2),(4*N+N*(N-1)/2))
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (4*N), (4*N))
      H12 = matrix(0, (4*N), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else if(type=="half-diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+3)),(nT+1)))
      
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(N+3)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),(N*(N+3)+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+3)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(N+3)), (N*(N+3)))
      H12 = matrix(0, (N*(N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((2*N*(N+1)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(2*N*(N+1)),(nT+1)))
      
      
      DvDTpl      =  array(0,dim=c((N*N),(2*N*(N+1)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),(2*N*(N+1)+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(2*N*(N+1)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2),nT))
      G      = matrix(0, (2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (2*N*(N+1)), (2*N*(N+1)))
      H12 = matrix(0, (2*N*(N+1)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }
  
  
  
  DTepsitDphi  = DTepsitDphi[,,(2:(nT+1))]
  DhtDTphi     = DhtDTphi[,,(2:(nT+1))]
  DhtDTlambda  = DhtDTlambda[,,(2:(nT+1))]
  DTepsitDpl   = DTepsitDpl[,,(2:(nT+1))]
  DhtDTpl      = DhtDTpl[,,(2:(nT+1))]
  
  
  invDRD  =  array(0,dim=c(N,N,nT))
  for (i in 1:nT){
    D              = diag(sq.h[i,])
    invD           = diag(sq.invh[i,])
    invH           = diag(invh[i,])
    z              = eta[i,]%o%eta[i,]
    invDRD[,,i]    = invD%*%invR%*%invD
    zeta           = matrix(1,N,1)-diag(eta[i,])%*%invR%*%eta[i,]
    DlDphi[,i]     = DTepsitDphi[,,i]%*%invDRD[,,i]%*%epsit[i,]+0.5*t(DhtDTphi[,,i])%*%invH%*%zeta
    DlDlambda[,i]  = 0.5*t(DhtDTlambda[,,i])%*%invH%*%zeta
    DlDdelta[,i]   = 0.5*DTGammaDdelta%*%as.vector(invR - invR%*%z%*%invR)
    
    
    
    H11           = H11 + w[i]*(DTepsitDpl[,,i]%*%invDRD[,,i]%*%t(DTepsitDpl[,,i]) + 0.25*t(DhtDTpl[,,i])%*%invH%*%Cn%*%invH%*%DhtDTpl[,,i])
    H12           = H12 + w[i]*(0.5*t(DhtDTpl[,,i])%*%invH%*%C1%*%P)
    H22           = H22 + w[i]*0.5*t(P)%*%P
    
    
    grad[,i]      = c(DlDphi[,i],DlDlambda[,i],DlDdelta[,i])
    grads[,,i]    = w[i]*w[i]*(grad[,i]%*%t(grad[,i]))
    G             = G + grads[,,i] 
    
    
    
    # for model checking 
    # partition matrix "DhtDTpl" into small matrices
    
    partitioned_matrices =  matrix(0, nrow = N, ncol = n1)
    matrix_original      =  DhtDTpl[,,i]
    for (ii in 1:N) {
      partitioned_matrices[ii, ] <- matrix_original[ii, ]
    }
    
    # construct quasi-diagonal matrix
    
    quasi_diagonal_matrix <- matrix(0, nrow = N, ncol = N * n1)
    for (ii in 1:N) {
      start_col <- (ii - 1) * n1 + 1
      end_col <- ii * n1
      quasi_diagonal_matrix[ii, start_col:end_col] <- partitioned_matrices[ii, ]
    }
    quasi_diagonal_matrix = 0.5*invD%*%quasi_diagonal_matrix
    
    #  vectorized the quasi-diagonal matrix 
    vectorized_matrix <- matrix(0, nrow = N * N, ncol = n1)
    for (ii in 1:(N * N)) {
      start_col <- ((ii - 1) %% N) * n1 + 1
      end_col <- ((ii - 1) %% N + 1) * n1
      vectorized_matrix[ii, ] <- quasi_diagonal_matrix[(ii - 1) %/% N + 1, start_col:end_col]
    }
    
    
    DvDTpl[,,i]      =  (diag(N)%x%(D%*%para.mat$Gamma))%*%vectorized_matrix+((D%*%para.mat$Gamma)%x%diag(N))%*%vectorized_matrix
    DvDTdelta[,,i]   =  (D%x%D)%*%t(DTGammaDdelta)
    DvDTtheta[,,i]   =  cbind(DvDTpl[,,i],DvDTdelta[,,i])
    
    
  }
  
  
  H   = rbind(cbind(H11,H12),cbind(t(H12),H22))         # analytical Hessian
  
  
  
  
  # ---------------------------------------------------------------------------
  # calculation of asymptotic variance 
  
  
  teta   =  matrix(0,1,nT)
  for (ii in 1:nT) {
    teta[ii] =w[ii]*(epsit[ii,,drop=FALSE]%*%invDRD[,,ii]%*%t(epsit[ii,,drop=FALSE])-N)
  }
  teta   = as.vector(teta)
  k0     = sum(teta^2)
  
  
  R       =  matrix(0,nrow=1,ncol=lag)
  for (i in 1:lag) {
    x11     =  teta[(i+1):nT]
    x22     =  teta[1:(nT-i)]
    R[i]    =  t(x11)%*%x22
    R[i]    =  R[i]/k0
    
    
  }
  
  
  X   = matrix(0,nrow=n,ncol=lag)
  xl  = matrix(0,nrow=n,ncol=1)
  
  for (i in 1:lag) {
    x111     =  w[(i+1):nT]
    x222     =  DvDTtheta[,,(i+1):nT]
    x333     =  invDRD[,,(i+1):nT]
    x444     =  teta[1:(nT-i)]
    for (ii in 1:(nT-i)) {
      xl     =  xl + x111[ii]*t(x222[,,ii])%*%as.vector(x333[,,ii]*x444[ii])
    }
    X[,i]    =  xl/(nT-i)
  }
  
  
  
  
  Z1 = matrix(0, nrow = nT-6, ncol = 1)
  Z2 = matrix(0, nrow = nT-6, ncol = 1)
  Z3 = matrix(0, nrow = nT-6, ncol = 1)
  Z4 = matrix(0, nrow = nT-6, ncol = 1)
  Z5 = matrix(0, nrow = nT-6, ncol = 1)
  Z6 = matrix(0, nrow = nT-6, ncol = 1)
  Z7 = matrix(0, nrow = nT-6, ncol = n)
  Z  = matrix(0, nrow = nT-6, ncol = 6+n)
  
  for (i in 7:nT) {
    Z1[i-6,] = teta[i]*teta[i-1]  
    Z2[i-6,] = teta[i]*teta[i-2] 
    Z3[i-6,] = teta[i]*teta[i-3] 
    Z4[i-6,] = teta[i]*teta[i-4]
    Z5[i-6,] = teta[i]*teta[i-5] 
    Z6[i-6,] = teta[i]*teta[i-6]
    Z7[i-6,] = -w[i]*t(grad[,i]) 
  }
  Z  = cbind(Z1,Z2,Z3,Z4,Z5,Z6,Z7)
  Z0 = (t(Z)%*%Z)/(nT-6)
  
  
  
  gradi      = G/nT                                    # quasi-information matrix calculated by the procudt of gradient
  invHe      = try(solve(H/nT),silent = TRUE)
  if(class(invHe)[1]=="try-error")
  {
    message("Caught an error 1")
    return(list(error=1, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0, n, n), QMM = 0))
  }else{
    Avar    = (invHe%*%gradi%*%invHe)/nT
    rob.se  = sqrt(diag(Avar))                                       # robust standard errors
    out.se  = try(sqrt(diag(solve(gradi)/nT)),silent = TRUE)         # standard errors based on the outer product of the gradient
    if(class(out.se)[1]=="try-error")
    {
      message("Caught an error 1")
      return(list(error=1, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0, n, n), QMM = 0))
    }else{
      H.se    = sqrt(diag(invHe/nT))               #  standard errors based on the inverted Hessian
      Avar1   = invHe%*%gradi%*%invHe              #  for wald test
      V       = cbind(diag(lag), -t(X)%*%invHe)
      MV      = (V%*%Z0%*%t(V))/(k0/nT)^2
      # QMM     = nT*R%*%solve(MV)%*%t(R)
      QMM     = (nT-lag-2*N-1-2+1)*R%*%solve(MV)%*%t(R)
      return(list(error=0,vari = cbind(H.se, out.se, rob.se), vari1 = Avar1, QMM = QMM))
    }
  }
}



#########################################    local QMLE     #################################################################################


"Localfit"<- function(data, par, par1, type, include.AR = F)
{
  #########################################################################################################
  ##                                                                                                     ##
  ##  Description:                                                                                       ##
  ##  Local QMLE for VARMA-GARCH model specified by LING & McAleer (2003)                                ##
  ##  Specified as AR(1)-GARCH(1,1) or mu-GARCH(1,1)                                                     ##
  ##  par --  self-weighted or local QMLE                                                                ##
  ##  data -- a nT-by-N matrix of data                                                                   ##
  ##  type -- only GARCH part is "diagonal", "half-diagonal", or "extended"                              ##
  ##                                                                                                     ##
  #########################################################################################################
  
  # -----------------------------------------------------------------------------------------
  # passin the arguments
  
  nT     = dim(data)[1]    # number of observations
  N      = dim(data)[2]    # number of dimensions
  n      = length(par)     # length of parameters 
  
  
  
  GH <- function(par, data){
    if (include.AR) {
      para.mat = p.match(par, include.AR, type, N)
      heps     = .Call("garch_ar", data, para.mat$mu, para.mat$phi, para.mat$W, para.mat$A, para.mat$B)
    }else{
      para.mat = p.match(par, include.AR, type, N)
      heps     = .Call("garch_mu", data, para.mat$mu, para.mat$W, para.mat$A, para.mat$B)
    }
    
    
    para.mat$Gamma <- make.pd(para.mat$Gamma)
    
    # constructing volatilities
    h       = heps[[1]]             # estimated volatility
    epsit   = heps[[2]]
    invh    = 1/h
    sq.h    = sqrt(h)
    sq.invh = 1/sqrt(h)
    eta     = heps[[2]]/sq.h     # eta or standardized residuals
    invR    = solve(para.mat$Gamma)
    Cn      = invR*para.mat$Gamma + diag(N)
    
    
    # -----------------------------------------------------------------------
    # the score function and information matrix( i.e. the first derivatives and second derivatives)
    
    
    # partial derivatives of Gamma w.r.t. correlation coefficients
    
    DTGammaDdelta = matrix(0,(N*(N-1)/2),N*N)
    
    for (i in 1:(N*(N-1)/2)){
      delta = matrix(0,(N*(N-1)/2),1)
      delta[i] = 1
      Gamma1 = matrix(0, N, N)
      Gamma1[lower.tri(Gamma1)] = delta
      Gamma1 = t(Gamma1)
      Gamma1[lower.tri(Gamma1)] = delta
      DTGammaDdelta[i,]=as.vector(Gamma1)
    }
    
    
    P = (diag(N)%x%invR)%*%t(DTGammaDdelta)
    C1=matrix(0,N,(N*N))
    for (i in 1:N) {
      C1[i,((i-1)*N+i)]=1
    }
    
    
    
    r.ini=matrix(0,1,N)
    h.ini=matrix(0,1,N)
    epsi.ini=matrix(0,1,N)
    rt=rbind(r.ini,data)
    hthat=rbind(h.ini,h)
    epsithat=rbind(epsi.ini,epsit)
    
    
    
    
    if(include.AR){
      
      DTepsitDphi  = array(0,dim=c((N*(N+1)),N,(nT+1)))
      DhtDTphi     = array(0,dim=c(N,(N*(N+1)),(nT+1)))
      
      if(type=="diagonal"){
        DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
        DTepsitDpl   = array(0,dim=c((N*(N+4)),N,(nT+1)))
        DhtDTpl      = array(0,dim=c(N,(N*(N+4)),(nT+1)))
        
        
        # the score functions ( i.e. the first derivatives )
        
        DlDphi       = matrix(0,(N*(N+1)),nT)
        DlDlambda    = matrix(0,(3*N),nT)
        DlDdelta     = matrix(0,(N*(N-1)/2),nT)
        
        
        grad   = matrix(0,(N*(N+4)+N*(N-1)/2),nT)
        G      = matrix(0, (N*(N+4)+N*(N-1)/2),1)
        
        
        # the Hessian matrix ( i.e. the second derivatives )
        H11 = matrix(0, (N*(N+4)), (N*(N+4)))
        H12 = matrix(0, (N*(N+4)), (N*(N-1)/2))
        H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
        
        
        
        for (i in 2:(nT+1)){
          DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
          DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
          DhtDTlambda[,,i]  = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
          DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
          DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
        }
        
      }else if(type=="half-diagonal"){
        DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
        DTepsitDpl   = array(0,dim=c((N*(2*N+3)),N,(nT+1)))
        DhtDTpl      = array(0,dim=c(N,(N*(2*N+3)),(nT+1)))
        
        
        
        # the score functions ( i.e. the first derivatives )
        DlDphi       = matrix(0,(N*(N+1)),nT)
        DlDlambda    = matrix(0,(N*(N+2)),nT)
        DlDdelta     = matrix(0,(N*(N-1)/2),nT)
        
        
        grad   = matrix(0,(N*(2*N+3)+N*(N-1)/2),nT)
        grads  = array(0, dim = c((N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2),nT))
        G      = matrix(0, (N*(2*N+3)+N*(N-1)/2),1)
        
        
        # the Hessian matrix ( i.e. the second derivatives )
        H11 = matrix(0, (N*(2*N+3)), (N*(2*N+3)))
        H12 = matrix(0, (N*(2*N+3)), (N*(N-1)/2))
        H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
        
        
        for (i in 2:(nT+1)){
          DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
          DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
          DhtDTlambda[,,i]  = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
          DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
          DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
        }
      }else{
        
        DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
        DTepsitDpl   = array(0,dim=c((N*(3*N+2)),N,(nT+1)))
        DhtDTpl      = array(0,dim=c(N,(N*(3*N+2)),(nT+1)))
        
        
        # the score functions ( i.e. the first derivatives )
        DlDphi       = matrix(0,(N*(N+1)),nT)
        DlDlambda    = matrix(0,(N*(2*N+1)),nT)
        DlDdelta     = matrix(0,(N*(N-1)/2),nT)
        
        
        grad   = matrix(0,(N*(3*N+2)+N*(N-1)/2),nT)
        G      = matrix(0, (N*(3*N+2)+N*(N-1)/2), 1)
        
        
        # the Hessian matrix ( i.e. the second derivatives )
        H11 = matrix(0, (N*(3*N+2)), (N*(3*N+2)))
        H12 = matrix(0, (N*(3*N+2)), (N*(N-1)/2))
        H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
        
        
        
        for (i in 2:(nT+1)){
          DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
          DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
          DhtDTlambda[,,i]  = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
          DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
          DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
        }
      }
    }else{
      if(type=="diagonal"){
        DTepsitDphi  = array(0,dim=c(N,N,(nT+1)))
        DhtDTphi     = array(0,dim=c(N,N,(nT+1)))
        DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
        DTepsitDpl   = array(0,dim=c((4*N),N,(nT+1)))
        DhtDTpl      = array(0,dim=c(N,(4*N),(nT+1)))
        
        
        # the score functions ( i.e. the first derivatives )
        DlDphi       = matrix(0,N,nT)
        DlDlambda    = matrix(0,(3*N),nT)
        DlDdelta     = matrix(0,(N*(N-1)/2),nT)
        
        
        grad   = matrix(0,(4*N+N*(N-1)/2),nT)
        G      = matrix(0, (4*N+N*(N-1)/2), 1)
        
        
        
        # the Hessian matrix ( i.e. the second derivatives )
        H11 = matrix(0, (4*N), (4*N))
        H12 = matrix(0, (4*N), (N*(N-1)/2))
        H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
        
        
        
        
        for (i in 2:(nT+1)){
          DTepsitDphi[,,i]   = -diag(N)
          DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
          DhtDTlambda[,,i]   = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
          DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
          DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
        }
      }else if(type=="half-diagonal"){
        
        DTepsitDphi  = array(0,dim=c(N,N,(nT+1)))
        DhtDTphi     = array(0,dim=c(N,N,(nT+1)))
        DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
        DTepsitDpl   = array(0,dim=c((N*(N+3)),N,(nT+1)))
        DhtDTpl      = array(0,dim=c(N,(N*(N+3)),(nT+1)))
        
        
        
        # the score functions ( i.e. the first derivatives )
        DlDphi       = matrix(0,N,nT)
        DlDlambda    = matrix(0,(N*(N+2)),nT)
        DlDdelta     = matrix(0,(N*(N-1)/2),nT)
        
        
        grad   = matrix(0,(N*(N+3)+N*(N-1)/2),nT)
        grads  = array(0, dim = c((N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2),nT))
        G      = matrix(0, (N*(N+3)+N*(N-1)/2),1)
        
        
        # the Hessian matrix ( i.e. the second derivatives )
        H11 = matrix(0, (N*(N+3)), (N*(N+3)))
        H12 = matrix(0, (N*(N+3)), (N*(N-1)/2))
        H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
        
        
        
        for (i in 2:(nT+1)){
          DTepsitDphi[,,i]   = -diag(N)
          DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
          DhtDTlambda[,,i]   = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
          DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
          DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
        }
      }else{
        
        DTepsitDphi  = array(0,dim=c(N,N,(nT+1)))
        DhtDTphi     = array(0,dim=c(N,N,(nT+1)))
        DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
        DTepsitDpl   = array(0,dim=c((2*N*(N+1)),N,(nT+1)))
        DhtDTpl      = array(0,dim=c(N,(2*N*(N+1)),(nT+1)))
        
        
        # the score functions ( i.e. the first derivatives )
        DlDphi       = matrix(0,N,nT)
        DlDlambda    = matrix(0,(N*(2*N+1)),nT)
        DlDdelta     = matrix(0,(N*(N-1)/2),nT)
        
        
        grad   = matrix(0,(2*N*(N+1)+N*(N-1)/2),nT)
        G      = matrix(0, (2*N*(N+1)+N*(N-1)/2),1)
        
        
        
        # the Hessian matrix ( i.e. the second derivatives )
        H11 = matrix(0, (2*N*(N+1)), (2*N*(N+1)))
        H12 = matrix(0, (2*N*(N+1)), (N*(N-1)/2))
        H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
        
        
        for (i in 2:(nT+1)){
          DTepsitDphi[,,i]   = -diag(N)
          DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
          DhtDTlambda[,,i]   = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
          DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
          DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
        }
      }
    }
    
    
    
    
    DTepsitDphi  = DTepsitDphi[,,(2:(nT+1))]
    DhtDTphi     = DhtDTphi[,,(2:(nT+1))]
    DhtDTlambda  = DhtDTlambda[,,(2:(nT+1))]
    DTepsitDpl   = DTepsitDpl[,,(2:(nT+1))]
    DhtDTpl      = DhtDTpl[,,(2:(nT+1))]
    
    
    for (i in 1:nT){
      D              = diag(sq.h[i,])
      invD           = diag(sq.invh[i,])
      invH           = diag(invh[i,])
      z              = eta[i,]%o%eta[i,]
      invDRD         = invD%*%invR%*%invD
      zeta           = matrix(1,N,1)-diag(eta[i,])%*%invR%*%eta[i,]
      DlDphi[,i]     = DTepsitDphi[,,i]%*%invDRD%*%epsit[i,]+0.5*t(DhtDTphi[,,i])%*%invH%*%zeta
      DlDlambda[,i]  = 0.5*t(DhtDTlambda[,,i])%*%invH%*%zeta
      DlDdelta[,i]   = 0.5*DTGammaDdelta%*%as.vector(invR - invR%*%z%*%invR)
      
      
      H11           = H11 + (DTepsitDpl[,,i]%*%invDRD%*%t(DTepsitDpl[,,i]) + 0.25*t(DhtDTpl[,,i])%*%invH%*%Cn%*%invH%*%DhtDTpl[,,i])
      H12           = H12 + (0.5*t(DhtDTpl[,,i])%*%invH%*%C1%*%P)
      H22           = H22 + 0.5*t(P)%*%P
      
      
      grad[,i]      = c(DlDphi[,i],DlDlambda[,i],DlDdelta[,i])
      G             = G + grad[,i]
      
    }
    
    
    H  =  rbind(cbind(H11,H12),cbind(t(H12),H22))             # analytical Hessian
    G  =  G/nT
    H  =  H/nT
    
    return(list(G=G,H=H))
  }
  
  GH       = GH(par, data)
  gradient = GH$G                             # gradients calculated by equations
  Hessian  = GH$H
  invHe    = try(solve(Hessian),silent = TRUE)
  
  if(class(invHe)[1]=="try-error"){
    message("Caught an error 1")
    return(list(error=1, fit = matrix(0,n,1)))
  }else{
    fit      = par - invHe%*%gradient
    para.mat = p.match(fit, include.AR, type, N)
    if(include.AR){
      if(max(abs(eigen(para.mat$B)$values))>=(1-1e-8)|max(abs(eigen(para.mat$phi)$values))>=(1-1e-8)){
        message("Caught an error 2")
        return(list(error=2, fit = matrix(0,n,1)))
      }else{
        if(all(fit[(N*(N+1)+1):(N*(N+2))]<=0.3)){
          if(type=="diagonal"){
            if(all(fit[(N*(N+1)+1):(N*(N+4))]>0)){
              return(list(error=0, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N*(N+1)+1):(N*(N+4))]<=0)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(N+4))]<=0)+(N*(N+1))]
              return(list(error=3, fit = fit))
            }
          }else if(type=="half-diagonal"){
            if(all(fit[(N*(N+1)+1):(N*(2*N+3))]>0)){
              return(list(error=0, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N*(N+1)+1):(N*(2*N+3))]<=0)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(2*N+3))]<=0)+(N*(N+1))]
              return(list(error=3, fit = fit))
            }
          }else{
            if(all(fit[(N*(N+1)+1):(N*(3*N+2))]>0)){
              return(list(error=0, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N*(N+1)+1):(N*(3*N+2))]<=0)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(3*N+2))]<=0)+(N*(N+1))]
              return(list(error=3, fit = fit))
            }
          }
        }else{
          message("Caught an error 3")
          fit[which(fit[(N*(N+1)+1):(N*(N+2))]>0.3)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(N+2))]>0.3)+(N*(N+1))]
          # fit[which(fit[(N*(N+1)+1):(N*(N+2))]<0.05)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(N+2))]<0.05)+(N*(N+1))]
          if(type=="diagonal"){
            if(all(fit[(N*(N+1)+1):(N*(N+4))]>0)){
              return(list(error=3, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N*(N+1)+1):(N*(N+4))]<=0)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(N+4))]<=0)+(N*(N+1))]
              return(list(error=3, fit = fit))
            }
          }else if(type=="half-diagonal"){
            if(all(fit[(N*(N+1)+1):(N*(2*N+3))]>0)){
              return(list(error=3, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N*(N+1)+1):(N*(2*N+3))]<=0)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(2*N+3))]<=0)+(N*(N+1))]
              return(list(error=3, fit = fit))
            }
          }else{
            if(all(fit[(N*(N+1)+1):(N*(3*N+2))]>0)){
              return(list(error=3, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N*(N+1)+1):(N*(3*N+2))]<=0)+(N*(N+1))]=par1[which(fit[(N*(N+1)+1):(N*(3*N+2))]<=0)+(N*(N+1))]
              return(list(error=3, fit = fit))
            }
          }
        }
      }
    }else{
      if(max(abs(eigen(para.mat$B)$values))>=(1-1e-8)){
        message("Caught an error 2")
        return(list(error=2, fit = matrix(0,n,1)))
      }else{
        if(all(fit[(N+1):(2*N)]<=0.3)){
          if(type=="diagonal"){
            if(all(fit[(N+1):(4*N)]>0)){
              return(list(error=0, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N+1):(4*N)]<=0)+N]=par1[which(fit[(N+1):(4*N)]<=0)+N]
              return(list(error=3, fit = fit))
            }
          }else if(type=="half-diagonal"){
            if(all(fit[(N+1):(N*(N+3))]>0)){
              return(list(error=0, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N+1):(N*(N+3))]<=0)+N]=par1[which(fit[(N+1):(N*(N+3))]<=0)+N]
              return(list(error=3, fit = fit))
            }
          }else{
            if(all(fit[(N+1):(2*N*(N+1))]>0)){
              return(list(error=0, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N+1):(2*N*(N+1))]<=0)+N]=par1[which(fit[(N+1):(2*N*(N+1))]<=0)+N]
              return(list(error=3, fit = fit))
            }
          }
        }else{
          message("Caught an error 3")
          fit[which(fit[(N+1):(2*N)]>0.3)+N]=par1[which(fit[(N+1):(2*N)]>0.3)+N]
          # fit[which(fit[(N+1):(2*N)]<0.05)+N]=par1[which(fit[(N+1):(2*N)]<0.05)+N]
          if(type=="diagonal"){
            if(all(fit[(N+1):(4*N)]>0)){
              return(list(error=3, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N+1):(4*N)]<=0)+N]=par1[which(fit[(N+1):(4*N)]<=0)+N]
              return(list(error=3, fit = fit))
            }
          }else if(type=="half-diagonal"){
            if(all(fit[(N+1):(N*(N+3))]>0)){
              return(list(error=3, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N+1):(N*(N+3))]<=0)+N]=par1[which(fit[(N+1):(N*(N+3))]<=0)+N]
              return(list(error=3, fit = fit))
            }
          }else{
            if(all(fit[(N+1):(2*N*(N+1))]>0)){
              return(list(error=3, fit = fit))
            }else{
              message("Caught an error 3")
              fit[which(fit[(N+1):(2*N*(N+1))]<=0)+N]=par1[which(fit[(N+1):(2*N*(N+1))]<=0)+N]
              return(list(error=3, fit = fit))
            }
          }
        }
      }
    }
  }
}


###################    local asymptotic variance    #########################################################################################


"LocalVA"<- function(data, par, type, include.AR = F)
{
  #########################################################################################################
  ##                                                                                                     ##
  ##  Description:                                                                                       ##
  ##  Local asymptotic variance for VARMA-GARCH model specified by LING & McAleer (2003)                 ##
  ##  Specified as AR(1)-GARCH(1,1) or mu-GARCH(1,1)                                                     ##
  ##  par --  self-weighted or local QMLE                                                                ##
  ##  data -- a nT-by-N matrix of data                                                                   ##
  ##  type -- only GARCH part is "diagonal", "half-diagonal", or "extended"                              ##
  ##                                                                                                     ##
  #########################################################################################################
  
  # -----------------------------------------------------------------------------------------
  # passin the arguments
  
  nT     = dim(data)[1]    # number of observations
  N      = dim(data)[2]    # number of dimensions
  n      = length(par)     # length of parameters 
  
  
  
  if (include.AR) {
    para.mat = p.match(par, include.AR, type, N)
    heps = .Call("garch_ar", data, para.mat$mu, para.mat$phi, para.mat$W, para.mat$A, para.mat$B)
  }else{
    para.mat = p.match(par, include.AR, type, N)
    heps = .Call("garch_mu", data, para.mat$mu, para.mat$W, para.mat$A, para.mat$B)
  }
  
  
  para.mat$Gamma <- make.pd(para.mat$Gamma)
  
  # constructing volatilities
  h       = heps[[1]]             # estimated volatility
  epsit   = heps[[2]]
  invh    = 1/h
  sq.h    = sqrt(h)
  sq.invh = 1/sqrt(h)
  eta     = heps[[2]]/sq.h     # eta or standardized residuals
  invR    = solve(para.mat$Gamma)
  Cn      = invR*para.mat$Gamma + diag(N)
  
  
  # -----------------------------------------------------------------------
  # the score function and information matrix( i.e. the first derivatives and second derivatives)
  
  
  # partial derivatives of Gamma w.r.t. correlation coefficients
  
  DTGammaDdelta = matrix(0,(N*(N-1)/2),N*N)
  
  for (i in 1:(N*(N-1)/2)){
    delta = matrix(0,(N*(N-1)/2),1)
    delta[i] = 1
    Gamma1 = matrix(0, N, N)
    Gamma1[lower.tri(Gamma1)] = delta
    Gamma1 = t(Gamma1)
    Gamma1[lower.tri(Gamma1)] = delta
    DTGammaDdelta[i,]=as.vector(Gamma1)
  }
  
  
  P = (diag(N)%x%invR)%*%t(DTGammaDdelta)
  C1=matrix(0,N,(N*N))
  for (i in 1:N) {
    C1[i,((i-1)*N+i)]=1
  }
  
  
  
  r.ini=matrix(0,1,N)
  h.ini=matrix(0,1,N)
  epsi.ini=matrix(0,1,N)
  rt=rbind(r.ini,data)
  hthat=rbind(h.ini,h)
  epsithat=rbind(epsi.ini,epsit)
  
  
  
  
  if(include.AR){
    
    DTepsitDphi  = array(0,dim=c((N*(N+1)),N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,(N*(N+1)),(nT+1)))
    
    if(type=="diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+4)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+4)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+4)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2),nT))
      G      = matrix(0,(N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(N+4)), (N*(N+4)))
      H12 = matrix(0, (N*(N+4)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
      
    }else if(type=="half-diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(2*N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(2*N+3)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(2*N+3)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(2*N+3)), (N*(2*N+3)))
      H12 = matrix(0, (N*(2*N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(3*N+2)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(3*N+2)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(3*N+2)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2))
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(3*N+2)), (N*(3*N+2)))
      H12 = matrix(0, (N*(3*N+2)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }else{
    
    DTepsitDphi  = array(0,dim=c(N,N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,N,(nT+1)))
    
    if(type=="diagonal"){
      
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((4*N),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(4*N),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(4*N+N*(N-1)/2),nT)
      grads  = array(0, dim = c((4*N+N*(N-1)/2),(4*N+N*(N-1)/2),nT))
      G      = matrix(0, (4*N+N*(N-1)/2),(4*N+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (4*N), (4*N))
      H12 = matrix(0, (4*N), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else if(type=="half-diagonal"){
      
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+3)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+3)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(N+3)), (N*(N+3)))
      H12 = matrix(0, (N*(N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((2*N*(N+1)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(2*N*(N+1)),(nT+1)))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(2*N*(N+1)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2),nT))
      G      = matrix(0, (2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (2*N*(N+1)), (2*N*(N+1)))
      H12 = matrix(0, (2*N*(N+1)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }
  
  
  
  
  DTepsitDphi  = DTepsitDphi[,,(2:(nT+1))]
  DhtDTphi     = DhtDTphi[,,(2:(nT+1))]
  DhtDTlambda  = DhtDTlambda[,,(2:(nT+1))]
  DTepsitDpl   = DTepsitDpl[,,(2:(nT+1))]
  DhtDTpl      = DhtDTpl[,,(2:(nT+1))]
  
  
  invDRD  =  array(0,dim=c(N,N,nT))
  for (i in 1:nT){
    D              = diag(sq.h[i,])
    invD           = diag(sq.invh[i,])
    invH           = diag(invh[i,])
    z              = eta[i,]%o%eta[i,]
    invDRD[,,i]    = invD%*%invR%*%invD
    zeta           = matrix(1,N,1)-diag(eta[i,])%*%invR%*%eta[i,]
    DlDphi[,i]     = DTepsitDphi[,,i]%*%invDRD[,,i]%*%epsit[i,]+0.5*t(DhtDTphi[,,i])%*%invH%*%zeta
    DlDlambda[,i]  = 0.5*t(DhtDTlambda[,,i])%*%invH%*%zeta
    DlDdelta[,i]   = 0.5*DTGammaDdelta%*%as.vector(invR - invR%*%z%*%invR)
    
    
    
    H11           = H11 + (DTepsitDpl[,,i]%*%invDRD[,,i]%*%t(DTepsitDpl[,,i]) + 0.25*t(DhtDTpl[,,i])%*%invH%*%Cn%*%invH%*%DhtDTpl[,,i])
    H12           = H12 + (0.5*t(DhtDTpl[,,i])%*%invH%*%C1%*%P)
    H22           = H22 + 0.5*t(P)%*%P
    
    grad[,i]      = c(DlDphi[,i],DlDlambda[,i],DlDdelta[,i])
    grads[,,i]    = grad[,i]%*%t(grad[,i])
    G             = G + grads[,,i]
    
    
    
  }
  
  
  
  H    = rbind(cbind(H11,H12),cbind(t(H12),H22))             # analytical Hessian
  
  
  
  # ---------------------------------------------------------------------------
  # calculation of asymptotic variance 
  
  
  
  gradi       = G/nT                                         # quasi-information matrix calculated by the procudt of gradient
  invHe       = try(solve(H/nT), silent = TRUE)
  if(class(invHe)[1]=="try-error")
  {
    message("Caught an error 4")
    return(list(error=4, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0,n,n)))
  }else{
    Avar    = (invHe%*%gradi%*%invHe)/nT
    rob.se  = sqrt(diag(Avar))                                      # robust standard errors
    out.se  = try(sqrt(diag(solve(gradi)/nT)),silent = TRUE)        # standard errors based on the outer product of the gradient
    if(class(out.se)[1]=="try-error")
    {
      message("Caught an error 4")
      return(list(error=4, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0,n,n)))
    }else{
      H.se   = sqrt(diag(invHe/nT))                                  # standard errors based on the inverted Hessian
      Avar1  = invHe%*%gradi%*%invHe                                 # for wald test
      return(list(error=0,vari = cbind(H.se, out.se, rob.se), vari1 = Avar1))  
    }
  }
}



###################    local asymptotic variance  and model checking  ######################################################################


"LocalVAMC"<- function(data, par, type, include.AR = F, lag = 6)
{
  #########################################################################################################
  ##                                                                                                     ##
  ##  Description:                                                                                       ##
  ##  Local asymptotic variance for VARMA-GARCH model specified by LING & McAleer (2003)                 ##
  ##  Specified as AR(1)-GARCH(1,1) or mu-GARCH(1,1)                                                     ##
  ##  par --  self-weighted or local QMLE                                                                ##
  ##  data -- a nT-by-N matrix of data                                                                   ##
  ##  type -- only GARCH part is "diagonal", "half-diagonal", or "extended"                              ##
  ##                                                                                                     ##
  #########################################################################################################
  
  # -----------------------------------------------------------------------------------------
  # passin the arguments
  
  nT     = dim(data)[1]    # number of observations
  N      = dim(data)[2]    # number of dimensions
  n      = length(par)     # length of parameters 
  n1     = n - (N*(N-1)/2)   # length of parameters except correlation coefficients
  
  
  
  if (include.AR) {
    para.mat = p.match(par, include.AR, type, N)
    heps = .Call("garch_ar", data, para.mat$mu, para.mat$phi, para.mat$W, para.mat$A, para.mat$B)
  }else{
    para.mat = p.match(par, include.AR, type, N)
    heps = .Call("garch_mu", data, para.mat$mu, para.mat$W, para.mat$A, para.mat$B)
  }
  
  
  para.mat$Gamma <- make.pd(para.mat$Gamma)
  
  # constructing volatilities
  h       = heps[[1]]             # estimated volatility
  epsit   = heps[[2]]
  invh    = 1/h
  sq.h    = sqrt(h)
  sq.invh = 1/sqrt(h)
  eta     = heps[[2]]/sq.h     # eta or standardized residuals
  invR    = solve(para.mat$Gamma)
  Cn      = solve(para.mat$Gamma)*para.mat$Gamma + diag(N)
  
  
  # -----------------------------------------------------------------------
  # the score function and information matrix( i.e. the first derivatives and second derivatives)
  
  
  # partial derivatives of Gamma w.r.t. correlation coefficients
  
  DTGammaDdelta = matrix(0,(N*(N-1)/2),N*N)
  
  for (i in 1:(N*(N-1)/2)){
    delta = matrix(0,(N*(N-1)/2),1)
    delta[i] = 1
    Gamma1 = matrix(0, N, N)
    Gamma1[lower.tri(Gamma1)] = delta
    Gamma1 = t(Gamma1)
    Gamma1[lower.tri(Gamma1)] = delta
    DTGammaDdelta[i,]=as.vector(Gamma1)
  }
  
  
  P = (diag(N)%x%solve(para.mat$Gamma))%*%t(DTGammaDdelta)
  C1=matrix(0,N,(N*N))
  for (i in 1:N) {
    C1[i,((i-1)*N+i)]=1
  }
  
  
  
  r.ini=matrix(0,1,N)
  h.ini=matrix(0,1,N)
  epsi.ini=matrix(0,1,N)
  rt=rbind(r.ini,data)
  hthat=rbind(h.ini,h)
  epsithat=rbind(epsi.ini,epsit)
  
  
  
  
  if(include.AR){
    
    DTepsitDphi  = array(0,dim=c((N*(N+1)),N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,(N*(N+1)),(nT+1)))
    
    if(type=="diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+4)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+4)),(nT+1)))
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(N+4)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),((N*(N+4))+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+4)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2),nT))
      G      = matrix(0,(N*(N+4)+N*(N-1)/2),(N*(N+4)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(N+4)), (N*(N+4)))
      H12 = matrix(0, (N*(N+4)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
      
    }else if(type=="half-diagonal"){
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(2*N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(2*N+3)),(nT+1)))
      
      
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(2*N+3)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),((N*(2*N+3))+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,((N*(2*N+3))+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(2*N+3)+N*(N-1)/2),(N*(2*N+3)+N*(N-1)/2))
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(2*N+3)), (N*(2*N+3)))
      H12 = matrix(0, (N*(2*N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
      
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(3*N+2)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(3*N+2)),(nT+1)))
      
      
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(3*N+2)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),((N*(3*N+2))+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,(N*(N+1)),nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(3*N+2)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(3*N+2)+N*(N-1)/2),(N*(3*N+2)+N*(N-1)/2))
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(3*N+2)), (N*(3*N+2)))
      H12 = matrix(0, (N*(3*N+2)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]  = -t(cbind((diag(N)-para.mat$phi),t(rt[i-1,]-para.mat$mu)%x%diag(N)))
        DhtDTphi[,,i]     = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]  = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]   = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]      = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }else{
    
    DTepsitDphi  = array(0,dim=c(N,N,(nT+1)))
    DhtDTphi     = array(0,dim=c(N,N,(nT+1)))
    
    if(type=="diagonal"){
      
      DhtDTlambda  = array(0,dim=c(N,(3*N),(nT+1)))
      DTepsitDpl   = array(0,dim=c((4*N),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(4*N),(nT+1)))
      
      
      DvDTpl      =  array(0,dim=c((N*N),(4*N),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),((4*N)+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(3*N),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(4*N+N*(N-1)/2),nT)
      grads  = array(0, dim = c((4*N+N*(N-1)/2),(4*N+N*(N-1)/2),nT))
      G      = matrix(0, (4*N+N*(N-1)/2),(4*N+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (4*N), (4*N))
      H12 = matrix(0, (4*N), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),diag(epsithat[i-1,]**2),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(3*N),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else if(type=="half-diagonal"){
      
      DhtDTlambda  = array(0,dim=c(N,(N*(N+2)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((N*(N+3)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(N*(N+3)),(nT+1)))
      
      
      
      DvDTpl      =  array(0,dim=c((N*N),(N*(N+3)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),(N*(N+3)+N*(N-1)/2),nT))
      
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(N+2)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(N*(N+3)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2),nT))
      G      = matrix(0, (N*(N+3)+N*(N-1)/2),(N*(N+3)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (N*(N+3)), (N*(N+3)))
      H12 = matrix(0, (N*(N+3)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),t(epsithat[i-1,]**2)%x%diag(N),diag(hthat[i-1,]))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(N+2)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }else{
      
      DhtDTlambda  = array(0,dim=c(N,(N*(2*N+1)),(nT+1)))
      DTepsitDpl   = array(0,dim=c((2*N*(N+1)),N,(nT+1)))
      DhtDTpl      = array(0,dim=c(N,(2*N*(N+1)),(nT+1)))
      
      
      DvDTpl      =  array(0,dim=c((N*N),(2*N*(N+1)),nT))
      DvDTdelta   =  array(0,dim=c((N*N),(N*(N-1)/2),nT))
      DvDTtheta   =  array(0,dim=c((N*N),(2*N*(N+1)+N*(N-1)/2),nT))
      
      
      # the score functions ( i.e. the first derivatives )
      DlDphi       = matrix(0,N,nT)
      DlDlambda    = matrix(0,(N*(2*N+1)),nT)
      DlDdelta     = matrix(0,(N*(N-1)/2),nT)
      
      
      grad   = matrix(0,(2*N*(N+1)+N*(N-1)/2),nT)
      grads  = array(0, dim = c((2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2),nT))
      G      = matrix(0, (2*N*(N+1)+N*(N-1)/2),(2*N*(N+1)+N*(N-1)/2))
      
      
      # the Hessian matrix ( i.e. the second derivatives )
      H11 = matrix(0, (2*N*(N+1)), (2*N*(N+1)))
      H12 = matrix(0, (2*N*(N+1)), (N*(N-1)/2))
      H22 = matrix(0, (N*(N-1)/2), (N*(N-1)/2))
      
      
      
      for (i in 2:(nT+1)){
        DTepsitDphi[,,i]   = -diag(N)
        DhtDTphi[,,i]      = 2*para.mat$A%*%(diag(epsithat[i-1,])%*%t(DTepsitDphi[,,i-1]))+para.mat$B%*%DhtDTphi[,,i-1]
        DhtDTlambda[,,i]   = cbind(diag(N),cbind(t(epsithat[i-1,]**2),t(hthat[i-1,]))%x%diag(N))+para.mat$B%*%DhtDTlambda[,,i-1]
        DTepsitDpl[,,i]    = rbind(DTepsitDphi[,,i],matrix(0,(N*(2*N+1)),N))
        DhtDTpl[,,i]       = cbind(DhtDTphi[,,i],DhtDTlambda[,,i])
      }
    }
  }
  
  
  
  
  DTepsitDphi  = DTepsitDphi[,,(2:(nT+1))]
  DhtDTphi     = DhtDTphi[,,(2:(nT+1))]
  DhtDTlambda  = DhtDTlambda[,,(2:(nT+1))]
  DTepsitDpl   = DTepsitDpl[,,(2:(nT+1))]
  DhtDTpl      = DhtDTpl[,,(2:(nT+1))]
  
  
  invDRD  =  array(0,dim=c(N,N,nT))
  for (i in 1:nT){
    D              = diag(sq.h[i,])
    invD           = diag(sq.invh[i,])
    invH           = diag(invh[i,])
    z              = eta[i,]%o%eta[i,]
    invDRD[,,i]    = invD%*%invR%*%invD
    zeta           = matrix(1,N,1)-diag(eta[i,])%*%invR%*%eta[i,]
    DlDphi[,i]     = DTepsitDphi[,,i]%*%invDRD[,,i]%*%epsit[i,]+0.5*t(DhtDTphi[,,i])%*%invH%*%zeta
    DlDlambda[,i]  = 0.5*t(DhtDTlambda[,,i])%*%invH%*%zeta
    DlDdelta[,i]   = 0.5*DTGammaDdelta%*%as.vector(invR - invR%*%z%*%invR)
    
    
    
    H11           = H11 + (DTepsitDpl[,,i]%*%invDRD[,,i]%*%t(DTepsitDpl[,,i]) + 0.25*t(DhtDTpl[,,i])%*%invH%*%Cn%*%invH%*%DhtDTpl[,,i])
    H12           = H12 + (0.5*t(DhtDTpl[,,i])%*%invH%*%C1%*%P)
    H22           = H22 + 0.5*t(P)%*%P
    
    grad[,i]      = c(DlDphi[,i],DlDlambda[,i],DlDdelta[,i])
    grads[,,i]    = grad[,i]%*%t(grad[,i])
    G             = G + grads[,,i]
    
    
    # for model checking 
    # partition matrix "DhtDTpl" into small matrices
    
    partitioned_matrices =  matrix(0, nrow = N, ncol = n1)
    matrix_original      =  DhtDTpl[,,i]
    for (ii in 1:N) {
      partitioned_matrices[ii, ] <- matrix_original[ii, ]
    }
    
    # construct quasi-diagonal matrix
    
    quasi_diagonal_matrix <- matrix(0, nrow = N, ncol = N * n1)
    for (ii in 1:N) {
      start_col <- (ii - 1) * n1 + 1
      end_col <- ii * n1
      quasi_diagonal_matrix[ii, start_col:end_col] <- partitioned_matrices[ii, ]
    }
    quasi_diagonal_matrix = 0.5*invD%*%quasi_diagonal_matrix
    
    #  vectorized the quasi-diagonal matrix 
    vectorized_matrix <- matrix(0, nrow = N * N, ncol = n1)
    for (ii in 1:(N * N)) {
      start_col <- ((ii - 1) %% N) * n1 + 1
      end_col <- ((ii - 1) %% N + 1) * n1
      vectorized_matrix[ii, ] <- quasi_diagonal_matrix[(ii - 1) %/% N + 1, start_col:end_col]
    }
    
    
    DvDTpl[,,i]      =  (diag(N)%x%(D%*%para.mat$Gamma))%*%vectorized_matrix+((D%*%para.mat$Gamma)%x%diag(N))%*%vectorized_matrix
    DvDTdelta[,,i]   =  (D%x%D)%*%t(DTGammaDdelta)
    DvDTtheta[,,i]   =  cbind(DvDTpl[,,i],DvDTdelta[,,i])
    
  }
  
  
  
  H    = rbind(cbind(H11,H12),cbind(t(H12),H22))             # analytical Hessian
  
  
  
  # ---------------------------------------------------------------------------
  # calculation of asymptotic variance 
  
  
  
  teta   =  matrix(0,1,nT)
  for (ii in 1:nT) {
    teta[ii] = t(epsit[ii, ])%*%invDRD[,,ii]%*%epsit[ii,]-N
  }
  teta   = as.vector(teta)
  k0     = sum(teta^2)
  
  
  R         =  matrix(0,nrow=1,ncol=lag)
  for (i in 1:lag) {
    x11     =  teta[(i+1):nT]
    x22     =  teta[1:(nT-i)]
    R[i]    =  t(x11)%*%x22
    R[i]    =  R[i]/k0
  }
  
  k   = k0/(nT*N)
  X   = matrix(0,nrow=n,ncol=lag)
  xl  = matrix(0,nrow=n,ncol=1)
  
  for (i in 1:lag) {
    x111     =  DvDTtheta[,,(i+1):nT]
    x222     =  invDRD[,,(i+1):nT]
    x333     =  teta[1:(nT-i)]
    for (ii in 1:(nT-i)) {
      xl     =  xl + t(x111[,,ii])%*%as.vector(x222[,,ii]*x333[ii])
    }
    X[,i]    =  xl/(nT-i)
  }
  
  
  
  Z1 = matrix(0, nrow = nT-6, ncol = 1)
  Z2 = matrix(0, nrow = nT-6, ncol = 1)
  Z3 = matrix(0, nrow = nT-6, ncol = 1)
  Z4 = matrix(0, nrow = nT-6, ncol = 1)
  Z5 = matrix(0, nrow = nT-6, ncol = 1)
  Z6 = matrix(0, nrow = nT-6, ncol = 1)
  Z7 = matrix(0, nrow = nT-6, ncol = n)
  Z  = matrix(0, nrow = nT-6, ncol = 6+n)
  
  for (i in 7:nT) {
    Z1[i-6,] = teta[i]*teta[i-1]
    Z2[i-6,] = teta[i]*teta[i-2]
    Z3[i-6,] = teta[i]*teta[i-3]
    Z4[i-6,] = teta[i]*teta[i-4]
    Z5[i-6,] = teta[i]*teta[i-5]
    Z6[i-6,] = teta[i]*teta[i-6]
    Z7[i-6,] = -t(grad[,i])
  }
  Z  = cbind(Z1,Z2,Z3,Z4,Z5,Z6,Z7)
  Z0 = (t(Z)%*%Z)/(nT-6)
  
  
  
  
  gradi       = G/nT                                         # quasi-information matrix calculated by the procudt of gradient
  invHe       = try(solve(H/nT), silent = TRUE)
  # Hessian     = hessian(LLT, par, method="Richardson")
  # invHe = try(solve(Hessian), silent = TRUE)
  if(class(invHe)[1]=="try-error")
  {
    message("Caught an error 4")
    return(list(error=4, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0,n,n), QMM = 0))
  }else{
    Avar    = (invHe%*%gradi%*%invHe)/nT
    rob.se  = sqrt(diag(Avar))                                      # robust standard errors
    out.se  = try(sqrt(diag(solve(gradi)/nT)),silent = TRUE)        # standard errors based on the outer product of the gradient
    if(class(out.se)[1]=="try-error")
    {
      message("Caught an error 4")
      return(list(error=4, vari = cbind(matrix(0,n,1), matrix(0,n,1), matrix(0,n,1)), vari1 = matrix(0,n,n), QMM = 0 ))
    }else{
      H.se   = sqrt(diag(invHe/nT))                                  # standard errors based on the inverted Hessian
      Avar1  = invHe%*%gradi%*%invHe                                 # for wald test
      V      = cbind(diag(lag), -t(X)%*%invHe)
      MV     = (V%*%Z0%*%t(V))/((k0/nT)^2)
      # MV     = diag(lag) + t(X)%*%(-k*invHe+Avar1)%*%X/((k0/nT)^2)
      # QMM    = nT*R%*%solve(MV)%*%t(R)
      QMM    = (nT-lag-2*N-1-2+1)*R%*%solve(MV)%*%t(R)
      return(list(error=0,vari = cbind(H.se, out.se, rob.se), vari1 = Avar1, QMM = QMM))  
    }
  }
}


############################################################################################################################################


