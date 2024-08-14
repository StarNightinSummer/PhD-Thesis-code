# constructing parameter vectors and matrices
# input: parameters, include.AR, type and dimension
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

