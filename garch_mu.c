#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemv and daxpy */

SEXP garch_mu(SEXP Data, SEXP Mean, SEXP Const, SEXP Arch, SEXP Garch)
{
  int i, j, nobs, ndim, ione = 1;
  double *rx, *rmu, *rw, *rA, *rB, *rh, *reps, *rtmp, *rtmp1, *rx1, *rx2, one = 1.0, zero = 0.0;
  SEXP x, mu, w, A, B, h, eps, tmp, tmp1, x1, x2, output;
  
  nobs = Rf_nrows(Data); ndim = Rf_ncols(Data);
  PROTECT(x = duplicate(Data));
  PROTECT(mu = duplicate(Mean));
  PROTECT(w = duplicate(Const));
  PROTECT(A = duplicate(Arch));
  PROTECT(B = duplicate(Garch));
  PROTECT(h = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(eps = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(tmp = allocVector(REALSXP, ndim));
  PROTECT(tmp1 = allocVector(REALSXP, ndim));
  PROTECT(x1 = allocVector(REALSXP, ndim));
  PROTECT(x2 = allocVector(REALSXP, ndim));
  PROTECT(output = allocVector(VECSXP, 2));
  
  rx = REAL(x);
  rmu = REAL(mu);
  rw = REAL(w);
  rA = REAL(A);
  rB = REAL(B);
  rh = REAL(h);
  reps = REAL(eps);
  rtmp = REAL(tmp); 
  rtmp1 = REAL(tmp1);
  rx1 = REAL(x1);   /* eps_{t-1}^2*/
  rx2 = REAL(x2);   /* h_{t-1} */
  
  /* Following Ling and McAleer(2003), initial values for conditional log likelihood, epsi_{0}, h_{0}, are zeros here*/
  for(j=0; j<ndim; j++){
    rx1[j] = 0.0;
    rx2[j] = 0.0;
  }
  
  /* recursion from T=1 to the end */
  for(i=0; i<nobs; i++){
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rA, &ndim, rx1, &ione, &zero, rtmp, &ione);  /* rtmp = A*eps_{t-1}^2 */
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rB, &ndim, rx2, &ione, &zero, rtmp1, &ione);  /* rtmp1 = B*h_{t-1} */
    
    for(j=0; j<ndim; j++){
      rh[i+j*nobs] = rw[j]+rtmp[j]+rtmp1[j];   /* h_{t} = w + A*eps_{t-1}^2 + B*h_{t-1} */
      reps[i+j*nobs] = rx[i+j*nobs]-rmu[j];   /* eps_{t} = r_{t} - mu */
      rx1[j] = (rx[i+j*nobs]-rmu[j])*(rx[i+j*nobs]-rmu[j]);    /* for the next loop */
      rx2[j] = rw[j]+rtmp[j]+rtmp1[j];         /* for the next loop */
    }
  }
  
  SET_VECTOR_ELT(output, 0, h);
  SET_VECTOR_ELT(output, 1, eps);
  
  UNPROTECT(12);
  return(output);
}
