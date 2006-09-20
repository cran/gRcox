/* File: matrixprodFun.c */
#include <Rinternals.h>
#include <R_ext/Applic.h> /* for dgemm */


static void matprod(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
		    x, &nrx, y, &nry, &zero, z, &nrx);
}


SEXP tr(SEXP x)
{
  int nrx, ncx, mode, i;
  SEXP tr, xdims;
  xdims = getAttrib(x, R_DimSymbol);
  mode = REALSXP;
  nrx = INTEGER(xdims)[0];
  ncx = INTEGER(xdims)[1];
  PROTECT(tr   = allocVector(mode,1));
  REAL(tr)[0] = 0;
  for (i=0; i< nrx; i++){
    REAL(tr)[0] = REAL(tr)[0] + REAL(x)[i*(nrx+1)];
  }
  UNPROTECT(1);
  return(tr);
}

SEXP trProd(SEXP x, SEXP y)
{
  int nrx, ncx, nry, ncy, mode;
  int i;
  SEXP tr;
  SEXP xdims, ydims, ans;
  xdims = getAttrib(x, R_DimSymbol);
  ydims = getAttrib(y, R_DimSymbol);
  mode = REALSXP;
  nrx = INTEGER(xdims)[0];
  ncx = INTEGER(xdims)[1];
  nry = INTEGER(ydims)[0];
  ncy = INTEGER(ydims)[1];
  PROTECT(ans = allocMatrix(mode, nrx, ncy));
  PROTECT(tr   = allocVector(mode,1));

  REAL(tr)[0] = 0;
  matprod(REAL(x), nrx, ncx, 
 	  REAL(y), nry, ncy, REAL(ans)); 
  
  for (i=0; i< nrx; i++){
    REAL(tr)[0] = REAL(tr)[0] + REAL(ans)[i*(nrx+1)];
  }
/*   Rprintf("%f", REAL(tr)[0]); */
  UNPROTECT(2);
  return(tr);
}

SEXP trProd2(SEXP x, SEXP y)
{
  int nrx, ncx, nry, ncy, mode, i;
  SEXP xdims, ydims, ans, ans2, tr;
  xdims = getAttrib(x, R_DimSymbol);
  ydims = getAttrib(y, R_DimSymbol);
  mode = REALSXP;
  nrx = INTEGER(xdims)[0];
  ncx = INTEGER(xdims)[1];
  nry = INTEGER(ydims)[0];
  ncy = INTEGER(ydims)[1];
  PROTECT(ans  = allocMatrix(mode, nrx, ncy));
  PROTECT(ans2 = allocMatrix(mode, nrx, ncy));
  PROTECT(tr   = allocVector(mode, 1));

  REAL(tr)[0] = 0;
  matprod(REAL(x), nrx, ncx, REAL(y), nry, ncy, REAL(ans)); 
  
  matprod(REAL(ans), nrx, ncy, REAL(ans), nrx, ncy, REAL(ans2)); 
  
  for (i=0; i< nrx; i++){
    REAL(tr)[0] = REAL(tr)[0] + REAL(ans2)[i*(nrx+1)];
  }

  UNPROTECT(3);
  return(tr);
}


SEXP trProd3(SEXP x, SEXP y, SEXP z)
{
  int nrx, ncx, nry, ncy, nrz, ncz, mode, i;
  SEXP xdims, ydims, zdims, ans, ans2, ans3, tr;
  xdims = getAttrib(x, R_DimSymbol);
  ydims = getAttrib(y, R_DimSymbol);
  zdims = getAttrib(z, R_DimSymbol);
  mode = REALSXP;
  nrx = INTEGER(xdims)[0];
  ncx = INTEGER(xdims)[1];
  nry = INTEGER(ydims)[0];
  ncy = INTEGER(ydims)[1];
  nrz = INTEGER(zdims)[0];
  ncz = INTEGER(zdims)[1];

  PROTECT(ans  = allocMatrix(mode, nrx, ncy));
  PROTECT(ans2 = allocMatrix(mode, nrz, ncy));
  PROTECT(ans3 = allocMatrix(mode, nrx, ncy));
  PROTECT(tr   = allocVector(mode, 1));

  REAL(tr)[0] = 0;

  /*   Calculate x * y; store result in ans */
  matprod(REAL(x), nrx, ncx, REAL(y), nry, ncy, REAL(ans)); 
  
  /*   Calculate z * y; store result in ans2  */
  matprod(REAL(z), nrz, ncz, REAL(y), nry, ncy, REAL(ans2)); 

  /*   Calculate ans * ans2; store result in ans3  */
  matprod(REAL(ans), nrx, ncy, REAL(ans2), nrz, ncy, REAL(ans3)); 

  /*   Calculate the trace */
  for (i=0; i< nrx; i++){
    REAL(tr)[0] = REAL(tr)[0] + REAL(ans3)[i*(nrx+1)];
  }

  UNPROTECT(4);
  return(tr);
}


SEXP trXYZ(SEXP x, SEXP y, SEXP z)
{
  int nrx, ncx, nry, ncy, nrz, ncz, mode, i;
  SEXP xdims, ydims, zdims, ans,  ans3, tr;
  xdims = getAttrib(x, R_DimSymbol);
  ydims = getAttrib(y, R_DimSymbol);
  zdims = getAttrib(z, R_DimSymbol);
  mode = REALSXP;
  nrx = INTEGER(xdims)[0];
  ncx = INTEGER(xdims)[1];
  nry = INTEGER(ydims)[0];
  ncy = INTEGER(ydims)[1];
  nrz = INTEGER(zdims)[0];
  ncz = INTEGER(zdims)[1];

  PROTECT(ans  = allocMatrix(mode, nrx, ncy));
  PROTECT(ans3 = allocMatrix(mode, nrx, ncy));
  PROTECT(tr   = allocVector(mode, 1));

  REAL(tr)[0] = 0;

  /*   Calculate x * y; store result in ans */
  matprod(REAL(x), nrx, ncx, REAL(y), nry, ncy, REAL(ans)); 
  
  /*   Calculate ans * z ; store result in ans3  */
  matprod(REAL(ans), nrx, ncy, REAL(z), nrz, ncz, REAL(ans3)); 

  /*   Calculate the trace */
  for (i=0; i< nrx; i++){
    REAL(tr)[0] = REAL(tr)[0] + REAL(ans3)[i*(nrx+1)];
  }

  UNPROTECT(3);
  return(tr);
}








