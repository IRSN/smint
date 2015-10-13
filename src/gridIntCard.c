#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#define NODEBUG 

/*============================================================================*
 * AUTHOR                                                                     *
 *                                                                            *
 *    Yves Deville <deville.yves@alpestat.com>                                *
 *                                                                            *
 * DESCRIPTION                                                                *
 *                                                                            *
 * Grid interpolation by looping over dimension using a Cardinal Basis hence  *
 * a linear interpolation method.
 *                                                                            *
 *============================================================================*/

SEXP grid_int_CB(SEXP F,          // vector array of response values dim 'dim'
		 SEXP nLev,       // vector of levels numbers
		 SEXP nStar,      // c(1, cumprod(nLev))
		 SEXP xLev,       // list of levels
		 SEXP Xout,       // matrix of output values nout * d
		 SEXP interpCB ,  // Cardinal Basis
		 SEXP rho) {      // An R environment
  
  SEXP dimXout, x, f, h, R_fcall, fout, foutprov, xout, Fprov;

  int  i, j, k, ell, d = length(nLev), nout,
    *inLev = INTEGER(nLev), *inStar = INTEGER(nStar); 

  double g, *rF = REAL(F);

  // check inputs
  if (!isVector(F)) error("'F' must be a vector");
  if (!isVector(nLev)) error("'nLev' must be a vector");
  if (!isVector(nStar)) error("'nStar' must be a vector");
  //if (!isList(xLev)) error("'xLev' must be a list");  
  if (!isMatrix(Xout)) error("'Xout' must be a matrix");
  if (!isFunction(interpCB)) error("'interpCB' must be a function");
  if (!isEnvironment(rho)) error("'rho' should be an environment");
 
  F = coerceVector(F, REALSXP);
  
  if (length(nStar) != d + 1) error("'nStar' must be of length 'd + 1'");
  if (length(xLev) != d) error("'xLev' must be of length 'd'");
  
  // find the number of output values 
  PROTECT(dimXout = getAttrib(Xout, R_DimSymbol));
  
  if (INTEGER(dimXout)[1] != d) {
    error("number of columns for 'Xout' must be equal to the length of 'nLev'");
  }
  nout = INTEGER(dimXout)[0]; 
  
#ifdef DEBUG 
  Rprintf("nout = %d\n", nout);
#endif 
  
  // prepare SEXP
  PROTECT(R_fcall = lang3(interpCB, x, xout));
  PROTECT(xout = allocVector(REALSXP, 1));
  PROTECT(fout = allocVector(REALSXP, nout)); 
  PROTECT(foutprov = allocVector(REALSXP, 1));
  PROTECT(Fprov = allocVector(REALSXP,inStar[d]));
 
  // allocate h to the max used number of levels. Yet the length
  // of the R object may be set to a smaller value in the loop
  ell = inLev[d - 1];
  for (j = d - 2; j >= 0; j--) {
    if (inLev[j] > ell) ell = inLev[j];
  }
  PROTECT(h = allocVector(REALSXP, ell));

  double *rFprov = REAL(Fprov);

#ifdef DEBUG 
  Rprintf("Start loop\n");
#endif 

  //===========================================================================
  // Perform 'd' 1D-interpolations starting from j = d - 1 to j = 0.
  //
  // o Each Cardinal Basis determination uses two formals 'x' and
  // 'xout' IN THAT ORDER!
  //  
  // o Each interpolation has 'xout' of length 1, hence its is a matrix with 
  // one row.
  //
  // o Within the j loop, 'Fprov' is an array with dimensions
  //  
  //      dim  = inLev[0] * inLev[1] * ... * inLev[j]
  //
  // 
  //===========================================================================
  
  for (k = 0; k < nout; k++) {
    
    // intialize temporary array 'Fprov' as a slice of 'F'
    for (i = 0; i < inStar[d]; i++) {
      rFprov[i] = rF[i];
    }

    for (j = d - 1; j >= 0; j--) {
    
      // pick 'xout' from within 'Xout' and set it as formal #2
      REAL(xout)[0] = REAL(Xout)[k + j * nout];
      SETCADDR(R_fcall, xout);

      // take xLev[[j]] as 'x' for call as the formal #1
      SETCADR(R_fcall, VECTOR_ELT(xLev, j));
      h = eval(R_fcall, rho);
      // SETLENGTH(h, inLev[j]); 
    
#ifdef DEBUG 
      Rprintf("xout : %6.3f\n", REAL(xout)[0]);
      Rprintf("nLev : %6d\n", inLev[j]);
      for (ell = 0; ell < inLev[j]; ell++) {
	Rprintf("%6.3f ", REAL(VECTOR_ELT(xLev, j))[ell]);
      } 
      Rprintf("\n");
#endif 

      // this is the 'apply' part 
      for (i = 0; i < inStar[j]; i++) {
	g = 0.0;
	for (ell = 0; ell < inLev[j]; ell++) {
	  g += REAL(h)[ell] * rFprov[i + inStar[j] * ell]; 
	}
	rFprov[i] = g;
      }

    }

    REAL(fout)[k] = rFprov[0];
    
  }
  
  UNPROTECT(7);
  return(fout);

}












