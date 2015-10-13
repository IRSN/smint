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
 * Grid interpolation by looping over dimension.                              *
 *                                                                            *
 * NOTE: some space could be saved by assuming that 'xLev' has length d - 1   *
 * instead of d and removing the d-th element in the list before calling.     *
 *                                                                            *
 *============================================================================*/

SEXP grid_int(SEXP F,          // vector array of response values dim 'dim'
	      SEXP nLev,       // vector of levels numbers
	      SEXP nStar,      // c(1, cumprod(nLev))
	      SEXP xLev,       // list of levels
	      SEXP Xout,       // matrix of output values nout * d
	      SEXP interpFun,  // interpolation fun with formals 'y', 'x', 'xout'
	      SEXP rho) {      // an R environment
  
  SEXP dimXout, fs, R_fcall, fout, foutprov, xout, Fprov;

  int  i, j, k, ell, d = length(nLev), nout,
    *inLev = INTEGER(nLev), *inStar = INTEGER(nStar); 

  double *rF = REAL(F);

  // check inputs
  if (!isVector(F)) error("'F' must be a vector");
  if (!isVector(nLev)) error("'nLev' must be a vector");
  if (!isVector(nStar)) error("'nStar' must be a vector");
  if (!isNewList(xLev)) error("'xLev' must be a list");  
  if (!isMatrix(Xout)) error("'Xout' must be a matrix");
  if (!isFunction(interpFun)) error("'interpFun' must be a function");
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
  PROTECT(xout = allocVector(REALSXP, 1));
  PROTECT(fout = allocVector(REALSXP, nout)); 
  PROTECT(foutprov = allocVector(REALSXP, 1));
  PROTECT(Fprov = allocVector(REALSXP, nout * inStar[d - 1]));
  
  // fs is a list of length d - 1 used temporarily to store 
  // the function values to be interpolated. 
  // would 'shallow_duplicate' be enough?
  PROTECT(fs = duplicate(xLev));

  PROTECT(R_fcall = lang4(interpFun, VECTOR_ELT(fs, 0), 
			  VECTOR_ELT(fs, 0), xout));

  double *rFprov = REAL(Fprov);

  //===========================================================================
  // Perform 'd - 1' 1D-interpolations starting from j = d - 2. The
  // interpolation for j = d - 1 was done before in the R caller and
  // the first dimension for the array F is used for the 'xout'
  // dimension
  //
  // o Each interpolation uses formals 'y', 'x', and 'xout' IN THAT ORDER! 
  //  
  // o Each interpolation has 'xout' of length 1, hence its result of
  // length 1.
  //
  // o Within the j loop, 'Fprov' is an array with dimensions
  //  
  //      dim  = nout * inLev[0] * inLev[1] * ... * inLev[j]
  //           = nout * inStar[j + 1]
  // 
  //===========================================================================
  
  for (k = 0; k < nout; k++) {
    
    // initialize temporary array 'Fprov' as a slice of 'F'
    for (i = 0; i < inStar[d - 1]; i++) {
      rFprov[i] = rF[ k + nout * i ];
    }

    for (j = d - 2; j >= 0; j--) {
    
      // pick 'xout' from within 'Xout' 
      REAL(xout)[0] = REAL(Xout)[k + j * nout];
      SETCADDDR(R_fcall, xout);

      // take xLev[[j]] as 'x' for call as the formal #2
      SETCADDR(R_fcall, VECTOR_ELT(xLev, j));
      
      // this is the 'apply' part 
      for (i = 0; i < inStar[j]; i++) {
	
	for (ell = 0; ell < inLev[j]; ell++) { 
	  REAL(VECTOR_ELT(fs, j))[ell] = rFprov[i + inStar[j] * ell];
	}
	
	SETCADR(R_fcall, VECTOR_ELT(fs, j));
	foutprov = eval(R_fcall, rho);
	rFprov[i] = REAL(foutprov)[0];
	
      }
      
    }

    REAL(fout)[k] = rFprov[0];

#ifdef DEBUG 
    Rprintf("Interpolation done for j = %d\n", j);
#endif 
    
  }
  
  UNPROTECT(7);
  return(fout);

}












