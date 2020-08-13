
/*
 * Code originally adapted from:
 *   http://systematicinvestor.github.io/Run-Correlation-Rcpp 
 *
 * version: 0.9
 */

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline NumericMatrix c_cov_helper_fast(const NumericMatrix& mat, const NumericVector& matref, const int rstart, const int rend,
                                     double means, double vars, double sds, const int cstart, const int cend) {
  //const int nc = matref.length();
  const int nperiod = rend - rstart;
  const int ksize = cend - cstart + 1;
  NumericMatrix routl(ksize, 3);
  
  // pre-compute outliers
    for (int col = 0; col < ksize; col++) {
      long double sXY = 0;
      long double aux;
      
      for (int r = rstart; r < rend; r++)
        sXY += matref(r) * mat(r, col + cstart);
      aux = sXY / (nperiod - 1);
      // moving average
      routl(col, 0) = aux / sds; //aggregated mean
      routl(col, 1) = aux * (means / vars); //not aggregating 
      routl(col, 2) = aux / vars;
    }
  
  return routl;
}

// [[Rcpp::export]]
NumericMatrix corCovAggFast(NumericMatrix mtx, NumericVector mtxref, double refmean, double refvar,
                          double refsds,  int colstart, int colend) {
  return c_cov_helper_fast(mtx, mtxref,  0, mtx.nrow(), refmean, refvar, refsds, colstart-1, colend-1);
}

