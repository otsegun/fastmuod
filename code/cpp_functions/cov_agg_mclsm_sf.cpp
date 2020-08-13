
/*
 * Code originally adapted from:
 *   http://systematicinvestor.github.io/Run-Correlation-Rcpp 
 *
 * version: 0.9
 */

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline NumericMatrix c_cov_helper_sf(const NumericMatrix& mat, const NumericMatrix& matref, const int rstart, const int rend,
                                  NumericVector means, NumericVector vars, NumericVector sds, const int cstart, const int cend) {
  const int nc = matref.ncol();
  const int nperiod = rend - rstart;
  const int ksize = cend - cstart + 1;
  NumericMatrix routl(ksize, 3);
  
  // pre-compute outliers
  for (int cur = 0; cur < nc; cur++) {
    for (int col = 0; col < ksize; col++) {
      long double sXY = 0;
      long double aux;
      
      for (int r = rstart; r < rend; r++)
        sXY += matref(r, cur) * mat(r, col + cstart);
      aux = sXY / (nperiod - 1);
      // moving average
      routl(col, 0) += (aux / sds[cur] - routl(col, 0)) / (cur+1); //aggregated mean
      routl(col, 1) += (aux * (means[cur] / vars[cur]) - routl(col, 1)) / (cur+1); //aggregated mean
      routl(col, 2) += (aux / vars[cur] - routl(col, 2)) / (cur+1);
    }
  }

  return routl;
}

// [[Rcpp::export]]
NumericMatrix corCovAggSf(NumericMatrix mtx, NumericMatrix mtxref, NumericVector refmean, NumericVector refvar,
                          NumericVector refsds,  int colstart, int colend) {
  return c_cov_helper_sf(mtx, mtxref,  0, mtxref.nrow(), refmean, refvar, refsds, colstart-1, colend-1);
}

