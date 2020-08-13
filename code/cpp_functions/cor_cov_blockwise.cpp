
/*
 * Code originally adapted from:
 *   http://systematicinvestor.github.io/Run-Correlation-Rcpp 
 *
 * version: 0.9
 */

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

inline NumericMatrix c_cov_helper(const NumericMatrix& mat, const int rstart, const int rend
  , NumericVector means, NumericVector vars, NumericVector sds
  , const int cstart, const int cend) {
  const int nc = mat.ncol();
  const int nperiod = rend - rstart;
  const int ksize = cend - cstart + 1;
  NumericMatrix routl(ksize, 3);

  // pre-compute outliers
  for (int cur = 0; cur < nc; cur++) {
    for (int col = 0; col < ksize; col++) {
      long double sXY = 0;
      long double aux;

      for (int r = rstart; r < rend; r++)
        sXY += mat(r, cur) * mat(r, col + cstart);
      aux = sXY / (nperiod - 1);
      // moving average
      routl(col, 0) += (aux / sds[cur] - routl(col, 0)) / (cur+1); //aggregated mean
      routl(col, 1) += (aux * (means[cur] / vars[cur]) - routl(col, 1)) / (cur+1); //aggregated mean
      routl(col, 2) += (aux / vars[cur] - routl(col, 2)) / (cur+1);
      /* classical average; step 1, accumulate
      // average (sum)
      routl(col, 0) += aux / sds[cur]; //aggregated mean
      routl(col, 1) += aux * (means[cur] / vars[cur]); //aggregated mean
      routl(col, 2) += aux / vars[cur];
      */
    }
  }
  // perform mean over routl
  /* classical average; step 2, divide by nc
  for (int col = 0; col < ksize; col++) {
    routl(col, 0) /= nc;
    routl(col, 1) /= nc;
    routl(col, 2) /= nc;
  }
  */
  return routl;
}

// [[Rcpp::export]]
NumericMatrix cor_cov_blockwise(NumericMatrix mat
  , NumericVector means, NumericVector vars, NumericVector sds
  , int c_ini, int c_end) {
  return c_cov_helper(mat, 0, mat.nrow(), means, vars, sds, c_ini-1, c_end-1);
}

