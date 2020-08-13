#include <RcppArmadillo.h>
using namespace Rcpp;



//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat covArmadillo(const arma::mat mtxi, const arma::mat mtxref){
  arma::mat out(mtxi.n_cols, mtxref.n_cols);
  
  // Degenerate cases 
  if (mtxi.n_cols == 0) {
    return out;
  } else if (mtxi.n_rows == 0 || mtxi.n_rows == 1) {
    out.fill(Rcpp::NumericVector::get_na());
  } else {
    out = arma::cov(mtxi, mtxref, 0);
  }
  
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
