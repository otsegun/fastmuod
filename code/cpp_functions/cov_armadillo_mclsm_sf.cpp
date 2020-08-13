#include <RcppArmadillo.h>
using namespace Rcpp;


//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat covArmadilloMclsmSf(const arma::mat mtxi, const arma::mat mtxref, int offset,
                               const arma::vec refmean, const arma::vec refvar, 
                               const arma::vec refsds, const arma::vec means, const arma::vec sds){
  arma::mat out = arma::cov( mtxref, mtxi);
  const int n = out.n_rows;
  const int p = out.n_cols;
  //declare matrices to hold result
  arma::vec temp(n);
  arma::mat preout(3,p);
  for(int i = 0; i<p ; ++i){
    temp = out.col(i)/refvar;
    preout(0,i) = arma::mean((out.col(i)/refsds))/sds(i + offset); //shape
    preout(1, i) = means(i + offset) - arma::mean(temp % refmean); //mag
    preout(2, i) = arma::mean(temp); //amplitude
  }
  
  return preout.t();
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
