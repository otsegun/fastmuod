#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

  Rcpp::List cov_armadillo_mclsm(const arma::uvec i, const arma::mat mtx_i, const arma::mat mtx_ref,
                          const double ref_mean, const double ref_var, 
                          const double ref_sds, const arma::vec means, const arma::vec sds){
    arma::mat out(mtx_i.n_cols, mtx_ref.n_cols);
    out = arma::cov(mtx_i, mtx_ref);
    const int n = out.n_rows;
    const int p = out.n_cols;
    //declare matrices to hold result
    arma::mat amplitude(n, p);
    arma::mat magnitude(n, p);
    arma::mat shape(n, p);
    // compute indices
    amplitude = out/ref_var;
    magnitude = arma::repelem(means.elem(i), 1, out.n_cols) - (amplitude * ref_mean);
    shape = (out/ref_sds)/arma::repelem(sds.elem(i), 1, out.n_cols);
    List answer = List::create(Named("amplitude") = amplitude, 
                               Named("magnitude") = magnitude,
                               Named("shape") = shape);
    return answer;
  }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
