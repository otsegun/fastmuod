#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;


//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]]

arma::mat covArmadilloMclsmNormal(arma::mat dti, arma::mat dt, int offset,
                                     arma::vec dtmeans, arma::vec dtvars,
                                     arma::vec dtsds){
  arma::mat out = arma::cov(dt, dti);

//  int n = out.n_rows;
//  int p = out.n_cols;
  //declare matrices to hold result
//  arma::vec temp(n);
//  arma::mat preout(3,p);
//  for(int i = 0; i<p ; ++i){
//    temp = out.col(i)/dtvars;
//    preout(0,i) = arma::mean((out.col(i)/dtsds))/dtsds(i + offset); //shape
//    preout(1, i) = dtmeans(i + offset) - arma::mean(temp % dtmeans); //mag
//    preout(2, i) = arma::mean(temp); //amplitude
//  }
  return out;
//  return preout.t();
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#getwd()

#sourceCpp("apply_mean_cpp.cpp")
#sourceCpp("apply_median_cpp.cpp")
# sourceCpp("apply_sd_cpp.cpp")
# sourceCpp("apply_var_cpp.cpp")
# source("../sim_models_w_plots.R")
# 
# n <- 100000
# p <- 100
# alpha <- .15
# # generate data
# data_ob <- model2(n = n, p = p, out.rate = alpha)
# data <- t(data_ob$data)
# 
# data_means <- colMeans(data, na.rm = T)
# data_vars <- apply_var_cpp(data, 2)
# data_sds <- apply_sd_cpp(data, 2)
# 
# i <- 1:n
# 
# pre_outl <- covArmadilloMclsmNormal(dti = data[,i], dt = data, offset = min(i) - 1L,
#                         dtmeans = data_means, dtvars = data_vars, dtsds = data_sds)

*/
