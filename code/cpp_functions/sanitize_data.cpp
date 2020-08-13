#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix sanitizeData(NumericMatrix x) {
  int xrow = x.nrow();
  for(int i = 0; i<xrow ; i++){
    NumericVector col = x(i, _);
//    NumericVector new_col = col[!is_na(col)];
    double impute = mean(col);
    col[is_na(col)] = impute;
    x(i, _) = col;
    }
   return x; 
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/