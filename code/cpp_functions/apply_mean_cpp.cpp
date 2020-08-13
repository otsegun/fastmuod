#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector apply_mean_cpp(NumericMatrix x, int dim){
  if(dim==1){
    NumericVector output(x.nrow());
    for(int i=0;i<x.nrow();i++){
      NumericVector temp = x(i,_);
      output[i] = mean(temp);
    }    
    return output ;
  }
  else if(dim==2){
    NumericVector output(x.ncol());
    for(int i=0; i<x.ncol(); i++){
      NumericVector temp = x(_, i);
      output[i] = mean(temp);
    }
    return output ;
  }
} 



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
