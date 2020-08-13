#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector apply_median_cpp(NumericMatrix x, int dim){
  if(dim==1){
    NumericVector output(x.nrow());
    for(int i=0;i<x.nrow();i++){
      NumericVector temp = x(i,_);

      output[i] = median(temp, true);
    }    
    return output ;
  }
  else if(dim==2){
    NumericVector output(x.ncol());
    for(int i=0; i<x.ncol(); i++){
      NumericVector temp = x(_, i);
      output[i] = median(temp, true);
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
