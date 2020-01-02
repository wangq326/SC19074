#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericMatrix rwme(double sigma,double xo,int N){
  NumericMatrix x(N,2);
  x(0,0) = xo; 
  x(1,0) = 1;
  NumericVector u = runif(N); 
  for (int i = 1; i < N ;i++){ 
    double y = as<double>(rnorm(1, x(i-1,0), sigma));
    double t = exp(-abs(y))/exp(-abs(x(i-1,0)));
    if (u[i] > t){
      x(i,0) = x(i-1,0); 
      x(i,1) = 0;
    }
    else{ 
      x(i,0) = y;
      x(i,1) = 1;} 
  };
  return x;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
