#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp
//' @param sigma the stardard deviation of the random walk
//' @param x_0 the start point
//' @param N the number of samples
//' @return a random sample of size \code{n} as well as the number of rejections
//' @examples
//' \dontrun{
//' rnC <- gibbsC(100,10)
//' par(mfrow=c(2,1));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
List rw_MetropolisC(double sigma,double x_0, int N) {
  NumericVector x(N);
  x[0] = x_0;
  NumericVector u = runif(N);
  double y;
  int k = 0;
  for(int i = 1; i < N; ++i){
    y = R::rnorm(x[i-1],sigma);
    if (u[i-1] <= exp(abs(x[i-1]) - abs(y))){
      x[i] = y;
    } else{
      x[i] = x[i-1];
      ++k;
    }
  }
  return List::create(
    _["x"] = x,
    _["k"] = k
  );
}
