#include <RcppArmadillo.h>
#include<iostream>
using namespace Rcpp;
using namespace std;
using namespace arma;
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "plinkfun.hpp"

// pass by object, which keeps the values in the original matrix
// [[Rcpp::export]]
RcppExport SEXP scaleC(arma::mat X) {
  long long n = X.n_rows, p = X.n_cols;
  // mat Xtmp(n,p);
  rowvec Xm = mean(X, 0);
  rowvec Xs = stddev(X, 0, 0); // see manual

  // X1tmp = (X1 - repmat(meanX1tmp, n1, 1))/ repmat(sdX1tmp, n1, 1) / sqrt(p1);
  X.each_row() -= Xm;
  X.each_row() /= Xs;
  X /= sqrt(p);

  List ret;
  ret["X"] = X;
  ret["Xs"] = Xs;
  ret["Xm"] = Xm;
  return ret;
}


// pass by reference, which changes the values in the original matrix
// [[Rcpp::export]]
RcppExport SEXP scaleC_ref(arma::mat& X) {
  long long n = X.n_rows, p = X.n_cols;
  // mat Xtmp(n,p);
  rowvec Xm = mean(X, 0);
  rowvec Xs = stddev(X, 0, 0); // see manual

  // X1tmp = (X1 - repmat(meanX1tmp, n1, 1))/ repmat(sdX1tmp, n1, 1) / sqrt(p1);
  X.each_row() -= Xm;
  X.each_row() /= Xs;
  X /= sqrt(p);

  List ret;
  ret["Xs"] = Xs;
  ret["Xm"] = Xm;
  return ret;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
