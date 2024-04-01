// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <stdio.h>
#include <typeinfo>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::DiagonalMatrix;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

struct lamRet {
  int a;
  // VectorXd<double> lambdaT;
  // vector<VectorXd> h_store;
};
typedef struct lamRet Struct;

// [[Rcpp::export]]
string timesTwo(VectorXd y, MatrixXd x) {
  VectorXd store(2);
  // Struct s;
  // cout << x;
  // s.a = 2;
  // s.h_store[0] = y;
  // VectorXd t(3);
  // store[1] = y;
  // for(int i=0;i<4;i++) {
  //   store[i] = y;
  //   // t = t + 2 * y;
  // }
  store = x.row(1);
  // store.diagonal() = y;
  // if (y > 3) {
  //   
  // }
  // store = y.asDiagonal();
  // store = x.inverse();
  // store.conservativeResize(store.rows() + 1, store.cols());
  // store.row(2) = y;
  // return store.colwise().sum();
  return typeid(store.transpose()).name();
}

// int test(VectorXd y, MatrixXd x) {
//   Struct res;
//   res = timesTwo(y, x);
//   return res.a;
// }
// [[Rcpp::export]]
vector<int> test() {
  vector<int> s(5);
  s[2] = 6;
  return s;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test()
timesTwo(c(9,3),matrix(c(42,2,6,7),nrow=2))
*/
