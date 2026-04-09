#include <RcppEigen.h>
#include <Eigen/Dense>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

//' Create a Kernel matrix from a design matrix G
//' 
//' @param G A numeric matrix.
//' @return A normalized kernel matrix K.
//' @export
// [[Rcpp::export]]
SEXP createK(const Eigen::Map<Eigen::MatrixXd> & G) {
  int n = G.cols();
  Eigen::VectorXd colSum = G.colwise().sum();
  
  for (int i = 0; i < n; ++i) {
    if (colSum(i) != 0.0) {
      colSum(i) = 1.0 / colSum(i);
    }
  }
  
  Eigen::MatrixXd betaG = colSum.asDiagonal();
  Eigen::MatrixXd K = G * betaG * G.adjoint();
  
  return Rcpp::wrap(K);
}