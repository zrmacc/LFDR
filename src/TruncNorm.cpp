// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

//' Newton-Raphson Update
//' 
//' 
//' @param Info Information matrix
//' @param S Score vector
//'   
// [[Rcpp::export]]
SEXP Delta(const Eigen::Map<Eigen::MatrixXd> Info, const Eigen::Map<Eigen::VectorXd> S){
  const Eigen::VectorXd d = Info.llt().solve(S);
  return Rcpp::wrap(d);
}