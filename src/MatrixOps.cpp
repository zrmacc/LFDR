// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//' Matrix Matrix Product
//'
//' Calculates the product \eqn{AB}. 
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix.
//' @export
// [[Rcpp::export]]

SEXP MMP(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
  const Eigen::MatrixXd C = A*B;
  return Rcpp::wrap(C);
}

//' Matrix Inner Product
//'
//' Calculates the inner product \eqn{A'B}.
//'
//' @param A Numeric matrix.
//' @param B Numeric matrix.
//' @return Numeric matrix. 
// [[Rcpp::export]]

SEXP matIP(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B){
  const Eigen::MatrixXd AtB = (A.transpose() * B);
  return Rcpp::wrap(AtB);
}

//' Matrix Inverse
//'
//' Calcualtes \eqn{A^{-1}}.
//' 
//' @param A Numeric matrix.
//' @return A numeric matrix. 
// [[Rcpp::export]]

SEXP matInv(const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd Ai = A.completeOrthogonalDecomposition().pseudoInverse();
  return Rcpp::wrap(Ai);
}

//' Matrix Determinant
//'
//' Calculates the determinant of matrix \eqn{A}. 
//'
//' @param A Numeric matrix.
//' @return Scalar. 
// [[Rcpp::export]]

SEXP det(const Eigen::Map<Eigen::MatrixXd> A){
  const double d = A.determinant();
  return Rcpp::wrap(d);
}

//' Matrix Quadratic Form
//' 
//' Calculates the quadratic form \eqn{X'AX}.
//' 
//' @param X Numeric matrix.
//' @param A Numeric matrix.
//' @return Numeric matrix
// [[Rcpp::export]]

SEXP matQF(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> A){
  const Eigen::MatrixXd q = X.transpose()*A*X;
  return Rcpp::wrap(q);
}

//' Schur complement
//'
//' Calculates the efficient information \eqn{I_{bb}-I_{ba}I_{aa}^{-1}I_{ab}}. 
//'
//' @param Ibb Information of target parameter
//' @param Iaa Information of nuisance parameter
//' @param Iba Cross information between target and nuisance parameters
//' @return Numeric matrix. 
// [[Rcpp::export]]
SEXP SchurC(const Eigen::Map<Eigen::MatrixXd> Ibb, const Eigen::Map<Eigen::MatrixXd> Iaa,
            const Eigen::Map<Eigen::MatrixXd> Iba){
  // Kernel matrix
  const Eigen::MatrixXd E = Ibb-(Iba*(Iaa.ldlt().solve(Iba.transpose())));
  return Rcpp::wrap(E);
}