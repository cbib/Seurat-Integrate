#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
int64_t n_zeros_dense_mat(const arma::mat &mat) {
  int64_t count = arma::accu(mat == 0);
  return count;
}

// [[Rcpp::export]]
int64_t n_zeros_sparse_mat(arma::sp_mat &mat) {
  int64_t count = mat.clean(0).n_elem - mat.n_nonzero;
  return count;
}

