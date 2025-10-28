#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Soft-thresholding function, returns scalar
// [[Rcpp::export]]
double soft_c(double a, double lambda){
  // Your function code goes here
  double abs_a = std::abs(a);
  if (abs_a <= lambda) return 0.0;
  return (a > 0.0) ? (abs_a - lambda) : -(abs_a - lambda);
}

// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Your function code goes here
  int n = Xtilde.n_rows;
  arma::colvec r = Ytilde - Xtilde * beta;
  double rss = arma::dot(r, r);
  double penalty = arma::sum(arma::abs(beta));
  double fval = 0.5 * rss / n + lambda * penalty;
  return fval;
}

// Lasso coordinate-descent on standardized data with one lamdba. Returns a vector beta.
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.001){
  // Your function code goes here
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
  arma::colvec beta = beta_start;
  arma::colvec r = Ytilde - Xtilde * beta;   // residual
  arma::colvec col_sq_norm_over_n = arma::sum(arma::square(Xtilde), 0).t() / n;
  
  double fmin = lasso_c(Xtilde, Ytilde, beta, lambda);
  
  while (true) {
    for (int j = 0; j < p; j++) {
      arma::colvec xj = Xtilde.col(j);
      double bj_old = beta[j];
      
      // zj = (1/n) * x_j^T (r + x_j * b_j)
      double zj = arma::dot(xj, r + xj * bj_old) / n;
  
      double denom = col_sq_norm_over_n[j];
      double bj_new = (denom > 0) ? soft_c(zj, lambda) / denom : 0.0;
    
      if (bj_new != bj_old) {
        double db = bj_new - bj_old;
        r -= xj * db;     // update residual cheaply
        beta[j] = bj_new;
      }
    }
  
    double fnew = lasso_c(Xtilde, Ytilde, beta, lambda);
    if (std::abs(fmin - fnew) < eps) {
      fmin = fnew;
      break;
    } 
    fmin = fnew;
  }
   
  return beta;
}  

// Lasso coordinate-descent on standardized data with supplied lambda_seq. 
// You can assume that the supplied lambda_seq is already sorted from largest to smallest, and has no negative values.
// Returns a matrix beta (p by number of lambdas in the sequence)
// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.001){
  // Your function code goes here
  int n = Xtilde.n_rows;
  int p = Xtilde.n_cols;
  int n_lambda = lambda_seq.n_elem;
  
  arma::mat beta_mat(p, n_lambda, arma::fill::zeros);
  arma::colvec beta_start(p, arma::fill::zeros);
   
  for (int i = 0; i < n_lambda; i++) {
    double lambda = lambda_seq[i];
    arma::colvec beta = fitLASSOstandardized_c(Xtilde, Ytilde, lambda, beta_start, eps);
    beta_mat.col(i) = beta;
    beta_start = beta; // warm start
  }

  return beta_mat;
}