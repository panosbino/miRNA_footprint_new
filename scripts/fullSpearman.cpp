#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Rank a single vector (no ties handled specially)
arma::vec rank_vec(const arma::vec &x) {
    arma::uvec idx = sort_index(x);
    arma::vec r(x.n_elem);

    for (size_t i = 0; i < x.n_elem; i++) {
        r(idx[i]) = i + 1;
    }
    return r;
}

// [[Rcpp::export]]
arma::mat spearman_full_cpp(const arma::mat &X) {
    // X: genes x samples (each row = gene)
    size_t n_genes = X.n_rows;

    // Allocate rank matrix
    arma::mat R(n_genes, X.n_cols);

    // Rank each gene (row)
    for (size_t i = 0; i < n_genes; i++) {
        R.row(i) = rank_vec(X.row(i).t()).t();
    }

    // Compute correlation on rank matrix
    arma::mat C = arma::cor(R.t());   // returns samples x samples, so transpose input!

    return C;
}
