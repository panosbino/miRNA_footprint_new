#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Rank a vector, giving tied values their average rank (matches R's rank(ties.method = "average"))
arma::vec rank_vec(const arma::vec &x) {
    const arma::uword n = x.n_elem;
    arma::uvec idx = arma::stable_sort_index(x);
    arma::vec r(n);

    arma::uword i = 0;
    while (i < n) {
        arma::uword j = i;
        while (j + 1 < n && x(idx(j + 1)) == x(idx(i))) j++;   // extend over the run of ties
        double avg_rank = (i + j) / 2.0 + 1.0;                 // average of ranks i+1 .. j+1
        for (arma::uword k = i; k <= j; k++) r(idx(k)) = avg_rank;
        i = j + 1;
    }
    return r;
}

// [[Rcpp::export]]
arma::mat spearman_full_cpp(const arma::mat &X) {
    // X: genes x samples (each row = gene)
    arma::mat R(X.n_rows, X.n_cols);
    for (arma::uword i = 0; i < X.n_rows; i++) {
        R.row(i) = rank_vec(X.row(i).t()).t();
    }
    // arma::cor correlates columns; R.t() is samples x genes,
    // so the result is genes x genes (gene-gene Spearman matrix)
    return arma::cor(R.t());
}