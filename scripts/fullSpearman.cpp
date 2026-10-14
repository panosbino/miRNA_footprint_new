#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Rank a single vector with AVERAGE-RANK tie handling (the standard Spearman
// convention: three equal values at sorted positions 1,2,3 all get rank 2).
// The original version gave tied values sequential, arbitrary ranks; on real
// scran-normalized data that produced correlations up to 0.063 away from R's
// cor(method = "spearman"). That matters here because the correlations are
// thresholded (>0.4) to build a graph, so a shift of that size can flip
// borderline gene pairs across the cutoff.
arma::vec rank_vec(const arma::vec &x) {
    arma::uvec idx = sort_index(x);
    arma::vec r(x.n_elem);
    const size_t n = x.n_elem;

    size_t i = 0;
    while (i < n) {
        size_t j = i;
        while (j + 1 < n && x(idx[j + 1]) == x(idx[i])) {
            j++;
        }
        const double avg_rank = (static_cast<double>(i) + static_cast<double>(j)) / 2.0 + 1.0;
        for (size_t k = i; k <= j; k++) {
            r(idx[k]) = avg_rank;
        }
        i = j + 1;
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
