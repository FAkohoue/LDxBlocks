// ld_core.cpp
// Core C++ routines for LDxBlocks
// Compiled via Rcpp + RcppArmadillo
//
// Functions exported to R:
//   compute_r2_cpp        -- standard r² LD matrix (no kinship)
//   compute_rV2_cpp       -- kinship-adjusted rV² LD matrix
//   maf_filter_cpp        -- fast MAF + monomorphic filter
//   build_adj_matrix_cpp  -- threshold adjacency from LD matrix
//   col_r2_cpp            -- r² of one column against all others (for boundary scan)
//   compute_r2_sparse_cpp -- sparse r² (only pairs within bp window)
//   boundary_scan_cpp     -- weak-LD cut position scan (subsegmentation)

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// ─────────────────────────────────────────────────────────────────────────────
// Internal: column-wise mean (handles NAs as 0 after imputation)
// ─────────────────────────────────────────────────────────────────────────────
static inline arma::vec col_means(const arma::mat& X) {
  return arma::mean(X, 0).t();
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. compute_r2_cpp
//    Standard r² LD matrix for a genotype window.
//    X : n × p matrix (individuals × SNPs), 0/1/2 coded, NAs allowed.
//    NA strategy: mean-impute each column before computing correlations.
//    Returns symmetric p × p matrix, diagonal = 0, values in [0,1].
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
arma::mat compute_r2_cpp(
    const arma::mat& X,
    int digits     = -1,
    int n_threads  = 1
) {
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;

  // Mean-impute NAs and standardise
  arma::mat Z(n, p);
  for (arma::uword j = 0; j < p; ++j) {
    arma::vec col = X.col(j);
    // replace NaN/NA with column mean
    double mu = 0.0; arma::uword cnt = 0;
    for (arma::uword i = 0; i < n; ++i) {
      if (std::isfinite(col(i))) { mu += col(i); cnt++; }
    }
    mu = (cnt > 0) ? mu / cnt : 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      Z(i, j) = std::isfinite(col(i)) ? col(i) - mu : 0.0;
    }
    double sd = arma::norm(Z.col(j), 2) / std::sqrt((double)(n - 1));
    if (sd > 1e-10) Z.col(j) /= sd;
    else            Z.col(j).zeros();
  }

  // r² = (Z'Z / (n-1))²  — but columns already unit-variance, so
  // Z'Z/(n-1) is already the correlation matrix
  arma::mat R(p, p, arma::fill::zeros);
#ifdef _OPENMP
  int n_thr = (n_threads == 0) ? omp_get_max_threads() : n_threads;
#else
  int n_thr = 1;  // OpenMP unavailable (e.g. Apple Clang)
  (void)n_threads;  // suppress unused-parameter warning
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(n_thr)
#endif
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = j + 1; k < p; ++k) {
      double r = arma::dot(Z.col(j), Z.col(k)) / (double)(n - 1);
      double r2 = r * r;
      // clamp numerical noise
      if (r2 > 1.0) r2 = 1.0;
      R(j, k) = r2;
      R(k, j) = r2;
    }
  }

  if (digits >= 0) {
    double fac = std::pow(10.0, (double)digits);
    R = arma::round(R * fac) / fac;
  }
  return R;
}


// ─────────────────────────────────────────────────────────────────────────────
// 2. compute_rV2_cpp
//    Kinship-adjusted rV² for a pre-whitened matrix X = V^{-1/2} G_centered.
//    The whitening is done in R (get_V_inv_sqrt); this function only
//    computes the squared correlations, matching compute_r2_cpp logic.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
arma::mat compute_rV2_cpp(
    const arma::mat& X,
    int digits    = -1,
    int n_threads = 1
) {
  // For rV², input is already whitened so we just standardise and correlate
  return compute_r2_cpp(X, digits, n_threads);
}


// ─────────────────────────────────────────────────────────────────────────────
// 3. maf_filter_cpp
//    Returns logical vector: TRUE = keep (MAF >= threshold, not monomorphic)
//    Faster than apply() in R for large matrices.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
LogicalVector maf_filter_cpp(const arma::mat& G, double maf_cut = 0.05) {
  arma::uword p = G.n_cols;
  arma::uword n = G.n_rows;
  LogicalVector keep(p);

  for (arma::uword j = 0; j < p; ++j) {
    double sum_g = 0.0; arma::uword cnt = 0;
    bool monomorphic = true;
    double first_val = -1.0;
    for (arma::uword i = 0; i < n; ++i) {
      double v = G(i, j);
      if (!std::isfinite(v)) continue;
      sum_g += v; cnt++;
      if (first_val < 0.0) first_val = v;
      else if (v != first_val) monomorphic = false;
    }
    if (monomorphic || cnt == 0) { keep[j] = false; continue; }
    double freq = sum_g / (2.0 * cnt);
    double maf  = (freq > 0.5) ? 1.0 - freq : freq;
    keep[j] = (maf >= maf_cut);
  }
  return keep;
}


// ─────────────────────────────────────────────────────────────────────────────
// 4. build_adj_matrix_cpp
//    Threshold an LD matrix into a 0/1 adjacency matrix.
//    More efficient than ifelse() in R because it writes in-place.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
arma::imat build_adj_matrix_cpp(const arma::mat& LD, double threshold) {
  arma::uword p = LD.n_rows;
  arma::imat A(p, p, arma::fill::zeros);
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = j + 1; k < p; ++k) {
      if (LD(j, k) >= threshold) {
        A(j, k) = 1; A(k, j) = 1;
      }
    }
  }
  return A;
}


// ─────────────────────────────────────────────────────────────────────────────
// 5. col_r2_cpp
//    r² between one query column and all other columns.
//    Used in cutsequence boundary scan — replaces many compute_r2 calls.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
arma::vec col_r2_cpp(const arma::mat& X, int query_col) {
  arma::uword p  = X.n_cols;
  arma::uword n  = X.n_rows;
  arma::uword qc = (arma::uword)(query_col - 1); // convert to 0-based
  arma::vec result(p, arma::fill::zeros);

  // Standardise query column
  arma::vec q = X.col(qc);
  double mu_q = arma::mean(q);
  q -= mu_q;
  double sd_q = arma::norm(q, 2) / std::sqrt((double)(n - 1));
  if (sd_q < 1e-10) return result; // constant column
  q /= sd_q;

  for (arma::uword k = 0; k < p; ++k) {
    if (k == qc) continue;
    arma::vec col = X.col(k);
    double mu_k = arma::mean(col);
    col -= mu_k;
    double sd_k = arma::norm(col, 2) / std::sqrt((double)(n - 1));
    if (sd_k < 1e-10) { result(k) = 0.0; continue; }
    col /= sd_k;
    double r = arma::dot(q, col) / (double)(n - 1);
    result(k) = std::min(r * r, 1.0);
  }
  return result;
}


// ─────────────────────────────────────────────────────────────────────────────
// 6. compute_r2_sparse_cpp
//    Only computes r² for pairs within a bp distance window.
//    Returns a sparse triplet representation (i, j, value).
//    Critical for large chromosomes — avoids O(p²) full matrix.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
List compute_r2_sparse_cpp(
    const arma::mat& X,
    const arma::ivec& bp,
    int    max_bp_dist,
    double threshold = 0.0
) {
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;

  // Pre-standardise all columns
  arma::mat Z(n, p, arma::fill::zeros);
  for (arma::uword j = 0; j < p; ++j) {
    arma::vec col = X.col(j);
    double mu = arma::mean(col);
    col -= mu;
    double sd = arma::norm(col, 2) / std::sqrt((double)(n - 1));
    if (sd > 1e-10) Z.col(j) = col / sd;
  }

  std::vector<int>    rows, cols;
  std::vector<double> vals;
  rows.reserve(p * 20); cols.reserve(p * 20); vals.reserve(p * 20);

  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = j + 1; k < p; ++k) {
      if (std::abs(bp(k) - bp(j)) > max_bp_dist) break; // sorted bp assumed
      double r  = arma::dot(Z.col(j), Z.col(k)) / (double)(n - 1);
      double r2 = std::min(r * r, 1.0);
      if (r2 >= threshold) {
        rows.push_back((int)j + 1); // 1-based for R
        cols.push_back((int)k + 1);
        vals.push_back(r2);
      }
    }
  }

  return List::create(
    Named("row") = IntegerVector(rows.begin(), rows.end()),
    Named("col") = IntegerVector(cols.begin(), cols.end()),
    Named("r2")  = NumericVector(vals.begin(), vals.end())
  );
}


// ─────────────────────────────────────────────────────────────────────────────
// 7. boundary_scan_cpp
//    The inner loop of cutsequence.modi — find weak-LD cut positions.
//    For each candidate cut point i, checks whether any pair (left j, right k)
//    has r² >= threshold. Returns integer vector of positions with weak LD.
// ─────────────────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
IntegerVector boundary_scan_cpp(
    const arma::mat& X,
    int start,
    int end,
    int half_w,
    double threshold
) {
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  int len = end - start + 1;
  IntegerVector result(len, 1); // default: valid cut

  // Pre-standardise
  arma::mat Z(n, p, arma::fill::zeros);
  for (arma::uword j = 0; j < p; ++j) {
    arma::vec col = X.col(j);
    double mu = arma::mean(col);
    col -= mu;
    double sd = arma::norm(col, 2) / std::sqrt((double)(n - 1));
    if (sd > 1e-10) Z.col(j) = col / sd;
  }

  for (int ci = 0; ci < len; ++ci) {
    int cut = start + ci - 1; // 0-based cut position
    int l_start = std::max(0, cut - half_w + 1);
    int l_end   = cut;
    int r_start = cut + 1;
    int r_end   = std::min((int)p - 1, cut + half_w);

    bool has_ld = false;
    for (int j = l_start; j <= l_end && !has_ld; ++j) {
      for (int k = r_start; k <= r_end && !has_ld; ++k) {
        double r  = arma::dot(Z.col(j), Z.col(k)) / (double)(n - 1);
        double r2 = std::min(r * r, 1.0);
        if (r2 >= threshold) has_ld = true;
      }
    }
    result[ci] = has_ld ? 0 : 1;
  }
  return result;
}
