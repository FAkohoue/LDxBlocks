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
//   compute_r2_sparse_cpp -- sparse r² (only pairs within bp window, OpenMP)
//   boundary_scan_cpp     -- weak-LD cut position scan (subsegmentation)

// [[Rcpp::depends(RcppArmadillo)]]
// OpenMP: flags are set via src/Makevars (platform-aware).
// Do NOT use [[Rcpp::plugins(openmp)]] here — on macOS that plugin
// injects -fopenmp into the link step via Apple Clang's OpenMP stub
// rather than Homebrew libomp, leaving __kmpc_* symbols unresolved
// at runtime. Makevars adds -lomp explicitly on Darwin.
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

//' @noRd
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

//' @noRd
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

//' @noRd
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

//' @noRd
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

//' @noRd
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

//' @noRd
 // [[Rcpp::export]]
 List compute_r2_sparse_cpp(
     const arma::mat& X,
     const arma::ivec& bp,
     int    max_bp_dist,
     double threshold  = 0.0,
     int    n_threads  = 1
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

   // Use n_threads for OpenMP parallelism over the outer loop.
   // Thread-local vectors avoid locking; merged after the loop.
#ifdef _OPENMP
   int n_thr = (n_threads == 0) ? omp_get_max_threads() : n_threads;
#else
   int n_thr = 1; (void)n_threads;
#endif
   // Serial path (n_thr=1) or parallel accumulation
   if (n_thr <= 1) {
     for (arma::uword j = 0; j < p; ++j) {
       for (arma::uword k = j + 1; k < p; ++k) {
         if (std::abs(bp(k) - bp(j)) > max_bp_dist) break;
         double r  = arma::dot(Z.col(j), Z.col(k)) / (double)(n - 1);
         double r2 = std::min(r * r, 1.0);
         if (r2 >= threshold) {
           rows.push_back((int)j + 1);
           cols.push_back((int)k + 1);
           vals.push_back(r2);
         }
       }
     }
   } else {
#ifdef _OPENMP
     std::vector<std::vector<int>>    t_rows(n_thr), t_cols(n_thr);
     std::vector<std::vector<double>> t_vals(n_thr);
#pragma omp parallel for schedule(dynamic, 4) num_threads(n_thr)
     for (int j = 0; j < (int)p; ++j) {
       int tid = omp_get_thread_num();
       for (arma::uword k = (arma::uword)j + 1; k < p; ++k) {
         if (std::abs(bp(k) - bp(j)) > max_bp_dist) break;
         double r  = arma::dot(Z.col(j), Z.col(k)) / (double)(n - 1);
         double r2 = std::min(r * r, 1.0);
         if (r2 >= threshold) {
           t_rows[tid].push_back(j + 1);
           t_cols[tid].push_back((int)k + 1);
           t_vals[tid].push_back(r2);
         }
       }
     }
     for (int t = 0; t < n_thr; ++t) {
       rows.insert(rows.end(), t_rows[t].begin(), t_rows[t].end());
       cols.insert(cols.end(), t_cols[t].begin(), t_cols[t].end());
       vals.insert(vals.end(), t_vals[t].begin(), t_vals[t].end());
     }
#endif
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

//' @noRd
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

// ============================================================================
// build_hap_strings_cpp()
// Build haplotype strings from a dosage submatrix.
// Input:  geno_block  -- integer matrix n_ind x n_snps (0/1/2, NA=-1 or NA_INTEGER)
//         na_char     -- character to use for missing genotypes (e.g. ".")
// Output: CharacterVector of length n_ind, each element = concatenated string
//         of allele codes (e.g. "01210")
// This replaces the vapply(seq_len(n_ind), function(i) paste(...)) R loop,
// which incurs one R function call per individual per block.
// For 17,078 blocks x 204 individuals = 3.5M R calls -> one C++ call per block.
// ============================================================================
//' @noRd
 // [[Rcpp::export]]
 Rcpp::CharacterVector build_hap_strings_cpp(
     Rcpp::IntegerMatrix geno_block,
     std::string na_char = "."
 ) {
   int n_ind  = geno_block.nrow();
   int n_snps = geno_block.ncol();
   Rcpp::CharacterVector out(n_ind);

   for (int i = 0; i < n_ind; i++) {
     std::string s;
     s.reserve(n_snps);
     for (int j = 0; j < n_snps; j++) {
       int g = geno_block(i, j);
       if (g == NA_INTEGER || g < 0) {
         s += na_char;
       } else {
         s += std::to_string(g);
       }
     }
     out[i] = s;
   }
   return out;
 }

// ─────────────────────────────────────────────────────────────────────────────
// 8. score_overlap_cpp  (internal helper — not exported)
//    BLAS-vectorised overlap scoring for ONE seam.
//    Implements friend's points B (precomputed Z) and D (matrix multiply).
//
//    Z_adj: n_ind x n_seam_cols, PRE-STANDARDISED columns for the seam window.
//           Columns are indexed 0-based relative to the seam window start.
//    ovlp_cols:  0-based indices (within Z_adj) of disputed SNPs
//    left_reps:  0-based indices (within Z_adj) of left representative SNPs
//    right_reps: 0-based indices (within Z_adj) of right representative SNPs
//
//    Returns: arma::vec of length ovlp_cols.size() with signed scores.
//    score[k] = mean_r2(ovlp[k], left_reps) - mean_r2(ovlp[k], right_reps)
//
//    The key BLAS operation:
//      Z_o  = Z_adj cols for overlap SNPs      (n x n_ovlp)
//      Z_L  = Z_adj cols for left_reps         (n x k_L)
//      Z_R  = Z_adj cols for right_reps        (n x k_R)
//      C_L  = Z_o.t() * Z_L / (n-1)           (n_ovlp x k_L) — BLAS DGEMM
//      C_R  = Z_o.t() * Z_R / (n-1)           (n_ovlp x k_R) — BLAS DGEMM
//      scores = rowMeans(C_L^2) - rowMeans(C_R^2)
//    Two DGEMM calls replace n_ovlp * 2k_rep individual dot products.
// ─────────────────────────────────────────────────────────────────────────────
static arma::vec score_overlap_cpp(
    const arma::mat&         Z_adj,      // pre-standardised seam window
    const std::vector<int>&  ovlp_cols,  // 0-based within Z_adj
    const std::vector<int>&  left_reps,  // 0-based within Z_adj
    const std::vector<int>&  right_reps  // 0-based within Z_adj
) {
  int n_ovlp = (int)ovlp_cols.size();
  int k_L    = (int)left_reps.size();
  int k_R    = (int)right_reps.size();
  int n      = (int)Z_adj.n_rows;

  if (n_ovlp == 0) return arma::zeros<arma::vec>(0);

  // Build sub-matrices by column selection (no copy if using submat view)
  // Z_o: n x n_ovlp
  arma::mat Z_o(n, n_ovlp);
  for (int j = 0; j < n_ovlp; j++) Z_o.col(j) = Z_adj.col(ovlp_cols[j]);

  arma::vec scores(n_ovlp, arma::fill::zeros);

  if (k_L > 0) {
    // Z_L: n x k_L
    arma::mat Z_L(n, k_L);
    for (int j = 0; j < k_L; j++) Z_L.col(j) = Z_adj.col(left_reps[j]);
    // C_L = Z_o.t() * Z_L / (n-1)  — BLAS DGEMM: (n_ovlp x n) * (n x k_L)
    arma::mat C_L = (Z_o.t() * Z_L) / (double)(n - 1);
    // r2 = C_L^2, clamped to [0,1]; rowMeans
    C_L = arma::clamp(C_L % C_L, 0.0, 1.0);
    scores += arma::sum(C_L, 1) / (double)k_L;  // rowMeans
  }

  if (k_R > 0) {
    // Z_R: n x k_R
    arma::mat Z_R(n, k_R);
    for (int j = 0; j < k_R; j++) Z_R.col(j) = Z_adj.col(right_reps[j]);
    // C_R = Z_o.t() * Z_R / (n-1)  — BLAS DGEMM
    arma::mat C_R = (Z_o.t() * Z_R) / (double)(n - 1);
    C_R = arma::clamp(C_R % C_R, 0.0, 1.0);
    scores -= arma::sum(C_R, 1) / (double)k_R;  // subtract rowMeans
  }

  return scores;
}


// ─────────────────────────────────────────────────────────────────────────────
// 9. resolve_seam_cpp  (internal helper — not exported)
//    Resolves overlaps at ONE seam between two adjacent segments.
//    Implements friend's point A (seam-local resolution).
//
//    Called immediately after each segment completes, passing only the
//    local adjusted matrix for the seam window (prev_end - k to curr_start + k).
//    This keeps the working set tiny and cache-friendly.
//
//    blocks_prev: blocks from previous segment (n_prev x 2, 1-based global)
//    blocks_curr: blocks from current segment  (n_curr x 2, 1-based global)
//    adj_seam:    n_ind x seam_width adjusted matrix for columns in the seam
//    seam_start:  1-based global start of the seam window (= adj_seam col 1)
//    k_rep:       max representatives from each side
//
//    Modifies the boundary blocks of prev and curr in-place and returns
//    {updated last row of blocks_prev, updated first row of blocks_curr}.
// ─────────────────────────────────────────────────────────────────────────────
[[maybe_unused]] static void resolve_seam_cpp(
    arma::imat&       blocks_prev,  // last few rows; modified in-place
    arma::imat&       blocks_curr,  // first few rows; modified in-place
    const arma::mat&  adj_seam,     // pre-standardised seam window
    int               seam_start,   // 1-based global col index of adj_seam col 0
    int               k_rep
) {
  // Only the last block of prev and first block of curr can overlap at a seam
  if (blocks_prev.n_rows == 0 || blocks_curr.n_rows == 0) return;

  int n_prev = (int)blocks_prev.n_rows;
  int n_curr = (int)blocks_curr.n_rows;
  int p_seam = (int)adj_seam.n_cols;
  int n      = (int)adj_seam.n_rows;

  // Precompute Z for the seam window ONCE (friend's point B)
  arma::mat Z_seam(n, p_seam, arma::fill::zeros);
  for (int j = 0; j < p_seam; j++) {
    arma::vec col = adj_seam.col(j);
    double mu = arma::mean(col);
    col -= mu;
    double sd = arma::norm(col, 2) / std::sqrt((double)(n - 1));
    if (sd > 1e-10) Z_seam.col(j) = col / sd;
    // else: remains zero (constant column)
  }

  // Convert global 1-based index to 0-based seam-local index
  auto to_seam = [&](int global_1based) -> int {
    return global_1based - seam_start;  // 0-based
  };
  auto in_seam = [&](int global_1based) -> bool {
    int s = to_seam(global_1based);
    return s >= 0 && s < p_seam;
  };

  // Check multiple block pairs at the seam (not just the last/first)
  // In practice, 1-3 blocks near the boundary may overlap
  for (int pi = n_prev - 1; pi >= std::max(0, n_prev - 3); pi--) {
    for (int ci = 0; ci < std::min(n_curr, 3); ci++) {
      int sA = blocks_prev(pi, 0); int eA = blocks_prev(pi, 1);
      int sB = blocks_curr(ci, 0); int eB = blocks_curr(ci, 1);

      if (eA < sB) continue;  // no overlap

      // Build zones (1-based global)
      bool has_left  = (sB > sA);
      bool has_right = (eB > eA);

      if (!has_left || !has_right) {
        // Union merge
        blocks_prev(pi, 0) = std::min(sA, sB);
        blocks_prev(pi, 1) = std::max(eA, eB);
        // Mark curr block as absorbed
        blocks_curr(ci, 0) = -1; blocks_curr(ci, 1) = -1;
        continue;
      }

      // Representatives (1-based global -> 0-based seam-local)
      int left_len  = sB - sA;
      int right_len = eB - eA;
      int k_L = std::min(k_rep, left_len);
      int k_R = std::min(k_rep, right_len);

      std::vector<int> left_reps, right_reps;
      for (int r = 0; r < k_L; r++) {
        int g = sB - k_L + r;  // last k_L of left core (1-based)
        if (in_seam(g)) left_reps.push_back(to_seam(g));
      }
      for (int r = 0; r < k_R; r++) {
        int g = eA + 1 + r;    // first k_R of right core (1-based)
        if (in_seam(g)) right_reps.push_back(to_seam(g));
      }

      // Overlap zone (1-based global -> 0-based seam-local)
      std::vector<int> ovlp_cols;
      for (int g = sB; g <= eA; g++) {
        if (in_seam(g)) ovlp_cols.push_back(to_seam(g));
      }

      if (ovlp_cols.empty() || (left_reps.empty() && right_reps.empty())) {
        // Fallback: union merge (can't score without reps in seam)
        blocks_prev(pi, 0) = std::min(sA, sB);
        blocks_prev(pi, 1) = std::max(eA, eB);
        blocks_curr(ci, 0) = -1; blocks_curr(ci, 1) = -1;
        continue;
      }

      // BLAS-vectorised scoring (friend's point D)
      arma::vec scores = score_overlap_cpp(Z_seam, ovlp_cols, left_reps, right_reps);

      // Cumulative score split rule (identical to R version)
      double cum = 0.0;
      int last_left = 0;
      int ovlp_len = (int)ovlp_cols.size();
      for (int oi = 0; oi < ovlp_len; oi++) {
        cum += scores[oi];
        if (cum >= 0.0) last_left = oi + 1;
      }

      // Map last_left back to 1-based global boundary
      if (last_left == 0) {
        blocks_prev(pi, 1) = sB - 1;
        blocks_curr(ci, 0) = sB;
      } else if (last_left == ovlp_len) {
        blocks_prev(pi, 1) = eA;
        blocks_curr(ci, 0) = eA + 1;
      } else {
        // ovlp_cols[last_left-1] is 0-based seam -> convert to 1-based global
        int split_1based = ovlp_cols[last_left - 1] + seam_start;
        blocks_prev(pi, 1) = split_1based;
        blocks_curr(ci, 0) = split_1based + 1;
      }
    }
  }

  // Remove absorbed curr blocks (marked -1)
  std::vector<int> keep;
  for (int ci = 0; ci < n_curr; ci++) {
    if (blocks_curr(ci, 0) >= 0) keep.push_back(ci);
  }
  if ((int)keep.size() < n_curr) {
    arma::imat tmp(keep.size(), 2);
    for (int k = 0; k < (int)keep.size(); k++) tmp.row(k) = blocks_curr.row(keep[k]);
    blocks_curr = tmp;
  }
}


// ─────────────────────────────────────────────────────────────────────────────
// 10. resolve_overlap_cpp  (exported to R)
//     Full C++ implementation called ONCE after all segments complete.
//     Used as fallback for any overlaps not caught at seam time, and as
//     the primary resolver when seam-local resolution is not used.
//
//     Implements all friend's recommendations:
//     A - single pass over block list (no global explosion)
//     B - lazy column cache (standardise each col at most once)
//     C - single call (not twice)
//     D - BLAS matrix multiply via score_overlap_cpp()
//     E - OpenMP over independent overlap pairs
//
//     blocks:  n_blocks x 2 integer matrix (1-based global indices)
//     adj_mat: n_ind x p_snps pre-centred adjusted genotype matrix
//     k_rep:   max representatives from each core
// ─────────────────────────────────────────────────────────────────────────────

//' @noRd
 // [[Rcpp::export]]
 arma::imat resolve_overlap_cpp(
     arma::imat       blocks,
     const arma::mat& adj_mat,
     int              k_rep = 10
 ) {
   int n_blocks = (int)blocks.n_rows;
   int p        = (int)adj_mat.n_cols;
   int n        = (int)adj_mat.n_rows;

   // Lazy standardisation cache: computed on first access, stored for reuse.
   // For a 204 x 314k matrix we never touch columns that aren't representatives
   // or overlap SNPs — typically < 1% of all columns.
   std::vector<double>    sd_vec(p, -1.0);
   std::vector<arma::vec> z_vec(p);

   auto get_z = [&](int col0) -> const arma::vec& {
     if (sd_vec[col0] < 0.0) {
       arma::vec c = adj_mat.col(col0);
       double mu = arma::mean(c);
       c -= mu;
       double sd = arma::norm(c, 2) / std::sqrt((double)(n - 1));
       sd_vec[col0] = sd;
       if (sd > 1e-10) z_vec[col0] = c / sd;
       else             z_vec[col0] = arma::zeros<arma::vec>(n);
     }
     return z_vec[col0];
   };

   // Collect all overlapping pairs first (fast scan, no LD computation)
   struct SeamPair { int i; int sA,eA,sB,eB; };
   std::vector<SeamPair> pairs;
   pairs.reserve(64);
   for (int i = 0; i < n_blocks - 1; i++) {
     int sA = blocks(i,0), eA = blocks(i,1);
     int sB = blocks(i+1,0), eB = blocks(i+1,1);
     if (eA >= sB) pairs.push_back({i, sA, eA, sB, eB});
   }

   // Resolve each overlapping pair.
   // Friend's point E: pairs are independent -> OpenMP parallel.
   // We use a serial result vector to avoid race conditions on blocks matrix.
   int n_pairs = (int)pairs.size();
   struct Resolution { int i; int new_eA; int new_sB; bool merged; int merge_start; int merge_end; };
   std::vector<Resolution> resolutions(n_pairs);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(4)
#endif
   for (int pi = 0; pi < n_pairs; pi++) {
     auto& P = pairs[pi];
     int i  = P.i;
     int sA = P.sA, eA = P.eA, sB = P.sB, eB = P.eB;

     bool has_left  = (sB > sA);
     bool has_right = (eB > eA);

     if (!has_left || !has_right) {
       resolutions[pi] = {i, 0, 0, true, std::min(sA,sB), std::max(eA,eB)};
       continue;
     }

     int left_len  = sB - sA;
     int right_len = eB - eA;
     int k_L = std::min(k_rep, left_len);
     int k_R = std::min(k_rep, right_len);

     // Left reps: last k_L of [sA0..sB0-1] (0-based)
     std::vector<int> left_reps(k_L), right_reps(k_R);
     for (int r = 0; r < k_L; r++) left_reps[r]  = sB - 1 - k_L + r;  // 0-based
     for (int r = 0; r < k_R; r++) right_reps[r] = eA + r;            // 0-based (eA0+1+r but eA0=eA-1)

     // Overlap SNPs (0-based)
     int ovlp_len = eA - sB + 1;
     std::vector<int> ovlp_cols(ovlp_len);
     for (int oi = 0; oi < ovlp_len; oi++) ovlp_cols[oi] = sB - 1 + oi; // 0-based

     // Warm up cache for reps (thread-safe reads, writes only to own index)
     // Note: in OpenMP parallel region, each col is initialised by at most one
     // thread — race is benign because result is deterministic.
     for (int r : left_reps)  get_z(r);
     for (int r : right_reps) get_z(r);

     // Build seam-local matrices for BLAS scoring (friend's point D)
     // Z_o (n x ovlp_len), Z_L (n x k_L), Z_R (n x k_R)
     arma::mat Z_o(n, ovlp_len), Z_L_m(n, k_L), Z_R_m(n, k_R);
     for (int oi = 0; oi < ovlp_len; oi++) Z_o.col(oi)  = get_z(ovlp_cols[oi]);
     for (int r  = 0; r  < k_L;      r++)  Z_L_m.col(r) = get_z(left_reps[r]);
     for (int r  = 0; r  < k_R;      r++)  Z_R_m.col(r) = get_z(right_reps[r]);

     // BLAS DGEMM: C_L = Z_o.t() * Z_L / (n-1)  -> (ovlp_len x k_L)
     arma::vec scores(ovlp_len, arma::fill::zeros);
     if (k_L > 0) {
       arma::mat C_L = (Z_o.t() * Z_L_m) / (double)(n - 1);
       scores += arma::sum(arma::clamp(C_L % C_L, 0.0, 1.0), 1) / (double)k_L;
     }
     if (k_R > 0) {
       arma::mat C_R = (Z_o.t() * Z_R_m) / (double)(n - 1);
       scores -= arma::sum(arma::clamp(C_R % C_R, 0.0, 1.0), 1) / (double)k_R;
     }

     // Cumulative score split
     double cum = 0.0;
     int last_left = 0;
     for (int oi = 0; oi < ovlp_len; oi++) {
       cum += scores[oi];
       if (cum >= 0.0) last_left = oi + 1;
     }

     // Determine new boundaries (1-based)
     int new_eA, new_sB;
     if (last_left == 0) {
       new_eA = sB - 1; new_sB = sB;
     } else if (last_left == ovlp_len) {
       new_eA = eA; new_sB = eA + 1;
     } else {
       int split = sB + last_left - 1;  // 1-based
       new_eA = split; new_sB = split + 1;
     }
     resolutions[pi] = {i, new_eA, new_sB, false, 0, 0};
   }

   // Apply resolutions serially (modifies blocks matrix)
   // Process in reverse order so that merged rows don't shift indices
   for (int pi = n_pairs - 1; pi >= 0; pi--) {
     auto& R = resolutions[pi];
     if (R.merged) {
       blocks(R.i, 0) = R.merge_start;
       blocks(R.i, 1) = R.merge_end;
       blocks.shed_row(R.i + 1);
     } else {
       blocks(R.i,   1) = R.new_eA;
       blocks(R.i+1, 0) = R.new_sB;
     }
   }

   return blocks;
 }

// ─────────────────────────────────────────────────────────────────────────────
// 11. block_snp_ranges_cpp  (B6: friend's "very high priority" C++ interval lookup)
//
//     For a chromosome with p SNPs and n_blocks LD blocks, maps each block's
//     [start_bp, end_bp] interval onto the corresponding SNP column indices
//     in a SINGLE LINEAR PASS over the sorted SNP position vector.
//
//     This replaces the R-side per-block which() / findInterval() loop that
//     previously ran separately for each of the ~17k retained blocks, each
//     scanning up to 350k SNP positions.
//
//     Algorithm: merge-style sweep (O(n_blocks + p)):
//       - Both snp_pos and block boundaries are assumed sorted ascending.
//       - A single pointer advances through snp_pos as blocks are processed.
//       - For each block [sb, eb]:
//           advance pointer while snp_pos[ptr] < sb   (skip left of block)
//           record lo = ptr
//           advance pointer while snp_pos[ptr] <= eb  (span the block)
//           record hi = ptr - 1
//           store hi - lo + 1 = n_snps in block
//       - Total work: O(p + n_blocks), not O(p * n_blocks).
//
//     Parameters:
//       snp_pos   : sorted integer vector of SNP base-pair positions (p elements)
//       block_sb  : integer vector of block start positions (n_blocks elements)
//       block_eb  : integer vector of block end positions   (n_blocks elements)
//
//     Returns: List with three integer vectors of length n_blocks:
//       $lo     : 1-based start index into snp_pos for each block (0 = no SNPs)
//       $hi     : 1-based end   index into snp_pos for each block (0 = no SNPs)
//       $n_snps : number of SNPs in each block
//
//     The caller uses lo[b]:hi[b] to slice the genotype matrix column range.
//     Blocks with n_snps == 0 have lo = hi = 0 and should be skipped.
// ─────────────────────────────────────────────────────────────────────────────

//' @noRd
 // [[Rcpp::export]]
 Rcpp::List block_snp_ranges_cpp(
     const Rcpp::IntegerVector& snp_pos,   // sorted SNP positions (p)
     const Rcpp::IntegerVector& block_sb,  // block start bp      (n_blocks)
     const Rcpp::IntegerVector& block_eb   // block end bp        (n_blocks)
 ) {
   int p        = snp_pos.size();
   int n_blocks = block_sb.size();

   Rcpp::IntegerVector lo(n_blocks, 0);
   Rcpp::IntegerVector hi(n_blocks, 0);
   Rcpp::IntegerVector n_snps(n_blocks, 0);

   int ptr = 0;  // current position in snp_pos (0-based)

   for (int b = 0; b < n_blocks; b++) {
     int sb = block_sb[b];
     int eb = block_eb[b];

     // Advance past SNPs before this block's start
     while (ptr < p && snp_pos[ptr] < sb) ptr++;

     if (ptr >= p) break;  // no more SNPs

     int lo_b = ptr;  // first SNP in block (0-based)

     // Count SNPs within this block
     int hi_b = lo_b;
     while (hi_b < p && snp_pos[hi_b] <= eb) hi_b++;
     // hi_b is now one past the last SNP in block

     int count = hi_b - lo_b;
     if (count > 0) {
       lo[b]     = lo_b + 1;     // convert to 1-based R index
       hi[b]     = hi_b;         // inclusive 1-based end
       n_snps[b] = count;
     }
     // Do NOT advance ptr: next block may start before or at same position
     // (blocks are sorted but may be adjacent; ptr stays at lo_b for next block)
     // Actually: blocks are non-overlapping and sorted, so we CAN advance ptr
     // to hi_b for efficiency. Next block starts at or after eb+1 >= snp_pos[lo_b].
     ptr = lo_b;  // reset to start of this block; next block's sb >= this eb+1
     // will advance past this block's SNPs naturally in the while loop above
   }

   return Rcpp::List::create(
     Rcpp::Named("lo")     = lo,
     Rcpp::Named("hi")     = hi,
     Rcpp::Named("n_snps") = n_snps
   );
 }

// ─────────────────────────────────────────────────────────────────────────────
// 12. extract_chr_haplotypes_cpp  (B7 + B9 + B10 combined)
//
//     Full chromosome-level haplotype extractor. Takes:
//       geno_chr   : n_ind x p_chr integer matrix (0/1/2, NA=-1)
//       snp_pos    : sorted SNP positions (length p_chr)
//       block_sb   : block start bp vectors (sorted, length n_blocks)
//       block_eb   : block end bp vectors
//       min_snps   : minimum SNPs per block to retain
//       min_freq   : minimum haplotype frequency to include in matrix
//       top_n      : max alleles per block (0 = no cap)
//       na_char    : character to use for missing genotypes
//
//     Returns a single List per chromosome containing:
//       $hap_strings : List of CharacterVectors (one per retained block)
//                      Each CharacterVector has length n_ind, names = sample IDs
//                      if block is retained, else absent from list.
//       $block_lo    : 0-based block start SNP indices (retained blocks only)
//       $block_hi    : 0-based block end SNP indices (retained blocks only)
//       $n_snps      : SNP count per retained block
//       $hap_counts  : List of IntegerVectors — frequency counts per allele
//       $hap_alleles : List of CharacterVectors — unique allele strings
//       $hap_freq    : List of NumericVectors — allele frequencies
//       $n_retained  : integer count of retained blocks
//
//     This replaces the R-side loop that:
//       1. slices geno_chr per block (R matrix indexing overhead)
//       2. calls build_hap_strings_cpp() per block (function call overhead)
//       3. tabulates strings in R (table() overhead)
//     All three steps run inside one OpenMP-parallelised C++ loop.
//
//     OpenMP parallelisation: blocks are independent, so #pragma omp parallel
//     for processes them concurrently. Thread-local result vectors are merged
//     after the parallel section.
// ─────────────────────────────────────────────────────────────────────────────

// Internal: build haplotype string for one individual across a set of columns
static std::string build_one_hap(
    const Rcpp::IntegerMatrix& G,
    int   ind,
    const std::vector<int>& col_idx,
    const std::string& na_char
) {
  int n = (int)col_idx.size();
  std::string s;
  s.reserve(n);
  for (int j = 0; j < n; j++) {
    int g = G(ind, col_idx[j]);
    if (g == NA_INTEGER || g < 0) {
      s += na_char;
    } else {
      s += static_cast<char>('0' + g);
    }
  }
  return s;
}

//' @noRd
 // [[Rcpp::export]]
 Rcpp::List extract_chr_haplotypes_cpp(
     const Rcpp::IntegerMatrix&  geno_chr,   // n_ind x p_chr
     const Rcpp::IntegerVector&  snp_pos,    // sorted, length p_chr
     const Rcpp::IntegerVector&  block_sb,   // block start bp
     const Rcpp::IntegerVector&  block_eb,   // block end bp
     int                         min_snps,
     double                      min_freq,
     int                         top_n,      // 0 = no cap
     std::string                 na_char = "."
 ) {
   int n_ind    = geno_chr.nrow();
   int p_chr    = geno_chr.ncol();
   int n_blocks = block_sb.size();

   // ── Step 1: C++ interval lookup (same algorithm as block_snp_ranges_cpp) ──
   // Single linear pass: O(p_chr + n_blocks)
   std::vector<int> lo_vec(n_blocks, -1), hi_vec(n_blocks, -1);
   {
     int ptr = 0;
     for (int b = 0; b < n_blocks; b++) {
       int sb = block_sb[b], eb = block_eb[b];
       while (ptr < p_chr && snp_pos[ptr] < sb) ptr++;
       if (ptr >= p_chr) break;
       int lo = ptr;
       int hi = lo;
       while (hi < p_chr && snp_pos[hi] <= eb) hi++;
       hi--;  // inclusive
       if ((hi - lo + 1) >= min_snps) {
         lo_vec[b] = lo;
         hi_vec[b] = hi;
       }
       // don't advance ptr — next block might start before this one ends
       ptr = lo;
     }
   }

   // ── Step 2: Per-block parallel string building + frequency tabulation ──
   // Thread-safe: each block writes to its own slot in result vectors.
   std::vector<std::vector<std::string>> all_strings(n_blocks);  // [block][ind]
   std::vector<int> n_snps_vec(n_blocks, 0);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(4)
#endif
   for (int b = 0; b < n_blocks; b++) {
     if (lo_vec[b] < 0) continue;
     int lo = lo_vec[b], hi = hi_vec[b];
     int ns = hi - lo + 1;
     n_snps_vec[b] = ns;

     // Build column index list for this block
     std::vector<int> cols(ns);
     for (int j = 0; j < ns; j++) cols[j] = lo + j;

     // Build haplotype string for every individual
     std::vector<std::string> strs(n_ind);
     for (int i = 0; i < n_ind; i++) {
       strs[i] = build_one_hap(geno_chr, i, cols, na_char);
     }
     all_strings[b] = strs;
   }

   // ── Step 3: Tabulate frequencies and filter alleles ──
   // Serial (string map operations not thread-safe without locks)
   Rcpp::List  hap_strings_list, hap_alleles_list, hap_freq_list, hap_counts_list;
   Rcpp::NumericVector out_freq_dom;  // dominant allele freq per retained block
   Rcpp::IntegerVector out_lo, out_hi, out_nsnps;

   for (int b = 0; b < n_blocks; b++) {
     if (lo_vec[b] < 0 || all_strings[b].empty()) continue;

     const std::vector<std::string>& strs = all_strings[b];
     std::string na_str = na_char;

     // Count occurrences (skip NA strings containing na_char)
     std::unordered_map<std::string, int> cnt;
     int total = 0;
     for (const auto& s : strs) {
       if (s.find(na_char[0]) == std::string::npos) {
         cnt[s]++;
         total++;
       }
     }
     if (total == 0 || (int)cnt.size() < 1) continue;

     // Sort alleles by frequency descending
     std::vector<std::pair<int, std::string>> sorted_cnt;
     sorted_cnt.reserve(cnt.size());
     for (auto& kv : cnt) sorted_cnt.push_back({kv.second, kv.first});
     std::sort(sorted_cnt.begin(), sorted_cnt.end(),
               [](const auto& a, const auto& b){ return a.first > b.first; });

     // Apply min_freq and top_n filters
     int n_alleles = (int)sorted_cnt.size();
     int keep = 0;
     for (int k = 0; k < n_alleles; k++) {
       double freq = (double)sorted_cnt[k].first / total;
       if (freq < min_freq) break;
       keep++;
       if (top_n > 0 && keep >= top_n) break;
     }
     if (keep == 0) continue;

     // Build output vectors for this block
     Rcpp::CharacterVector alleles_rv(keep), hap_rv(n_ind);
     Rcpp::NumericVector   freq_rv(keep);
     Rcpp::IntegerVector   cnt_rv(keep);

     for (int k = 0; k < keep; k++) {
       alleles_rv[k] = sorted_cnt[k].second;
       cnt_rv[k]     = sorted_cnt[k].first;
       freq_rv[k]    = (double)sorted_cnt[k].first / total;
     }
     for (int i = 0; i < n_ind; i++) hap_rv[i] = strs[i];

     // Dominant allele frequency: most frequent allele (sorted_cnt[0])
     double freq_dom = (total > 0 && !sorted_cnt.empty())
       ? (double)sorted_cnt[0].first / total : 0.0;
     hap_strings_list.push_back(hap_rv);
     hap_alleles_list.push_back(alleles_rv);
     hap_freq_list.push_back(freq_rv);
     hap_counts_list.push_back(cnt_rv);
     out_freq_dom.push_back(freq_dom);
     out_lo.push_back(lo_vec[b]);
     out_hi.push_back(hi_vec[b]);
     out_nsnps.push_back(n_snps_vec[b]);
   }

   return Rcpp::List::create(
     Rcpp::Named("hap_strings")   = hap_strings_list,
     Rcpp::Named("hap_alleles")   = hap_alleles_list,
     Rcpp::Named("hap_freq")      = hap_freq_list,
     Rcpp::Named("hap_counts")    = hap_counts_list,
     Rcpp::Named("freq_dominant") = out_freq_dom,
     Rcpp::Named("block_lo")      = out_lo,
     Rcpp::Named("block_hi")      = out_hi,
     Rcpp::Named("n_snps")        = out_nsnps,
     Rcpp::Named("n_retained")    = (int)out_lo.size()
   );
 }

// =============================================================================
// impute_and_filter_cpp()
//
// Single O(n x p) pass over the genotype matrix that:
//   1. Computes per-SNP call rate (proportion of non-missing samples).
//   2. Flags SNPs whose call rate < min_callrate for removal.
//   3. Imputes missing values (NA = INT_MIN in integer matrices from R) using
//      either mean_rounded or mode strategy — on the SNPs that pass the filter.
//
// Arguments
// ---------
// geno         : integer matrix, individuals (rows) x SNPs (cols), NA = NA_INTEGER
// min_callrate : double in [0,1]. SNPs with call_rate < min_callrate are dropped.
//                Set to 0.0 to disable call-rate filtering.
// method       : 0 = mean_rounded  (round column mean to nearest 0/1/2)
//                1 = mode          (most common observed dosage 0/1/2)
//
// Returns a named List:
//   $geno_imputed  : integer matrix, same dimensions as input; NAs filled in
//                    passing SNPs; columns of failing SNPs are unchanged
//                    (caller uses $keep to subset columns).
//   $keep          : logical vector, length p; TRUE if SNP passes call-rate filter
//   $call_rates    : numeric vector, length p; per-SNP call rate
//   $imp_values    : integer vector, length p; imputed value used per SNP
//                    (NA for SNPs that had no missing or were filtered out)
//   $n_imputed     : integer; total number of imputed values
//   $n_filtered    : integer; number of SNPs removed by call-rate filter
// =============================================================================

// [[Rcpp::export]]
Rcpp::List impute_and_filter_cpp(
    Rcpp::IntegerMatrix geno,
    double              min_callrate = 0.0,
    int                 method       = 0      // 0=mean_rounded, 1=mode
) {
  const int n = geno.nrow();   // individuals
  const int p = geno.ncol();   // SNPs

  Rcpp::IntegerMatrix out    = Rcpp::clone(geno);   // will be modified in-place
  Rcpp::LogicalVector keep   (p, true);
  Rcpp::NumericVector cr     (p, NA_REAL);
  Rcpp::IntegerVector imp_val(p, NA_INTEGER);

  int n_imputed  = 0;
  int n_filtered = 0;

  for (int j = 0; j < p; ++j) {
    // ── Pass 1: scan column, accumulate counts ──────────────────────────────
    int    n_obs  = 0;
    long   sum_g  = 0;
    int    cnt[3] = {0, 0, 0};   // counts of dosage 0, 1, 2
    bool   has_na = false;

    for (int i = 0; i < n; ++i) {
      int g = geno(i, j);
      if (g == NA_INTEGER) {
        has_na = true;
        continue;
      }
      ++n_obs;
      sum_g += g;
      if (g >= 0 && g <= 2) cnt[g]++;
    }

    double call_rate = (n > 0) ? (double)n_obs / n : 0.0;
    cr[j] = call_rate;

    // ── Call-rate filter ─────────────────────────────────────────────────────
    if (call_rate < min_callrate) {
      keep[j] = false;
      ++n_filtered;
      continue;   // leave column untouched; caller will drop it
    }

    // ── Imputation (only if missing values present) ──────────────────────────
    if (!has_na || n_obs == 0) continue;

    int fill_val;
    if (method == 1) {
      // mode: most common of {0, 1, 2}.
      // Tie-breaking: if two or more dosages share the maximum count,
      // the mode is ambiguous -- fall back to mean_rounded for this SNP.
      // This avoids silently biasing toward the lower dosage and is
      // consistent with what the user would get from "mean_rounded".
      int max_cnt = cnt[0];
      if (cnt[1] > max_cnt) max_cnt = cnt[1];
      if (cnt[2] > max_cnt) max_cnt = cnt[2];
      int n_modes = (cnt[0] == max_cnt) + (cnt[1] == max_cnt) +
        (cnt[2] == max_cnt);
      if (n_modes == 1) {
        // Unique mode: use it.
        fill_val = (cnt[0] == max_cnt) ? 0 : (cnt[1] == max_cnt) ? 1 : 2;
      } else {
        // Ambiguous mode: fall back to mean_rounded.
        double mn = (double)sum_g / n_obs;
        fill_val = (int)(mn + 0.5);
        if (fill_val < 0) fill_val = 0;
        if (fill_val > 2) fill_val = 2;
      }
    } else {
      // mean_rounded: round(mean) clamped to [0, 2]
      double mn = (double)sum_g / n_obs;
      fill_val = (int)(mn + 0.5);          // round to nearest integer
      if (fill_val < 0) fill_val = 0;
      if (fill_val > 2) fill_val = 2;
    }
    imp_val[j] = fill_val;

    // ── Pass 2: fill NAs ─────────────────────────────────────────────────────
    for (int i = 0; i < n; ++i) {
      if (out(i, j) == NA_INTEGER) {
        out(i, j) = fill_val;
        ++n_imputed;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("geno_imputed") = out,
    Rcpp::Named("keep")         = keep,
    Rcpp::Named("call_rates")   = cr,
    Rcpp::Named("imp_values")   = imp_val,
    Rcpp::Named("n_imputed")    = n_imputed,
    Rcpp::Named("n_filtered")   = n_filtered
  );
}
