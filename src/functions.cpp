// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/* ------------------------------------------------------------------
 * Build H_{ik} = ∑_{t_{k-1} ∈ [s_k , y_i)}  Δt / (1-G(t-))
 *
 * Inputs
 *   Y      : n  (synthetic event times  y_i  or  Y*_i)
 *   knots  : m  censoring jump times   s_k   (exclude 0)
 *   dt     : m  Δt_k  =  knots[k] - knots[k-1]
 *   Sc     : m  S_C(t_k) = 1 - G(t_k)           (left-continuous)
 *
 * We pre-compute a cumulative vector
 *   Ccum[k] = Σ_{j ≤ k}  Δt_j / S_C(t_j)
 * then H_{ik} = Ccum[idx_up(i,k)] - Ccum[k-1].
 * ------------------------------------------------------------------*/

// [[Rcpp::export]]
arma::mat H_mat_cpp(const arma::vec& Y,          // n
                    const arma::vec& knots,      // m
                    const arma::vec& dt,         // m
                    const arma::vec& Sc) {       // m
  int n = Y.n_elem, m = knots.n_elem;

  /* cumulative integral  C_k = Σ_{j≤k} Δt_j / S_C(t_j) */
  arma::vec K = dt / arma::clamp(Sc, 1e-12, arma::datum::inf);
  arma::vec Ccum = arma::cumsum(K);              // length m

  arma::mat H(n, m, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    double yi = Y(i);
    for (int k = 0; k < m; ++k) {
      double upper = yi;            // integrate to y_i
      if (upper <= knots(k)) break; // remainder of row is zero
      // largest index  j  with  knots[j] < upper
      int idx_up = std::upper_bound(knots.begin(), knots.end(), upper)
        - knots.begin() - 1;
      double lowerC = (k == 0) ? 0.0 : Ccum(k - 1);
      H(i, k) = Ccum(idx_up) - lowerC;
    }
  }
  return H;
}



/*  Compute add_i
 *  Inputs:
 *    yhat   : n   (ŷ_i = x_i'β)
 *    tgrid  : m   knot times  t_k   (exclude 0)
 *    cumK   : m   cumulative KΔt up to t_k
 *    Pk_dt  : m   vector  (1-F_n(t_k)) * Δt_k
 *    Kvals  : m+1   KΔt up to t_k
 *  Output: n-vector  add_i
 */
// [[Rcpp::export]]
arma::vec add_term_cpp2(const arma::vec& yhat,
                        const arma::vec& tgrid,
                        const arma::vec& cumK,
                        const arma::vec& Pk_dt,
                        const arma::vec& Kvals,
                        const bool tail_ok)
{
  int n = yhat.n_elem, m = tgrid.n_elem;
  arma::vec out(n, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    double gi = yhat(i);
    double acc = 0.0;

    for (int k = 0; k < m; ++k) {
      double upper = gi + tgrid(k);
      // if u <= first knot, the integral from -inf to u is approximated as 0
      // (because we have no grid mass below tgrid(0)). Skip.
      if (upper <= tgrid(0)) continue;

      // last index j with tgrid[j] <= upper
      int idx_up = std::upper_bound(tgrid.begin(), tgrid.end(), upper)
        - tgrid.begin() - 1;
      if (idx_up < 0) continue;              // safety; shouldn't happen due to check above
      if (idx_up >= m) idx_up = m - 1;       // safety clamp

      // integral up to tgrid[idx_up]
      double inner = cumK(idx_up);

      // add the partial leftover on (tgrid[idx_up], upper]
      double left = tgrid(idx_up);
      double Klast = Kvals(idx_up + 1);          // use the K of the interval containing 'upper'
      //inner += (upper - left) * Klast;
      if(tail_ok){
        inner = (upper - left) * Klast;
      }
      acc += Pk_dt(k) * inner;
    }
    out(i) = acc;
  }
  return out;
}
