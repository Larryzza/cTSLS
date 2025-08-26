#' Sandwich Variance Estimator for Two-Stage IV AFT
#'
#' Computes the asymptotic variance components \eqn{A} and \eqn{B} for the two-stage
#' instrumental variable accelerated failure time (IV AFT) model.
#'
#' Under Theorem 1, the estimator satisfies the asymptotic normality:
#' \deqn{\sqrt{n}(\hat{\theta} - \theta) \to N(0, A^{-1} B A^{-T}).}
#' This function computes \eqn{\hat{A}} and the naive estimate \eqn{\hat{B}} (excluding
#' the contribution from \eqn{\Psi^*}), which is currently left as a placeholder.
#'
#' @param fs A list returned by \code{first_stage()}, containing:
#'   \code{beta} (i.e., \eqn{\alpha}), \code{fitted}, and \code{resid1}.
#' @param ss A list returned by \code{second_stage_wls()}, containing:
#'   \code{theta} (i.e., \eqn{\beta}), \code{fitted}, \code{resid2}, and \code{omega2}.
#' @param Cov A numeric vector of covariates used in stage 1.
#' @param Cov2 A numeric vector of covariates (excluding the endogenous variable) used in stage 2.
#' @param Naive Logical; whether to use the naive sandwich estimator (ignoring \eqn{\Psi^*}).
#'
#' @return A list with components:
#' \describe{
#'   \item{A}{A \eqn{4 \times 4} matrix estimating \eqn{\hat{A}} (block derivatives).}
#'   \item{B_naive}{A \eqn{4 \times 4} matrix estimating naive \eqn{\hat{B} = \frac{1}{n} \sum_j \Psi_j \Psi_j^T}.}
#'   \item{B}{A \eqn{4 \times 4} matrix estimating full \eqn{\hat{B}} (currently equals \code{B_naive}).}
#' }
#' @export
compute_AB <- function(fs, ss, Cov, Cov2, Naive = FALSE) {

  n <- length(fs$fitted)
  # unpack parameters
  alpha <- fs$beta           # c(alpha0, alpha1)
  beta  <- ss$theta          # c(beta0, beta1)
  xhat  <- fs$fitted         # n×1
  resid1 <- fs$resid
  resid2 <- ss$resid         # Ystar - fitted second stage
  omega2 <- ss$w        # variance weights input

  #--- Estimate A_hat ----------------------------------------
  # dPsi1/dalpha = -n^{-1} sum [1;Z] %*% t([1;Z])
  M1 <- cbind(1, Cov)          # n×2
  d11 <- - crossprod(M1, M1) / n
  # dPsi2/dbeta  = -n^{-1} sum omega2 * [1; alpha0+alpha1 Z] %*% t([1; alpha0+alpha1 Z])
  M22 <- cbind(1, xhat, Cov2)        # n×2
  d22 <- - crossprod(M22 * sqrt(omega2), M22 * sqrt(omega2)) / n
  # dPsi2/dalpha = n^{-1} sum omega2 * [ -beta1; resid2 - beta1*xhat ] %*% t([1;Z])
  m211 <- t(omega2 * resid2) %*% M1
  m212 <- - beta[2] * t(omega2 * M22) %*% M1
  m212[2,] <- m212[2,] + m211
  d21 <- m212 / n
  #temp_matrix <- beta[2]*cbind(- rep(1, n), resid2 - xhat, - Cov2)
  #d21 <- - t(sqrt(omega2) * temp_matrix) %*% (sqrt(omega2) * M1) / n
  # d12 <- m122 / n                # 2×2
  # dPsi1/dbeta = 0
  d12 <- matrix(0, dim(d11)[1], dim(d22)[1])
  # assemble A_hat
  A_hat <- rbind(cbind(d11, d12), cbind(d21, d22))

  #--- Estimate B_naive -------------------------------------
  # psi1_j = resid1_j * [1; Z_j]
  Psi1 <- resid1 * M1
  # psi2_j = omega2_j * resid2_j * [1; xhat_j]
  Psi2 <- omega2 * resid2 * M22
  # full Psi = [Psi1 | Psi2]
  Psi   <- cbind(Psi1, Psi2)     # n×4
  B_naive <- crossprod(Psi) / n

  if(Naive == FALSE){
    ## ----  Ψ* contribution -------------------------------------------------
    ## --- 1. censoring KM
    G_surv <- ss$syn$fitG            # survfit object for S_C(t)
    jump_t <- G_surv$time
    Sjump  <- G_surv$surv[-length(G_surv$surv)]
    dt_g   <- diff(G_surv$time)           # interval widths Δt_k

    ## --- 2. observed (or synthetic) times
    Y_obs  <- ss$syn$time            # ≡ \tilde y_i  or Ystar

    ## --- 3. at-risk counts at each jump time s > 0
    at_risk <- vapply(jump_t[-1], function(s) sum(Y_obs >= s), numeric(1))
    at_risk[at_risk == 0] <- 1

    ## --- 4. indicator matrices
    risk_ind  <- outer(Y_obs, jump_t[-1], `>=`)

    ## Censoring increments with tolerance (only for censored subjects)
    m <- length(jump_t[-1])
    DeltaN <- matrix(0L, n, m)
    rows   <- which(status == 0L)            # rows that were censored
    idx    <- nearest_knot_idx(Y_obs[rows], jump_t[-1])  # which knot each censored time maps to
    keep   <- idx > 0L
    DeltaN[cbind(rows[keep], idx[keep])] <- 1L

    ## Column totals (pooled censoring jumps)
    DeltaN_tot <- colSums(DeltaN)

    ## --- 5. H_{ij}(s) matrix
    m <- length(jump_t) - 1          # number of time points  sₖ
    H_mat <- t(vapply(Y_obs, function(y) {
      vapply(jump_t[-1], function(s) {
        if (s < y) {
          idx <- which(jump_t[-1] >= s & jump_t[-1] <= y)
          if (length(idx) == 0) 0 else
            sum(dt_g[idx] / pmax(Sjump[idx], .Machine$double.eps))
        } else 0
      }, numeric(1))
    }, numeric(m)))
    # H_mat is n × m  (rows = subjects, cols = jump times)


    p_beta <- ncol(M22)            # number of β–components

    ## 1) make an n × p_beta matrix of gradients
    gprime_mat <- M22              # each row = (1, x̂_i, W_i)

    ## 2) ω_i g′(θ) H_{ij}(s)
    #   element-wise multiply each column of H_mat by ω_i, then sweep by gprime
    W_H_arr <- array(0, dim = c(n, p_beta, length(jump_t) - 1))
    for (k in seq_len(p_beta)) {
      W_H_arr[, k, ] <- (omega2 * gprime_mat[, k]) * H_mat
    }

    ## 3) martingale increments  (n × m)
    Mart <- (DeltaN - risk_ind * rep(DeltaN_tot / at_risk, each = n)) /
      rep(at_risk, each = n)

    ## 4) Ψ2*_j matrix, n × p_beta
    Psi2_star <- matrix(0, n, p_beta)
    for (k in seq_len(p_beta)) {
      Psi2_star[, k] <- rowSums(W_H_arr[, k, ] * Mart)
    }

    ## 5) plug into B_full
    B_full <- crossprod(cbind(Psi1, Psi2 + Psi2_star)) / n
  }else{
    B_full <- B_naive
  }

  list(A = solve(A_hat), B_naive = B_naive, B = B_full)
}



#' Sandwich Variance Estimator for Two-Stage IV AFT
#'
#' Computes the asymptotic variance components \eqn{A} and \eqn{B} for the two-stage
#' instrumental variable accelerated failure time (IV AFT) model.
#'
#' Under Theorem 1, the estimator satisfies the asymptotic normality:
#' \deqn{\sqrt{n}(\hat{\theta} - \theta) \to N(0, A^{-1} B A^{-T}).}
#' This function computes \eqn{\hat{A}} and the naive estimate \eqn{\hat{B}} (excluding
#' the contribution from \eqn{\Psi^*}), which is currently left as a placeholder.
#'
#' @param fs A list returned by \code{first_stage()}, containing:
#'   \code{beta} (i.e., \eqn{\alpha}), \code{fitted}, and \code{resid1}.
#' @param ss A list returned by \code{second_stage_wls()}, containing:
#'   \code{theta} (i.e., \eqn{\beta}), \code{fitted}, \code{resid2}, and \code{omega2}.
#' @param Cov A numeric vector of covariates used in stage 1.
#' @param Cov2 A numeric vector of covariates (excluding the endogenous variable) used in stage 2.
#' @param Naive Logical; whether to use the naive sandwich estimator (ignoring \eqn{\Psi^*}).
#'
#' @return A list with components:
#' \describe{
#'   \item{A}{A \eqn{4 \times 4} matrix estimating \eqn{\hat{A}} (block derivatives).}
#'   \item{B_naive}{A \eqn{4 \times 4} matrix estimating naive \eqn{\hat{B} = \frac{1}{n} \sum_j \Psi_j \Psi_j^T}.}
#'   \item{B}{A \eqn{4 \times 4} matrix estimating full \eqn{\hat{B}} (currently equals \code{B_naive}).}
#' }
#' @export
compute_AB2 <- function(fs, ss, Cov, Cov2, Naive = FALSE) {

  n <- length(fs$fitted)
  # unpack parameters
  alpha <- fs$beta           # c(alpha0, alpha1)
  beta  <- ss$theta          # c(beta0, beta1)
  xhat  <- fs$fitted         # n×1
  #x     <- fs$fitted+fs$resid
  resid1 <- fs$resid
  resid2 <- ss$resid         # Ystar - fitted second stage
  omega2 <- ss$w        # variance weights input

  #--- Estimate A_hat ----------------------------------------
  # dPsi1/dalpha = -n^{-1} sum [1;Z] %*% t([1;Z])
  M1 <- cbind(1, Cov)          # n×2
  d11 <- - crossprod(M1, M1) / n
  # dPsi2/dbeta  = -n^{-1} sum omega2 * [1; alpha0+alpha1 Z] %*% t([1; alpha0+alpha1 Z])
  M22 <- cbind(1, xhat, Cov2)        # n×2
  d22 <- - crossprod(M22 * sqrt(omega2), M22 * sqrt(omega2)) / n
  # dPsi2/dalpha = n^{-1} sum omega2 * [ -beta1; resid2 - beta1*xhat ] %*% t([1;Z])
  m211 <- t(omega2 * resid2) %*% M1
  m212 <- - beta[2] * t(omega2 * M22) %*% M1
  m212[2,] <- m212[2,] + m211
  d21 <- m212 / n
  # dPsi1/dbeta = 0
  d12 <- matrix(0, dim(d11)[1], dim(d22)[1])
  # assemble A_hat
  A_hat <- rbind(cbind(d11, d12), cbind(d21, d22))

  #--- Estimate B_naive -------------------------------------
  # psi1_j = resid1_j * [1; Z_j]
  Psi1 <- resid1 * M1
  # psi2_j = omega2_j * resid2_j * [1; xhat_j]
  Psi2 <- omega2 * resid2 * M22
  # full Psi = [Psi1 | Psi2]
  Psi   <- cbind(Psi1, Psi2)     # n×4
  B_naive <- crossprod(Psi) / n

  if(Naive == FALSE){
    ## ----  Ψ* contribution -------------------------------------------------
    ## --- 1. censoring KM
    G_surv <- ss$syn$fitG            # survfit object for S_C(t)
    jump_t <- G_surv$time
    #Sjump  <- G_surv$surv[-length(G_surv$surv)]
    #dt_g   <- diff(G_surv$time)           # interval widths Δt_k
    Sjump  <- G_surv$surv      # left-continuous S_C(t-)
    dt_g   <- c(diff(jump_t), 0)
    ## --- 2. observed (or synthetic) times
    Y_obs  <- ss$syn$time            # ≡ \tilde y_i  or Ystar

    ## --- 3. at-risk counts at each jump time s > 0
    at_risk <- vapply(jump_t[-1], function(s) sum(Y_obs >= s), numeric(1))
    at_risk[at_risk == 0] <- 1

    ## --- 4. indicator matrices
    risk_ind  <- outer(Y_obs, jump_t[-1], `>=`) # n * n-1

    ## Censoring increments with tolerance (only for censored subjects)
    m <- length(jump_t[-1])
    DeltaN <- matrix(0L, n, m)
    rows   <- which(ss$syn$status == 0L)            # rows that were censored
    idx    <- nearest_knot_idx(Y_obs[rows], jump_t[-1])  # which knot each censored time maps to
    keep   <- idx > 0L
    DeltaN[cbind(rows[keep], idx[keep])] <- 1L

    ## Column totals (pooled censoring jumps)
    DeltaN_tot <- colSums(DeltaN)

    ## --- 5. H_{ij}(s) matrix

    H_mat <- H_mat_cpp(Y_obs, jump_t[-1], dt_g, Sjump)

    # H_mat is n × m  (rows = subjects, cols = jump times)


    p_beta <- ncol(M22)            # number of β–components
    gprime_mat <- M22                # v_i = ∂g_i/∂β row-wise
    Psi2_star  <- matrix(0, n, p_beta)

    ## Build martingale increments (n × m) without recycling big vectors
    # Mart = (DeltaN - risk_ind * (DeltaN_tot / at_risk)) / at_risk
    haz   <- DeltaN_tot / at_risk                    # length m
    denom <- at_risk                                 # length m
    Mart  <- sweep(DeltaN,  2, denom, "/") -
      sweep(risk_ind, 2, haz / (denom), "*")  # equivalent to dividing by denom after

    # Precompute once to avoid repeating work inside the apply
    WH0 <- H_mat * Mart                    # n × m
    Psi2_star <- sapply(seq_len(p_beta), function(k) {
      rowSums(WH0 * (omega2 * gprime_mat[, k]))  # row-scale by omega2*v_{ik}, then sum over m
    })

    # Ensure matrix shape when p_beta = 1
    Psi2_star <- as.matrix(Psi2_star)           # n × p_beta
    ## 5) plug into B_full
    B_full <- crossprod(cbind(Psi1, Psi2 + Psi2_star)) / n
  }else{
    B_full <- B_naive
  }

  list(A = solve(A_hat), B_naive = B_naive, B = B_full)
}


#' Sandwich Variance Estimator for Two-Stage IV AFT
#'
#' Computes the asymptotic variance components \eqn{A} and \eqn{B} for the two-stage
#' instrumental variable accelerated failure time (IV AFT) model.
#'
#' Under Theorem 1, the estimator satisfies the asymptotic normality:
#' \deqn{\sqrt{n}(\hat{\theta} - \theta) \to N(0, A^{-1} B A^{-T}).}
#' This function computes \eqn{\hat{A}} and the naive estimate \eqn{\hat{B}} (excluding
#' the contribution from \eqn{\Psi^*}), which is currently left as a placeholder.
#'
#' @param fs A list returned by \code{first_stage()}, containing:
#'   \code{beta} (i.e., \eqn{\alpha}), \code{fitted}, and \code{resid1}.
#' @param ss A list returned by \code{second_stage_wls()}, containing:
#'   \code{theta} (i.e., \eqn{\beta}), \code{fitted}, \code{resid2}, and \code{omega2}.
#' @param Cov A numeric vector of covariates used in stage 1.
#' @param Cov2 A numeric vector of covariates (excluding the endogenous variable) used in stage 2.
#' @param Naive Logical; whether to use the naive sandwich estimator (ignoring \eqn{\Psi^*}).
#'
#' @return A list with components:
#' \describe{
#'   \item{A}{A \eqn{4 \times 4} matrix estimating \eqn{\hat{A}} (block derivatives).}
#'   \item{B_naive}{A \eqn{4 \times 4} matrix estimating naive \eqn{\hat{B} = \frac{1}{n} \sum_j \Psi_j \Psi_j^T}.}
#'   \item{B}{A \eqn{4 \times 4} matrix estimating full \eqn{\hat{B}} (currently equals \code{B_naive}).}
#' }
#' @export
compute_AB3 <- function(fs, ss, Cov, Cov2, Naive = FALSE) {

  n <- length(fs$fitted)
  # unpack parameters
  alpha <- fs$beta           # c(alpha0, alpha1)
  beta  <- ss$theta          # c(beta0, beta1)
  xhat  <- fs$fitted         # n×1
  #x     <- fs$fitted+fs$resid
  resid1 <- fs$resid
  resid2 <- ss$resid         # Ystar - fitted second stage
  omega2 <- ss$w        # variance weights input

  #--- Estimate A_hat ----------------------------------------
  # dPsi1/dalpha = -n^{-1} sum [1;Z] %*% t([1;Z])
  M1 <- cbind(1, Cov)          # n×2
  d11 <- crossprod(M1, M1) / n
  # dPsi2/dbeta  = -n^{-1} sum omega2 * [1; alpha0+alpha1 Z] %*% t([1; alpha0+alpha1 Z])
  M22 <- cbind(1, xhat, Cov2)        # n×2
  d22 <- crossprod(M22 * sqrt(omega2), M22 * sqrt(omega2)) / n
  # dPsi2/dalpha = n^{-1} sum omega2 * [ -beta1; resid2 - beta1*xhat ] %*% t([1;Z])
  m211 <- t(omega2 * resid2) %*% M1
  m212 <- beta[2] * t(omega2 * M22) %*% M1
  m212[2,] <- m212[2,] - m211
  d21 <- m212 / n
  # dPsi1/dbeta = 0
  d12 <- matrix(0, dim(d11)[1], dim(d22)[1])
  # assemble A_hat
  A_hat <- rbind(cbind(d11, d12), cbind(d21, d22))

  #--- Estimate B_naive -------------------------------------
  # psi1_j = resid1_j * [1; Z_j]
  Psi1 <- resid1 * M1
  # psi2_j = omega2_j * resid2_j * [1; xhat_j]
  Psi2 <- omega2 * resid2 * M22
  # full Psi = [Psi1 | Psi2]
  Psi   <- cbind(Psi1, Psi2)     # n×4
  B_naive <- crossprod(Psi) / n

  if(Naive == FALSE){
    ## ----  Ψ* contribution -------------------------------------------------
    ## --- 1. censoring KM
    G_surv <- ss$syn$fitG            # survfit object for S_C(t)
    jump_t <- G_surv$time
    #Sjump  <- G_surv$surv[-length(G_surv$surv)]
    #dt_g   <- diff(G_surv$time)           # interval widths Δt_k
    Sjump  <- G_surv$surv      # left-continuous S_C(t-)
    dt_g   <- c(diff(jump_t), 0)
    ## --- 2. observed (or synthetic) times
    Y_obs  <- ss$syn$time            # ≡ \tilde y_i  or Ystar

    ## --- 3. at-risk counts at each jump time s > 0
    at_risk <- vapply(jump_t[-1], function(s) sum(Y_obs >= s), numeric(1))
    at_risk[at_risk == 0] <- 1

    ## --- 4. indicator matrices
    risk_ind  <- outer(Y_obs, jump_t[-1], `>=`) # n * n-1

    ## Censoring increments with tolerance (only for censored subjects)
    m <- length(jump_t[-1])
    DeltaN <- matrix(0L, n, m)
    rows   <- which(ss$syn$status == 0L)            # rows that were censored
    idx    <- nearest_knot_idx(Y_obs[rows], jump_t[-1])  # which knot each censored time maps to
    keep   <- idx > 0L
    DeltaN[cbind(rows[keep], idx[keep])] <- 1L

    ## Column totals (pooled censoring jumps)
    DeltaN_tot <- colSums(DeltaN)

    ## --- 5. H_{ij}(s) matrix

    H_mat <- H_mat_cpp(Y_obs, jump_t[-1], dt_g, Sjump)

    # H_mat is n × m  (rows = subjects, cols = jump times)

    ## Build martingale increments (n × m) without recycling big vectors
    # Mart = (DeltaN - risk_ind * (DeltaN_tot / at_risk)) / at_risk
    haz   <- DeltaN_tot / at_risk                    # length m
    ## 1) martingale increments for subject j (no extra /Y here):
    Mart <- DeltaN - sweep(risk_ind, 2, haz, "*")     # n × m, entries = ΔM_j(s_k)

    ## 2) timewise sums over i for each β-component, then divide by Y(s_k):
    ##    S_mat (m × p_beta) with column r:  S_k^{(r)} / Y(s_k)
    S_mat <- t(H_mat) %*% (omega2 * M22)             # m × p_beta, sum_i ω_i v_{i,r} H_{i}(s_k)
    S_mat <- sweep(S_mat, 1, at_risk, "/")           # divide each row k by Y(s_k)

    ## 3) assemble Ψ2* for all subjects j and all β-components r:
    Psi2_star <- Mart %*% S_mat                      # (n × m) %*% (m × p_beta) -> n × p_beta

    # Ensure matrix shape when p_beta = 1
    Psi2_star <- as.matrix(Psi2_star)           # n × p_beta
    ## 5) plug into B_full
    B_full <- crossprod(cbind(Psi1, Psi2 + Psi2_star)) / n
  }else{
    B_full <- B_naive
  }

  list(A = solve(A_hat), B_naive = B_naive, B = B_full)
}



#' Map times to nearest KM knot (with tolerance)
#'
#' Returns the index (1..m) of the nearest knot in \code{knots} for each
#' element of \code{x}, or 0L if no knot is within tolerance.
#'
#' @param x Numeric vector of times to map.
#' @param knots Sorted (or unsorted) numeric vector of KM knot times.
#' @param tol Optional numeric tolerance. If \code{NULL}, a safe default based
#'   on knot spacing and machine precision is used.
#' @return An integer vector the same length as \code{x} with values in
#'   \code{0,1,2,...,length(knots)}; 0 means “no matching knot”.
#' @keywords internal
#' @noRd
nearest_knot_idx <- function(x, knots, tol = NULL) {
  knots <- sort(unique(knots))
  m <- length(knots)
  if (m == 0L) return(integer(length(x)))

  # default tol: half of the minimum spacing, bounded below by machine eps scale
  if (is.null(tol)) {
    d <- diff(knots)
    mind <- if (length(d)) min(d[d > 0]) else Inf
    tol <- max(1e-8, 100 * .Machine$double.eps * max(1, abs(knots)), 0.5 * mind)
  }

  # Voronoi midpoints define bins for nearest-neighbor via findInterval
  mids <- c(-Inf, (head(knots, -1) + tail(knots, -1)) / 2, Inf)
  idx  <- findInterval(x, mids)               # 1..m (bin index = nearest knot)

  # keep only if actually within tol of that nearest knot
  near <- abs(x - knots[idx]) <= tol
  idx[!near] <- 0L                            # 0L means "no matching knot"
  idx
}
