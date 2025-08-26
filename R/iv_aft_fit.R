#' @useDynLib cTSLS, .registration = TRUE
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
NULL

#' Fit a two‐stage IV accelerated‐failure‐time model
#'
#' @param formula_stage2 A second‐stage formula, e.g. \code{Surv(time, status) ~ W + X1 + X2}
#' @param formula_stage1 A first‐stage formula, e.g. \code{W ~ Z1 + Z2}
#' @param data           A data.frame containing all variables in the formulas
#' @param maxit          Integer; maximum number of iterations.
#' @param Naive          Whether use the naive sandwich estimator of variance
#' @return An object of class \code{cTSLS} with components:
#'   \describe{
#'     \item{theta}{Second‐stage coefficients (including intercept).}
#'     \item{var}{Placeholder for sandwich variance matrix (currently \code{NULL}).}
#'     \item{first_stage}{Output list from \code{first_stage()}.}
#'     \item{second_stage}{Output list from \code{second_stage_wls()}.}
#'     \item{call}{The matched call.}
#'   }
#' @export
#' @rdname cTSLS
iv_aft_fit <- function(formula_stage2, formula_stage1, data,
                       maxit = 10, Naive = FALSE) {
  # --- Stage 1: fit exposure model
  mf1 <- model.frame(formula_stage1, data)
  W   <- model.response(mf1)
  Cov   <- model.matrix(attr(mf1, "terms"), mf1)[, -1, drop = FALSE]
  fs  <- first_stage(W, Cov)

  # --- Synthetic outcome (Leurgans)
  endog   <- all.vars(formula_stage1[[2]])
  covs   <- setdiff(attr(terms(formula_stage2), "term.labels"), endog)
  formula_stage2 <- reformulate(covs, deparse(formula_stage2[[2]]))
  mf2    <- model.frame(formula_stage2, data)
  time   <- model.response(mf2)[,1]
  status <- model.response(mf2)[,2]
  syn    <- synthetic_leurgans(time, status)

  # --- Stage 2: weighted LS on synthetic outcome
  #    drop intercept column from covariates matrix

  Cov2 <- model.matrix(attr(mf2, "terms"), mf2)[, -1, drop = FALSE]
  ss <- second_stage_wls2(syn = syn, endogenous = fs$fitted,
                          covars = Cov2, maxit = maxit)
  AB <- compute_AB3(fs, ss, Cov, Cov2, Naive)


  n        <- nrow(data)
  p_alpha  <- 1 + NCOL(Cov)           # (Intercept) + Z’s
  p_beta   <- length(ss$theta)        # (Intercept) + endogenous + covars

  idx_alpha <- seq_len(p_alpha)
  idx_beta  <- p_alpha + seq_len(p_beta)

  V_full       <- AB$A %*% AB$B       %*% t(AB$A) / n
  V_full_naive <- AB$A %*% AB$B_naive %*% t(AB$A) / n

  ## --- Name the coefficients
  alpha_names <- c("(Intercept)", colnames(Cov))
  beta_names  <- c("(Intercept)", endog, colnames(Cov2))
  names(fs$beta) <- alpha_names
  names(ss$theta) <- beta_names

  ## --- Pull variances and SEs by block
  alpha_var        <- diag(V_full)[idx_alpha]
  alpha_var_naive  <- diag(V_full_naive)[idx_alpha]
  beta_var         <- diag(V_full)[idx_beta]
  beta_var_naive   <- diag(V_full_naive)[idx_beta]

  alpha_se <- sqrt(pmax(alpha_var, 0))
  beta_se  <- sqrt(pmax(beta_var, 0))

  structure(
    list(
      ## First stage
      alpha        = fs$beta,
      alpha_se     = alpha_se,
      alpha_var    = alpha_var,
      alpha_var_naive = alpha_var_naive,

      ## Second stage (keeps your original field names for backward compat)
      theta        = ss$theta,
      se           = beta_se,
      var          = beta_var,
      var_naive    = beta_var_naive,

      ## (Optional) expose full vcovs if you want them
      vcov_full        = V_full,
      vcov_full_naive  = V_full_naive,

      ## Raw stage outputs
      first_stage  = fs,
      second_stage = ss,
      converge     = ss$converge,
      call         = match.call()
    ),
    class = "cTSLS"
  )
}



#' Bootstrap standard errors for `cTSLS` coefficients
#'
#' Draws non–parametric bootstrap samples of the original data, re-fits
#' `iv_aft_fit()` on each resample, and returns bootstrap means, standard
#' errors, and percentile confidence intervals for the coefficient vector
#' \eqn{\theta = (\beta_0,\beta_1,\dots)}.
#'
#' @param formula_stage2 Second-stage formula (see `iv_aft_fit`).
#' @param formula_stage1 First-stage formula.
#' @param data           Original data frame.
#' @param R              Integer. Number of bootstrap replicates (default 500).
#' @param seed           Optional seed for reproducibility.
#' @param maxit          Integer; maximum number of iterations.
#' @return A list with components
#'   \itemize{
#'     \item \code{coef}: Vector of point estimates from the original sample.
#'     \item \code{boot}: \code{R × p} matrix of bootstrap estimates.
#'     \item \code{se}:   Bootstrap standard errors (column-wise sd).
#'     \item \code{ci}:   2.5\% and 97.5\% percentile intervals.
#'   }
#' @export
bootstrap_iv_aft <- function(formula_stage2, formula_stage1, data,
                             R = 500, seed = NULL, maxit = 10) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)

  # helper: fit on a resample index vector
  boot_fit <- function(idx) {
    d <- data[idx, , drop = FALSE]
    fit <- iv_aft_fit(formula_stage2, formula_stage1, d,
                      maxit = maxit, Naive = TRUE)
    if(is.null(fit$converge)) fit$converge <- TRUE
    return(c(fit$theta,fit$converge))
  }

  # run R bootstrap replicates sequentially
  boot_mat <- sapply(1:R, function(x)boot_fit(sample.int(n, n, replace = TRUE)))
  boot_mat <- replicate(R, boot_fit(sample.int(n, n, replace = TRUE)), simplify = "matrix")
  boot_mat <- boot_mat[-dim(boot_mat)[1], which(boot_mat[dim(boot_mat)[1],]==TRUE)]
  fit <- iv_aft_fit(formula_stage2, formula_stage1, data,
                    maxit = maxit, Naive = TRUE)
  point_est <- fit$theta
  if(is.null(fit$converge)==FALSE){
    if(isFALSE(fit$converge)) fit$theta <- NA
  }
  se_boot   <- apply(boot_mat, 1, stats::sd, na.rm = TRUE)
  ci_boot   <- t(apply(boot_mat, 1, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  colnames(ci_boot) <- c("2.5%", "97.5%")

  list(coef = point_est,
       boot = t(boot_mat),
       se   = se_boot,
       ci   = ci_boot)
}
