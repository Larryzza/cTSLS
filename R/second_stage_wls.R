## R/zzz.R  (or at the top of R/second_stage.R)


#' Second‐stage weighted least squares on synthetic outcome (via lm)
#'
#' @param syn results from \code{synthetic_leurgans()}.
#' @param endogenous  Numeric fitted exposure from \code{first_stage()}.
#' @param covars Optional matrix/data.frame of extra covariates.
#' @param tol Numeric; convergence tolerance on change in coefficients.
#' @param maxit Integer; maximum number of iterations.
#' @return A list with \code{theta}, \code{resid}, \code{qr2}, \code{WX}.
#' @export
second_stage_wls <- function(syn, endogenous, covars = NULL,
                             tol = 1e-3, maxit = 50) {
  Ystar <- syn$Ystar
  fitG <- syn$fitG
  w    <- rep(1/stats::var(Ystar), length(Ystar))
  df   <- as.data.frame(cbind(Ystar = Ystar, endogenous = endogenous, covars))
  # build formula string
  rhs  <- if(is.null(covars) | dim(covars)[2] == 0){
    "endogenous"
  }else{
    paste("endogenous +", paste(colnames(covars), collapse = "+"))
  }
  form <- stats::as.formula(paste("Ystar ~", rhs))
  # fit weighted lm
  ss <- stats::lm(form, data = df, weights = w)

  for (iter in seq_len(maxit)) {
    if(maxit == 1) break
    theta_old <- ss$coefficients
    resid     <- ss$residuals

    # 4b) estimate F_i(s)
    fitF <- survival::survfit(survival::Surv(resid, syn$status) ~ 1)

    # Build the variance function once
    var_fun <- make_var_fun(Ystar, fitF, fitG)
    # Compute all subject‐specific variances in one go
    Vnew <- vapply(ss$fitted.values, var_fun, numeric(1))

    # guard against nonpositive/NA
    Vnew[is.na(Vnew) | Vnew <= 0] <- stats::var(Ystar)
    w <- 1 / Vnew
    # 4d) re‐fit weighted LS
    ss <- stats::lm(form, data = df, weights = w)

    # 4e) check convergence
    if (max(abs(ss$coefficients - theta_old)) < tol) {
      message("Converged in ", iter, " iterations")
      break
    }
  }
  # extract pieces
  theta <- stats::coef(ss)
  resid <- stats::residuals(ss)
  #qr2   <- ss$qr
  #WX    <- stats::model.matrix(ss) * sqrt(w)
  list(theta = theta, resid = resid, w = w, syn = syn)
}

#' Second‐stage weighted least squares on synthetic outcome (via lm)
#'
#' @param syn results from \code{synthetic_leurgans()}.
#' @param endogenous  Numeric fitted exposure from \code{first_stage()}.
#' @param covars Optional matrix/data.frame of extra covariates.
#' @param tol Numeric; convergence tolerance on change in coefficients.
#' @param maxit Integer; maximum number of iterations.
#' @return A list with \code{theta}, \code{resid}, \code{qr2}, \code{WX}.
#' @export
second_stage_wls2 <- function(syn, endogenous, covars = NULL,
                              tol = 1e-3, maxit = 50) {

  Ystar <- syn$Ystar
  fitG  <- syn$fitG
  w     <- rep(1 / var(Ystar), length(Ystar))

  ## ---- static pieces ----------------------------------------------------
  jump_t <- fitG$time          # knots  t₀,…,t_m  (include 0)
  S_C    <- fitG$surv          # S_C(tₖ) = 1-G(tₖ)
  tail_ok <- S_C[length(S_C)] > 0

  dt_g   <- diff(jump_t)                              # Δt_k
  Kvec   <- (1 - S_C) / pmax(S_C, 1e-12)      #  G / (1-G) at left of int.

  cumK   <- cumsum(Kvec[-length(S_C)] * dt_g)      # C_l = Σ_{j≤l} K_j Δt_j   (length m)

  ## store for later
  pre_grid <- list(t     = jump_t[-1],   # length m
                   dt    = dt_g,         # length m
                   cumK  = cumK)

  ## ---- model frame ------------------------------------------------------
  df  <- data.frame(Ystar = Ystar, endogenous, covars)
  rhs  <- if(is.null(covars) | dim(covars)[2] == 0){
    "endogenous"
  }else{
    paste("endogenous +", paste(colnames(covars), collapse = "+"))
  }
  form <- stats::as.formula(paste("Ystar ~", rhs))
  ss   <- stats::lm(form, data = df, weights = w)
  converge <- NULL
  for (iter in seq_len(maxit)) {
    if(maxit == 1) break
    if(mean(syn$status) == 1) break

    theta_old <- coef(ss)
    resid     <- residuals(ss)

    ## 1. update empirical residual CDF
    fitF <- survival::survfit(survival::Surv(resid, syn$status) ~ 1)
    Pk   <- stepfun(fitF$time, c(1, fitF$surv))(pre_grid$t)
    Pk_dt <- Pk * pre_grid$dt                          # vector length m

    ## 2. subject-specific Var(Y*)  (C++ or R)

    add <- add_term_cpp2(ss$fitted.values,
                         pre_grid$t,
                         pre_grid$cumK,
                         Pk_dt,
                         Kvec,
                         tail_ok)
    Vnew <- var(Ystar) + 2 * as.numeric(add)
    Vnew[Vnew <= 0 | is.na(Vnew)] <- var(Ystar)   # guard
    w <- 1 / Vnew
    ss <- stats::lm(form, data = df, weights = w)

    if (max(abs(coef(ss) - theta_old)) < tol) {
      #message("Converged in ", iter, " iterations")
      converge <- TRUE
      break
    }else{
      converge <- FALSE
    }
  }

  list(theta = coef(ss),
       resid = residuals(ss),
       w     = w,
       syn   = syn,
       converge = converge)
}



#' Build a per-subject variance function for y*_i
#'
#' @param Ystar  Vector of synthetic outcomes.
#' @param fitF   survfit object for Surv(resid, status) — the KM for residuals.
#' @param fitG   survfit object for Surv(time, 1-status) — the KM for censoring.
#' @return        A function f(xi_beta) that returns Var(y*_i) for xi_beta.
#' @keywords internal
make_var_fun <- function(Ystar, fitF, fitG) {
  # 1) global var(Ystar)
  varY   <- stats::var(Ystar)
  # 3) extract censor‐KM knots
  knot_t <- fitG$time
  S_c    <- fitG$surv      # S_censor(t) = P(C > t)
  F_c    <- 1 - S_c              # F_C(t) = P(C <= t)
  dt_g   <- diff(knot_t)
  # 2) extract residual‐KM knots
  knot_s <- fitF$time[which(fitF$time<max(knot_t))]
  S_r    <- fitF$surv[which(fitF$time<max(knot_t))]      # S_resid(s)
  ds_f   <- diff(c(knot_s, max(knot_t)))

  # 4) return the actual variance function
  function(xi_beta) {
    # inner integral at each s-knot
    inner_vals <- vapply(knot_s, function(s) {
      upper <- xi_beta + s
      idx   <- which(knot_t <= upper)
      if (length(idx) == 0) return(0)
      # adjust the last bin width to exactly hit `upper`
      dt_sel <- dt_g[idx]
      last   <- length(idx)
      # new width = (upper − knot_t[last])
      dt_sel[last] <- upper - knot_t[last]
      # integrand = F_c(t)/S_c(t)
      ratios <- F_c[idx] / pmax(S_c[idx], 1e-5)
      sum(ratios * dt_sel)
    }, numeric(1))

    # outer integral approximation via Riemann sum
    # ∫[1−F(s)]·inner(s) ds ≈ sum((S_r)*inner_vals*ds_f)
    outer_val <- sum(S_r * inner_vals * ds_f)

    # final Var(y*_i)
    varY + 2 * outer_val
  }
}


