#' @rdname cTSLS
#' @export
coef.cTSLS <- function(object, stage = c("second", "first", "both"), ...) {
  stage <- match.arg(stage)
  if (stage == "second") return(object$theta)
  if (stage == "first")  return(object$alpha)
  # both
  out <- c(object$alpha, object$theta)
  names(out) <- c(paste0("alpha:", names(object$alpha)),
                  paste0("beta:",  names(object$theta)))
  out
}

#' @rdname cTSLS
#' @export
print.cTSLS <- function(x, ...) {
  cat("Two-stage AFT-IV model\n")
  if (!is.null(x$converge)) {
    cat("Second-stage WLS convergence:", if (isTRUE(x$converge)) "converged", " \n")
  }
  cat("\nFirst stage (alpha):\n")
  if (!is.null(x$alpha)) {
    se1 <- if (!is.null(x$alpha_se)) x$alpha_se else
      rep(NA_real_, length(x$alpha))
    tab1 <- cbind(Estimate = x$alpha, `Std. Error` = se1)
    print(tab1)
  } else {
    cat("(not available)\n")
  }
  cat("\nSecond stage (beta):\n")
  se2 <- if (!is.null(x$se)) x$se else {
    # fall back if only variances were stored
    if (!is.null(x$var)) sqrt(pmax(x$var, 0)) else rep(NA_real_, length(x$theta))
  }
  tab2 <- cbind(Estimate = x$theta, `Std. Error` = se2)
  print(tab2)
  invisible(x)
}

#' @export
print.cTSLS_table <- function(x, ...) {
  # print as a plain matrix/data.frame without the class footer
  if (is.matrix(x)) {
    print(unclass(x), ...)   # drops the class
  } else if (is.data.frame(x)) {
    print(x, ...)            # data.frame never shows the class footer
  } else {
    NextMethod()
  }
  invisible(x)
}

#' @export
print.cTSLS_table <- function(x, ...) {
  # print as a plain matrix/data.frame without the class footer
  if (is.matrix(x)) {
    print(unclass(x), ...)   # drops the class
  } else if (is.data.frame(x)) {
    print(x, ...)            # data.frame never shows the class footer
  } else {
    NextMethod()
  }
  invisible(x)
}

#' @rdname cTSLS
#' @export
summary.cTSLS <- function(object, ...) {
  ## First-stage table
  alpha <- object$alpha
  se1 <- if (!is.null(object$alpha_se)) object$alpha_se else {
    # try to extract from full vcov if present
    if (!is.null(object$vcov_full) && !is.null(alpha)) {
      p1 <- length(alpha)
      sqrt(pmax(diag(object$vcov_full)[seq_len(p1)], 0))
    } else rep(NA_real_, length(alpha))
  }
  z1 <- if (length(alpha)) alpha / se1 else numeric()
  p1 <- if (length(alpha)) 2 * stats::pnorm(-abs(z1)) else numeric()
  tab1 <- if (length(alpha)) {
    out <- cbind(Estimate = alpha, `Std. Error` = se1, `z value` = z1, `Pr(>|z|)` = p1)
    class(out) <- "cTSLS_table"
    out
  } else NULL

  ## Second-stage table
  beta <- object$theta
  se2  <- if (!is.null(object$se)) object$se else {
    if (!is.null(object$vcov_full) && length(beta)) {
      p1 <- length(object$alpha %||% numeric(0))
      idx_beta <- p1 + seq_along(beta)
      sqrt(pmax(diag(object$vcov_full)[idx_beta], 0))
    } else if (!is.null(object$var)) sqrt(pmax(object$var, 0)) else rep(NA_real_, length(beta))
  }
  z2 <- if (length(beta)) beta / se2 else numeric()
  p2 <- if (length(beta)) 2 * stats::pnorm(-abs(z2)) else numeric()
  tab2 <- if (length(beta)) {
    out <- cbind(Estimate = beta, `Std. Error` = se2, `z value` = z2, `Pr(>|z|)` = p2)
    class(out) <- "cTSLS_table"
    out
  } else NULL

  res <- list(first_stage  = tab1,
              second_stage = tab2,
              converge     = object$converge,
              call         = object$call)
  class(res) <- "summary.cTSLS"

  ## pretty-print
  cat("Two-stage AFT-IV model\n")
  if (!is.null(object$converge)) {
    cat("Second-stage WLS convergence:", if (isTRUE(object$converge)) "converged", " \n\n")
  }
  cat("First stage (alpha):\n")
  if (is.null(tab1)) cat("(not available)\n") else print(tab1)
  cat("\nSecond stage (beta):\n")
  if (is.null(tab2)) cat("(not available)\n") else print(tab2)
  invisible(res)
}

#' @rdname cTSLS
#' @param stage Which covariance to return: "second" (default), "first", or "full".
#' @param naive Logical; if TRUE, return the naive sandwich VCOV (if stored).
#' @export
vcov.cTSLS <- function(object, stage = c("second", "first", "full"), naive = FALSE, ...) {
  stage <- match.arg(stage)
  full <- if (!is.null(object$vcov_full)) {
    if (naive && !is.null(object$vcov_full_naive)) object$vcov_full_naive else object$vcov_full
  } else {
    # Fall back to diagonal-only if that's what we have
    p1 <- length(object$alpha %||% numeric(0))
    p2 <- length(object$theta %||% numeric(0))
    V  <- matrix(0, p1 + p2, p1 + p2)
    if (p1 && !is.null(object$alpha_var)) diag(V)[seq_len(p1)] <- object$alpha_var
    if (p2 && !is.null(object$var))       diag(V)[p1 + seq_len(p2)] <- object$var
    V
  }

  p1 <- length(object$alpha %||% numeric(0))
  p2 <- length(object$theta %||% numeric(0))
  if (stage == "full")  return(full)
  if (stage == "first") return(full[seq_len(p1), seq_len(p1), drop = FALSE])
  # second
  full[p1 + seq_len(p2), p1 + seq_len(p2), drop = FALSE]
}

#' @rdname cTSLS
#' @export
predict.cTSLS <- function(object, newdata = NULL, type = c("response","link"), ...) {
  type <- match.arg(type)
  stop("predict.cTSLS() not yet implemented for this object; ",
       "please construct the model matrix for the second stage and use %*% coefficients.")
}

## tiny helper (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b
