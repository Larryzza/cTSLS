# R/first_stage.R

#' First‐stage OLS for the endogenous exposure (via lm)
#'
#' @param W Numeric vector of the endogenous variable.
#' @param Cov Numeric matrix or data.frame of instrument(s).
#' @param add_intercept Logical; include intercept? Default TRUE.
#' @return A list with \code{beta}, \code{fitted}, \code{resid}, and the lm \code{model}.
#' @export
first_stage <- function(W, Cov, add_intercept = TRUE) {
  # Build a data frame for lm()

  df <- as.data.frame(Cov)
  df$W <- W

  # Create formula: W ~ ...  (or without intercept if requested)
  if (add_intercept) {
    form <- stats::as.formula(paste("W ~", paste(names(df)[-ncol(df)], collapse = " + ")))
  } else {
    form <- stats::as.formula(paste("W ~", paste(names(df)[-ncol(df)], collapse = " + "), "- 1"))
  }

  model  <- stats::lm(form, data = df)
  fitted <- as.vector(stats::fitted(model))
  resid  <- as.vector(stats::residuals(model))

  list(
    beta   = stats::coef(model),
    fitted = fitted,
    resid  = resid,
    model  = model
  )
}
