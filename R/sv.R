#' Leurgans synthetic outcome for right‐censoring
#'
#' @param time  Numeric vector of observed times.
#' @param status Integer censoring indicator (1=event, 0=censored).
#' @return A list with \code{Ystar}, \code{sigma2}, \code{Ghat}.
#' @export
synthetic_leurgans <- function(time, status) {
  fitG <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  Ghat <- stepfun(fitG$time, c(1, fitG$surv))
  ord  <- order(time)
  dt   <- diff(c(0, time[ord]))
  invG <- 1 / Ghat(time[ord] - 1e-5)
  Ystar <- cumsum(dt * invG)[order(ord)]
  list(Ystar = Ystar, time = time, status = status, fitG = fitG)
}
