#' Generate bivariate error terms for simulation
#'
#' @description
#' Simulate bivariate error terms under six different scenarios for two-stage IV simulations.
#'
#' @param num Integer. Number of observations to generate.
#' @param scenario Integer in 1:6. Which error scenario to simulate:
#'   1. Single Gaussian cluster
#'   2. Five-component Gaussian (shifted means)
#'
#' @return Numeric matrix with two columns: \code{eps1}, \code{eps2}.
#'
#' @details generate bivariate error data
#'
#' @importFrom MASS mvrnorm
#' @importFrom MDBED rBED
#' @export
generate_error_data <- function(num = 100, scenario = 1) {
  if (!is.numeric(num) || num <= 0) stop("`num` must be a positive integer.")
  scenario <- as.integer(scenario)
  error_data <- switch(
    as.character(scenario),
    # 1. single Gaussian
    '1' = mvrnorm(n = num,
                  mu    = c(0.0, 0.0),
                  Sigma = matrix(c(0.5, -0.3,
                                   -0.3, 1), 2, 2)),
    # 2. Three-component heterogeneous Gaussian
    '2' = {
      mus <- list(
        c(5, 4), c(5, 1), c(5, 5)
      )
      Sigmas <- list(
        matrix(c(0.2,  0.70*sqrt(0.2*1.0),  0.70*sqrt(0.2*1.0),  1.0), 2, 2),  # sd=0.4, corr=0.3
        matrix(c(0.4,  0.50*sqrt(0.4*0.5),  0.50*sqrt(0.4*0.5),  0.5), 2, 2),  # sd=(0.6,0.4), corr=-0.5
        matrix(c(0.3, -0.90*sqrt(0.3*2.0), -0.90*sqrt(0.3*2.0),  2.0), 2, 2)   # sd=0.5, corr=-0.2
      )
      props <- c(0.5, 0.3, 0.2)
      do.call(rbind, lapply(1:3, function(i) {
        mvrnorm(n = round(num * props[i]), mu = mus[[i]], Sigma = Sigmas[[i]])
      }))
    },
    stop("`scenario` must be 1 - 2")
  )
  colnames(error_data) <- c("eps1", "eps2")
  return(error_data)
}

#' Generate censoring times to achieve a desired censoring rate
#'
#' @description
#' For true event times \code{Y}, find a normal-based censoring distribution that
#' yields approximately \code{c_rate} proportion censored.
#'
#' @param Y Numeric vector of true event times.
#' @param mu  mean censoring time
#' @param sd_c SD of the normal censoring noise.
#' @return A data.frame with columns:
#'   \describe{
#'     \item{time}{Observed time: \code{pmin(censor_time, Y)}.}
#'     \item{status}{Event indicator (1=event observed, 0=censored).}
#'   }
#' @importFrom stats sd rnorm
#' @export
generate_censoring_fixed <- function(Y, mu, sd_c) {
  C  <- stats::rnorm(length(Y), mu, sd_c)
  time   <- pmin(Y, C)
  status <- as.integer(C >= Y)      # 1 = event observed
  data.frame(time = time, status = status)
}


#' Simulate two‐stage IV survival data
#'
#' @description
#' Generate a right-censored dataset with the following structure:
#' \itemize{
#'   \item Two instruments \eqn{G_1, G_2 \sim N(0, 0.8^2)}
#'   \item Two observed confounders \eqn{U_1, U_2 \sim N(0, 1)}
#'   \item Bivariate error \eqn{(e_1, e_2)} drawn from one of six specified scenarios
#'   \item Exposure: \eqn{X = a_0 + G a_1 + U a_2 + e_1}
#'   \item True event time: \eqn{Y = b_0 + b_1 X + U b_2 + e_2}
#'   \item Right-censoring time generated from a normal distribution, calibrated to match the desired censoring rate (\code{c_rate})
#' }
#'
#' @param n Integer; sample size.
#' @param scenario Integer in 1:6; error scenario (passed to \code{generate_error_data}).
#' @param c_rate Numeric in (0,1); target censoring proportion.
#' @param c_scale Numeric; multiplies the SD of the censoring noise.
#' @param a0 Numeric intercept in exposure model.
#' @param a1 Numeric vector length-2; coefficients for \(G_1,G_2\).
#' @param a2 Numeric vector length-2; coefficients for \(U_1,U_2\).
#' @param b0 Numeric intercept in outcome model.
#' @param b1 Numeric; coefficient on \(X\) in outcome model.
#' @param b2 Numeric vector length-2; coefficients for \(U_1,U_2\) in outcome.
#' @param seed Optional integer; if set, passed to \code{set.seed()}.
#' @return A tibble with columns
#'   \item{time}{Observed (possibly censored) time}
#'   \item{status}{Event indicator (1=event, 0=censored)}
#'   \item{X}{Endogenous exposure}
#'   \item{G1, G2}{Instruments}
#'   \item{U1, U2}{Observed confounders}
#' @importFrom tibble tibble
#' @export
simulate_iv_data <- function(n,
                             scenario,
                             c_rate,
                             c_scale = 1,
                             mu_cen = NULL, sd_cen = NULL,
                             a0 = 0, a1 = c(0.5, 0.5), a2 = c(0.3, 0.3),
                             b0 = 0, b1 = 1, b2 = c(0.5, 0.5),
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # 1. instruments & confounders
  G1 <- rnorm(n,   sd = 0.8)
  G2 <- rnorm(n,   sd = 0.8)
  U1 <- rnorm(n)
  U2 <- rnorm(n)
  G  <- cbind(G1, G2)
  U  <- cbind(U1, U2)

  # 2. error terms
  eps <- generate_error_data(num = n, scenario = scenario)
  e1  <- eps[,1]
  e2  <- eps[,2]

  # 3. exposure model
  X <- as.numeric(a0 + G %*% a1 + U %*% a2 + e1)

  # 4. true event time
  Y_true <- as.numeric(b0 + b1 * X + U %*% b2 + e2)

  if(c_rate == 0){
    data.out <- tibble::tibble(
      time   = Y_true,
      status = 1,
      X      = X,
      G1     = G1,
      G2     = G2,
      U1     = U1,
      U2     = U2
    )
    return(list(data.out = data.out,
                mu_cen = NA,
                sd_cen = NA))
  }
  # 5. apply censoring
  if (is.null(mu_cen) || is.null(sd_cen)) {
    temp <- calibrate_censor_mu(Y_true, c_rate, c_scale = 1)
    mu_cen <- temp$mu
    sd_cen <- temp$sd_c
    cens <- generate_censoring_fixed(Y_true, mu_cen, sd_cen)
  } else {
    cens <- generate_censoring_fixed(Y_true, mu_cen, sd_cen)
  }

  # 6. assemble tibble
  data.out <- tibble::tibble(
    time   = cens$time,
    status = cens$status,
    X      = X,
    G1     = G1,
    G2     = G2,
    U1     = U1,
    U2     = U2
  )
  return(list(data.out = data.out,
              mu_cen = mu_cen,
              sd_cen = sd_cen))
}

#' Calibrate a normal censoring shift for a target rate
#'
#' @param Y_ref Numeric vector drawn from the *same distribution* you will use
#'   for truth-time generation in the simulation (a pilot sample, or a very
#'   large draw).
#' @param c_rate Desired censoring proportion (0–1).
#' @param c_scale Multiplier for sd(Y_ref) to define σ of the censoring noise.
#' @return A list with elements \code{mu} and \code{sd_c}; plug these into
#'   \code{generate_censoring_fixed()} below.
#' @export
calibrate_censor_mu <- function(Y_ref, c_rate, c_scale = 1) {
  if (!(c_rate > 0 && c_rate < 1)) stop("c_rate must be in (0,1)")
  sd_c <- stats::sd(Y_ref) * c_scale
  f <- function(mu) stats::pnorm((Y_ref - mu) / sd_c) |> mean() - c_rate
  mu <- stats::uniroot(f, c(min(Y_ref) - 5 * sd_c,
                            max(Y_ref) + 5 * sd_c))$root
  list(mu = mu, sd_c = sd_c)
}
