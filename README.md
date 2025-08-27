cTSLS: Two-Stage IV Estimation for AFT with Right-Censored Data
================

`cTSLS` package implements a censored two-stage least squares (cTSLS)
estimator for semiparametric accelerated failure time (AFT) models with
right-censored outcomes. It combines:

- Leurgans’ synthetic outcome to handle censoring

- Subject-specific weights for heteroskedastic Var($`Y_i^*`$)

- A sandwich variance estimator that adjusts for replacing G with
  $`\hat G`$

A nonparametric bootstrap for SEs is also included.

# Installation

``` r
# from source
# setwd("path/to/ivAFT")
# install.packages(".", repos = NULL, type = "source")

# or from GitHub (replace user/repo)
# install.packages("devtools")
# devtools::install_github("user/ivAFT")
```

# Quick start

``` r
library(cTSLS)
set.seed(123)

# simulate a toy dataset
sim <- simulate_iv_data(
  n        = 500,
  scenario = 1,   # 1: single Gaussian; 2: Gaussian mixture
  c_rate   = 0.25
)
df <- sim$data.out
str(df)
#> tibble [500 × 7] (S3: tbl_df/tbl/data.frame)
#>  $ time  : num [1:500] -1.64 -1.903 -0.715 0.422 -2.223 ...
#>  $ status: int [1:500] 1 1 0 0 1 0 1 1 0 1 ...
#>  $ X     : num [1:500] -0.41 -1.313 1.386 0.215 -1.48 ...
#>  $ G1    : num [1:500] -0.4484 -0.1841 1.247 0.0564 0.1034 ...
#>  $ G2    : num [1:500] -0.482 -0.795 0.821 0.601 -1.207 ...
#>  $ U1    : num [1:500] -0.996 -1.04 -0.018 -0.132 -2.549 ...
#>  $ U2    : num [1:500] -0.821 -0.307 -0.902 0.627 1.12 ...
```

``` r
# columns: time, status (1=event, 0=censored), X, G1, G2, U1, U2

# formulas
f2 <- survival::Surv(time, status) ~ X + U1 + U2   # stage 2
f1 <- X ~ G1 + G2 + U1 + U2                        # stage 1

# fit cTSLS
fit <- iv_aft_fit(
  formula_stage2 = f2,
  formula_stage1 = f1,
  data  = df,
  maxit = 10,        # WLS iterations for weights
  Naive = FALSE      # variance that accounts for estimating G
)

print(fit)
#> Two-stage AFT-IV model
#> Second-stage WLS convergence: converged  
#> 
#> First stage (alpha):
#>               Estimate Std. Error
#> (Intercept) 0.01788932 0.02989381
#> G1          0.55407997 0.03936244
#> G2          0.55130540 0.03725466
#> U1          0.24463179 0.02857051
#> U2          0.27000487 0.02810334
#> 
#> Second stage (beta):
#>                Estimate Std. Error
#> (Intercept) -0.07728071 0.04681048
#> X            0.93440881 0.09614824
#> U1           0.47061731 0.06512804
#> U2           0.45950350 0.06264622
```

``` r
summary(fit)   # first-stage (alpha) and second-stage (beta) tables
#> Two-stage AFT-IV model
#> Second-stage WLS convergence: converged  
#> 
#> First stage (alpha):
#>               Estimate Std. Error    z value     Pr(>|z|)
#> (Intercept) 0.01788932 0.02989381  0.5984289 5.495538e-01
#> G1          0.55407997 0.03936244 14.0763625 5.307205e-45
#> G2          0.55130540 0.03725466 14.7982939 1.502453e-49
#> U1          0.24463179 0.02857051  8.5623877 1.105537e-17
#> U2          0.27000487 0.02810334  9.6075717 7.428009e-22
#> 
#> Second stage (beta):
#>                Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept) -0.07728071 0.04681048 -1.650927 9.875338e-02
#> X            0.93440881 0.09614824  9.718419 2.516592e-22
#> U1           0.47061731 0.06512804  7.226032 4.973089e-13
#> U2           0.45950350 0.06264622  7.334896 2.218939e-13
```

``` r
vcov(fit)      # sandwich variance for beta block
#>                                                        
#>  2.191221e-03  1.189208e-05  0.0003427349  0.0002541799
#>  1.189208e-05  9.244483e-03 -0.0030547441 -0.0030552835
#>  3.427349e-04 -3.054744e-03  0.0042416612  0.0017973047
#>  2.541799e-04 -3.055284e-03  0.0017973047  0.0039245489
```

``` r
coef(fit)      # beta coefficients
#> (Intercept)           X          U1          U2 
#> -0.07728071  0.93440881  0.47061731  0.45950350
```

`Tip`: set Naive = TRUE to see the variance that ignores uncertainty in
G.
