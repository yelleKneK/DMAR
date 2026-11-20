# Simulate Data From a Multivariate Normal Multiple-Regression Model

Generates random data \\(Y, X_1, \ldots, X_p)\\ jointly multivariate
normal with user-specified marginal means, marginal SDs, and full
correlation structure. Useful for sensitivity analyses, Monte Carlo
studies of regression sample size methods, and teaching demonstrations.

## Usage

``` r
simulate_regression_data(
  N,
  p,
  true_cor_YX,
  true_cor_XX = NULL,
  true_mean_Y = 0,
  true_mean_X = 0,
  true_sd_Y = 1,
  true_sd_X = 1,
  seed = NULL,
  column_names = NULL
)
```

## Arguments

- N:

  The total sample size (a positive integer \\\ge p + 2\\).

- p:

  The number of predictor variables.

- true_cor_YX:

  A numeric vector of length `p` giving the population correlations
  between \\Y\\ and each predictor \\X_j\\.

- true_cor_XX:

  A \\p \times p\\ symmetric correlation matrix for the predictors.
  Defaults to the identity matrix (orthogonal predictors).

- true_mean_Y:

  The population mean of \\Y\\ (default `0`).

- true_mean_X:

  A numeric vector of length `p` giving the population mean of each
  predictor (default `0`, recycled across predictors).

- true_sd_Y:

  The population standard deviation of \\Y\\ (default `1`, in which case
  the simulated \\Y\\ is on the standardized scale).

- true_sd_X:

  A single number or a numeric vector of length `p` giving the
  population standard deviation of each predictor (default `1`,
  standardized).

- seed:

  Optional integer random seed for reproducibility (default `NULL`).

- column_names:

  Optional character vector of length `p + 1` giving column names for
  the returned `data.frame`; defaults to
  `c("y", "x1", "x2", ..., "xp")`.

## Value

A `data.frame` with `N` rows and `p + 1` columns: the outcome \\Y\\
(first column) followed by predictors \\X_1, \ldots, X_p\\.

## Details

Internally the joint correlation matrix is assembled as \$\$R =
\begin{pmatrix} 1 & \rho\_{YX}^\top \\ \rho\_{YX} & R\_{XX}
\end{pmatrix},\$\$ where \\\rho\_{YX}\\ is `true_cor_YX` and \\R\_{XX}\\
is `true_cor_XX`, converted to a covariance matrix via the supplied SDs,
and `N` draws are taken using
[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html). The resulting
\\Y\\ and predictors satisfy the requested marginal means and standard
deviations and (in expectation) the requested correlation structure.

## See also

[`simulate_ancova_data`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
[`simulate_anova_data`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)

Other data simulators:
[`simulate_ancova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
[`simulate_ancova_factorial_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_factorial_data.md),
[`simulate_anova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
[`simulate_longitudinal_gompertz()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_gompertz.md),
[`simulate_longitudinal_logistic()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_logistic.md),
[`simulate_longitudinal_negative_exponential()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_negative_exponential.md),
[`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md),
[`simulate_longitudinal_richards()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five orthogonal predictors, each correlating .30 with Y.
set.seed(113)
d <- simulate_regression_data(
  N      = 200,
  p      = 5,
  true_cor_YX = rep(0.30, 5)
)
summary(lm(y ~ ., data = d))$r.squared   # about 0.45: five predictors, each 0.30^2
#> [1] 0.5074721

# Predictors with shared structure (exchangeable correlation matrix).
true_cor_XX <- matrix(0.5, nrow = 5, ncol = 5); diag(true_cor_XX) <- 1
simulate_regression_data(
  N      = 300,
  p      = 5,
  true_cor_YX = c(.50, .40, .30, .20, .10),
  true_cor_XX = true_cor_XX,
  seed   = 113
)[1:3, ]
#>            y         x1         x2         x3         x4         x5
#> 1 -0.6456770 -0.3331409  0.3744711 -0.2785978  0.2520220 -0.1200495
#> 2 -0.7330446 -0.2591847 -1.7987493 -0.3018566 -1.6813382 -1.2684455
#> 3 -0.5521862 -1.0963196 -0.7098246  0.1520311 -0.2132198 -0.8667458
```
