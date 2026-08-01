# Simulate Data From a Multivariate Normal Multiple-Regression Model

Generates random data \\(Y, X_1, \ldots, X_p)\\ jointly multivariate
normal with user-specified marginal means, marginal SDs, and full
correlation structure. Useful as a backbone for sensitivity analyses,
Monte Carlo studies of regression sample size methods, and pedagogical
demonstrations.

## Usage

``` r
simulate_regression_data(
  N,
  p,
  rho_YX,
  rho_XX = NULL,
  mu_Y = 0,
  mu_X = 0,
  sigma_Y = 1,
  sigma_X = 1,
  seed = NULL,
  column_names = NULL
)
```

## Arguments

- N:

  The total sample size (a positive integer \\\ge p + 2\\).

- p:

  The number of predictor variables.

- rho_YX:

  A numeric vector of length `p` giving the population correlations
  between \\Y\\ and each predictor \\X_j\\.

- rho_XX:

  A \\p \times p\\ symmetric correlation matrix for the predictors.
  Defaults to the identity matrix (orthogonal predictors).

- mu_Y:

  The population mean of \\Y\\ (default `0`).

- mu_X:

  A numeric vector of length `p` giving the population mean of each
  predictor (default `0`, recycled across predictors).

- sigma_Y:

  The population standard deviation of \\Y\\ (default `1`, in which case
  the simulated \\Y\\ is on the standardized scale).

- sigma_X:

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
\end{pmatrix},\$\$ converted to a covariance matrix via the supplied
SDs, and `N` draws are taken using
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
[`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five orthogonal predictors, each correlating .30 with Y.
set.seed(113)
d <- simulate_regression_data(
  N      = 200,
  p      = 5,
  rho_YX = rep(0.30, 5)
)
summary(lm(y ~ ., data = d))$r.squared   # ~ 5 * 0.30^2 = 0.45
#> [1] 0.5074721

# Predictors with shared structure (exchangeable correlation matrix).
rho_XX <- matrix(0.5, nrow = 5, ncol = 5); diag(rho_XX) <- 1
simulate_regression_data(
  N      = 300,
  p      = 5,
  rho_YX = c(.50, .40, .30, .20, .10),
  rho_XX = rho_XX,
  seed   = 113
)[1:3, ]
#>            y         x1         x2         x3         x4         x5
#> 1 -0.3960700 -0.6224146  0.2159206 -0.3064253  0.3549176  0.1135693
#> 2 -0.1355309 -0.9516532 -2.1782904 -0.3684705 -1.4350247 -0.7092047
#> 3 -0.6852586 -0.9420999 -0.6252969  0.1668666 -0.2680763 -0.9912944
```
