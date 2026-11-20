# Simulate Data From a Factorial ANCOVA Design (up to Four Factors, Any Number of Covariates)

Generates random data appropriate for an analysis of covariance with up
to four crossed fixed factors (\\A, B, C, D\\) and one or more
continuous covariates. Within every cell the outcome \\Y\\ and the
covariates \\X_1, \ldots, X_q\\ are jointly multivariate normal with a
common (homogeneous) within-cell covariance structure, the standard
assumption underlying classical ANCOVA. Per-cell sample sizes may be
equal or unequal.

## Usage

``` r
simulate_ancova_factorial_data(
  a,
  b = 1,
  c = 1,
  d = 1,
  n_covariates = 1,
  true_mean_Y,
  true_mean_X,
  true_sd_Y,
  true_sd_X,
  true_cor_YX,
  true_cor_XX = NULL,
  n_per_cell,
  randomized = TRUE,
  seed = NULL
)
```

## Arguments

- a:

  Number of levels of the first factor (must be at least `2`).

- b:

  Number of levels of the second factor (default `1`; use `1` if
  absent).

- c:

  Number of levels of the third factor (default `1`; use `1` if absent).

- d:

  Number of levels of the fourth factor (default `1`; use `1` if
  absent). So `a = 3, b = 2, c = 1, d = 1` specifies a 3 \\\times\\ 2
  design.

- n_covariates:

  Integer \\\ge 1\\: the number of continuous covariates. Default `1`.

- true_mean_Y:

  Numeric vector of length \\a \cdot b \cdot c \cdot d\\ giving the
  population cell means of \\Y\\, in the same row order as
  [`expand.grid`](https://rdrr.io/r/base/expand.grid.html)`(A, B, C, D)`
  (which varies factor `A` fastest, then `B`, then `C`, then `D`).

- true_mean_X:

  Numeric matrix of dimension \\(a \cdot b \cdot c \cdot d) \times q\\
  giving the population cell means of each covariate (rows = cells in
  the same order as `true_mean_Y`; columns = covariates). When
  `randomized = TRUE` every column must be constant across cells.

- true_sd_Y:

  Within-cell standard deviation of \\Y\\ (a single positive number,
  common across cells).

- true_sd_X:

  Within-cell standard deviations of the covariates. Either a single
  number (recycled to all `n_covariates`) or a numeric vector of length
  `n_covariates`.

- true_cor_YX:

  Within-cell correlations between \\Y\\ and each covariate. Either a
  single number (recycled to all `n_covariates`) or a numeric vector of
  length `n_covariates`.

- true_cor_XX:

  Within-cell correlation matrix among the covariates, \\q \times q\\.
  Default `NULL`, interpreted as the `n_covariates`-dimensional
  identity.

- n_per_cell:

  A single number (equal sample size per cell) or a numeric vector of
  length \\a \cdot b \cdot c \cdot d\\ giving per-cell sample sizes (in
  the same row order as `true_mean_Y`).

- randomized:

  Logical. `TRUE` (the default) for a randomized design, each cell is
  sampled from the same population covariate distribution, so
  `true_mean_X` must be constant across cells. `FALSE` for
  non-randomized / preexisting-groups designs in which the cells'
  population covariate means may differ.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A long-format `data.frame` with one row per simulated subject and the
following columns:

- `A`, `B`, `C`, `D`:

  Factor columns for each present design factor (omitted when the factor
  is absent, i.e., its corresponding `a/b/c/d` argument is `1`).

- `x1`, `x2`, ..., `x<q>`:

  Numeric simulated covariates.

- `y`:

  Numeric simulated outcome.

## Details

This is the factorial generalization of
[`simulate_ancova_data`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md)
(which is the special case `b = c = d = 1, n_covariates = 1`). All cells
share the same within-cell covariance structure (homogeneity of
regression, the classical ANCOVA assumption); the difference between
randomized and non-randomized designs lies entirely in whether the cell
covariate means are constrained to be equal.

**Why `true_mean_X` must be constant across cells when
`randomized = TRUE`.** Random assignment forms each cell as an
exchangeable random sample from the same joint distribution of
covariates and outcome. Cell-specific covariate means would silently
break that interpretation. The function checks the constraint and stops
with an informative error if it is violated. (The same logic that powers
[`simulate_ancova_data`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md).)

**Cell ordering.** The function uses
[`expand.grid`](https://rdrr.io/r/base/expand.grid.html)'s convention,
factor `A` varies fastest, then `B`, then `C`, then `D`. So for a \\2
\times 3\\ design, the six cells of `true_mean_Y` are \\(A_1 B_1), (A_2
B_1), (A_1 B_2), (A_2 B_2), (A_1 B_3), (A_2 B_3)\\. If you build the
cell specification by passing the factor levels to `expand.grid` in the
same order, the row indexing automatically matches.

## See also

[`simulate_ancova_data`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md)
(one factor, one covariate special case),
[`simulate_anova_data`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
[`simulate_regression_data`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md)

Other data simulators:
[`simulate_ancova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
[`simulate_anova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
[`simulate_longitudinal_gompertz()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_gompertz.md),
[`simulate_longitudinal_logistic()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_logistic.md),
[`simulate_longitudinal_negative_exponential()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_negative_exponential.md),
[`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md),
[`simulate_longitudinal_richards()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md),
[`simulate_regression_data()`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. 2 x 2 randomized design, single covariate.
set.seed(113)
design_2x2 <- expand.grid(A = factor(1:2), B = factor(1:2))
design_2x2$true_mean_Y <- c(50, 60, 55, 65)   # cell means in expand.grid order

d1 <- simulate_ancova_factorial_data(
  a            = 2, b = 2,
  true_mean_Y  = design_2x2$true_mean_Y,
  true_mean_X  = matrix(10, nrow = 4, ncol = 1),  # constant covariate mean
  true_sd_Y      = 8,
  true_sd_X      = 3,
  true_cor_YX  = 0.40,
  n_per_cell   = 30
)
aggregate(y ~ A + B, data = d1, FUN = mean)
#>   A B        y
#> 1 1 1 48.65943
#> 2 2 1 60.48404
#> 3 1 2 58.24399
#> 4 2 2 65.81071

# 2. 3 x 2 nonrandomized design with two covariates and unequal n.
set.seed(113)
a <- 3; b <- 2; q <- 2
n_cells <- a * b
d2 <- simulate_ancova_factorial_data(
  a            = a, b = b,
  n_covariates = q,
  true_mean_Y  = c(50, 55, 60,  52, 58, 64),
  true_mean_X  = matrix(c(10, 11, 12,  9, 10, 11,
                           5,  6,  7,  4,  5,  6),
                         nrow = n_cells, ncol = q),
  true_sd_Y      = 8,
  true_sd_X      = c(3, 2),
  true_cor_YX  = c(0.40, 0.25),
  true_cor_XX  = matrix(c(1,   0.3,
                          0.3, 1),
                        nrow = 2),
  n_per_cell   = c(30, 25, 20,  35, 30, 25),
  randomized   = FALSE
)
head(d2)
#>   A B        x1         x2        y
#> 1 1 1 13.254525  8.5360960 50.29890
#> 2 1 1 11.325283  4.5916342 61.15584
#> 3 1 1  9.805804  6.4857128 56.11123
#> 4 1 1  8.257290 -0.4105572 39.99210
#> 5 1 1 10.668032  5.1076307 45.26514
#> 6 1 1  4.298724  3.7604657 36.75722

# 3. 2 x 2 x 2 randomized design with one covariate.
set.seed(113)
d3 <- simulate_ancova_factorial_data(
  a       = 2, b = 2, c = 2,
  true_mean_Y = c(50, 55, 52, 57,  53, 58, 55, 60),  # 8 cells in A-fastest order
  true_mean_X = matrix(10, nrow = 8, ncol = 1),
  true_sd_Y = 8,
  true_sd_X = 3,
  true_cor_YX = 0.40,
  n_per_cell = 25
)
table(d3$A, d3$B, d3$C)  # 25 per cell, 200 total
#> , ,  = 1
#> 
#>    
#>      1  2
#>   1 25 25
#>   2 25 25
#> 
#> , ,  = 2
#> 
#>    
#>      1  2
#>   1 25 25
#>   2 25 25
#> 
```
