# Simulate Data From a One-Covariate ANCOVA Model

Generates random data appropriate for an analysis of covariance with one
continuous outcome (\\Y\\) and one continuous covariate (\\X\\) crossed
with \\a\\ fixed groups. The covariate is treated as a random variable;
\\Y\\ and \\X\\ are jointly multivariate normal within each group. Both
the randomized-design case (population covariate mean common across
groups) and the non-randomized / preexisting-groups case (population
covariate means differ across groups; \\Y\\-on-\\X\\ correlation may
also differ) are supported. Per-group sample sizes may be equal or
unequal.

## Usage

``` r
simulate_ancova_data(
  mu_y,
  mu_x,
  sigma_y,
  sigma_x,
  rho,
  a,
  n,
  randomized = TRUE
)
```

## Arguments

- mu_y:

  A numeric vector of length `a` giving the population mean of \\Y\\ in
  each group.

- mu_x:

  When `randomized = TRUE`, a single number giving the common population
  mean of \\X\\ across all groups. When `randomized = FALSE`, a numeric
  vector of length `a` giving the population covariate mean in each
  group.

- sigma_y:

  The population standard deviation of \\Y\\ (assumed common across
  groups).

- sigma_x:

  The population standard deviation of \\X\\ (assumed common across
  groups).

- rho:

  The population correlation between \\Y\\ and \\X\\. When
  `randomized = TRUE`, must be a single number (see Details). When
  `randomized = FALSE`, may be a single number (recycled to all `a`
  groups) or a numeric vector of length `a` giving a distinct
  correlation in each group.

- a:

  The number of fixed levels of the grouping factor (i.e., the number of
  conditions in a fixed-effects ANCOVA design). Use this argument when
  groups are the design levels of interest. (For sample-selected groups,
  e.g., classrooms or schools randomly drawn from a population, the
  convention in this package is to use `J` instead.)

- n:

  A single number (equal sample size per group) or a numeric vector of
  length `a` giving the sample size in each group.

- randomized:

  Logical. `TRUE` (the default) for a randomized design (random
  assignment of subjects to groups, so that the population covariate
  mean is the same in every group). `FALSE` for a non-randomized /
  preexisting-groups design.

## Value

A long-format `data.frame` with one row per simulated subject and three
columns:

- `group`:

  A factor with `a` levels (`"1"`, ..., `as.character(a)`) identifying
  each subject's group.

- `y`:

  Numeric simulated outcome.

- `x`:

  Numeric simulated covariate.

This format is directly usable with
[`aov()`](https://rdrr.io/r/stats/aov.html),
[`lm()`](https://rdrr.io/r/stats/lm.html), and other model-fitting
functions.

## Details

Each group's \\(Y, X)\\ pairs are drawn from a bivariate normal
distribution with mean \\(\mu\_{Y,j}, \mu\_{X,j})\\ and covariance
\$\$\Sigma_j = \begin{pmatrix} \sigma_Y^2 & \rho_j\\\sigma_Y\\\sigma_X
\\ \rho_j\\\sigma_Y\\\sigma_X & \sigma_X^2 \end{pmatrix}.\$\$

**Why `rho` must be a single number when `randomized = TRUE`.** Random
assignment forms each group as an exchangeable random sample from the
same population. The bivariate distribution of \\(Y, X)\\ is therefore
the same in every group, including the correlation. Allowing `rho` to
differ across groups would silently break that interpretation and
produce data that no randomized design could plausibly have generated.
The function therefore stops with an error in that case; use
`randomized = FALSE` if you genuinely want group-specific correlations.

**Convention on group labels.** The argument `a` is used here (and
throughout DMAR's experimental-design functions) for the number of
*fixed* levels of a designed factor, the levels you intend to compare.
The letter `J` is reserved for the number of *sample- selected* groups,
e.g., when classrooms or schools are randomly sampled from a population
(a random-effects context).

## See also

[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html),
[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ss_aipe_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova.md)

Other data simulators:
[`simulate_ancova_factorial_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_factorial_data.md),
[`simulate_anova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
[`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md),
[`simulate_regression_data()`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Randomized design, two groups, equal n.
set.seed(113)
simple <- simulate_ancova_data(
  mu_y       = c(3, 5),
  mu_x       = 10,
  sigma_y    = 1,
  sigma_x    = 2,
  rho        = 0.8,
  a          = 2,
  n          = 20
)
head(simple)
#>   group         y         x
#> 1     1 2.9708554 10.327468
#> 2     1 3.2960080 13.118574
#> 3     1 3.0362070 11.751984
#> 4     1 2.0565073  7.354049
#> 5     1 2.0626997  9.086845
#> 6     1 0.6000121  6.949545

# 2. Four preexisting groups with different correlations and unequal n.
#    The first two groups share rho = 0.30; the second two share a larger
#    rho = 0.60. The four groups are not used in the data-generation
#    machinery beyond their per-group means and correlations -- they are
#    just four distinct populations being sampled. In a downstream
#    analysis these four groups could be cross-classified as a 2 x 2
#    factorial design (e.g., the first factor distinguishing groups 1-2
#    from groups 3-4, and the second factor distinguishing groups 1, 3
#    from groups 2, 4) and analyzed via factorial ANCOVA.
set.seed(113)
preexisting <- simulate_ancova_data(
  mu_y       = c(50, 55, 60, 65),
  mu_x       = c(10, 12, 11, 13),
  sigma_y    = 8,
  sigma_x    = 3,
  rho        = c(0.30, 0.30, 0.60, 0.60),
  a          = 4,
  n          = c(40, 35, 45, 30),
  randomized = FALSE
)
aggregate(cbind(y, x) ~ group, data = preexisting,
          FUN = function(z) round(c(mean = mean(z), sd = sd(z)), 2))
#>   group y.mean  y.sd x.mean  x.sd
#> 1     1  48.81  8.25   9.60  3.36
#> 2     2  56.16  7.65  13.30  2.55
#> 3     3  58.98  7.42  10.97  3.12
#> 4     4  65.39  7.44  13.37  3.36
```
