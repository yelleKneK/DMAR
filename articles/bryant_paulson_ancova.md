# Simultaneous Comparison of Adjusted Means in ANCOVA With DMAR

This vignette shows the end-to-end workflow for comparing
covariate-adjusted means in the analysis of covariance (ANCOVA): fit the
model, read off the adjusted means and the error term, and then place
**simultaneous** confidence intervals on the differences with
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
which uses the Bryant–Paulson generalized studentized range
([`qbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)).
The companion vignette, *Bryant-Paulson Simultaneous Intervals: A
Simulation Study*, shows **why** the ordinary Tukey distribution is the
wrong reference when the covariate is random. Here we focus on **how**
to do the analysis.

## 1. The Worked Example from Bryant & Bruvold (1980)

The `test_market` data ship with DMAR. A company compared six marketing
strategies (“panels”) for a brand across four blocks of retail outlets,
with the *remaining category movement* in each outlet as a random
covariate.

``` r

data(test_market)
descriptives(test_market[c("brand_movement", "category_movement")])
#> $descriptives
#>            variable    type  n n_missing prop_missing      mean median
#> 1    brand_movement numeric 24         0            0  4.220833   4.25
#> 2 category_movement numeric 24         0            0 12.196667  12.66
#>          sd  min   max    q25    q75     skewness   kurtosis
#> 1 0.6999995 2.98  5.61 3.8750  4.575  0.020064949 -0.3790459
#> 2 2.8304196 7.88 16.51 9.5925 14.235 -0.005730893 -1.5356582
#> 
#> $correlations
#> NULL
```

### Fit the ANCOVA

This is a randomized-block ANCOVA: panel (the treatment of interest),
block, and the covariate. We fit it with
[`lm()`](https://rdrr.io/r/stats/lm.html) and confirm it reproduces the
published quantities (covariate slope 0.4079, error mean square 0.01326
on 14 degrees of freedom).

``` r

fit <- lm(brand_movement ~ panel + block + category_movement, data = test_market)

s_ancova <- summary(fit)$sigma           # ANCOVA error SD = sqrt(MS error)
nu       <- fit$df.residual              # error degrees of freedom
c(slope = unname(coef(fit)["category_movement"]),
  error_MS = s_ancova^2, df = nu)
#>       slope    error_MS          df 
#>  0.40788407  0.01325859 14.00000000
```

### Adjusted Panel Means

The adjusted mean for each panel is the model’s prediction at the
covariate grand mean, averaged over blocks.

``` r

xbar <- mean(test_market$category_movement)
adj_means <- vapply(levels(test_market$panel), function(p) {
  nd <- data.frame(panel = factor(p, levels = levels(test_market$panel)),
                   block = factor(1:4, levels = levels(test_market$block)),
                   category_movement = xbar)
  mean(predict(fit, nd))
}, numeric(1))
round(adj_means, 3)     # 3.595 3.619 4.102 4.515 4.618 4.876
#>     1     2     3     4     5     6 
#> 3.595 3.619 4.102 4.515 4.618 4.876
```

### Simultaneous Bryant–Paulson Intervals

With six panels there are 15 pairwise comparisons. We want all of them
to hold at a *familywise* 95% level. Because the design has `s = 4`
blocks, the standard error of a single adjusted mean is built from
`n = 4`, and the error degrees of freedom (14) must be supplied directly
(the blocked design is not the one-way default `N - k - p`).

``` r

bp <- ci_c_ancova_bp(adj_means = adj_means, s_ancova = s_ancova,
                     n = 4, num_covariates = 1, df = nu)
bp
```

| contrast          | estimate | lower_limit | upper_limit |
|:------------------|:---------|:------------|:------------|
| group_1 - group_2 | -0.0243  | -0.302      | 0.254       |
| group_1 - group_3 | -0.507   | -0.785      | -0.228      |
| group_1 - group_4 | -0.92    | -1.2        | -0.642      |
| group_1 - group_5 | -1.02    | -1.3        | -0.745      |
| group_1 - group_6 | -1.28    | -1.56       | -1          |
| group_2 - group_3 | -0.482   | -0.76       | -0.204      |
| group_2 - group_4 | -0.896   | -1.17       | -0.618      |
| group_2 - group_5 | -0.998   | -1.28       | -0.72       |
| group_2 - group_6 | -1.26    | -1.53       | -0.978      |
| group_3 - group_4 | -0.413   | -0.691      | -0.135      |
| group_3 - group_5 | -0.516   | -0.794      | -0.238      |
| group_3 - group_6 | -0.774   | -1.05       | -0.496      |
| group_4 - group_5 | -0.103   | -0.381      | 0.175       |
| group_4 - group_6 | -0.361   | -0.639      | -0.0827     |
| group_5 - group_6 | -0.258   | -0.536      | 0.02        |

Confidence level: 95%

Every pairwise interval has the same half-width, and the critical
difference is exactly the 0.278 reported in the paper:

``` r

attr(bp, "critical_value")                      # q_.05;1,6,14 = 4.83
#> [1] 4.829856
unique(round((bp$upper_limit - bp$lower_limit) / 2, 3))   # 0.278
#> [1] 0.278
```

Reading the table: panels 1 and 2 are statistically indistinguishable
(their interval covers 0), but panel 1 differs from panels 3 through 6,
and so on, the same conclusions Bryant and Bruvold reached. In
substantive terms, the best panel (6) outsold the weakest (1) by 1.28
hundred cases (95% simultaneous CI \[1, 1.56\]), a gap the random
covariate adjustment was sharp enough to resolve.

## 2. Simultaneous vs. Per-Comparison Intervals

DMAR already provides
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md)
for a **single** contrast of adjusted means. It uses a *t* critical
value and includes the per-pair covariate term in the standard error,
correct for one pre-planned comparison, but it does **not** protect the
familywise error rate across many comparisons, and it treats the
covariate adjustment as fixed.
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
is its familywise, random-covariate counterpart.

``` r

# Per-comparison interval for panel 1 vs. panel 4 (needs the covariate means
# and the within-group SS of the covariate).
cov_means <- tapply(test_market$category_movement, test_market$panel, mean)
SSwx <- sum((test_market$category_movement -
             ave(test_market$category_movement, test_market$panel))^2)

per_comparison <- ci_c_ancova(
  adj_means = adj_means, s_ancova = s_ancova,
  c_weights = c(1, 0, 0, -1, 0, 0), n = 4,
  cov_means = cov_means, SSwithin_x = SSwx)
per_comparison
```

| term        | value  |
|:------------|:-------|
| lower_limit | -1.09  |
| psi         | -0.92  |
| upper_limit | -0.748 |

Confidence level: 95%

``` r


# The same contrast, but as one member of the simultaneous family:
ci_c_ancova_bp(adj_means = adj_means, s_ancova = s_ancova, n = 4,
               num_covariates = 1, df = nu,
               c_weights = c(1, 0, 0, -1, 0, 0))
```

| contrast   | estimate | lower_limit | upper_limit |
|:-----------|:---------|:------------|:------------|
| contrast_1 | -0.92    | -1.2        | -0.642      |

Confidence level: 95%

The simultaneous interval is wider; that is the price of protecting all
15 comparisons at once rather than just this one.

## 3. A One-Way ANCOVA With `ancova()`

For the common one-way case (a single treatment factor with several
levels and one or more covariates), DMAR’s
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md) gives
a tidy fit (omnibus *F*, effect size CIs, adjusted means, and a
homogeneity-of-regression check) and feeds directly into
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md).
Here is a four-group example.

``` r

set.seed(113)
k <- 4; n <- 25
group <- factor(rep(c("control", "low", "medium", "high"), each = n),
                levels = c("control", "low", "medium", "high"))
x <- rnorm(k * n, 50, 10)                          # random covariate
mu <- c(control = 0, low = 1.5, medium = 3, high = 3.2)   # adjusted effects
y <- 0.5 * (x - 50) + mu[as.integer(group)] + rnorm(k * n, 0, 5)
dat <- data.frame(group, x, y)

fit_tidy <- ancova(dat, outcome = "y", treatment = "group", covariates = "x")
fit_tidy
```

| term                         | value     |
|:-----------------------------|:----------|
| F_value                      | 9.74      |
| df_1                         | 3         |
| df_2                         | 95        |
| p_value                      | \< 0.0001 |
| sum_of_squares_type          | 3         |
| eta_squared_partial          | 0.235     |
| eta_squared_partial_lower    | 0.0843    |
| eta_squared_partial_upper    | 0.35      |
| omega_squared_partial        | 0.208     |
| omega_squared_partial_lower  | 0.0843    |
| omega_squared_partial_upper  | 0.35      |
| adjusted_mean\[control\]     | -0.695    |
| adjusted_mean\[low\]         | -0.521    |
| adjusted_mean\[medium\]      | 3.52      |
| adjusted_mean\[high\]        | 4.52      |
| se_adjusted_mean\[control\]  | 0.86      |
| se_adjusted_mean\[low\]      | 0.862     |
| se_adjusted_mean\[medium\]   | 0.863     |
| se_adjusted_mean\[high\]     | 0.859     |
| F_homogeneity_of_regression  | 0.621     |
| df_homogeneity_of_regression | 3         |
| p_homogeneity_of_regression  | 0.6030    |

Confidence level: 95%

Pull the adjusted means and the ANCOVA error SD out of the tidy table
(or any fitted model) and pass them to
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md).
With one covariate and `N - k - p` error degrees of freedom, the default
`df` is correct, so it need not be supplied.

``` r

adj <- fit_tidy$value[grepl("^adjusted_mean", fit_tidy$term)]
names(adj) <- levels(group)

# ANCOVA error SD from the fitted model.
s_yx <- summary(lm(y ~ group + x, data = dat))$sigma

ci_c_ancova_bp(adj_means = adj, s_ancova = s_yx, n = n, num_covariates = 1)
```

| contrast          | estimate | lower_limit | upper_limit |
|:------------------|:---------|:------------|:------------|
| group_1 - group_2 | -0.174   | -3.37       | 3.02        |
| group_1 - group_3 | -4.22    | -7.41       | -1.03       |
| group_1 - group_4 | -5.21    | -8.4        | -2.02       |
| group_2 - group_3 | -4.05    | -7.24       | -0.854      |
| group_2 - group_4 | -5.04    | -8.23       | -1.85       |
| group_3 - group_4 | -0.993   | -4.18       | 2.2         |

Confidence level: 95%

The intervals that exclude zero identify the groups whose adjusted means
differ, with familywise 95% protection and, through the Bryant–Paulson
critical value, an accounting for the extra sampling variability the
random covariate introduces.

## 4. Critical Values on Their Own

If you only need the critical value (for a table, a power calculation,
or a hand computation), call
[`qbryant_paulson()`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
directly. It is the ANCOVA analogue of
[`qtukey()`](https://rdrr.io/r/stats/Tukey.html).

``` r

# 95% critical value for k = 5 groups, p = 2 covariates, nu = 40 error df:
qbryant_paulson(0.95, num_covariates = 2, num_groups = 5, df = 40)
#> [1] 4.145129

# How much larger than the (incorrect) Tukey value?
qbryant_paulson(0.95, 2, 5, 40) / qtukey(0.95, nmeans = 5, df = 40)
#> [1] 1.026245
```

## References

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey’s method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*, 631–638.

Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures
in the analysis of covariance. *Journal of the American Statistical
Association, 75*(372), 874–880.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.
