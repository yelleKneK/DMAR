# Unbiased and Adjusted Estimators of the Population Squared Multiple Correlation

Estimates the population squared multiple correlation coefficient
\\\rho^2\\ from an observed sample \\R^2\\, correcting the well-known
positive (upward) bias of \\R^2\\. Two estimators are available: the
(essentially) unbiased Olkin and Pratt (1958) estimator (the default),
and the classic Ezekiel (1930) adjusted-\\R^2\\ shrinkage formula
reported by [`summary.lm`](https://rdrr.io/r/stats/summary.lm.html) as
`adj.r.squared`. This is the inverse direction of
[`expected_R2`](https://yelleknek.github.io/DMAR/reference/expected_R2.md),
which gives the forward expectation \\E\[R^2 \mid \rho^2\]\\.

## Usage

``` r
unbiased_R2(R2, N, p, method = c("olkin_pratt", "ezekiel"))
```

## Arguments

- R2:

  Observed sample squared multiple correlation coefficient, in \\\[0,
  1\]\\.

- N:

  Sample size.

- p:

  Number of predictor variables.

- method:

  Which estimator to compute: `"olkin_pratt"` (default) for the
  (essentially) unbiased Olkin-Pratt (1958) estimator, or `"ezekiel"`
  for the Ezekiel (1930) adjusted \\R^2\\ (the value
  [`summary.lm`](https://rdrr.io/r/stats/summary.lm.html) reports as
  `adj.r.squared`).

## Value

A 1-row `data.frame` (class `dmar_tbl`) with columns `term` and `value`.
The `term` is `"unbiased_population_R2"` when `method = "olkin_pratt"`
and `"adjusted_population_R2"` when `method = "ezekiel"`; `value` is the
corresponding estimate of \\\rho^2\\.

## Details

The sample \\R^2\\ overestimates \\\rho^2\\; the bias is larger for
smaller samples and for more predictors. Two corrections are offered.

The Ezekiel (1930) adjusted estimator is
\$\$\hat\rho^2\_{\mathrm{Ezekiel}} = 1 - \frac{N - 1}{N - p - 1}\\(1 -
R^2).\$\$ It reduces the bias but is not unbiased; it is exactly the
quantity `summary(lm(...))$adj.r.squared` reports.

The Olkin and Pratt (1958) estimator is (essentially) unbiased:
\$\$\hat\rho^2\_{\mathrm{OP}} = 1 - \frac{N - 3}{N - p - 1}\\(1 - R^2)\\
{}\_2F_1\\\left(1, 1; \frac{N - p + 1}{2}; 1 - R^2\right),\$\$ where
\\{}\_2F_1\\ is the Gaussian hypergeometric function (the same function
[`expected_R2`](https://yelleknek.github.io/DMAR/reference/expected_R2.md)
uses for the forward direction; see Stuart, Ord, & Arnold, 1999, section
28). Both estimators can fall below 0 for very small \\R^2\\; that is
expected behavior for a bias-corrected estimator and is not truncated
here (matching `adj.r.squared`, which is also allowed to be negative).

## References

Ezekiel, M. (1930). *Methods of correlation analysis*. Wiley.

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain
correlation coefficients. *The Annals of Mathematical Statistics,
29*(1), 201–211.

Stuart, A., Ord, J. K., & Arnold, S. (1999). *Kendall's advanced theory
of statistics, volume 2A: Classical inference and the linear model* (6th
ed.). Arnold.

## See also

[`expected_R2`](https://yelleknek.github.io/DMAR/reference/expected_R2.md)
for the forward expectation \\E\[R^2 \mid \rho^2\]\\, and
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`var_R2`](https://yelleknek.github.io/DMAR/reference/var_R2.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)
for interval, variance, and planning tools on the same effect size.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# An observed R^2 = .50 with N = 50 and p = 5 predictors overstates rho^2.
# The Olkin-Pratt (essentially unbiased) estimate is the default:
unbiased_R2(R2 = .50, N = 50, p = 5)
#>  term                   value
#>  unbiased_population_R2 0.454

# The Ezekiel adjusted R^2 (what summary(lm) reports) over-shrinks slightly,
# so it typically sits a little below the Olkin-Pratt value:
unbiased_R2(R2 = .50, N = 50, p = 5, method = "ezekiel")
#>  term                   value
#>  adjusted_population_R2 0.443

# The Ezekiel option reproduces summary(lm)$adj.r.squared exactly.
set.seed(113)
d   <- as.data.frame(matrix(rnorm(50 * 6), 50, 6))
fit <- lm(V1 ~ ., data = d)
s   <- summary(fit)
unbiased_R2(R2 = s$r.squared, N = 50, p = 5, method = "ezekiel")$value
#> [1] -0.09483852
s$adj.r.squared
#> [1] -0.09483852

# The bias (and so the correction) shrinks as N grows for fixed R^2 and p.
unbiased_R2(.50, 50, 5)
#>  term                   value
#>  unbiased_population_R2 0.454
unbiased_R2(.50, 500, 5)
#>  term                   value
#>  unbiased_population_R2 0.496
```
