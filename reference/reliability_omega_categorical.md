# Categorical Omega for Ordered-Categorical Items, With a Confidence Interval

Estimates categorical omega (\\\omega_C\\; Green & Yang, 2009; Kelley &
Pornprasertmanit, 2016) for a homogeneous composite of
ordered-categorical items and returns a bootstrap confidence interval.

## Usage

``` r
reliability_omega_categorical(
  data,
  ci_method = c("bca", "percentile", "bootstrap_se", "bootstrap_se_logistic", "none"),
  conf_level = 0.95,
  B = 10000,
  seed = NULL
)

reliability_omega_c(
  data,
  ci_method = c("bca", "percentile", "bootstrap_se", "bootstrap_se_logistic", "none"),
  conf_level = 0.95,
  B = 10000,
  seed = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of ordered-categorical item scores
  (integer codes for the categories). Rows with any missing values are
  listwise-deleted.

- ci_method:

  Method for constructing the confidence interval. See *Details*. When
  not supplied, no interval is computed: every interval for categorical
  omega is bootstrap based, and a bootstrap is never run unless
  requested. Ask for `"bca"` (the recommended method) or `"percentile"`
  to obtain one.

- conf_level:

  Confidence level for the interval. Defaults to `0.95`.

- B:

  Number of bootstrap replications. Defaults to `10000`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, which leaves the user's current RNG state intact; supply an
  integer for reproducibility.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample \\\omega_C\\), `"se"` (bootstrap standard deviation across
replications, `NA` for `ci_method = "none"`), `"lower_limit"` and
`"upper_limit"` (clamped to \[0, 1\]), `"conf_level"`, `"N"`,
`"N_complete"` (the complete cases; equal to `"N"` here, and carried so
the whole reliability family returns one shape), and `"J"`. Attributes
`coefficient` (`"omega_categorical"`), `ci_method`, and `B` record the
computation.

## Details

Categorical omega is designed for items measured on an ordered
categorical scale (e.g., Likert items). It uses a probit-link
single-factor model in which each observed item \\X_j\\ is modeled as a
categorization of an underlying continuous response variable
\\X_j^{\*}\\ via thresholds \\t\_{j,c}\\ (Muthén, 1984; Millsap &
Yun-Tein, 2004), fit by diagonally weighted least squares with mean- and
variance-adjusted test statistic (WLSMV), the standard estimator for
ordered categorical items. With the delta parameterization
(\\Var(X_j^{\*}) = 1\\), the population categorical omega is
\$\$\omega_C = \frac{\sum\_{j=1}^{J} \sum\_{j'=1}^{J}
\sigma\_{jj'}\\\left(\lambda\_{j}\lambda\_{j'}\right)}{\sum\_{j=1}^{J}
\sum\_{j'=1}^{J} \sigma\_{jj'}\\\left(\rho\_{X\_{j}^{\*}
X\_{j'}^{\*}}\right)},\$\$ where \\\sigma\_{jj'}(r)\\ is the model
implied covariance of \\(X\_{j}, X\_{j'})\\ computed from a bivariate
normal CDF over pairs of category thresholds and a correlation \\r\\
(Green & Yang, 2009, Eq. 13–14; Kelley & Pornprasertmanit, 2016, Eq.
17–18). The numerator uses model implied polychoric correlations
(\\\lambda_j \lambda\_{j'}\\), while the denominator uses observed
polychoric correlations estimated from the data via a saturated
bivariate model.

Kelley and Pornprasertmanit (2016) found in extensive Monte Carlo
simulation that the bias-corrected and accelerated (BCa) bootstrap
confidence interval for categorical omega achieved acceptable coverage
across a wide variety of threshold patterns, sample sizes, item counts,
and population reliability values. They specifically recommend BCa for
categorical omega. Because no bootstrap runs in DMAR unless the user
requests one, the default output is the point estimate with a message
naming the call that produces the recommended interval; request
`ci_method = "bca"` to obtain it.

**When to use.** Use `reliability_omega_categorical` when items are
ordered-categorical, especially when (a) the number of categories is
small (e.g., two to five), (b) item distributions are skewed, or (c)
threshold patterns differ markedly across items. In Kelley and
Pornprasertmanit's (2016) Study 3, treating ordered items as continuous
and using
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
with the observed total variance in the denominator achieved acceptable
coverage only when threshold patterns were similar across items; in
their experience that condition is rare in practice.

Available confidence interval methods (set via `ci_method`). Every
interval here is bootstrap based: the rows of `data` are resampled with
replacement `B` times (10000 by default) and categorical omega is
recomputed, with the full WLSMV model refit, on each replication (Efron
& Tibshirani, 1993). Replications whose refit fails or does not
converge, most common with small samples and sparse response categories,
are dropped, and the interval is computed from the replications that
return a value. Bootstrap results vary from run to run; supply `seed`
for reproducibility.

- `"bca"`:

  The bias-corrected and accelerated bootstrap, the default and the
  specific recommendation of Kelley and Pornprasertmanit (2016). Where
  the percentile interval reads its limits directly off the empirical
  quantiles of the bootstrap estimates, BCa adjusts the two quantile
  positions for median bias (estimated from the bootstrap distribution)
  and for the rate at which the estimator's variance changes with the
  parameter (the acceleration, estimated by the jackknife), making it
  second-order accurate where the percentile interval is first-order
  accurate (DiCiccio & Efron, 1996). The adjusted quantile positions sit
  farther into the tails than the percentile interval uses, which is why
  the default `B = 10000` is larger than the customary 2000; reduce `B`
  for exploration, not for a reported analysis.

- `"percentile"`:

  Percentile bootstrap: the interval limits are the empirical quantiles
  of the bootstrap estimates.

- `"bootstrap_se"`, `"bootstrap_se_logistic"`:

  Wald intervals using the bootstrap standard deviation as a standard
  error, built on the logit scale for the `_logistic` variant so the
  endpoints respect \[0, 1\].

- `"none"`:

  Return only the point estimate.

**Comparison with other packages.** The psych package's
[`omega`](https://rdrr.io/pkg/psych/man/omega.html) fits a Schmid-Leiman
hierarchical factor model on continuous (or treated-as-continuous) items
and does not implement categorical omega in the sense of Green and Yang
(2009). For ordered-categorical items `reliability_omega_categorical` is
the appropriate choice;
[`polychoric`](https://rdrr.io/pkg/psych/man/tetrachor.html) provides
polychoric correlation estimation as a separate tool.

## References

DiCiccio, T. J., & Efron, B. (1996). Bootstrap confidence intervals.
*Statistical Science, 11*(3), 189–228.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*, 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance in
ordered-categorical measures. *Multivariate Behavioral Research, 39*(3),
479–515.
[doi:10.1207/s15327906mbr3903_4](https://doi.org/10.1207/s15327906mbr3903_4)

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper),
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
(use for continuous items),
[`reliability_kr20`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md)
(for dichotomous items),
[`omega`](https://rdrr.io/pkg/psych/man/omega.html),
[`polychoric`](https://rdrr.io/pkg/psych/man/tetrachor.html).

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_kr20()`](https://yelleknek.github.io/DMAR/reference/reliability_kr20.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
set.seed(113)
# Six 5-category items with a single latent factor.
N <- 500
J <- 6
loadings <- rep(0.7, J)
eta <- rnorm(N)
latent <- outer(eta, loadings) +
          matrix(rnorm(N * J), N, J) %*% diag(sqrt(1 - loadings^2))
items <- apply(latent, 2, function(x)
  as.integer(cut(x, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
                 labels = FALSE)))
colnames(items) <- paste0("y", seq_len(J))

# Default: point estimate only, with a message naming the call that
# produces the recommended interval.
reliability_omega_categorical(data = items)
#> Categorical omega is reported without a confidence interval by default because its interval is bootstrap based. Request it with ci_method = "bca" (the recommended method) or "percentile"; B = 10000 replications is the default when you do.
#>  term        value
#>  estimate    0.83 
#>  se          <NA> 
#>  lower_limit <NA> 
#>  upper_limit <NA> 
#>  conf_level  0.95 
#>  N           500  
#>  N_complete  500  
#>  J           6    

# The recommended BCa bootstrap, requested explicitly (reduced B for
# a fast example).
reliability_omega_categorical(data = items, ci_method = "bca", B = 200)
#>  term        value 
#>  estimate    0.83  
#>  se          0.0117
#>  lower_limit 0.8   
#>  upper_limit 0.85  
#>  conf_level  0.95  
#>  N           500   
#>  N_complete  500   
#>  J           6     
# }
```
