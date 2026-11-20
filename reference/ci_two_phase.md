# Confidence Interval for a Mean Under Two-Phase Sampling With a Proxy Measurement

The confidence interval for the population mean of a target measurement
that is costly to obtain and so is observed on only a random subsample
of units, when an inexpensive proxy measurement is observed on every
unit of a larger pool. The design is two-phase, or double, sampling
(Neyman, 1938): the first phase draws the pool and measures the proxy on
all of it, and the second draws a simple random subsample of the pool,
the labeled units, and measures the target on them. The regression
estimator of the target mean, the row `mean_Y` of the table,
\$\$\hat{\mu}\_Y = \bar y_n + \hat b\\(\bar x_M - \bar x_n),\$\$ adjusts
the labeled units' target mean \\\bar y_n\\ by how far their proxy mean
\\\bar x_n\\ sits from the proxy mean of the whole pool \\\bar x_M\\,
through the least squares slope \\\hat b\\ of the target on the proxy in
the labeled units. When the two measurements are strongly correlated,
the interval is much narrower than the one the labeled units give alone,
which is the point of the design.

## Usage

``` r
ci_two_phase(
  y,
  x,
  x_unlabeled,
  conf_level = 0.95,
  method = c("regression", "prediction_powered"),
  cluster = NULL,
  cluster_unlabeled = NULL
)
```

## Arguments

- y:

  The target measurement on the \\n\\ labeled units (numeric, or logical
  for a binary measurement, taken as 0 and 1).

- x:

  The proxy measurement on the same \\n\\ units, aligned with `y`.

- x_unlabeled:

  The proxy measurement on the units of the pool that have no target
  value, \\N = M - n\\ of them, where \\M\\ is the size of the pool.

- conf_level:

  Confidence level (default `0.95`).

- method:

  `"regression"` (default), the regression estimator above; or
  `"prediction_powered"`, the prediction-powered estimate of
  Angelopoulos et al. (2023), which gives the proxy the fixed weight
  \\\lambda = 1\\, the row `proxy_weight` (see Details).

- cluster, cluster_unlabeled:

  Optional cluster identifiers for the labeled and the unlabeled units
  (for example the respondent who contributed several units, or the
  prompt that generated several items), supplied together or not at all.
  When given, the standard error is computed from the cluster totals of
  each unit's contribution to the estimate, and the limits use the *t*
  distribution on \\Q_n - 1\\ degrees of freedom, where \\Q_n\\ is the
  number of clusters among the labeled units (see Details).

## Value

A `data.frame` (class `dmar_tbl`), its `method` attribute naming the
method, with rows

- mean_Y:

  the estimate of \\\mu_Y\\, the population mean of the target.

- se, df:

  its standard error and the degrees of freedom of the *t* reference
  distribution.

- ci_lower, ci_upper:

  the confidence limits.

- slope:

  the least squares slope \\\hat b\\ of the target on the proxy in the
  labeled units. The regression estimate applies it; the
  prediction-powered estimate does not, and the row is reported for
  comparison (`NA` when the proxy is constant on the labeled units).

- proxy_weight:

  the weight \\\lambda\\ the estimate gives the proxy mean of the
  unlabeled units: \\\hat b N / M\\ for the regression estimate and 1
  for the prediction-powered estimate.

- correlation:

  the correlation of target and proxy in the labeled units.

- proxy_bias:

  the proxy mean minus the target mean in the labeled units, the proxy's
  bias when the two share a scale. It does not enter the regression
  estimate.

- se_ratio_target_only:

  the standard error divided by that of the labeled units' target mean
  alone (\\s_Y/\sqrt n\\, or its cluster counterpart), so values below 1
  are the precision the proxy buys.

- n_labeled, n_pool:

  the numbers of labeled units and of units in the pool, \\n\\ and
  \\M\\.

## Details

**The estimand.** The interval is for \\\mu_Y\\, the mean of the target
measurement in the population from which the pool was drawn at random,
the labeled units being a simple random subsample of the pool. The proxy
enters only as a covariate. It need not be on the target's scale,
unbiased for it, or even measured in its units, and whatever its bias,
that bias does not enter the estimate; what the proxy contributes is its
correlation \\\rho\\ with the target. In the evaluation of a system
whose outputs an automated judge labels, for example, the target is a
human assessor's label of an output, the proxy is the judge's label of
the same output, and the pool is the set of outputs the judge labels.

**Standard error.** To first order, \$\$\mathrm{Var}(\hat{\mu}\_Y) =
\frac{\sigma_Y^2 (1 - \rho^2)}{n} + \frac{\rho^2 \sigma_Y^2}{M},\$\$ the
first term the residual variation the regression leaves among the \\n\\
labeled units and the second the sampling variation of the pool itself.
The reported standard error estimates \\\sigma_e^2/n +
\beta^2\sigma_X^2/M\\, with the residual mean square of the target on
the proxy in the labeled units (on \\n - 2\\ degrees of freedom) for
\\\sigma_e^2\\ and \\\hat b^2\\ times the proxy's variance over the
whole pool for \\\beta^2 \sigma_X^2\\, and the limits use the *t*
distribution on the \\n - 2\\ degrees of freedom of the first term,
which dominates when the pool is large. If the pool is itself the whole
population of interest, a fixed collection of units every one of which
carries the proxy, the estimand is the pool's own target mean and the
reported standard error is conservative for it: the second term is then
not needed, and the first can be multiplied by \\1 - n/M\\.

**Relation to prediction-powered inference.** Writing \\\bar x_N\\ for
the proxy mean of the unlabeled units, \\\hat{\mu}\_Y = \lambda \bar
x_N + (\bar y_n - \lambda \bar x_n)\\ with \\\lambda = \hat b N / M\\.
With \\\lambda = 1\\ in its place this is the prediction-powered
estimate of Angelopoulos et al. (2023): the proxy mean of the unlabeled
units, corrected by the mean target-minus-proxy difference of the
labeled units. That estimate is unbiased for any proxy, but it uses the
proxy efficiently only when the proxy is on the target's scale with a
slope near one, while the regression estimator learns the weight from
the data; the `proxy_weight` row reports \\\lambda\\, the weight each
method applied. `method = "prediction_powered"` reports the interval
with standard error \\\\s\_{Y - X}^2/n + s\_{X,N}^2/N\\^{1/2}\\ and the
*t* distribution on \\n - 1\\ degrees of freedom; Angelopoulos et al.
state it with the normal quantile, which the *t* interval approaches as
\\n\\ grows.

**Clusters.** With `cluster` and `cluster_unlabeled`, each unit's
contribution to the linearized estimate is summed within its cluster,
and the standard error is \\\\Q_n/(Q_n - 1) \sum_q T_q^2\\^{1/2}\\, the
sum running over every cluster of the pool, \\T_q\\ the total of cluster
\\q\\, and \\Q_n\\ the number of clusters among the labeled units. For
the regression estimate a labeled unit contributes \\e_i/n + \hat
b(x_i - \bar x_M)/M\\ and an unlabeled unit \\\hat b(x_j - \bar
x_M)/M\\, where \\e_i\\ is its residual; for the prediction-powered
estimate a labeled unit contributes its target-minus-proxy difference
less the mean difference, over \\n\\, and an unlabeled unit its proxy
less \\\bar x_N\\, over \\N\\. The term that dominates the variance, the
residual term, comes from the labeled clusters alone, so they, and not
all the clusters of the pool, set both the degrees of freedom, \\Q_n -
1\\, and the small-sample factor \\Q_n/(Q_n - 1)\\. When the labeled
units come from a few clusters this matters: the interval then rests on
those few, however many clusters the pool has. With one unit per cluster
the labeled part of the prediction-powered standard error reduces to
\\s\_{Y - X}^2/n\\, its form without clusters. At least two labeled
clusters are required.

**A constant proxy or target.** The regression estimate needs both
measurements to vary among the labeled units. A proxy that takes one
value there leaves the slope undefined, and a target that takes one
value (a binary target that is 1 on every labeled unit, for example)
makes the residual variation and the slope both zero, so the interval
would have zero width. Each is refused with an error rather than
reported as an interval; more labeled units are the remedy. The
prediction-powered estimate does not rest on the regression and is
defined in both cases.

**Missing values.** A labeled unit enters only when its target, its
proxy, and (if given) its cluster are all present, and an unlabeled unit
only when its proxy (and cluster) is present; the others are set aside,
which leaves the design intact when values are missing completely at
random. The `n_labeled` and `n_pool` rows count the units used.

## References

Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., & Zrnic,
T. (2023). Prediction-powered inference. *Science, 382*(6671), 669–674.
[doi:10.1126/science.adi6000](https://doi.org/10.1126/science.adi6000)

Neyman, J. (1938). Contribution to the theory of sampling human
populations. *Journal of the American Statistical Association, 33*(201),
101–116.
[doi:10.1080/01621459.1938.10503378](https://doi.org/10.1080/01621459.1938.10503378)

## See also

[`ss_aipe_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md)
and
[`ss_aipe_two_phase_cost`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md)
to plan the numbers of labeled and pool units,
[`ss_aipe_two_phase_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)
to evaluate a plan, and
[`ss_seq_two_phase`](https://yelleknek.github.io/DMAR/reference/ss_seq_two_phase.md)
for the sequential rule that needs no planning values.

Other confidence intervals:
[`ci_crossed()`](https://yelleknek.github.io/DMAR/reference/ci_crossed.md),
[`ci_proportion()`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md),
[`ci_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/ci_proportion_paired.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A brief screening scale (the proxy) is scored for all 2,000
# respondents of a survey; a structured interview (the target) is given
# to a random 150 of them.
set.seed(113)
screen    <- rnorm(2000, mean = 50, sd = 10)
interview <- 20 + 0.6 * screen + rnorm(2000, sd = 5)
labeled   <- sample(2000, 150)
ci_two_phase(y = interview[labeled], x = screen[labeled],
             x_unlabeled = screen[-labeled])
#>  term                 value 
#>  mean_Y               51.045
#>  se                   0.405 
#>  df                   148   
#>  ci_lower             50.245
#>  ci_upper             51.845
#>  slope                0.608 
#>  proxy_weight         0.562 
#>  correlation          0.801 
#>  proxy_bias           -0.055
#>  se_ratio_target_only 0.637 
#>  n_labeled            150   
#>  n_pool               2000  
#> 
#> Confidence level: 95%

# The screener is not on the interview's scale: the interview rises
# about 0.6 points per screener point. The prediction-powered estimate
# gives the screener a weight of 1 regardless, so it stays unbiased but
# its interval is wider than the regression estimator's.
ci_two_phase(y = interview[labeled], x = screen[labeled],
             x_unlabeled = screen[-labeled],
             method = "prediction_powered")
#>  term                 value 
#>  mean_Y               50.387
#>  se                   0.552 
#>  df                   149   
#>  ci_lower             49.295
#>  ci_upper             51.478
#>  slope                0.608 
#>  proxy_weight         1     
#>  correlation          0.801 
#>  proxy_bias           -0.055
#>  se_ratio_target_only 0.869 
#>  n_labeled            150   
#>  n_pool               2000  
#> 
#> Confidence level: 95%
```
