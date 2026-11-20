# Confidence Interval for the Ratio of Two Independent Means

Computes a confidence set for the ratio of two population means
estimated from independent groups, \\\rho = \mu_1 / \mu_2\\, or, with a
third, baseline group, for the ratio of the two groups' advantages over
the baseline, \\\rho = (\mu_1 - \mu_b) / (\mu_2 - \mu_b)\\. A ratio
states an effect in relative terms (a new procedure's mean is 90 percent
of the standard procedure's), which is the natural scale when the
outcome has a meaningful zero. The second form is the fraction of a
reference treatment's advantage over a placebo or a trivial baseline
that a new treatment retains, the retention of effect of the three-arm
noninferiority trial (Pigeot et al., 2003). The default interval is
Fieller's, which inverts the *t* test of \\\mu_1 - \rho \mu_2 = 0\\ and
is reported as what it is: a bounded interval, the complement of an
interval, or the whole line. The log-scale interval of the response
ratio and percentile and BCa bootstrap intervals are offered beside it.
For two means measured on the same units, see
[`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md).

## Usage

``` r
ci_ratio_of_means(
  x = NULL,
  y = NULL,
  baseline = NULL,
  mean_1 = NULL,
  sd_1 = NULL,
  n_1 = NULL,
  mean_2 = NULL,
  sd_2 = NULL,
  n_2 = NULL,
  mean_baseline = NULL,
  sd_baseline = NULL,
  n_baseline = NULL,
  var_equal = FALSE,
  ci_method = c("fieller", "log", "percentile", "bca"),
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  B = 10000L,
  seed = NULL
)
```

## Arguments

- x, y:

  Numeric vectors of observations from the numerator group and the
  denominator group. The groups are independent and need not be the same
  size; `NA`s are removed from each vector separately.

- baseline:

  An optional numeric vector of observations from a third, independent
  baseline group (a placebo arm, a trivial comparison condition). When
  it is supplied the ratio is that of the advantages, \\(\mu_1 - \mu_b)
  / (\mu_2 - \mu_b)\\.

- mean_1, mean_2:

  Group sample means, for analysis from summary statistics in place of
  `x` and `y`.

- sd_1, sd_2:

  Group sample standard deviations (with \\n - 1\\ in the denominator of
  the variance).

- n_1, n_2:

  Group sample sizes.

- mean_baseline, sd_baseline, n_baseline:

  The baseline group's sample mean, standard deviation, and size, for
  the ratio of the two groups' advantages over it from summary
  statistics (the three arms a noninferiority trial reports). Supply all
  three or none.

- var_equal:

  Logical. `FALSE` (the default) lets the group variances differ and
  inverts the Welch-type test, whose degrees of freedom change with the
  hypothesized ratio; `TRUE` pools the variances, which gives Fieller's
  exact interval under normality and equal variances. Used by
  `"fieller"` and `"log"`; the bootstrap resamples each group from
  itself and needs neither.

- ci_method:

  The interval: `"fieller"` (the default), `"log"` (the response ratio
  interval, formed on the log scale and defined only when both means, or
  both advantages, are positive), `"percentile"`, or `"bca"` (bootstrap
  intervals, which need the raw data).

- conf_level:

  Confidence level (default `0.95`).

- alpha_lower:

  The tail area at the lower limit. An omitted confidence level yields
  to supplied tail areas; an explicit level must agree with them. Each
  tail area is at most .5.

- alpha_upper:

  The tail area at the upper limit, supplied together with
  `alpha_lower`. The table records the coverage,
  `1 - alpha_lower - alpha_upper`, as its level and, when the two tail
  areas differ, the areas themselves, which the line beneath the printed
  table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound,
  unbounded on that side.

- B:

  Number of bootstrap replications for `"percentile"` and `"bca"`, a
  single whole number (ignored otherwise). Default `10000`: a
  replication of a ratio of means costs little, so a reported interval
  deserves the default, and the BCa adjustment moves the working
  quantiles into the tails, where a small `B` is least stable. The
  fewest accepted with 0.025 in each tail is 400 for the percentile
  interval and 2000 for BCa (200 and 1000 at a 0.90 level, 2000 and
  10000 at 0.99); with unequal tail areas the smaller one decides. A
  smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`). See *Bootstrap replications* under *Details*
  and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- seed:

  Optional integer seed for the bootstrap. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  for the call and the prior state restored on exit.

## Value

A long `data.frame` (classes `dmar_ci_ratio_of_means`, `dmar_ci_long`,
`dmar_tbl`) with columns `term` and `value`, whose rows are `ci_lower`,
`ratio`, and `ci_upper` (the set's outer ends, `-Inf` or `Inf` when it
is unbounded on that side); `excluded_lower` and `excluded_upper` (the
ends of the excluded gap of a complement or a union of two intervals,
`NA` otherwise); `mean_1`, `mean_2`, and `mean_baseline` (the group
means; `mean_baseline` is `NA` without a baseline); `se_1`, `se_2`, and
`se_baseline` (the standard errors of the means, from the pooled
variance when `var_equal = TRUE`); `t_denominator` (the denominator,
\\\bar x_2\\ or \\\bar x_2 - \bar x_b\\, over its standard error); `df`
(the degrees of freedom of that *t*: pooled, or Satterthwaite's for the
denominator when `var_equal = FALSE`); and `n_1`, `n_2`, and
`n_baseline`. The rows are the same for every call. Attributes:
`confidence_set` (`"interval"`, `"complement of an interval"`,
`"whole line"`, or, with unequal tail areas only,
`"union of two intervals"`; `NA` when a bootstrap interval was withheld
and both limits are `NA`), `ci_method`, `var_equal`, `conf_level`, and
the tail areas when they differ; `df_limits` when `var_equal = FALSE`
and the interval is Fieller's or the log ratio's, the degrees of freedom
of the test at each end of the set, a vector with elements `lower`,
`upper`, `excluded_lower`, and `excluded_upper` that is `NA` at an end
no test decides (an infinite limit, the open end of a one-sided bound,
or a gap the set does not have); `B`, `B_used`, `n_dropped`,
`B_required`, and `bootstrap_status` for a bootstrap interval (see
*Bootstrap replications* under Details); and `fieller_set` (the form of
the Fieller set) when the interval is a log or bootstrap interval and
the Fieller set at the same level is unbounded. `tidy()` returns the
ratio, its limits, the ends of any excluded gap, and the form of the set
in one row (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

**Fieller's confidence set.** For a candidate ratio \\\rho\\ the
contrast \\\bar x_1 - \rho \bar x_2\\ (with a baseline, \\(\bar x_1 -
\bar x_b) - \rho (\bar x_2 - \bar x_b)\\) estimates zero exactly when
\\\rho\\ is the population ratio, and its *t* statistic, \$\$T(\rho) =
\frac{\hat\theta_1 - \rho \hat\theta_2} {\sqrt{v\_{11} - 2 \rho
v\_{12} + \rho^2 v\_{22}}},\$\$ where \\\hat\theta_1\\ and
\\\hat\theta_2\\ are the numerator and denominator and the \\v\\ terms
their estimated variances and covariance, is a pivot: it involves the
data and \\\rho\\ alone (Fieller, 1941, 1954). The 1941 reference is
Fieller's study of the biological standardization of insulin, which
Fieller (1954) cites with the year 1940. The confidence set is every
\\\rho\\ that the test does not reject, the solution of a quadratic
inequality in \\\rho\\. Its form is decided by the denominator's *t*,
\\\hat\theta_2 / \sqrt{v\_{22}}\\, reported in the `t_denominator` row
(Fieller, 1954, p. 177). The row keeps the sign of the denominator, and
only its absolute value matters: when the absolute value of the
denominator's *t* exceeds the critical value (with unequal tail areas,
each side of the estimate is bounded when that absolute value exceeds
the side's own critical value) the set is a bounded interval. When it
does not, the set contains arbitrarily large values of the ratio, as it
should when the data cannot distinguish the denominator from zero: it is
then either the complement of an interval, every value outside a gap,
or, when the data are consistent with every ratio at the chosen level,
the whole line. The table reports all three as they are. The `ci_lower`
and `ci_upper` rows are the outer ends of the set (`-Inf` and `Inf` when
it is unbounded), and the `excluded_lower` and `excluded_upper` rows are
the ends of the gap of a complement (`NA` otherwise): the set is every
value from `ci_lower` to `ci_upper` except those strictly between
`excluded_lower` and `excluded_upper`. The `confidence_set` attribute
names the form, and so does the line printed under the table;
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
writes the set in that form, for example as every value except those
between the ends of the gap. With unequal tail areas one side of the set
can be bounded while the other is not, which can also give an interval
unbounded on one side or the union of two intervals. von Luxburg and
Franz (2009) derive the same sets geometrically, from confidence
intervals for the means of linear combinations of the two variables.

**Equal or unequal variances.** With `var_equal = TRUE` the variances
are pooled over the groups, on \\n_1 + n_2 - 2\\ degrees of freedom
(\\n_1 + n_2 + n_b - 3\\ with a baseline), and the set is Fieller's
exact one under normality and equal variances. With the default
`var_equal = FALSE` each group keeps its own variance and the contrast
is tested as Welch's *t* tests a difference, with Satterthwaite degrees
of freedom; those degrees of freedom depend on \\\rho\\, so the set has
no closed form and is found by solving for the values of \\\rho\\ at
which the test's *p*-value equals its tail area. The interval therefore
agrees with the test: without a baseline and with equal tail areas, each
finite end \\L\\ of the set (a limit or an end of the gap) is a value at
which `t.test(x, L * y)` has *p*-value `1 - conf_level`; with a baseline
the test is that of the contrast \\\bar x_1 - L \bar x_2 - (1 - L) \bar
x_b\\, with the Satterthwaite degrees of freedom of its three variances.
The degrees of freedom of the test at each finite end of the set are in
the `df_limits` attribute, and the `df` row holds those of the
denominator's *t*. This inversion is the one the mratios function
`ttestratio()` performs with `iterativeCI = TRUE`, which its
documentation describes as searching iteratively for the limits by
inverting the test; its default instead puts degrees of freedom computed
at the estimate into Fieller's formula, and its documentation notes that
those limits need not agree with the test. Hasler, Vonk, and Hothorn
(2008) treat the three-arm design when the group variances differ.

**The log ratio.** `ci_method = "log"` forms the interval for \\\log
\rho\\ from the delta method standard error,
\\\sqrt{v\_{11}/\hat\theta_1^2 + v\_{22}/\hat\theta_2^2 - 2
v\_{12}/(\hat\theta_1 \hat\theta_2)}\\, and maps it back by the
exponential. This is the response ratio of meta-analysis (Hedges,
Gurevitch, & Curtis, 1999), offered so a ratio reported here can be
compared with one entered in a meta-analysis. It is defined only when
both terms of the ratio are positive, and is refused otherwise with a
condition of class `dmar_log_ratio_undefined` (which inherits
`dmar_refusal`). It is always bounded; when the denominator is well
determined it is close to Fieller's interval, and the two part company
as the denominator's *t* falls. With unequal variances its degrees of
freedom are the Satterthwaite degrees of freedom of the linearized log
ratio.

**The bootstrap.** `"percentile"` and `"bca"` resample each group from
itself (the baseline group too), as boot, the software of Davison and
Hinkley (1997), does for several samples given as strata, and compute
the ratio in every replication. The BCa interval reads its limits at
quantile positions moved from those of the percentile interval by a bias
correction, the median bias of the replications in normal units, and by
an acceleration, the rate at which the standard error of the estimate
changes with the parameter on a normalized scale (Efron & Tibshirani,
1993); its acceleration is estimated by the jackknife, deleting one
observation at a time within its group, with each group's influence
values scaled by its size. A bootstrap interval for a ratio is always
bounded, whatever the denominator's *t*, so it can undercover badly when
the denominator is near zero, exactly where Fieller's set becomes
unbounded. When the Fieller set at the same level is unbounded, the line
under the table says so. Pigeot et al. (2003) compare the percentile
interval with Fieller's in the three-arm design. Summary statistics
cannot be resampled, so a bootstrap from them is refused with a
condition of class `dmar_bootstrap_needs_data`.

**Bootstrap replications.** The limits are quantiles of the `B`
replicated ratios. A limit that cuts off a tail area \\\alpha\\ is
located by the replications beyond it, so the percentile interval is run
only when `B` puts at least 10 replications beyond each limit asked for
(400 with 0.025 in each tail; with unequal tail areas the smaller one
decides, and a one-sided bound needs 10 beyond its one limit), and the
BCa interval only with 50 beyond each nominal limit (2000 with 0.025 in
each tail). A replication in which the ratio is undefined (a denominator
of exactly zero, which discrete data can produce) is dropped and
counted, never replaced, and the interval is computed from the rest. The
result records `B` (requested), `B_used` (the number the interval rests
on), `n_dropped`, and `B_required` (the fewest replications a request
may ask for; that count is held to `B`, not to `B_used`) as attributes,
and the line under the table states them, with `bootstrap_status`
(`"complete"`, `"dropped"`, `"exploratory"`, `"insufficient_count"`, or
`"insufficient_share"`). One warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the 2.5 percent of `B` that lies beyond one limit at the 95
percent level, because they could then decide by themselves where that
limit falls. When more than half were dropped, the limits asked for are
`NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`. For the BCa interval the count of 10 is
applied again where the bias correction and the acceleration place the
limits, and when fewer than 10 replications lie beyond one of them no
BCa interval is reported (limits `NA`, the same class, with the two
positions as `positions` and a suggested `B` as `B_suggest`). A BCa
interval whose bias correction is infinite (every replicate on one side
of the estimate) or whose jackknife is undefined (deleting one
observation leaves a denominator of zero) stops with an error that names
`ci_method = "percentile"`. `options(dmar.bootstrap_exploratory = TRUE)`
lets a smaller `B` (2 or more) run for a first look or a demonstration;
every such call gives a warning of class `dmar_bootstrap_exploratory`
and the result is not one to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

**Reading the interval for noninferiority.** With a baseline, a new
treatment is shown to retain more than a fraction \\f\\ of the reference
treatment's advantage, at one-sided level `alpha_lower`, exactly when
\\f\\ is below the lower limit of a set that is a bounded interval (or,
for a one-sided bound, an interval unbounded only above) and whose
denominator is positive. The one-sided bound is
`alpha_lower = 0.05, alpha_upper = 0`, and a two-sided 95 percent
interval gives the bound at one-sided level .025. Reporting the limit
lets each reader apply his or her own fraction; the function itself
tests no fraction and states no verdict. When the denominator is
negative, or the set is unbounded below, the lower limit does not carry
that reading.

## References

Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap methods and their
application*. Cambridge University Press.
[doi:10.1017/CBO9780511802843](https://doi.org/10.1017/CBO9780511802843)

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Fieller, E. C. (1941). The biological standardization of insulin.
*Supplement to the Journal of the Royal Statistical Society, 7*(1),
1–54. [doi:10.2307/2983630](https://doi.org/10.2307/2983630)

Fieller, E. C. (1954). Some problems in interval estimation. *Journal of
the Royal Statistical Society, Series B, 16*(2), 175–185.
[doi:10.1111/j.2517-6161.1954.tb00159.x](https://doi.org/10.1111/j.2517-6161.1954.tb00159.x)

Hasler, M., Vonk, R., & Hothorn, L. A. (2008). Assessing non-inferiority
of a new treatment in a three-arm trial in the presence of
heteroscedasticity. *Statistics in Medicine, 27*(4), 490–503.
[doi:10.1002/sim.3052](https://doi.org/10.1002/sim.3052)

Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis
of response ratios in experimental ecology. *Ecology, 80*(4), 1150–1156.
[doi:10.1890/0012-9658(1999)080\[1150:TMAORR\]2.0.CO;2](https://doi.org/10.1890/0012-9658%281999%29080%5B1150%3ATMAORR%5D2.0.CO%3B2)

Pigeot, I., Schäfer, J., Röhmel, J., & Hauschke, D. (2003). Assessing
non-inferiority of a new treatment in a three-arm clinical trial
including a placebo. *Statistics in Medicine, 22*(6), 883–899.
[doi:10.1002/sim.1450](https://doi.org/10.1002/sim.1450)

von Luxburg, U., & Franz, V. H. (2009). A geometric approach to
confidence sets for ratios: Fieller's theorem, generalizations and
bootstrap. *Statistica Sinica, 19*(3), 1095–1117.

## See also

[`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md)
for two means measured on the same units;
[`summary_t_test`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md)
and [`welch_t`](https://yelleknek.github.io/DMAR/reference/welch_t.md)
for the difference between two means.

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_ratio_of_means_paired()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Drinks per week among the homeless alcohol-dependent adults of
# drinks_trial: the mean under the community reinforcement approach (CRA)
# relative to the mean under standard treatment. The group variances are
# not assumed equal.
drinks <- drinks_trial$drinks_per_week
cra <- drinks[drinks_trial$treatment == "CRA"]
standard <- drinks[drinks_trial$treatment == "Standard"]
ci_ratio_of_means(cra, standard)
#>  term           value 
#>  ci_lower       0.116 
#>  ratio          0.521 
#>  ci_upper       1.865 
#>  excluded_lower <NA>  
#>  excluded_upper <NA>  
#>  mean_1         26.636
#>  mean_2         51.110
#>  mean_baseline  <NA>  
#>  se_1           9.982 
#>  se_2           17.456
#>  se_baseline    <NA>  
#>  t_denominator  2.928 
#>  df             36    
#>  n_1            32    
#>  n_2            37    
#>  n_baseline     <NA>  
#> 
#> Confidence level: 95%
#> 
#> Fieller confidence set: a bounded interval.

# The drinking counts are strongly skewed, so the response ratio interval
# and a BCa bootstrap interval are worth setting beside Fieller's;
# B = 2000 keeps the example quick, and a reported interval deserves the
# default B = 10000.
ci_ratio_of_means(cra, standard, ci_method = "log")
#>  term           value 
#>  ci_lower       0.189 
#>  ratio          0.521 
#>  ci_upper       1.435 
#>  excluded_lower <NA>  
#>  excluded_upper <NA>  
#>  mean_1         26.636
#>  mean_2         51.110
#>  mean_baseline  <NA>  
#>  se_1           9.982 
#>  se_2           17.456
#>  se_baseline    <NA>  
#>  t_denominator  2.928 
#>  df             36    
#>  n_1            32    
#>  n_2            37    
#>  n_baseline     <NA>  
#> 
#> Confidence level: 95%
#> 
#> Log ratio interval: a bounded interval.
ci_ratio_of_means(cra, standard, ci_method = "bca", B = 2000, seed = 113)
#>  term           value 
#>  ci_lower       0.156 
#>  ratio          0.521 
#>  ci_upper       1.276 
#>  excluded_lower <NA>  
#>  excluded_upper <NA>  
#>  mean_1         26.636
#>  mean_2         51.110
#>  mean_baseline  <NA>  
#>  se_1           9.982 
#>  se_2           17.456
#>  se_baseline    <NA>  
#>  t_denominator  2.928 
#>  df             36    
#>  n_1            32    
#>  n_2            37    
#>  n_baseline     <NA>  
#> 
#> Confidence level: 95%
#> 
#> BCa bootstrap interval: a bounded interval.
#> 
#> Bootstrap: 2000 replications, 2000 used (0 dropped; at least 2000 must be requested)

# From the summary statistics a paper reports, with pooled variances.
ci_ratio_of_means(mean_1 = 104, sd_1 = 14, n_1 = 40,
                  mean_2 = 96, sd_2 = 15, n_2 = 42, var_equal = TRUE)
#>  term           value  
#>  ci_lower       1.016  
#>  ratio          1.083  
#>  ci_upper       1.155  
#>  excluded_lower <NA>   
#>  excluded_upper <NA>   
#>  mean_1         104.000
#>  mean_2         96.000 
#>  mean_baseline  <NA>   
#>  se_1           2.296  
#>  se_2           2.241  
#>  se_baseline    <NA>   
#>  t_denominator  42.845 
#>  df             80     
#>  n_1            40     
#>  n_2            42     
#>  n_baseline     <NA>   
#> 
#> Confidence level: 95%
#> 
#> Fieller confidence set: a bounded interval.

# A three-arm trial reported as summary statistics: the fraction of the
# reference arm's advantage over placebo that the new arm retains.
ci_ratio_of_means(mean_1 = 60.1, sd_1 = 12.2, n_1 = 60,
                  mean_2 = 62.4, sd_2 = 10.1, n_2 = 60,
                  mean_baseline = 50.3, sd_baseline = 9.8,
                  n_baseline = 60)
#>  term           value  
#>  ci_lower       0.512  
#>  ratio          0.810  
#>  ci_upper       1.171  
#>  excluded_lower <NA>   
#>  excluded_upper <NA>   
#>  mean_1         60.100 
#>  mean_2         62.400 
#>  mean_baseline  50.300 
#>  se_1           1.575  
#>  se_2           1.304  
#>  se_baseline    1.265  
#>  t_denominator  6.660  
#>  df             117.893
#>  n_1            60     
#>  n_2            60     
#>  n_baseline     60     
#> 
#> Confidence level: 95%
#> 
#> Fieller confidence set: a bounded interval.

# A three-arm trial: the fraction of the reference treatment's advantage
# over placebo that a new treatment retains, with a one-sided 95 percent
# lower bound. Noninferiority at a retained fraction f holds for every f
# below the lower limit.
set.seed(113)
placebo <- rnorm(60, mean = 50, sd = 10)
reference <- rnorm(60, mean = 62, sd = 10)
new_drug <- rnorm(60, mean = 60, sd = 12)
ci_ratio_of_means(new_drug, reference, baseline = placebo,
                  alpha_lower = 0.05, alpha_upper = 0)
#>  term           value  
#>  ci_lower       0.357  
#>  ratio          0.678  
#>  ci_upper       Inf    
#>  excluded_lower <NA>   
#>  excluded_upper <NA>   
#>  mean_1         57.699 
#>  mean_2         60.465 
#>  mean_baseline  51.864 
#>  se_1           1.353  
#>  se_2           1.359  
#>  se_baseline    1.312  
#>  t_denominator  4.553  
#>  df             117.855
#>  n_1            60     
#>  n_2            60     
#>  n_baseline     60     
#> 
#> Confidence level: 95% (one-sided: the lower limit is the 95% lower confidence bound, with tail area 0.05 below it and none above the upper limit)
#> 
#> Fieller confidence set: an interval unbounded above.
```
