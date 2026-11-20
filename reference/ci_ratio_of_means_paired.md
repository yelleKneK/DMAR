# Confidence Interval for the Ratio of Two Means From Paired Observations

Computes a confidence set for the ratio of two population means, \\\rho
= \mu_x / \mu_y\\, when both variables are measured on the same units
(the same participants before and after treatment, two instruments
applied to the same cases, two procedures scored on the same items), or,
with a third aligned variable, for the ratio of the two variables'
advantages over it, \\\rho = (\mu_x - \mu_z) / (\mu_y - \mu_z)\\. The
second form is the fraction of a reference procedure's advantage over a
trivial baseline that a candidate procedure retains, the retention of
effect of three-arm noninferiority designs (Pigeot et al., 2003) when
all three are evaluated on the same units. The default interval is
Fieller's, which inverts the one-sample *t* test of the per-unit
contrast \\x - \rho y\\ and is reported as what it is: a bounded
interval, the complement of an interval, or the whole line. When the
units are grouped (several items from each of a set of topics, several
patients from each clinic), `cluster` gives the cluster-robust version.
The log-scale interval of the response ratio and percentile and BCa
bootstrap intervals are offered beside it. For independent groups, see
[`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md).

## Usage

``` r
ci_ratio_of_means_paired(
  x = NULL,
  y = NULL,
  baseline = NULL,
  cluster = NULL,
  mean_1 = NULL,
  sd_1 = NULL,
  mean_2 = NULL,
  sd_2 = NULL,
  r = NULL,
  N = NULL,
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

  Numeric vectors of the numerator and denominator variables, aligned by
  unit. A unit with a missing value on any of `x`, `y`, `baseline`, or
  `cluster` is dropped.

- baseline:

  An optional numeric vector of a third variable aligned with `x` and
  `y` (a trivial or placebo condition measured on the same units). When
  it is supplied the ratio is that of the per-unit advantages, \\x - z\\
  over \\y - z\\.

- cluster:

  An optional vector (numeric, character, or factor) naming the cluster
  of each unit. When it is supplied the covariance matrix of the two
  means is the cluster-robust one, formed from the cluster sums, on
  \\Q - 1\\ degrees of freedom for \\Q\\ clusters, and the bootstrap
  resamples whole clusters.

- mean_1, mean_2:

  Means of the numerator and denominator variables, for analysis from
  summary statistics in place of `x` and `y`. With a baseline, supply
  the means, standard deviations, and correlation of the two advantages
  (\\x - z\\ and \\y - z\\).

- sd_1, sd_2:

  Standard deviations of the two variables (with \\N - 1\\ in the
  denominator of the variance).

- r:

  The correlation between the two variables over the units.

- N:

  The number of units, the total sample size of the paired design.

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
`NA` otherwise); `mean_1`, `mean_2`, and `mean_baseline` (the means of
`x`, `y`, and `baseline`; `mean_baseline` is `NA` without a baseline);
`se_numerator` and `se_denominator` (the standard errors of the
numerator and denominator, the means of `x` and `y` or of the two
advantages); `r` (the correlation between the numerator and denominator
estimates, which without clusters is the correlation of the two
variables over the units); `t_denominator` (the denominator over its
standard error); `df`; `N` (the number of units); and `n_clusters` (`NA`
without `cluster`). The rows are the same for every call. Attributes:
`confidence_set` (`"interval"`, `"complement of an interval"`,
`"whole line"`, or, with unequal tail areas only,
`"union of two intervals"`; `NA` when a bootstrap interval was withheld
and both limits are `NA`), `ci_method`, `conf_level`, and the tail areas
when they differ; `B`, `B_used`, `n_dropped`, `B_required`, and
`bootstrap_status` for a bootstrap interval (see *Bootstrap
replications* under Details); and `fieller_set` (the form of the Fieller
set) when the interval is a log or bootstrap interval and the Fieller
set at the same level is unbounded. `tidy()` returns the ratio, its
limits, the ends of any excluded gap, and the form of the set in one row
(see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

**Fieller's confidence set.** For a candidate ratio \\\rho\\ the
per-unit contrast \\d_i = x_i - \rho y_i\\ (with a baseline, \\(x_i -
z_i) - \rho (y_i - z_i)\\) has mean zero exactly when \\\rho\\ is the
population ratio, so the one-sample *t* test of \\d\\ on \\N - 1\\
degrees of freedom tests each candidate, and the confidence set is every
\\\rho\\ it does not reject (Fieller, 1941, 1954). The 1941 reference is
Fieller's study of the biological standardization of insulin, which
Fieller (1954) cites with the year 1940. Written through the means and
their covariance matrix \\V\\, the test statistic is \$\$T(\rho) =
\frac{\hat\theta_1 - \rho \hat\theta_2} {\sqrt{v\_{11} - 2 \rho
v\_{12} + \rho^2 v\_{22}}},\$\$ with \\V\\ the sample covariance matrix
of the two variables divided by \\N\\, and the set is the solution of a
quadratic inequality in \\\rho\\. With equal tail areas and no clusters,
each finite end \\L\\ of the set (a limit or an end of the gap) is a
value at which `t.test(x - L * y)` has *p*-value `1 - conf_level` (with
a baseline, `t.test((x - baseline) - L * (y - baseline))`), and under
bivariate normality the coverage is exact. The form of the set is
decided by the denominator's *t*, \\\hat\theta_2 / \sqrt{v\_{22}}\\, in
the `t_denominator` row (Fieller, 1954, p. 177). The row keeps the sign
of the denominator, and only its absolute value matters: when the
absolute value of the denominator's *t* exceeds the critical value (with
unequal tail areas, each side of the estimate is bounded when that
absolute value exceeds the side's own critical value) the set is a
bounded interval. When it does not, the set contains arbitrarily large
values of the ratio, as it should when the data cannot distinguish the
denominator from zero: it is then either the complement of an interval,
every value outside a gap, or, when the data are consistent with every
ratio at the chosen level, the whole line. The `ci_lower` and `ci_upper`
rows are the outer ends of the set (`-Inf` and `Inf` when it is
unbounded), and the `excluded_lower` and `excluded_upper` rows are the
ends of the gap of a complement (`NA` otherwise): the set is every value
from `ci_lower` to `ci_upper` except those strictly between
`excluded_lower` and `excluded_upper`. The `confidence_set` attribute
names the form, and so does the line printed under the table;
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
writes the set in that form, for example as every value except those
between the ends of the gap. With unequal tail areas one side of the set
can be bounded while the other is not, which can also give an interval
unbounded on one side or the union of two intervals. Fieller's own
example in the examples below shows the three forms from one data set,
at three confidence levels. von Luxburg and Franz (2009) derive the same
sets geometrically, from confidence intervals for the means of linear
combinations of the two variables.

**Clustered units.** When units come in clusters whose members are
correlated, the per-unit standard error understates the uncertainty.
With `cluster`, \\V\\ is the cluster-robust covariance matrix formed
from the cluster sums of the centered variables, \$\$V = \frac{Q}{Q - 1}
\frac{1}{N^2} \sum\_{q = 1}^{Q} u_q u_q^\top,\$\$ where \\u_q\\ holds
the sums over cluster \\q\\ of the two variables' deviations from their
means; with \\Q/(Q - 1)\\ as its small-sample factor it is the variance
of an intercept that sandwich's `vcovCL(type = "HC1")` gives. The
reference distribution is *t* on \\Q - 1\\ degrees of freedom, and the
number of clusters is the `n_clusters` row.

**The log ratio.** `ci_method = "log"` forms the interval for \\\log
\rho\\ from the delta method standard error,
\\\sqrt{v\_{11}/\hat\theta_1^2 + v\_{22}/\hat\theta_2^2 - 2
v\_{12}/(\hat\theta_1 \hat\theta_2)}\\, on the same degrees of freedom,
and maps it back by the exponential. This is the response ratio of
meta-analysis (Hedges, Gurevitch, & Curtis, 1999) with the covariance
term that correlated means need (Lajeunesse, 2011). It is defined only
when both terms of the ratio are positive, and is refused otherwise with
a condition of class `dmar_log_ratio_undefined` (which inherits
`dmar_refusal`). It is always bounded; when the denominator is well
determined it is close to Fieller's interval, and the two part company
as the denominator's *t* falls.

**The bootstrap.** `"percentile"` and `"bca"` resample units (pairs, or
triples with a baseline), or whole clusters when `cluster` is supplied,
and compute the ratio in every replication (Davison & Hinkley, 1997).
The BCa interval reads its limits at quantile positions moved from those
of the percentile interval by a bias correction, the median bias of the
replications in normal units, and by an acceleration, the rate at which
the standard error of the estimate changes with the parameter on a
normalized scale (Efron & Tibshirani, 1993); its acceleration is
estimated by the jackknife, deleting one unit (or one cluster) at a
time. A bootstrap interval for a ratio is always bounded, whatever the
denominator's *t*, so it can undercover badly when the denominator is
near zero, exactly where Fieller's set becomes unbounded. When the
Fieller set at the same level is unbounded, the line under the table
says so. Summary statistics cannot be resampled, so a bootstrap from
them is refused with a condition of class `dmar_bootstrap_needs_data`.

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

**Reading the interval for noninferiority.** With a baseline, the
candidate is shown to retain more than a fraction \\f\\ of the
reference's advantage, at one-sided level `alpha_lower`, exactly when
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

Hedges, L. V., Gurevitch, J., & Curtis, P. S. (1999). The meta-analysis
of response ratios in experimental ecology. *Ecology, 80*(4), 1150–1156.
[doi:10.1890/0012-9658(1999)080\[1150:TMAORR\]2.0.CO;2](https://doi.org/10.1890/0012-9658%281999%29080%5B1150%3ATMAORR%5D2.0.CO%3B2)

Lajeunesse, M. J. (2011). On the meta-analysis of response ratios for
studies with correlated and multi-group designs. *Ecology, 92*(11),
2049–2055. [doi:10.1890/11-0423.1](https://doi.org/10.1890/11-0423.1)

Pigeot, I., Schäfer, J., Röhmel, J., & Hauschke, D. (2003). Assessing
non-inferiority of a new treatment in a three-arm clinical trial
including a placebo. *Statistics in Medicine, 22*(6), 883–899.
[doi:10.1002/sim.1450](https://doi.org/10.1002/sim.1450)

von Luxburg, U., & Franz, V. H. (2009). A geometric approach to
confidence sets for ratios: Fieller's theorem, generalizations and
bootstrap. *Statistica Sinica, 19*(3), 1095–1117.

## See also

[`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md)
for independent groups;
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md)
for another effect size of paired observations.

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
[`ci_ratio_of_means()`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md),
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
# Depressive symptoms at the end of treatment relative to the start,
# for the ten participants given the SSRI: the ratio of the posttest
# mean to the pretest mean of the Beck Depression Inventory.
ssri <- depression_bdi[depression_bdi$condition == "ssri", ]
ci_ratio_of_means_paired(ssri$bdi_post, ssri$bdi_pre)
#>  term           value 
#>  ci_lower       0.218 
#>  ratio          0.429 
#>  ci_upper       0.632 
#>  excluded_lower <NA>  
#>  excluded_upper <NA>  
#>  mean_1         7.300 
#>  mean_2         17.000
#>  mean_baseline  <NA>  
#>  se_numerator   1.693 
#>  se_denominator 1.291 
#>  r              0.442 
#>  t_denominator  13.168
#>  df             9     
#>  N              10    
#>  n_clusters     <NA>  
#> 
#> Confidence level: 95%
#> 
#> Fieller confidence set: a bounded interval.

# Fieller's (1954, p. 180) example, the ratio of the mean effects of two
# hypnotic drugs given to the same ten patients, from his summary
# statistics (variances of the means 0.3201 and 0.4009, covariance
# 0.2848). Its denominator's t is 3.68, so the set is a bounded interval
# at 95 percent, the complement of an interval at 99.5 percent, and the
# whole line at 99.9 percent, as Fieller reported.
drug_a_sd <- sqrt(10 * 0.3201)
drug_b_sd <- sqrt(10 * 0.4009)
r_ab <- 0.2848 / sqrt(0.3201 * 0.4009)
ci_ratio_of_means_paired(mean_1 = 0.75, sd_1 = drug_a_sd, mean_2 = 2.33,
                         sd_2 = drug_b_sd, r = r_ab, N = 10)
#>  term           value 
#>  ci_lower       -0.485
#>  ratio          0.322 
#>  ci_upper       0.657 
#>  excluded_lower <NA>  
#>  excluded_upper <NA>  
#>  mean_1         0.750 
#>  mean_2         2.330 
#>  mean_baseline  <NA>  
#>  se_numerator   0.566 
#>  se_denominator 0.633 
#>  r              0.795 
#>  t_denominator  3.680 
#>  df             9     
#>  N              10    
#>  n_clusters     <NA>  
#> 
#> Confidence level: 95%
#> 
#> Fieller confidence set: a bounded interval.
ci_ratio_of_means_paired(mean_1 = 0.75, sd_1 = drug_a_sd, mean_2 = 2.33,
                         sd_2 = drug_b_sd, r = r_ab, N = 10,
                         conf_level = 0.995)
#>  term           value  
#>  ci_lower       -Inf   
#>  ratio          0.322  
#>  ci_upper       Inf    
#>  excluded_lower 0.896  
#>  excluded_upper 147.018
#>  mean_1         0.750  
#>  mean_2         2.330  
#>  mean_baseline  <NA>   
#>  se_numerator   0.566  
#>  se_denominator 0.633  
#>  r              0.795  
#>  t_denominator  3.680  
#>  df             9      
#>  N              10     
#>  n_clusters     <NA>   
#> 
#> Confidence level: 99.5%
#> 
#> Fieller confidence set: the complement of an interval, every value except those strictly between excluded_lower and excluded_upper. The set is unbounded because |t_denominator| does not exceed the critical value.
ci_ratio_of_means_paired(mean_1 = 0.75, sd_1 = drug_a_sd, mean_2 = 2.33,
                         sd_2 = drug_b_sd, r = r_ab, N = 10,
                         conf_level = 0.999)
#>  term           value
#>  ci_lower       -Inf 
#>  ratio          0.322
#>  ci_upper       Inf  
#>  excluded_lower <NA> 
#>  excluded_upper <NA> 
#>  mean_1         0.750
#>  mean_2         2.330
#>  mean_baseline  <NA> 
#>  se_numerator   0.566
#>  se_denominator 0.633
#>  r              0.795
#>  t_denominator  3.680
#>  df             9    
#>  N              10   
#>  n_clusters     <NA> 
#> 
#> Confidence level: 99.9%
#> 
#> Fieller confidence set: the whole line; no value of the ratio is excluded at this level. The set is unbounded because |t_denominator| does not exceed the critical value.

# A retained fraction: a candidate procedure, a reference procedure, and
# a trivial baseline scored on the same 300 items, which come from 60
# topics of five items each. The ratio is the share of the reference's
# gain over the baseline that the candidate keeps, with a one-sided
# 95 percent lower bound and topics as clusters.
set.seed(113)
topic <- rep(1:60, each = 5)
topic_effect <- rnorm(60, sd = 0.08)[topic]
trivial <- 0.20 + topic_effect + rnorm(300, sd = 0.10)
reference <- trivial + 0.15 + rnorm(300, sd = 0.10)
candidate <- trivial + 0.13 + rnorm(300, sd = 0.10)
ci_ratio_of_means_paired(candidate, reference, baseline = trivial,
                         cluster = topic, alpha_lower = 0.05,
                         alpha_upper = 0)
#>  term           value 
#>  ci_lower       0.725 
#>  ratio          0.800 
#>  ci_upper       Inf   
#>  excluded_lower <NA>  
#>  excluded_upper <NA>  
#>  mean_1         0.338 
#>  mean_2         0.371 
#>  mean_baseline  0.208 
#>  se_numerator   0.007 
#>  se_denominator 0.006 
#>  r              0.109 
#>  t_denominator  28.365
#>  df             59    
#>  N              300   
#>  n_clusters     60    
#> 
#> Confidence level: 95% (one-sided: the lower limit is the 95% lower confidence bound, with tail area 0.05 below it and none above the upper limit)
#> 
#> Fieller confidence set: an interval unbounded above.

# A percentile bootstrap interval for the SSRI ratio; B = 2000 keeps the
# example quick, and a reported interval deserves the default B = 10000.
ci_ratio_of_means_paired(ssri$bdi_post, ssri$bdi_pre,
                         ci_method = "percentile", B = 2000, seed = 113)
#>  term           value 
#>  ci_lower       0.253 
#>  ratio          0.429 
#>  ci_upper       0.588 
#>  excluded_lower <NA>  
#>  excluded_upper <NA>  
#>  mean_1         7.300 
#>  mean_2         17.000
#>  mean_baseline  <NA>  
#>  se_numerator   1.693 
#>  se_denominator 1.291 
#>  r              0.442 
#>  t_denominator  13.168
#>  df             9     
#>  N              10    
#>  n_clusters     <NA>  
#> 
#> Confidence level: 95%
#> 
#> Percentile bootstrap interval: a bounded interval.
#> 
#> Bootstrap: 2000 replications, 2000 used (0 dropped; at least 400 must be requested)
```
