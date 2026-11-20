# Confidence Interval for the Square Root of the Signal-to-Noise Ratio

Computes the exact confidence interval for the square root of the
signal-to-noise ratio, the standard deviation of the group means
relative to the error standard deviation. On this root scale the
quantity is an effect size in standard deviation units, the multi-group
analog of the standardized mean difference.

## Usage

``` r
ci_srsnr(
  F_value = NULL,
  df_1 = NULL,
  df_2 = NULL,
  N = NULL,
  means = NULL,
  var_within = NULL,
  n_per_group = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- F_value:

  Observed *F*-value from the analysis of variance. Use this argument
  when re-analyzing existing data.

- df_1:

  Numerator degrees of freedom

- df_2:

  Denominator degrees of freedom

- N:

  Sample size

- means:

  Numeric vector of group means, one per group and at least two, given
  with `var_within` and `n_per_group` in place of `F_value`: the
  observed sample means of a one-way design, from which the function
  computes the observed *F*, or population means for a hypothetical
  interval (see Details).

- var_within:

  The within-group variance that goes with `means`, a single positive
  number: the pooled within-group mean square (the mean square error of
  the one-way analysis of variance), or a population within-group
  variance for a hypothetical interval.

- n_per_group:

  The group sizes that go with `means`: a single size for equal groups,
  or one size per group, each a whole number of at least 2.

- conf_level:

  Confidence interval coverage (i.e., 1 - Type I error rate); default is
  .95

- alpha_lower:

  Type I error for the lower confidence limit

- alpha_upper:

  Type I error for the upper confidence limit. The table records the
  coverage, `1 - alpha_lower - alpha_upper`, as its level and, when the
  two tail areas differ, the areas themselves, which the line beneath
  the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound.

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 2-row `data.frame` with columns `term` and `value`. The `term` values
are `"ci_lower"` and `"ci_upper"`, giving the square roots of the
corresponding signal-to-noise-ratio confidence limits.

## Details

The confidence level is given in one of two ways: as the coverage,
`conf_level`, or as the two tail areas, `alpha_lower` and `alpha_upper`.

The square root of the signal-to-noise ratio is the standard deviation
of the group means, each weighted by its share of the total sample size
\\n_j / N\\, divided by the within-group standard deviation; it is
\\\sqrt{\lambda / N}\\ for the noncentrality parameter \\\lambda\\ of
the *F* test, and it is Cohen's *f* with those weights (see
[`cohen_f`](https://yelleknek.github.io/DMAR/reference/cohen_f.md)).
This function uses the confidence interval transformation principle
(Steiger, 2004) to transform the confidence limits for the noncentrality
parameter to the confidence limits for the square root of the
signal-to-noise ratio. The confidence interval for the noncentrality
parameter itself is given by
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md).

The interval has two routes. The first takes the observed *F*, its
degrees of freedom, and *N*. The second takes the group means
\\\bar{y}\_j\\, the pooled within-group mean square \\MS\_{within}\\,
and the group sizes \\n_j\\ of a one-way design with *k* groups, and
computes \$\$F = \frac{\sum_j n_j (\bar{y}\_j - \bar{y})^2 / (k -
1)}{MS\_{within}}, \qquad \bar{y} = \frac{1}{N} \sum_j n_j
\bar{y}\_j,\$\$ with \\df_1 = k - 1\\, \\df_2 = N - k\\, and \\N =
\sum_j n_j\\. With the observed sample means and mean square, this is
the observed *F*, and the interval is the one the first route gives for
the same data. Values of `df_1`, `df_2`, or `N` given with `means` must
agree with these. Group means that are all equal give \\F = 0\\, which
is refused, as it is on the first route: at \\F = 0\\ the upper
confidence limit is undefined.

Population means and a population within-group variance can be given in
the same places. The result is then a hypothetical interval: the
interval a study with these group sizes would report if its sample means
and its pooled within-group mean square equaled the values given. It is
not an interval from observed data, and it is not the expected interval
of such a study: the limits of the interval vary from sample to sample,
and the interval computed at the population values need not equal their
average.

## References

Bargmann, R. E. (1970). Interpretation and use of a generalized
discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), *Essays in probability
and statistics* (pp. 35–60). University of North Carolina Press.

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.
[doi:10.1177/001316448004000309](https://doi.org/10.1177/001316448004000309)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

Venables, W. (1975). Calculation of confidence intervals for
noncentrality parameters. *Journal of the Royal Statistical Society,
Series B, 37*(3), 406–412.
[doi:10.1111/j.2517-6161.1975.tb01554.x](https://doi.org/10.1111/j.2517-6161.1975.tb01554.x)

## See also

[`ci_snr`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md)

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
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
## Bargmann (1970) gave an example of a one-way ANOVA with an observed
## F value of 11.2213 on 4 and 50 degrees of freedom, which Venables (1975),
## Fleishman (1980), and Steiger (2004) used in turn. Steiger describes it
## as five groups of 11 and reports F as 11.221, the value used here.
## The exact confidence interval for the square root of the signal-to-noise ratio of that example:

ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#>  term     value
#>  ci_lower 0.541
#>  ci_upper 1.193
#> 
#> Confidence level: 95%

ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#>  term     value
#>  ci_lower 0.594
#>  ci_upper 1.141
#> 
#> Confidence level: 90%

ci_srsnr(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = .02, alpha_upper = .03)
#>  term     value
#>  ci_lower 0.525
#>  ci_upper 1.180
#> 
#> Confidence level: 95% (unequal tail areas: 0.02 below the lower limit, 0.03 above the upper limit)

# The interval from the summaries of an observed one-way design: the
# group means, the pooled within-group mean square, and the group sizes
# of the three arms of drinks_trial. The F route gives the same interval
# for the same data.
fit <- stats::lm(log_drinks ~ treatment, data = drinks_trial)
ci_srsnr(means = tapply(drinks_trial$log_drinks, drinks_trial$treatment, mean),
         var_within = sum(stats::residuals(fit)^2) / stats::df.residual(fit),
         n_per_group = as.vector(table(drinks_trial$treatment)))
#>  term     value
#>  ci_lower 0.050
#>  ci_upper 0.510
#> 
#> Confidence level: 95%
ci_srsnr(F_value = stats::anova(fit)[1, "F value"], df_1 = 2, df_2 = 85, N = 88)
#>  term     value
#>  ci_lower 0.050
#>  ci_upper 0.510
#> 
#> Confidence level: 95%

# A hypothetical interval: with population means of 94, 91, 92, and 83, a
# within-group variance of 67.375, and 20 per group, the interval a study
# would report if its summaries equaled these values. It is not the
# expected interval of that study.
ci_srsnr(means = c(94, 91, 92, 83), var_within = 67.375, n_per_group = 20)
#>  term     value
#>  ci_lower 0.238
#>  ci_upper 0.722
#> 
#> Confidence level: 95%
```
