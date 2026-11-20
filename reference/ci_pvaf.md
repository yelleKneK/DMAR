# Confidence Interval for the Proportion of Variance Accounted for (in the Dependent Variable by Knowing the Levels of the Factor)

Computes the exact confidence limits for the proportion of variance in
the dependent variable accounted for by knowing the levels of the factor
(group status in a single factor design) in a fixed effects analysis of
variance, so an omnibus *F*-test is accompanied by an effect size with a
statement of its precision.

## Usage

``` r
ci_pvaf(
  F_value = NULL,
  df_1 = NULL,
  df_2 = NULL,
  N = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  ...
)
```

## Arguments

- F_value:

  Observed *F*-value from fixed effects analysis of variance

- df_1:

  Numerator degrees of freedom

- df_2:

  Denominator degrees of freedom

- N:

  Sample size

- conf_level:

  Confidence interval coverage (i.e., 1-Type I error rate); default is
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

A 4-row `data.frame` with columns `term`, `value`, `prob_less`, and
`prob_greater`. The `term` values are `"ci_lower"` (the lower confidence
limit on the proportion of variance accounted for, on the \[0, 1\]
scale), `"pvaf"` (the sample proportion of variance accounted for,
`df_1 * F_value / (df_1 * F_value + df_2)`, the same value that eta
squared reports, so the point estimate sits between its confidence
limits), `"ci_upper"` (the upper confidence limit), and
`"actual_coverage"` (the coverage probability of the interval,
`1 - alpha_lower - alpha_upper`; `NA` when the upper limit is
undefined). The `prob_less` and `prob_greater` columns hold, on the row
of each limit, the two tail areas at the observed `F_value` of the
noncentral *F* distribution whose noncentrality parameter corresponds to
that limit: `prob_less` is the area at or above `F_value` and
`prob_greater` the area at or below it (the reverse of the columns of
the same names in
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
where `prob_less` is the area at or below `F_value`), so `prob_less` is
`alpha_lower` on the `"ci_lower"` row and `prob_greater` is
`alpha_upper` on the `"ci_upper"` row. With a tail area of 0 no limit is
computed on that side (the lower limit is 0, or the upper limit is 1),
and the row of that limit holds the requested area, 0, and its
complement, 1. When `alpha_lower` is larger than 0 and the observed *F*
is below the `alpha_lower` critical value, so that the lower limit is 0,
the `"ci_lower"` row holds the tail areas of the central *F*
distribution, whose area at or above `F_value` is the *p*-value of the
*F*-test and is larger than `alpha_lower`; those areas describe the
central distribution and not the coverage of the interval. Both columns
are `NA` on the `"pvaf"` and `"actual_coverage"` rows, and on the
`"ci_upper"` row when the upper limit is undefined.

## Details

The confidence level must be specified in one of following two ways:
using confidence interval coverage (`conf_level`), or lower and upper
confidence limits (`alpha_lower` and `alpha_upper`).

This function uses the confidence interval transformation principle
(Steiger, 2004) to transform the confidence limits for the noncentrality
parameter to the confidence limits for the population proportion of
variance accounted for by knowing the group status. The confidence
interval for the noncentral *F* parameter can be obtained from the
function
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
which is used within this function.

**A lower limit of 0.** When the observed *F* is below the `alpha_lower`
critical value of the central *F* distribution (an *F*-test that is not
significant at `alpha_lower`), no positive noncentrality parameter
leaves as little as `alpha_lower` above the observed *F*, and the lower
limit is 0, with a warning. The interval from 0 to the upper limit has
the stated coverage, `1 - alpha_lower - alpha_upper`. A lower limit
misses only when it lies above the population value, which happens with
probability `alpha_lower` whatever that value is, and a lower limit of 0
never lies above it, so a limit of 0 gives up no coverage. In 20,000
simulated analyses of variance of 5 groups of 11 at each of five
population values from 0 to .20, the 95% interval contained the
population value in 95.0% to 95.1% of them (Monte Carlo standard error
0.15 percentage points), although its lower limit was 0 in 27% to 98% of
them. An analysis whose *F* was so small that the central *F*
distribution leaves `alpha_upper` or less below it has no upper limit
(the `"ci_upper"` and `"actual_coverage"` rows are `NA`) and was counted
as not containing the population value: 2.5% of the analyses at a
population value of 0, and fewer at each larger value.

## Note

This function can be used for single or factorial ANOVA designs.

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

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the Analysis of Variance and
Contrast Analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

Venables, W. (1975). Calculation of confidence intervals for
noncentrality parameters. *Journal of the Royal Statistical Society,
Series B, 37*(3), 406–412.
[doi:10.1111/j.2517-6161.1975.tb01554.x](https://doi.org/10.1111/j.2517-6161.1975.tb01554.x)

## See also

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
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
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
## The exact confidence interval for the proportion of variance accounted for of that example:
ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55)
#>  term            value prob_less prob_greater
#>  ci_lower        0.226 0.025     0.975       
#>  pvaf            0.473 <NA>      <NA>        
#>  ci_upper        0.587 0.975     0.025       
#>  actual_coverage 0.950 <NA>      <NA>        
#> 
#> Confidence level: 95%

ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, conf_level = .90)
#>  term            value prob_less prob_greater
#>  ci_lower        0.261 0.050     0.950       
#>  pvaf            0.473 <NA>      <NA>        
#>  ci_upper        0.565 0.950     0.050       
#>  actual_coverage 0.900 <NA>      <NA>        
#> 
#> Confidence level: 90%

ci_pvaf(F_value = 11.221, df_1 = 4, df_2 = 50, N = 55, alpha_lower = 0, alpha_upper = .05)
#>  term            value prob_less prob_greater
#>  ci_lower        0.000 0.000     1.000       
#>  pvaf            0.473 <NA>      <NA>        
#>  ci_upper        0.565 0.950     0.050       
#>  actual_coverage 0.950 <NA>      <NA>        
#> 
#> Confidence level: 95% (one-sided: the upper limit is the 95% upper confidence bound, with tail area 0.05 above it and none below the lower limit)
```
