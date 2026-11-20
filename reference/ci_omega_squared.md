# Confidence Interval for Omega Squared (Effect Size for ANOVA)

Computes the point estimate of omega squared (\\\omega^2\\) of each
effect of an analysis of variance and an exact, noncentrality-based
confidence interval for its population value, for the *partial* and the
*total* estimand (see
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md)).
With `estimand = "all"`, the default, the function reports both; the
interval from an effect's *F* statistic is an interval for the partial
value, and the total has one only where it equals the partial value, in
a one-way design (see Details). `"total"` or `"partial"` reports one
estimand. The function takes a fitted model
([`aov`](https://rdrr.io/r/stats/aov.html),
[`lm`](https://rdrr.io/r/stats/lm.html), or an `aovlist` from
[`aov()`](https://rdrr.io/r/stats/aov.html) with an `Error()` term) or
an *F* statistic with its degrees of freedom and *N*.

## Usage

``` r
ci_omega_squared(
  object = NULL,
  estimand = c("all", "total", "partial"),
  F_value = NULL,
  df_effect = NULL,
  df_error = NULL,
  N = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL
)
```

## Arguments

- object:

  Optional. A fitted [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` object. When
  supplied, the function returns one row per effect, with the stratum
  named for an `aovlist`. Give either `object` or the *F* interface, not
  both.

- estimand:

  Which estimand to report: `"all"` (the default, both, with `NA` where
  the inputs determine no estimate or no interval), `"total"`, or
  `"partial"`. An estimand asked for by name for which the inputs
  determine no interval is refused (see Details).

- F_value:

  Observed *F*-value from the fixed-effects ANOVA, or a vector of the
  *F* statistics of several effects tested against the same error term.

- df_effect:

  Numerator degrees of freedom for the effect (one value for each
  effect, or one for all).

- df_error:

  Error (residual) degrees of freedom.

- N:

  Total sample size, the number of observations
  ([`nobs`](https://rdrr.io/r/stats/nobs.html)`(object)` is used for a
  fitted model). For an *F* from a design with an `Error()` stratum,
  give the size of the effect's stratum as `N`, as
  [`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md)
  and
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
  use it for such a fit: the number of subjects for an effect tested
  between subjects, and for an effect tested within subjects the
  stratum's total degrees of freedom, which is `df_effect + df_error`
  when the effect is the only one tested in its stratum. `N` must be at
  least `df_effect + df_error`.

- conf_level:

  Desired confidence coverage; default `0.95`. Used only when
  `alpha_lower` and `alpha_upper` are both `NULL`.

- alpha_lower, alpha_upper:

  Optional Type I error on the lower and upper side. If both are `NULL`,
  a symmetric interval at `conf_level` is used. If both are supplied,
  `conf_level` is recomputed as `1 - alpha_lower - alpha_upper`. The
  table records the coverage, `1 - alpha_lower - alpha_upper`, as its
  level and, when the two tail areas differ, the areas themselves, which
  the line beneath the printed table,
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md),
  and the `alpha_lower` and `alpha_upper` columns of `tidy()` state;
  with a tail area of 0 the interval is a one-sided confidence bound.

## Value

A `data.frame` (class `dmar_ci_anova`, printed through
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md))
with one row per effect and the columns `effect`; for each estimand
reported, its estimate and limits, `omega_squared_total`,
`omega_squared_total_ci_lower`, `omega_squared_total_ci_upper`,
`omega_squared_partial`, `omega_squared_partial_ci_lower`, and
`omega_squared_partial_ci_upper`; `stratum` (for an `aovlist` fit); and
`F_value`, `df_effect`, `df_error`, and `N`. With the *F* interface
`effect` is `"overall"` (the names of `F_value`, when it has them, or
`effect_1`, `effect_2`, and so on, for several). `N` is the number of
observations; the estimate and limits of an effect of an `aovlist` fit
are computed with the size of its stratum, as in
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md).
The table records the confidence level of its intervals in its
`conf_level` attribute (`conf_level`, or `1 - alpha_lower - alpha_upper`
when the tail areas are supplied), which the console prints beneath it
and
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
reports; the `estimand` attribute records the estimand asked for, and
the `estimand_notes` attribute names each estimate or interval that is
`NA` with the reason, which the console also prints.
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives one row
per effect and estimand, with the limits in `ci_lower` and `ci_upper`
(see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

**Point estimates.** The estimates are those of
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md):
the partial value, which for a fixed effect with *F* on
\\df\_{\text{effect}}\\ degrees of freedom and *N* observations is
\$\$\hat{\omega}^2\_{\text{partial}} = \frac{df\_{\text{effect}} (F -
1)}{df\_{\text{effect}} (F - 1) + N},\$\$ and the total value of a
fitted model with one error stratum, \$\$\hat{\omega}^2\_{\text{total}}
= \frac{\mathit{SS}\_{\text{effect}} - df\_{\text{effect}}
\mathit{MS}\_{\text{error}}}{\mathit{SS}\_{\text{total}} +
\mathit{MS}\_{\text{error}}}\$\$ (Hays, 1994; Keppel, 1991; Olejnik &
Algina, 2000, Table 10). In a one-way design the two are the same
number. Values below zero are reported as zero. An *F* statistic does
not determine the total value, which the *F* interface reports as `NA`,
and a fit with an `Error()` term has no total \\\omega^2\\ (see
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md)).

**Confidence interval.** The CI is constructed by Steiger's (2004,
Proposition 1) confidence interval transformation principle: a CI for
the noncentrality parameter \\\lambda\\ of the *F* distribution is
obtained (via
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md)) and
then mapped through \$\$\omega^2\_{\text{bound}} =
\frac{\lambda\_{\text{bound}}}{\lambda\_{\text{bound}} + N},\$\$ with
the size of its stratum in place of *N* for an effect of an `aovlist`
fit. When the lower CI on \\\lambda\\ is not identified (i.e., the
observed *F* is below the one-sided critical value), the lower limit on
\\\omega^2\\ is set to 0, matching the convention used in
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md). The
population \\\omega^2\\ and \\\eta^2\\ are the same quantity, so the
interval is the one
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
gives (and, in a one-way design,
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md)); the
functions differ only in their *sample* estimators (an implication of
the confidence interval transformation principle of Steiger, 2004).

**Which estimand the interval is for.** For an effect tested between
subjects, \\\lambda/(\lambda + N)\\ is the population *partial* value,
the effect's share of the variance due to the effect and the error
(Steiger, 2004, Eqs. 13 and 16), so the interval is for the partial
\\\omega^2\\. In a one-way design that is also the total, and the
interval of a fitted one-way model is reported for both estimands. For
an effect of a design with other effects, the total depends on their
noncentrality parameters too, which one *F* statistic does not inform,
and no interval is formed for it: with `estimand = "all"` its estimate
is reported with `NA` limits, and `estimand = "total"` is refused with
an error of class `dmar_estimand_undetermined`. The Details of
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
give the sources and the within-subjects case.

**Sums of squares in factorial designs.** When a fitted model is
supplied, the function reads the *F*-values from
[`anova()`](https://rdrr.io/r/stats/anova.html), which in base R uses
Type I (sequential) sums of squares. For balanced designs, Types I, II,
and III give identical *F*-values; for unbalanced designs they differ.
If Type II or III *F*-values are required, compute them with e.g.\\
`car::Anova(object, type = 3)` and pass the relevant *F* / df into the
*F* interface. See the Details of
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)
on the coding of factors that Type III tests depend on.
[`factorial_anova`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md)
applies sum-to-zero contrasts itself and reports partial \\\eta^2\\ and
partial \\\omega^2\\ with their confidence intervals for
between-subjects factorial designs.

## References

Bargmann, R. E. (1970). Interpretation and use of a generalized
discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), *Essays in probability
and statistics* (pp. 35–60). University of North Carolina Press.

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.
[doi:10.1177/001316448004000309](https://doi.org/10.1177/001316448004000309)

Hays, W. L. (1994). *Statistics* (5th ed.). Fort Worth, TX: Harcourt
Brace College Publishers.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Keppel, G. (1991). *Design and analysis: A researcher's handbook* (3rd
ed.). Prentice Hall.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\omega^2\\ and its interval, and
Chapter 7 on factorial designs.)

Olejnik, S., & Algina, J. (2000). Measures of effect size for
comparative studies: Applications, interpretations, and limitations.
*Contemporary Educational Psychology, 25*(3), 241–286.
[doi:10.1006/ceps.2000.1040](https://doi.org/10.1006/ceps.2000.1040)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

Venables, W. (1975). Calculation of confidence intervals for
noncentrality parameters. *Journal of the Royal Statistical Society,
Series B, 37*(3), 406–412.
[doi:10.1111/j.2517-6161.1975.tb01554.x](https://doi.org/10.1111/j.2517-6161.1975.tb01554.x)

## See also

[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
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
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. The F interface. Bargmann's (1970) example, used in turn by
#    Venables (1975), Fleishman (1980), and Steiger (2004): a one-way
#    ANOVA with five groups of 11 and an observed F of 11.221 on 4 and
#    50 df, at 95% and at 90%. F determines the partial value, which in
#    a one-way design is the total.
ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                 estimand = "partial")
#>  effect  omega_squared_partial omega_squared_partial_ci_lower
#>  overall 0.426                 0.226                         
#>  omega_squared_partial_ci_upper F_value df_effect df_error N 
#>  0.587                          11.221  4         50       55
#> 
#> Confidence level: 95%
ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
                 estimand = "partial", conf_level = 0.90)
#>  effect  omega_squared_partial omega_squared_partial_ci_lower
#>  overall 0.426                 0.261                         
#>  omega_squared_partial_ci_upper F_value df_effect df_error N 
#>  0.565                          11.221  4         50       55
#> 
#> Confidence level: 90%

# 2. One-way ANOVA from a fitted model: mean IQ gain differs across
#    the six grades of the pygmalion data (N = 310).
fit_one <- aov(iq_gain ~ factor(grade), data = pygmalion)
ci_omega_squared(fit_one)
#>  effect        omega_squared_total omega_squared_total_ci_lower
#>  factor(grade) 0.104               0.048                       
#>  omega_squared_total_ci_upper omega_squared_partial
#>  0.176                        0.104                
#>  omega_squared_partial_ci_lower omega_squared_partial_ci_upper F_value
#>  0.048                          0.176                          8.159  
#>  df_effect df_error N  
#>  5         304      310
#> 
#> Confidence level: 95%

# 3. Two-factor ANOVA for the manipulated expectancy treatment and the
#    measured grade classification (pygmalion data, N = 310). The
#    treatment by grade interaction is weak here (F = 1.19), so the
#    additive model is used.
fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
ci_omega_squared(fit_additive)
#>  effect        omega_squared_total omega_squared_total_ci_lower
#>  treatment     0.017               <NA>                        
#>  factor(grade) 0.028               <NA>                        
#>  omega_squared_total_ci_upper omega_squared_partial
#>  <NA>                         0.017                
#>  <NA>                         0.028                
#>  omega_squared_partial_ci_lower omega_squared_partial_ci_upper F_value
#>  0.001                          0.062                          6.519  
#>  0.001                          0.082                          2.787  
#>  df_effect df_error N  
#>  1         303      310
#>  5         303      310
#> 
#> Confidence level: 95%
#> 
#> NA omega_squared_total limits: No confidence interval is formed for the total omega squared of an effect in a design with other effects or with an Error() term: the noncentral F interval from the effect's F is an interval for the partial value, and the total depends as well on the other effects of the design (see the Details of ?ci_omega_squared).
```
