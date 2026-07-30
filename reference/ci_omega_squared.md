# Confidence Interval for Omega Squared (Effect Size for ANOVA)

Computes the point estimate and an exact, noncentrality-based confidence
interval for the population omega squared (\\\omega^2\\), the proportion
of variance in the dependent variable accounted for by a fixed effect.
Accepts either the raw ANOVA summary (*F*, effect df, error df, total
*N*) or a fitted [`aov`](https://rdrr.io/r/stats/aov.html) /
[`lm`](https://rdrr.io/r/stats/lm.html) object, in which case the
function returns a row per effect (partial \\\omega^2\\ in factorial
designs).

## Usage

``` r
ci_omega_squared(
  object = NULL,
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

  Optional. A fitted [`aov`](https://rdrr.io/r/stats/aov.html) or
  [`lm`](https://rdrr.io/r/stats/lm.html) object. When supplied, the
  function loops over the non-`Residuals` rows of
  [`anova`](https://rdrr.io/r/stats/anova.html)`(object)` and returns
  one row per effect.

- F_value:

  Observed *F*-value from the fixed-effects ANOVA (ignored if `object`
  is supplied).

- df_effect:

  Numerator degrees of freedom for the effect (ignored if `object` is
  supplied).

- df_error:

  Error (residual) degrees of freedom (ignored if `object` is supplied).

- N:

  Total sample size (ignored if `object` is supplied;
  [`nobs`](https://rdrr.io/r/stats/nobs.html)`(object)` is used
  instead).

- conf_level:

  Desired confidence coverage; default `0.95`. Used only when
  `alpha_lower` and `alpha_upper` are both `NULL`.

- alpha_lower, alpha_upper:

  Optional Type I error on the lower and upper side. If both are `NULL`,
  a symmetric interval at `conf_level` is used. If both are supplied,
  `conf_level` is recomputed as `1 - alpha_lower - alpha_upper`.

## Value

A `data.frame` with one row per effect and the columns `effect`,
`omega_squared` (point estimate), `lower_limit`, `upper_limit`,
`F_value`, `df_effect`, `df_error`, and `N`. When the raw-argument
interface is used, `effect` is `"overall"`.

## Details

**Point estimate.** The function reports the usual sample omega squared,
which for a one-way design can be written as \$\$\hat{\omega}^2 =
\frac{\mathit{SS}\_{\text{effect}} - df\_{\text{effect}} \cdot
\mathit{MS}\_{\text{error}}}{\mathit{SS}\_{\text{total}} +
\mathit{MS}\_{\text{error}}} = \frac{df\_{\text{effect}} (F -
1)}{df\_{\text{effect}} (F - 1) + N}\$\$ (Hays, 1994; Keppel, 1991). For
factorial designs the same formula applied per effect yields *partial*
omega squared (Olejnik & Algina, 2003); values below zero are truncated
to zero.

**Confidence interval.** The CI is constructed by Steiger's (2004, Eq.
17) confidence interval transformation principle: a CI for the
noncentrality parameter \\\lambda\\ of the *F* distribution is obtained
(via
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md))
and then mapped through \$\$\omega^2\_{\text{bound}} =
\frac{\lambda\_{\text{bound}}}{\lambda\_{\text{bound}} + N}.\$\$ When
the lower CI on \\\lambda\\ is not identified (i.e., the observed *F* is
below the one-sided critical value), the lower limit on \\\omega^2\\ is
set to 0, matching the convention used in
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md). In a
one-way design, the interval produced here is identical to the CI for
\\\eta^2\\ from
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md); the
two estimands coincide in the population and differ only in their
*sample* estimators (Steiger, 2004).

**Sums of squares in factorial designs.** When a fitted model is
supplied, the function reads the *F*-values from
[`anova()`](https://rdrr.io/r/stats/anova.html), which in base R uses
Type I (sequential) sums of squares. For balanced designs, Types I, II,
and III give identical *F*-values; for unbalanced designs they differ.
If Type II or III *F*-values are required, compute them with e.g.\\
`car::Anova(object, type = 3)` and pass the relevant *F* / df into the
raw-argument interface.

## References

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.

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
ed.). Routledge. (See Chapter 3 on \\\eta^2\\, Chapter 7 on factorial
designs, and Chapter 11 on generalized \\\eta^2\\ for within-subjects
designs.)

Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: Measures of effect size for some common research designs.
*Psychological Methods, 8*(4), 434–447.
[doi:10.1037/1082-989X.8.4.434](https://doi.org/10.1037/1082-989X.8.4.434)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
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
# 1. Raw-argument interface. Bargman's (1970) example, also used in
#        Venables (1975), Fleishman (1980), and Steiger (2004): a 5-group
#        one-way ANOVA with 11 subjects per group, observed F = 11.221.
ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#>  effect  omega_squared lower_limit upper_limit F_value df_effect df_error N 
#>  overall 0.426         0.226       0.587       11.2    4         50       55

# Same example with a 90% confidence interval.
ci_omega_squared(
  F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
  conf_level = 0.90
)
#>  effect  omega_squared lower_limit upper_limit F_value df_effect df_error N 
#>  overall 0.426         0.261       0.565       11.2    4         50       55

# 2. One way ANOVA from a fitted model (PlantGrowth, 3 groups, N = 30).
fit_one <- aov(weight ~ group, data = PlantGrowth)
ci_omega_squared(fit_one)
#>  effect omega_squared lower_limit upper_limit F_value df_effect df_error N 
#>  group  0.204         0.0099      0.464       4.85    2         27       30

# 3. Factorial ANOVA: partial omega squared per effect (balanced
#        warpbreaks data: 2 wool types x 3 tensions, 9 per cell, N = 54).
fit_factorial <- aov(breaks ~ wool * tension, data = warpbreaks)
ci_omega_squared(fit_factorial)
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution; the lower noncentrality limit has been clamped to 0 and the reported 'prob_greater' on the lower_limit row reflects the actual upper-tail probability at lambda = 0.
#>  effect       omega_squared lower_limit upper_limit F_value df_effect df_error
#>  wool         0.0487        0           0.222       3.77    1         48      
#>  tension      0.217         0.0558      0.411       8.5     2         48      
#>  wool:tension 0.106         0.00191     0.298       4.19    2         48      
#>  N 
#>  54
#>  54
#>  54
```
