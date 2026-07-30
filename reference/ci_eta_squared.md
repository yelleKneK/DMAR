# Confidence Interval for Eta Squared (Effect Size for ANOVA)

Computes the point estimate and an exact, noncentrality-based confidence
interval for the population eta squared (\\\eta^2\\), the proportion of
variance in the dependent variable accounted for by a fixed effect.
Accepts either the raw ANOVA summary (*F*, effect df, error df, total
*N*) or a fitted model object. Supports both between-subjects designs
([`aov`](https://rdrr.io/r/stats/aov.html) /
[`lm`](https://rdrr.io/r/stats/lm.html)) and within-subjects / mixed
designs (`aovlist` fits with an `Error()` term in the formula). For
factorial and within-subjects designs the function returns one row per
effect with the CI for *partial* \\\eta^2\\ computed against that
effect's own error stratum.

## Usage

``` r
ci_eta_squared(
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

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (multi-stratum
  aov fit, e.g.\\ `aov(y ~ A + Error(subject/A), data = d)`). For
  multi-stratum fits the function walks every error stratum and returns
  one row per non-`Residuals` effect, identifying the stratum used.

- F_value:

  Observed *F*-value (ignored if `object` is supplied).

- df_effect:

  Numerator degrees of freedom for the effect (ignored if `object` is
  supplied).

- df_error:

  Error (residual) degrees of freedom (ignored if `object` is supplied).

- N:

  Total sample size (ignored if `object` is supplied; derived
  automatically from a fitted model, via
  [`nobs`](https://rdrr.io/r/stats/nobs.html)`(object)` for
  single-stratum fits, or by summing the degrees of freedom across error
  strata for `aovlist` fits).

- conf_level:

  Desired confidence coverage; default `0.95`. Used only when
  `alpha_lower` and `alpha_upper` are both `NULL`.

- alpha_lower, alpha_upper:

  Optional Type I error on the lower and upper side. If both are `NULL`,
  a symmetric interval at `conf_level` is used. If both are supplied,
  `conf_level` is recomputed as `1 - alpha_lower - alpha_upper`.

## Value

A `data.frame` with one row per effect. Single-stratum fits and the raw
interface return columns `effect`, `eta_squared`, `lower_limit`,
`upper_limit`, `F_value`, `df_effect`, `df_error`, `N`. `aovlist`
(within-subjects / mixed) fits additionally include a `stratum` column.
With the raw-argument interface `effect` is `"overall"`.

## Details

**Point estimate.** \\\hat{\eta}^2 = df\_{\text{effect}} \cdot F /
(df\_{\text{effect}} \cdot F + df\_{\text{error}})\\, which equals
\\\mathit{SS}\_{\text{effect}}/(\mathit{SS}\_{\text{effect}} +
\mathit{SS}\_{\text{error}})\\. In a one-way ANOVA this is also
\\\mathit{SS}\_{\text{effect}}/\mathit{SS}\_{\text{total}}\\. In a
factorial design the same expression gives the per-effect *partial*
\\\eta^2\\.

**Confidence interval.** The CI is constructed by Steiger's (2004)
confidence interval transformation principle: a CI for the noncentrality
parameter \\\lambda\\ of the *F* distribution is obtained (via
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md))
and then mapped through \$\$\eta^2\_{\text{bound}} =
\frac{\lambda\_{\text{bound}}}{\lambda\_{\text{bound}} + N}.\$\$ This is
the same transformation used by
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md) and
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md);
the three functions share CI machinery and differ only in their sample
point estimators. When the lower CI on \\\lambda\\ is not identified
(i.e., the observed *F* is below the one-sided critical value), the
lower limit on \\\eta^2\\ is set to 0.

**Designs supported.**

- *Between-subjects ANOVA*: fitted `aov`/`lm`, or raw *F*/df/N.

- *Within-subjects or mixed ANOVA*: fitted `aovlist`. Each effect's CI
  is built from its own stratum's *F* and residual *df*; `N` is the
  total number of observations across all strata. The reported `stratum`
  column identifies which error term each row used.

**Sums of squares in factorial designs.** When a fitted model is
supplied, *F*-values are read from
[`anova()`](https://rdrr.io/r/stats/anova.html) (single-stratum) or
[`summary()`](https://rdrr.io/r/base/summary.html) (multi-stratum), both
of which use Type I (sequential) sums of squares in base R. For balanced
designs Types I, II, and III agree; for unbalanced designs they differ.
If Type II or III *F*-values are required, compute them with e.g.\\
`car::Anova(object, type = 3)` and pass the relevant *F*, degrees of
freedom, and *N* into the raw-argument interface.

## References

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\eta^2\\, Chapter 7 on factorial
designs, and Chapter 11 on generalized \\\eta^2\\ for within-subjects
designs.)

Smithson, M. (2001). Correct confidence intervals for various regression
effect sizes and parameters: The importance of noncentral distributions
in computing intervals. *Educational and Psychological Measurement, 61*,
605–632.
[doi:10.1177/00131640121971392](https://doi.org/10.1177/00131640121971392)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
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
# 1. Raw-argument interface. Bargman's (1970) example.
ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#>  effect  eta_squared lower_limit upper_limit F_value df_effect df_error N 
#>  overall 0.473       0.226       0.587       11.2    4         50       55

# Same example with a 90% confidence interval.
ci_eta_squared(
  F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
  conf_level = 0.90
)
#>  effect  eta_squared lower_limit upper_limit F_value df_effect df_error N 
#>  overall 0.473       0.261       0.565       11.2    4         50       55

# 2. One way ANOVA from a fitted model.
fit_one <- aov(weight ~ group, data = PlantGrowth)
ci_eta_squared(fit_one)
#>  effect eta_squared lower_limit upper_limit F_value df_effect df_error N 
#>  group  0.264       0.0099      0.464       4.85    2         27       30

# 3. Factorial ANOVA: partial eta squared per effect.
fit_factorial <- aov(breaks ~ wool * tension, data = warpbreaks)
ci_eta_squared(fit_factorial)
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution; the lower noncentrality limit has been clamped to 0 and the reported 'prob_greater' on the lower_limit row reflects the actual upper-tail probability at lambda = 0.
#>  effect       eta_squared lower_limit upper_limit F_value df_effect df_error N 
#>  wool         0.0727      0           0.222       3.77    1         48       54
#>  tension      0.261       0.0558      0.411       8.5     2         48       54
#>  wool:tension 0.149       0.00191     0.298       4.19    2         48       54

# 4. Within-subjects ANOVA. CI computed against the within-subjects
#        error stratum (the row reports which one via 'stratum').
set.seed(113)
n <- 20
rm_data <- data.frame(
  subject = factor(rep(seq_len(n), each = 3)),
  time    = factor(rep(c("Pre", "Mid", "Post"), n),
                   levels = c("Pre", "Mid", "Post")),
  y       = rnorm(n, sd = 1.5)[rep(seq_len(n), each = 3)] +
            0.7 * rep(1:3, n) + rnorm(n * 3, sd = 1.2)
)
fit_rm <- aov(y ~ time + Error(subject/time), data = rm_data)
ci_eta_squared(fit_rm)
#>  effect eta_squared lower_limit upper_limit stratum      F_value df_effect
#>  time   0.231       0.0152      0.322       subject:time 5.7     2        
#>  df_error N 
#>  38       60
```
