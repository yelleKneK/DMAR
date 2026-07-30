# Omega Squared (Effect Size for ANOVA)

Computes the sample omega squared (\\\omega^2\\), Hays' (1994)
bias-corrected estimator of the proportion of variance in the dependent
variable accounted for by a fixed effect. Accepts either the raw ANOVA
summary (*F*, the effect and error degrees of freedom, and total *N*) or
a fitted [`aov`](https://rdrr.io/r/stats/aov.html) or
[`lm`](https://rdrr.io/r/stats/lm.html) object, in which case the
function returns one row per effect (partial \\\omega^2\\ in factorial
designs).

## Usage

``` r
omega_squared(
  object = NULL,
  F_value = NULL,
  df_effect = NULL,
  df_error = NULL,
  N = NULL
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

## Value

A `data.frame` with one row per effect. The columns are `effect`,
`omega_squared` (point estimate), `F_value`, `df_effect`, `df_error`,
and `N`. When the raw-argument interface is used, `effect` is
`"overall"`.

## Details

The confidence interval is provided by the separate
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
paralleling the existing
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md)/[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)
and
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)/[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
pairings.

**Point estimate.** The reported value is Hays' (1994) sample omega
squared, which for a one-way design is \$\$\hat{\omega}^2 =
\frac{df\_{\text{effect}} (F - 1)}{df\_{\text{effect}} (F - 1) + N}.\$\$
For factorial designs the same formula applied per effect yields
*partial* omega squared (Olejnik & Algina, 2003); negative values are
truncated to zero. This is the same point-estimate convention used by
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
so the two functions agree on the point estimate row by row.

**Hand-in-hand with
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).**
Pair this function with
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
when reporting effect sizes: `omega_squared()` returns the point
estimate(s), and
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
returns the same point estimate plus its noncentrality-based confidence
limits (Steiger, 2004). The columns shared by the two functions
(`effect`, `omega_squared`, `F_value`, `df_effect`, `df_error`, `N`) are
aligned so the outputs compose cleanly with
[`merge()`](https://rdrr.io/r/base/merge.html) or a join.

**Sums of squares in factorial designs.**
[`anova()`](https://rdrr.io/r/stats/anova.html) on an `aov`/`lm` uses
Type I (sequential) sums of squares. For balanced designs all three
types agree; for unbalanced designs they differ. If Type II or III
*F*-values are required, compute them with e.g.\\
`car::Anova(object, type = 3)` and pass the relevant *F* and degrees of
freedom into the raw-argument interface.

## References

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

[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Raw-argument interface. Bargman's (1970) 5-group one-way ANOVA,
#        also used in Venables (1975), Fleishman (1980), and Steiger (2004):
#        11 subjects per group, observed F = 11.221.
omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#>    effect omega_squared F_value df_effect df_error  N
#> 1 overall     0.4263902  11.221         4       50 55

# 2. One way ANOVA from a fitted model (PlantGrowth, 3 groups, N = 30).
fit_one <- aov(weight ~ group, data = PlantGrowth)
omega_squared(fit_one)
#>   effect omega_squared  F_value df_effect df_error  N
#> 1  group     0.2040788 4.846088         2       27 30

# 3. Factorial ANOVA: partial omega squared per effect (balanced
#        warpbreaks data: 2 wool types x 3 tensions, 9 per cell, N = 54).
fit_factorial <- aov(breaks ~ wool * tension, data = warpbreaks)
omega_squared(fit_factorial)
#>         effect omega_squared  F_value df_effect df_error  N
#> 1         wool    0.04871442 3.765288         1       48 54
#> 2      tension    0.21734699 8.498047         2       48 54
#> 3 wool:tension    0.10563655 4.189069         2       48 54

# 4. omega_squared() and ci_omega_squared() compose: the point
#        estimates agree row-by-row.
pt  <- omega_squared(fit_factorial)
ci  <- ci_omega_squared(fit_factorial)
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution; the lower noncentrality limit has been clamped to 0 and the reported 'prob_greater' on the lower_limit row reflects the actual upper-tail probability at lambda = 0.
merge(pt, ci, by = c("effect", "omega_squared",
                     "F_value", "df_effect", "df_error", "N"))
#>         effect omega_squared  F_value df_effect df_error  N lower_limit
#> 1      tension    0.21734699 8.498047         2       48 54 0.055751214
#> 2         wool    0.04871442 3.765288         1       48 54 0.000000000
#> 3 wool:tension    0.10563655 4.189069         2       48 54 0.001906446
#>   upper_limit
#> 1   0.4106053
#> 2   0.2222776
#> 3   0.2983809
```
