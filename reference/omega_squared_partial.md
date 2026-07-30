# Partial Omega Squared (Effect Size for ANOVA)

Computes the sample *partial* omega squared (\\\omega^2_p\\), Hays'
(1994) bias-corrected estimator of the proportion of population variance
in the dependent variable accounted for by a fixed effect after the
variance attributable to the other effects in the model has been
removed: \$\$\hat{\omega}^2_p \\=\\ \frac{\mathit{SS}\_{\text{effect}} -
df\_{\text{effect}} \cdot \mathit{MS}\_{\text{error}}}
{\mathit{SS}\_{\text{effect}} + (N - df\_{\text{effect}}) \cdot
\mathit{MS}\_{\text{error}}} \\=\\ \frac{df\_{\text{effect}} (F -
1)}{df\_{\text{effect}} (F - 1) + N}.\$\$ Accepts either the raw ANOVA
summary (*F*, effect df, error df, total *N*) or a fitted
`aov`/`lm`/`aovlist` object, in which case the function returns one row
per effect (with stratum identification for within-subjects fits).

## Usage

``` r
omega_squared_partial(
  object = NULL,
  F_value = NULL,
  df_effect = NULL,
  df_error = NULL,
  N = NULL
)
```

## Arguments

- object:

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (multi-stratum
  aov fit, e.g.\\ `aov(y ~ A + Error(subject/A), data = d)`). When
  supplied, the function loops over the non-`Residuals` rows and returns
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
`omega_squared_partial` (point estimate), `F_value`, `df_effect`,
`df_error`, and `N`. When the raw-argument interface is used, `effect`
is `"overall"`. Negative point estimates (which occur whenever *F* \< 1)
are truncated to zero, matching the convention used by
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md)
and
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).

## Details

This function is the explicitly-named counterpart of
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md).
The two share the same point-estimate formula , in a one-way ANOVA they
coincide with the total \\\omega^2\\; in a factorial ANOVA both return
the per-effect *partial* value computed against the model's residual
mean square. `omega_squared_partial` is provided so that user code that
explicitly intends partial \\\omega^2\\ carries that meaning in its
name, parallel to the
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)
/
[`eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md)
pair.

**Why partial omega squared and not total.** In a one-way ANOVA, partial
\\\omega^2\\ reduces to total \\\omega^2\\; in a factorial ANOVA the two
diverge. Total \\\omega^2\\ for an effect divides its variance
contribution by the *total* population variance of \\Y\\, so adding
orthogonal factors to a study mechanically shrinks each effect's total
\\\omega^2\\. Partial \\\omega^2\\ divides instead by the variance that
is left after the other effects in the model have been partialled out,
so a given fixed effect's partial \\\omega^2\\ is (approximately)
invariant to whether additional orthogonal factors are present (Olejnik
& Algina, 2003; Maxwell, Delaney, & Kelley, 2027, Sections 7.4.4 and
8.4). For that reason, partial \\\omega^2\\ is the chapter's preferred
effect size index when off-factors are "extrinsic" (i.e., would not vary
in a hypothetical full replication of the population setup).

**Bias correction vs.\\ partial eta squared.** \\\hat{\eta}^2_p\\ , the
sample partial eta squared, is the proportion of *sample* variance
accounted for and is upward-biased as an estimator of the population
\\\eta^2_p\\. \\\hat{\omega}^2_p\\ subtracts \\df\_{\text{effect}} \cdot
\mathit{MS}\_{\text{error}}\\ from the effect's sum of squares and
rescales, yielding an estimator of the *population* variance proportion
with substantially smaller bias (Hays, 1994; Olejnik & Algina, 2000;
Kelley, 2007). Truncation at zero is conventional when the unbiased
estimator goes negative because \\\omega^2 \ge 0\\ by definition.

**Hand-in-hand with
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).**
Pair this function with
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
when reporting effect sizes: `omega_squared_partial()` returns the point
estimate(s) and
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
returns the same point estimate plus its noncentral *F* confidence
limits (Steiger, 2004; Kelley, 2007). The columns shared by the two
functions are aligned so the outputs compose cleanly with
[`merge()`](https://rdrr.io/r/base/merge.html) or a join.

**Sums of squares in unbalanced factorial designs.**
[`anova()`](https://rdrr.io/r/stats/anova.html) on an `aov`/`lm` uses
Type I (sequential) sums of squares. For balanced designs all three SS
types agree; for unbalanced designs they differ. If Type II or III
*F*-values are required, compute them with e.g.\\
`car::Anova(object, type = 3)` and pass the relevant *F* and degrees of
freedom into the raw-argument interface.

## References

Cohen, J. (1973). Eta-squared and partial eta-squared in fixed factor
ANOVA designs. *Educational and Psychological Measurement, 33*(1),
107–112.

Hays, W. L. (1994). *Statistics* (5th ed.). Fort Worth, TX: Harcourt
Brace College Publishers.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Keppel, G., & Wickens, T. D. (2004). *Design and analysis: A
researcher's handbook* (4th ed.). Pearson Prentice Hall.

Keren, G., & Lewis, C. (1979). Partial omega squared for ANOVA designs.
*Educational and Psychological Measurement, 39*(1), 119–128.

Maxwell, S. E., Camp, C. J., & Arvey, R. D. (1981). Measures of strength
of association: A comparative examination. *Journal of Applied
Psychology, 66*(5), 525–534.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\eta^2\\, Chapter 7 on factorial
designs, and Chapter 11 on generalized \\\eta^2\\ for within-subjects
designs.)

Olejnik, S., & Algina, J. (2000). Measures of effect size for
comparative studies: Applications, interpretations, and limitations.
*Contemporary Educational Psychology, 25*(3), 241–286.
[doi:10.1006/ceps.2000.1040](https://doi.org/10.1006/ceps.2000.1040)

Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: Measures of effect size for some common research designs.
*Psychological Methods, 8*(4), 434–447.
[doi:10.1037/1082-989X.8.4.434](https://doi.org/10.1037/1082-989X.8.4.434)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md)

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
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Raw-argument interface (Bargman 1970 / Steiger 2004 example):
#        five groups of 11, observed F = 11.221.
omega_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
#>    effect omega_squared_partial F_value df_effect df_error  N
#> 1 overall             0.4263902  11.221         4       50 55

# 2. Factorial ANOVA: partial omega squared per effect on the
#        balanced warpbreaks data (2 wool x 3 tension, 9/cell, N = 54).
fit_factorial <- aov(breaks ~ wool * tension, data = warpbreaks)
omega_squared_partial(fit_factorial)
#>         effect omega_squared_partial  F_value df_effect df_error  N
#> 1         wool            0.04871442 3.765288         1       48 54
#> 2      tension            0.21734699 8.498047         2       48 54
#> 3 wool:tension            0.10563655 4.189069         2       48 54

# 3. omega_squared_partial() and ci_omega_squared() agree on the
#        point estimate row-by-row.
pt  <- omega_squared_partial(fit_factorial)
ci  <- ci_omega_squared(fit_factorial)
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution; the lower noncentrality limit has been clamped to 0 and the reported 'prob_greater' on the lower_limit row reflects the actual upper-tail probability at lambda = 0.
pt$omega_squared_partial
#> [1] 0.04871442 0.21734699 0.10563655
ci$omega_squared
#> [1] 0.04871442 0.21734699 0.10563655

# 4. The named pair: omega_squared() and omega_squared_partial()
#        report identical numbers in this design; the only difference
#        is the name of the value column, which makes the user's
#        intent (partial) explicit.
omega_squared(fit_factorial)$omega_squared
#> [1] 0.04871442 0.21734699 0.10563655
omega_squared_partial(fit_factorial)$omega_squared_partial
#> [1] 0.04871442 0.21734699 0.10563655
```
