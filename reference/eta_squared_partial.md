# Partial Eta Squared (Effect Size for ANOVA)

Computes the sample *partial* eta squared (\\\eta^2_p\\), the proportion
of variance accounted for by a fixed effect after the variance
attributable to the other effects in the model has been removed:
\$\$\hat{\eta}^2_p =
\frac{\mathit{SS}\_{\text{effect}}}{\mathit{SS}\_{\text{effect}} +
\mathit{SS}\_{\text{error}}} = \frac{df\_{\text{effect}} \cdot
F}{df\_{\text{effect}} \cdot F + df\_{\text{error}}}.\$\$ Accepts either
the raw ANOVA summary (*F*, effect df, error df) or a fitted
`aov`/`lm`/`aovlist` object, in which case the function returns one row
per effect (with stratum identification for within-subjects fits).

## Usage

``` r
eta_squared_partial(
  object = NULL,
  F_value = NULL,
  df_effect = NULL,
  df_error = NULL
)
```

## Arguments

- object:

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (multi-stratum
  aov fit, e.g.\\ `aov(y ~ A + Error(subject/A), data = d)`).

- F_value:

  Observed *F*-value (ignored if `object` is supplied).

- df_effect:

  Numerator degrees of freedom for the effect (ignored if `object` is
  supplied).

- df_error:

  Error (residual) degrees of freedom (ignored if `object` is supplied).

## Value

A `data.frame` with one row per effect. Single-stratum fits and the raw
interface return columns `effect`, `eta_squared_partial`, `F_value`,
`df_effect`, `df_error`. `aovlist` (within-subjects / mixed) fits
additionally include a `stratum` column identifying which error term
each effect's *F* test came from. With the raw-argument interface
`effect` is `"overall"`.

## Details

This function is the explicitly-named counterpart of
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md).
The two share the same point-estimate formula, in a one-way ANOVA they
coincide with *total* \\\eta^2\\; in a factorial or within-subjects
ANOVA both functions return the per-effect *partial* value computed
against that effect's own error stratum. `eta_squared_partial` is
provided so that user code that explicitly intends partial \\\eta^2\\
carries that meaning in its name.

**Designs supported.** Single-stratum `aov`/`lm` fits and multi-stratum
`aovlist` fits (within-subjects and mixed designs) are both handled by
the model interface. For multi-stratum fits, each effect uses its own
stratum's residual *df*, so a within-subjects factor's partial
\\\eta^2\\ is computed against the within-subjects error and a
between-subjects factor's is computed against the between-subjects
error.

## References

Cohen, J. (1973). Eta-squared and partial eta-squared in fixed factor
ANOVA designs. *Educational and Psychological Measurement, 33*(1),
107–112.

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

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ci_eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_r()`](https://yelleknek.github.io/DMAR/reference/expected_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Raw-argument interface.
eta_squared_partial(F_value = 11.221, df_effect = 4, df_error = 50)
#>  effect  eta_squared_partial F_value df_effect df_error
#>  overall 0.473               11.2    4         50      

# Factorial ANOVA: partial eta squared per effect (pygmalion data:
# expectancy treatment x grade, 2 x 6 with unequal cell sizes,
# N = 310). The treatment is manipulated; grade is a measured
# classification of the pupils.
fit <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
eta_squared_partial(fit)
#>  effect                  eta_squared_partial F_value df_effect df_error
#>  treatment               0.0215              6.54    1         298     
#>  factor(grade)           0.0448              2.8     5         298     
#>  treatment:factor(grade) 0.0196              1.19    5         298     

# Within-subjects ANOVA: per-effect partial eta squared with stratum.
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
eta_squared_partial(fit_rm)
#>  effect eta_squared_partial stratum      F_value df_effect df_error
#>  time   0.231               subject:time 5.7     2         38      
```
