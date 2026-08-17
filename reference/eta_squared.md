# Eta Squared (Effect Size for ANOVA)

Computes the sample eta squared (\\\eta^2\\), the proportion of variance
in the dependent variable accounted for by a fixed effect. Accepts
either the raw ANOVA summary (*F* and the effect and error degrees of
freedom) or a fitted model object. Supports both between-subjects
designs (single-stratum [`aov`](https://rdrr.io/r/stats/aov.html) or
[`lm`](https://rdrr.io/r/stats/lm.html) fits) and within-subjects /
mixed designs (`aovlist` fits produced by
[`aov`](https://rdrr.io/r/stats/aov.html) with an `Error()` term in the
formula). For factorial and within-subjects designs the function returns
*partial* \\\eta^2\\ per effect (one row per non-`Residuals` effect
across all strata), each computed against its own stratum's error term.

## Usage

``` r
eta_squared(object = NULL, F_value = NULL, df_effect = NULL, df_error = NULL)
```

## Arguments

- object:

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (a multi-stratum
  aov fit such as `aov(y ~ A + Error(subject/A), data = d)`). When
  supplied, the function loops over the non-`Residuals` effects across
  all error strata and returns one row per effect, with the stratum
  reported.

- F_value:

  Observed *F*-value from the fixed-effects ANOVA (ignored if `object`
  is supplied).

- df_effect:

  Numerator degrees of freedom for the effect (ignored if `object` is
  supplied).

- df_error:

  Error (residual) degrees of freedom (ignored if `object` is supplied).

## Value

A `data.frame` with one row per effect. For single-stratum (`aov`/`lm`)
fits and the raw-argument interface the columns are `effect`,
`eta_squared`, `F_value`, `df_effect`, `df_error`. For multi-stratum
(`aovlist`) fits an additional `stratum` column reports which error
stratum each effect's *F* test came from. When the raw-argument
interface is used, `effect` is `"overall"`.

## Details

The confidence interval is provided by the separate
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
paralleling the existing
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md)/[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md)
pairing.

**Point estimate.** The function uses the algebraically equivalent
*F*-and-df form \$\$\hat{\eta}^2 = \frac{df\_{\text{effect}} \cdot
F}{df\_{\text{effect}} \cdot F + df\_{\text{error}}},\$\$ which equals
\\\mathit{SS}\_{\text{effect}} / (\mathit{SS}\_{\text{effect}} +
\mathit{SS}\_{\text{error}})\\. In a one-way ANOVA this is also
\\\mathit{SS}\_{\text{effect}} / \mathit{SS}\_{\text{total}}\\, the
conventional *total* \\\eta^2\\. In a factorial design the same
expression yields *partial* \\\eta^2\\ for each effect, because
\\\mathit{SS}\_{\text{error}}\\ appears in the denominator instead of
\\\mathit{SS}\_{\text{total}}\\; this matches the convention used by
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).

**Designs supported.**

- *Between-subjects ANOVA* (one-way or factorial): supply either a
  fitted `aov`/`lm` model or the raw *F* and degrees of freedom for a
  single effect.

- *Within-subjects or mixed ANOVA*: supply a fitted `aovlist` model
  produced with an `Error()` term, e.g.\\
  `aov(y ~ time + Error(subject/time), data = d)`. The function walks
  every error stratum returned by `summary(object)` and reports each
  effect with the stratum's residual *df*, so the \\\eta^2\\ value uses
  the stratum's specific error term. The reported `stratum` column tells
  you which one.

For more advanced model classes (`lmerMod`, `lme`, etc.) the
fitted-model interface is not yet supported; supply the relevant *F* and
degrees of freedom via the raw interface.

**Sums of squares in factorial designs.**
[`anova()`](https://rdrr.io/r/stats/anova.html) on an `aov`/`lm` uses
Type I (sequential) sums of squares. For balanced designs all three
types agree; for unbalanced designs they differ. If Type II or III
*F*-values are required, compute them with e.g.\\
`car::Anova(object, type = 3)` and pass the relevant *F* and degrees of
freedom into the raw-argument interface.

**Generalized eta squared, comparable across designs.** The basic
\\\eta^2\\ (and partial \\\eta^2\\) returned here are not comparable
across studies that differ in factor structure. See
[`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)
for a comparable alternative (Olejnik & Algina, 2003; Bakeman, 2005).

## References

Bakeman, R. (2005). Recommended effect size statistics for repeated
measures designs. *Behavior Research Methods, 37*(3), 379–384.
[doi:10.3758/BF03192707](https://doi.org/10.3758/BF03192707)

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

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

Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: Measures of effect size for some common research designs.
*Psychological Methods, 8*(4), 434–447.
[doi:10.1037/1082-989X.8.4.434](https://doi.org/10.1037/1082-989X.8.4.434)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
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
# 1. Raw-argument interface. Bargman's (1970) 5-group one-way ANOVA,
#        also used in Venables (1975), Fleishman (1980), and Steiger (2004):
#        11 subjects per group, observed F = 11.221.
eta_squared(F_value = 11.221, df_effect = 4, df_error = 50)
#>  effect  eta_squared F_value df_effect df_error
#>  overall 0.473       11.2    4         50      

# 2. One way ANOVA from a fitted model (depression_bdi: three
#        treatment arms, 10 per arm, N = 30).
fit_one <- aov(bdi_post ~ condition, data = depression_bdi)
eta_squared(fit_one)
#>  effect    eta_squared F_value df_effect df_error
#>  condition 0.184       3.03    2         27      

# 3. Factorial ANOVA: partial eta squared per effect (pygmalion
#        data: expectancy treatment x grade, 2 x 6 with unequal cell
#        sizes, N = 310). The treatment is manipulated; grade is a
#        measured classification of the pupils.
fit_factorial <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
eta_squared(fit_factorial)
#>  effect                  eta_squared F_value df_effect df_error
#>  treatment               0.0215      6.54    1         298     
#>  factor(grade)           0.0448      2.8     5         298     
#>  treatment:factor(grade) 0.0196      1.19    5         298     

# 4. Within-subjects (repeated measures) ANOVA. Simulated 20-subject
#        x 3-time design; each subject is measured at Pre, Mid, Post.
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
eta_squared(fit_rm)   # 'stratum' column identifies the within-subjects error
#>  effect eta_squared stratum      F_value df_effect df_error
#>  time   0.231       subject:time 5.7     2         38      
```
