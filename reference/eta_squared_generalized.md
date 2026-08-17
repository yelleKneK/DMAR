# Generalized Eta Squared (Effect Size for ANOVA, Comparable Across Designs)

Computes the sample generalized eta squared (\\\eta^2_G\\; Olejnik &
Algina, 2003; Bakeman, 2005), the proportion of variance in the
dependent variable accounted for by a fixed effect after the variance
attributable to *measured* (observed) factors, but not other
*manipulated* factors, has been left in the denominator. This makes
\\\eta^2_G\\ comparable across designs that differ in which factors are
present, which regular \\\eta^2\\ and partial \\\eta^2\\ are not.

## Usage

``` r
eta_squared_generalized(
  object = NULL,
  observed = NULL,
  SS_effect = NULL,
  SS_observed = NULL,
  SS_error = NULL,
  F_effect = NULL,
  df_effect = NULL,
  F_observed = NULL,
  df_observed = NULL,
  df_error = NULL
)
```

## Arguments

- object:

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (multi-stratum
  aov fit, e.g.\\ `aov(y ~ A + Error(subject/A), data = d)`).

- observed:

  Character vector naming the *measured* (rather than *manipulated*)
  factors. Their sums of squares, and the SS of every interaction
  containing a listed factor, are kept in the denominator of
  \\\eta^2_G\\ per Olejnik and Algina (2003, Eq. 5). Manipulated effects
  are excluded from the denominator when they are not the focal effect.
  For `aovlist` fits the names must match the effect labels visible in
  `summary(object)`.

- SS_effect:

  Sum of squares for the focal effect (option 2).

- SS_observed:

  Sums of squares for the measured factors. Scalar or numeric vector
  (option 2).

- SS_error:

  Error (residual) sum of squares (option 2).

- F_effect:

  Observed *F*-value for the focal effect (option 3).

- df_effect:

  Numerator degrees of freedom for the focal effect (option 3).

- F_observed:

  Vector of *F*-values for the measured factors (option 3).

- df_observed:

  Numerator degrees of freedom for the measured factors, aligned with
  `F_observed` (option 3).

- df_error:

  Error degrees of freedom (option 3).

## Value

A `data.frame` with one row per focal effect. With a single-stratum fit
and the raw interfaces the columns are `effect` and
`eta_squared_generalized`; `effect` is `"overall"` for the raw
interfaces. With an `aovlist` fit a `stratum` column is added,
identifying which error stratum each effect came from.

## Details

The function accepts *one* of three input interfaces:

1.  a fitted model object ([`aov`](https://rdrr.io/r/stats/aov.html),
    [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` for
    within-subjects / mixed designs) together with an `observed` vector
    listing which factors are measured (rather than manipulated);

2.  raw sums of squares: `SS_effect`, `SS_observed` (one value per
    measured factor, or a scalar), and `SS_error`; or

3.  raw *F*-values and degrees of freedom: `F_effect`, `df_effect`,
    `F_observed` (vector aligned with `df_observed`), and `df_error`.

If the user supplies both the SS interface (option 2) and the F/df
interface (option 3), the function computes \\\eta^2_G\\ from each and
compares the results to within a 1e-6 tolerance. When the two interfaces
agree, the SS value is returned. When they disagree, the function stops
with a detailed message reporting both values.

**Formula.** \$\$\hat{\eta}^2_G =
\frac{\mathit{SS}\_{\text{effect}}}{\mathit{SS}\_{\text{effect}} +
\sum\_\text{obs} \mathit{SS}\_{\text{measured}} +
\mathit{SS}\_{\text{error}}}.\$\$ For the F/df interface the equivalent
ratio form is used, dividing through by \\\mathit{SS}\_{\text{error}}\\
so that no total-N argument is required: \\\mathit{SS}\_i /
\mathit{SS}\_{\text{error}} = F_i \cdot df_i / df\_{\text{error}}\\.

**Designs supported.**

- *Between-subjects ANOVA (single stratum).* The function reads
  `anova(object)`. For each focal effect, the denominator is
  \\\mathit{SS}\_{\text{focal}} + \sum \mathit{SS}\_{\text{measured
  (others)}} + \mathit{SS}\_{\text{error}}\\. Manipulated factors that
  are not the focal effect contribute nothing to the denominator.

- *Within-subjects and mixed ANOVA (`aovlist`, multi-stratum).* The
  function reads `summary(object)` and walks every error stratum. For
  each focal effect, the denominator is \\\mathit{SS}\_{\text{focal}} +
  \sum \mathit{SS}\_{\text{measured (others)}} + \sum\_{s}
  \mathit{SS}\_{\text{error}(s)}\\, where the last sum runs over every
  error stratum (both the between-subjects "subjects" stratum and any
  within-subjects error strata). This is the Olejnik & Algina (2003) /
  Bakeman (2005) rule that makes the subject-level variance act as an
  "always-measured" contributor in repeated measures designs.

**Focal-effect self-exclusion.** If the focal effect itself is listed in
`observed`, the function excludes it from the observed-sum component
(the focal effect's *own* SS already appears in the numerator and the
leading term of the denominator). Practically this means listing every
effect as `observed` reduces to total \\\eta^2\\ for between-subjects
designs.

**Higher-order interactions.** An effect is a measured source of
variance when *any* factor in its term is measured (Olejnik & Algina,
2003, Eq. 5; Bakeman, 2005), so listing a factor in `observed` also
places every interaction containing that factor in the denominator
automatically. In a design with manipulated `A` and measured `c`,
`observed = "c"` therefore puts `c` and `A:c` in the denominator, which
is what the cited papers' worked examples do. An explicit interaction
label in `observed` is honored as given, declaring that one term
measured without marking its constituent factors.

**Covariates.** Under Olejnik and Algina's Eq. 5, a covariate is a
measured source whose SS always belongs in the denominator, so in an
ANCOVA list the covariate in `observed`. Because
[`anova()`](https://rdrr.io/r/stats/anova.html) on an `lm`/`aov` fit
uses sequential sums of squares, enter the covariate before the
treatment factors in the model formula so its SS is adjusted the way the
ANCOVA decomposition intends.

**Confidence intervals.** See
[`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md)
for the corresponding CI function; both available CI methods are
approximate and require independent evaluation.

## References

Bakeman, R. (2005). Recommended effect size statistics for repeated
measures designs. *Behavior Research Methods, 37*(3), 379–384.
[doi:10.3758/BF03192707](https://doi.org/10.3758/BF03192707)

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

## See also

[`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_partial`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
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
# 1. Fitted model with `observed`. In the pygmalion expectancy
#        experiment, treatment is manipulated while grade is a measured
#        classification, so grade's variance stays in the denominator.
pyg <- pygmalion
pyg$grade <- factor(pyg$grade)
fit <- aov(iq_8 ~ treatment * grade, data = pyg)
eta_squared_generalized(fit, observed = "grade")
#>            effect eta_squared_generalized
#> 1       treatment              0.02015391
#> 2           grade              0.04396181
#> 3 treatment:grade              0.01872564

# 2. Raw sums of squares.
eta_squared_generalized(SS_effect = 100, SS_observed = c(40, 30),
                        SS_error = 200)
#>    effect eta_squared_generalized
#> 1 overall               0.2702703

# 3. Raw F-values and degrees of freedom.
eta_squared_generalized(F_effect = 6.0, df_effect = 2,
                        F_observed = c(2.5, 1.8),
                        df_observed = c(1, 2), df_error = 50)
#>    effect eta_squared_generalized
#> 1 overall               0.1762115

# 4. Within-subjects (repeated measures) ANOVA. The denominator
#        automatically includes the subject-level variance plus the
#        within-subjects error, following Bakeman (2005).
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
eta_squared_generalized(fit_rm)
#>   effect eta_squared_generalized      stratum
#> 1   time              0.09652627 subject:time

# 5. Mixed design with a measured between-subjects factor. Treat
#        'group' as observed; its SS stays in the denominator for
#        'time' and the 'group:time' interaction.
set.seed(113)
n_per_group <- 10
n <- n_per_group * 2
mixed_data <- data.frame(
  subject = factor(rep(seq_len(n), each = 3)),
  group   = factor(rep(c("Treatment", "Control"), each = 3 * n_per_group)),
  time    = factor(rep(c("Pre", "Mid", "Post"), n),
                   levels = c("Pre", "Mid", "Post")),
  y       = rnorm(n, sd = 1)[rep(seq_len(n), each = 3)] +
            0.5 * rep(1:3, n) + rnorm(n * 3, sd = 1)
)
fit_mixed <- aov(y ~ group * time + Error(subject/time), data = mixed_data)
eta_squared_generalized(fit_mixed, observed = "group")
#>       effect eta_squared_generalized      stratum
#> 1      group              0.01108528      subject
#> 2       time              0.09101617 subject:time
#> 3 group:time              0.01120063 subject:time
```
