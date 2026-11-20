# Eta Squared (Effect Size for ANOVA)

Computes the sample eta squared (\\\eta^2\\) of each effect of an
analysis of variance, the proportion of the variance of the dependent
variable that the effect accounts for, as two estimands: the *total*
\\\eta^2\\, the effect's sum of squares as a proportion of the corrected
total sum of squares, and the *partial* \\\eta^2\\, its sum of squares
as a proportion of its own sum of squares and the error sum of squares.
With `estimand = "all"`, the default, the function reports every
estimand its inputs determine; `"total"` or `"partial"` reports one. In
a one-way design the two are the same number.

## Usage

``` r
eta_squared(
  object = NULL,
  estimand = c("all", "total", "partial"),
  SS_effect = NULL,
  SS_error = NULL,
  SS_total = NULL,
  F_value = NULL,
  df_effect = NULL,
  df_error = NULL
)
```

## Arguments

- object:

  Optional. A fitted model of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (a multi-stratum
  fit such as `aov(y ~ A + Error(subject/A), data = d)`). The function
  returns one row per effect (every row of `anova(object)`, or of each
  stratum of `summary(object)`, other than `Residuals`), with the
  stratum named for an `aovlist`. Give either `object` or summary
  statistics, not both.

- estimand:

  Which estimand to report: `"all"` (the default, every estimand the
  inputs determine, `NA` for the others), `"total"`, or `"partial"`. An
  estimand asked for by name that the inputs do not determine is refused
  (see Details).

- SS_effect:

  The sum of squares of the effect, or a vector of the sums of squares
  of several effects of one analysis (named for the effects, which label
  the rows).

- SS_error:

  The error (residual) sum of squares. With `SS_effect`, it determines
  the partial \\\eta^2\\.

- SS_total:

  The corrected total sum of squares, \\\sum (y - \bar{y})^2\\ (see
  Details). With `SS_effect`, it determines the total \\\eta^2\\.

- F_value:

  The observed *F* statistic of the effect, or a vector of the *F*
  statistics of several effects tested against the same error term. With
  `df_effect` and `df_error`, it determines the partial \\\eta^2\\ only.

- df_effect:

  The numerator degrees of freedom of the effect (one value for each
  effect, or one for all). Required with `F_value`; with sums of squares
  it is optional and, with `df_error` and `SS_error`, fills the
  `F_value` column.

- df_error:

  The error (residual) degrees of freedom. Required with `F_value`;
  optional with sums of squares.

## Value

A `data.frame` (class `dmar_anova_es`, printed through
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md))
with one row per effect and the columns `effect`, `eta_squared_total`
and `eta_squared_partial` (only the one asked for when `estimand` is
`"total"` or `"partial"`), `stratum` (for an `aovlist` fit), `F_value`,
`df_effect`, and `df_error`. An effect given through the summary
interfaces without names is `"overall"` (or `effect_1`, `effect_2`, and
so on, for several), and a quantity the inputs do not give (the *F*
statistic of sums of squares given without degrees of freedom) is `NA`.
The `estimand` attribute records the estimand asked for, and the
`estimand_notes` attribute names each estimand column that is `NA` with
the reason, which the line beneath the printed table states.
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives one row
per effect and estimand (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

The function takes a fitted model
([`aov`](https://rdrr.io/r/stats/aov.html),
[`lm`](https://rdrr.io/r/stats/lm.html), or an `aovlist` from
[`aov()`](https://rdrr.io/r/stats/aov.html) with an `Error()` term), the
sums of squares of a published ANOVA table, or an *F* statistic with its
degrees of freedom. The confidence interval is computed by
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md).

**Total and partial.** For an effect with sum of squares
\\\mathit{SS}\_{\text{effect}}\\, error sum of squares
\\\mathit{SS}\_{\text{error}}\\, and corrected total sum of squares
\\\mathit{SS}\_{\text{total}}\\, \$\$\hat{\eta}^2\_{\text{total}} =
\frac{\mathit{SS}\_{\text{effect}}}{\mathit{SS}\_{\text{total}}}, \qquad
\hat{\eta}^2\_{\text{partial}} =
\frac{\mathit{SS}\_{\text{effect}}}{\mathit{SS}\_{\text{effect}} +
\mathit{SS}\_{\text{error}}} = \frac{df\_{\text{effect}}
F}{df\_{\text{effect}} F + df\_{\text{error}}}\$\$ (Olejnik & Algina,
2000, Table 10, p. 268; Maxwell, Delaney, & Kelley, 2027, Eqs. 7.40 to
7.42). The total value is a proportion of all the variation in the data,
and so depends on the other factors of the design: adding a factor with
a large effect lowers it. The partial value leaves the other effects out
of its denominator, so the partial values of the effects of one design
do not share a denominator and can sum to more than 1 (Olejnik & Algina,
2000, pp. 268–269). Maxwell, Delaney, and Kelley (2027, Chapter 7) frame
the choice by the other factors: when another factor is extrinsic to the
question (a study of the effect could as well hold it constant, as a
study of drugs could hold biofeedback constant) the partial value
answers it, and when the factor is intrinsic (it varies wherever the
effect is studied, as gender does) the total value does. In a one-way
design \\\mathit{SS}\_{\text{total}} = \mathit{SS}\_{\text{effect}} +
\mathit{SS}\_{\text{error}}\\, and the two estimands are the same
number.

**What the inputs determine.** A fitted model determines both estimands
of every effect, and so do sums of squares that include both `SS_error`
and `SS_total`. `SS_effect` with `SS_total` alone determines the total,
and with `SS_error` alone the partial. An *F* statistic and its degrees
of freedom determine the partial value only: the total of an effect also
depends on the sums of squares of the design's other effects (from the
*F* statistics of every effect it can be written as
\\df\_{\text{effect}} F / (\sum df F + df\_{\text{error}})\\, Maxwell,
Delaney, & Kelley, 2027, Eq. 7.43, which a fitted model provides). The
function cannot tell from *F* whether the design has other effects, so
the *F* interface reports the total as `NA` even for a one-way design,
where it equals the partial value. With `estimand = "all"` an estimand
the inputs do not determine is `NA`, with the reason in the
`estimand_notes` attribute and beneath the printed table; asked for by
name, it is refused with an error of class `dmar_estimand_undetermined`
(see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).
Sums of squares and *F* are two ways of giving the same effect, and a
call that gives both, or a fitted model with summary statistics, is
refused, as is a set of sums of squares in which `SS_effect + SS_error`
exceeds `SS_total`.

**The corrected total.** `SS_total` is the corrected total sum of
squares, \\\sum (y - \bar{y})^2\\, the sum of every sum of squares of an
ANOVA table with an intercept, the error included. Some programs also
print an uncorrected total, \\\sum y^2\\, which includes the grand mean:
Kim (2016, Table 3) lists both, a “Total” of 53758.357 and a “Corrected
total” of 5863.715, and the second is `SS_total`. A fitted model without
an intercept does not partition the corrected total, and its total
\\\eta^2\\ is `NA`.

**Sums of squares of a fitted model.** The sums of squares and *F*
statistics of a fitted model are those of
[`anova()`](https://rdrr.io/r/stats/anova.html) (one stratum) or
[`summary()`](https://rdrr.io/r/base/summary.html) (an `aovlist`), which
are Type I (sequential): each effect is adjusted for the effects before
it in the formula. Sequential sums of squares add up to the corrected
total, so the total \\\eta^2\\ of every effect and the error's share sum
to 1. In a balanced design Types I, II, and III agree; in an unbalanced
one they differ, and the order of the terms changes the Type I values.
For Type II or III sums of squares, compute them with, for example,
`car::Anova(object, type = 3)`, and give the effect's sum of squares
with the error and corrected total sums of squares to the
sums-of-squares interface (Type III sums of squares need not add up to
the total in an unbalanced design, so the total values of the effects
then need not sum with the error's share to 1). Type III tests depend on
how the factors are coded. The Type III main-effect tests compare
unweighted marginal means only when the factors are coded with contrasts
that sum to zero over their levels, such as
[`contr.sum`](https://rdrr.io/r/stats/contrast.html) (`contr.helmert`
and `contr.poly` also qualify). Under R's default treatment contrasts,
the Type III test of a factor that enters an interaction compares its
levels only at the reference levels of the factors it interacts with,
which is seldom the hypothesis intended. Set the coding before the model
is fitted, either in the call
(`lm(y ~ A * B, data = d, contrasts = list(A = contr.sum, B = contr.sum))`)
or for the session
(`options(contrasts = c("contr.sum", "contr.poly"))`). The help page of
[`car::Anova()`](https://rdrr.io/pkg/car/man/Anova.html) gives the same
warning, and notes that when a model has numeric covariates and
interactions, the Type III test of a factor compares the groups where
the covariates are 0; Type II tests do not depend on the coding.
[`factorial_anova`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md)
applies sum-to-zero contrasts itself and reports partial \\\eta^2\\ and
partial \\\omega^2\\ with their confidence intervals for
between-subjects factorial designs.

**Interactions and covariates.** An interaction is an effect like any
other, with its own sum of squares and degrees of freedom, and so is a
covariate of an `lm` fit: each is a row, its partial value computed
against the error, its total against the corrected total.

**Within-subjects designs.** For an `aovlist` fit, each effect's partial
\\\eta^2\\ is computed against the error of its own stratum, so a
within-subjects factor's partial value leaves the differences among
subjects out of its denominator. Its total value divides by the
corrected total over all strata (the sum of every sum of squares that
`summary(object)` reports), which includes the differences among
subjects: in a design with one within-subjects factor the total value is
the factor's sum of squares over the sum of those for subjects, the
factor, and their interaction (Bakeman, 2005, p. 380), and is smaller
than the partial value. Neither is comparable with the same effect in a
between-subjects design; see
[`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md)
for an index that is (Olejnik & Algina, 2003; Bakeman, 2005). Other
model classes (`lmerMod`, `lme`) are not taken; give their sums of
squares or *F* statistics to the summary interfaces.

## References

Bakeman, R. (2005). Recommended effect size statistics for repeated
measures designs. *Behavior Research Methods, 37*(3), 379–384.
[doi:10.3758/BF03192707](https://doi.org/10.3758/BF03192707)

Bargmann, R. E. (1970). Interpretation and use of a generalized
discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), *Essays in probability
and statistics* (pp. 35–60). University of North Carolina Press.

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.
[doi:10.1177/001316448004000309](https://doi.org/10.1177/001316448004000309)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Kim, H.-Y. (2016). Statistical notes for clinical researchers: Sample
size calculation 3. Comparison of several means using one-way ANOVA.
*Restorative Dentistry & Endodontics, 41*(3), 231–234.
[doi:10.5395/rde.2016.41.3.231](https://doi.org/10.5395/rde.2016.41.3.231)

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

Venables, W. (1975). Calculation of confidence intervals for
noncentrality parameters. *Journal of the Royal Statistical Society,
Series B, 37*(3), 406–412.
[doi:10.1111/j.2517-6161.1975.tb01554.x](https://doi.org/10.1111/j.2517-6161.1975.tb01554.x)

## See also

[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
for the confidence interval,
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
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
# 1. The sums of squares of a published two-way table: Kim (2016,
#    Table 3), light by material, 2 by 4 with 10 per cell (N = 80).
#    SS_total is the corrected total, 5863.715.
eta_squared(SS_effect = c(LIGHT = 34.716, MATERIAL = 1996.998,
                          "LIGHT:MATERIAL" = 1570.389),
            SS_error = 2261.612, SS_total = 5863.715,
            df_effect = c(1, 3, 3), df_error = 72)
#>  effect         eta_squared_total eta_squared_partial F_value df_effect
#>  LIGHT          0.006             0.015               1.105   1        
#>  MATERIAL       0.341             0.469               21.192  3        
#>  LIGHT:MATERIAL 0.268             0.410               16.665  3        
#>  df_error
#>  72      
#>  72      
#>  72      

# 2. An F statistic and its degrees of freedom determine the partial
#    value only. Bargmann's (1970) example, used in turn by Venables
#    (1975), Fleishman (1980), and Steiger (2004): a one-way ANOVA with
#    five groups of 11 and an observed F of 11.221 on 4 and 50 df. In a
#    one-way design the partial value is the total.
eta_squared(F_value = 11.221, df_effect = 4, df_error = 50,
            estimand = "partial")
#>  effect  eta_squared_partial F_value df_effect df_error
#>  overall 0.473               11.221  4         50      

# 3. One-way ANOVA from a fitted model (depression_bdi: three
#    treatment arms, 10 per arm, N = 30); the two estimands agree.
fit_one <- aov(bdi_post ~ condition, data = depression_bdi)
eta_squared(fit_one)
#>  effect    eta_squared_total eta_squared_partial F_value df_effect df_error
#>  condition 0.184             0.184               3.035   2         27      

# 4. Factorial ANOVA (pygmalion: expectancy treatment by grade, 2 by 6
#    with unequal cell sizes, N = 310). The treatment is manipulated;
#    grade is a measured classification of the pupils.
fit_factorial <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
eta_squared(fit_factorial)
#>  effect                  eta_squared_total eta_squared_partial F_value
#>  treatment               0.020             0.021               6.539  
#>  factor(grade)           0.043             0.045               2.795  
#>  treatment:factor(grade) 0.018             0.020               1.191  
#>  df_effect df_error
#>  1         298     
#>  5         298     
#>  5         298     

# 5. Within-subjects ANOVA: 20 simulated subjects, each measured at Pre,
#    Mid, and Post. The total value includes the differences among
#    subjects in its denominator; the partial value does not.
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
eta_squared(fit_rm)
#>  effect eta_squared_total eta_squared_partial stratum      F_value df_effect
#>  time   0.097             0.231               subject:time 5.701   2        
#>  df_error
#>  38      
```
