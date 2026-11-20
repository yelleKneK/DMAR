# Omega Squared (Effect Size for ANOVA)

Computes the sample omega squared (\\\omega^2\\), Hays' (1994)
bias-corrected estimator of the proportion of the population variance of
the dependent variable that a fixed effect accounts for, as two
estimands: the *total* \\\omega^2\\, a proportion of the total variance,
and the *partial* \\\omega^2\\, a proportion of the variance due to the
effect and error alone. With `estimand = "all"`, the default, the
function reports every estimand its inputs determine; `"total"` or
`"partial"` reports one. In a one-way design the two are the same
number.

## Usage

``` r
omega_squared(
  object = NULL,
  estimand = c("all", "total", "partial"),
  SS_effect = NULL,
  SS_error = NULL,
  SS_total = NULL,
  MS_error = NULL,
  F_value = NULL,
  df_effect = NULL,
  df_error = NULL,
  N = NULL
)
```

## Arguments

- object:

  Optional. A fitted model of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (a multi-stratum
  fit such as `aov(y ~ A + Error(subject/A), data = d)`). The function
  returns one row per effect, with the stratum named for an `aovlist`.
  Give either `object` or summary statistics, not both.

- estimand:

  Which estimand to report: `"all"` (the default, every estimand the
  inputs determine, `NA` for the others), `"total"`, or `"partial"`. An
  estimand asked for by name that the inputs do not determine is refused
  (see Details).

- SS_effect:

  The sum of squares of the effect, or a vector of the sums of squares
  of several effects of one analysis (named for the effects, which label
  the rows). Given with `df_effect` and the error mean square.

- SS_error:

  The error (residual) sum of squares, given with `df_error` as an
  alternative to `MS_error`.

- SS_total:

  The corrected total sum of squares, \\\sum (y - \bar{y})^2\\ (see the
  Details of
  [`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)).
  With the sums of squares, it determines the total \\\omega^2\\.

- MS_error:

  The error (residual) mean square. Give either `MS_error` or `SS_error`
  with `df_error`.

- F_value:

  The observed *F* statistic of the effect, or a vector of the *F*
  statistics of several effects tested against the same error term. With
  `df_effect`, `df_error`, and `N`, it determines the partial
  \\\omega^2\\ only.

- df_effect:

  The numerator degrees of freedom of the effect (one value for each
  effect, or one for all). Required with `F_value` and with the sums of
  squares.

- df_error:

  The error (residual) degrees of freedom. Required with `F_value` and
  with `SS_error`.

- N:

  The total sample size, the number of observations. With the sums of
  squares or `F_value`, it determines the partial \\\omega^2\\. For an
  *F* from a design with an `Error()` stratum, give the size of the
  effect's stratum instead (see Details): the number of subjects for an
  effect tested between subjects, and for an effect tested within
  subjects the stratum's total degrees of freedom, which is
  `df_effect + df_error` when the effect is the only one tested in its
  stratum. `N` must be at least `df_effect + df_error`.

## Value

A `data.frame` (class `dmar_anova_es`, printed through
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md))
with one row per effect and the columns `effect`, `omega_squared_total`
and `omega_squared_partial` (only the one asked for when `estimand` is
`"total"` or `"partial"`), `stratum` (for an `aovlist` fit), `F_value`,
`df_effect`, `df_error`, and `N`. An effect given through the summary
interfaces without names is `"overall"` (or `effect_1`, `effect_2`, and
so on, for several), and a quantity the inputs do not give is `NA`. `N`
is the number of observations (for the summary interfaces, `N` as
given); the partial value of an effect of an `aovlist` fit is computed
with the size of its own stratum (see Details), so that its row's *F*,
degrees of freedom, and that size, given to the *F* interface, reproduce
it. The `estimand` attribute records the estimand asked for, and the
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
degrees of freedom and *N*. The confidence interval is computed by
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).

**Total and partial.** For a fixed effect with sum of squares
\\\mathit{SS}\_{\text{effect}}\\ on \\df\_{\text{effect}}\\ degrees of
freedom, error mean square \\\mathit{MS}\_{\text{error}}\\, corrected
total sum of squares \\\mathit{SS}\_{\text{total}}\\, and *N*
observations, \$\$\hat{\omega}^2\_{\text{total}} =
\frac{\mathit{SS}\_{\text{effect}} - df\_{\text{effect}}
\mathit{MS}\_{\text{error}}}{\mathit{SS}\_{\text{total}} +
\mathit{MS}\_{\text{error}}},\$\$ \$\$\hat{\omega}^2\_{\text{partial}} =
\frac{\mathit{SS}\_{\text{effect}} - df\_{\text{effect}}
\mathit{MS}\_{\text{error}}}{\mathit{SS}\_{\text{effect}} + (N -
df\_{\text{effect}}) \mathit{MS}\_{\text{error}}} =
\frac{df\_{\text{effect}} (F - 1)}{df\_{\text{effect}} (F - 1) + N}\$\$
(Olejnik & Algina, 2000, p. 266 and Table 10, p. 268; Maxwell, Delaney,
& Kelley, 2027, Eqs. 3.97, 3.98, and 7.44 to 7.46). The total value
estimates the effect's share of the total population variance, all other
effects and the error included; the partial value estimates its share of
the variance due to the effect and the error alone, the variance that
remains when the other effects are held constant. Which answers the
question depends on whether the design's other factors would vary where
the effect is studied (see the Details of
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)).
In a one-way design \\N - df\_{\text{effect}} = df\_{\text{error}} +
1\\, so the two denominators agree and the estimands are the same
number. Some sources print the total's denominator as
\\\mathit{SS}\_{\text{total}} - \mathit{MS}\_{\text{error}}\\ (Kim,
2016, Table 4, which gives .326 for an effect whose total
\\\hat{\omega}^2\\ is .323); the denominator of Hays' estimator adds the
error mean square, as here.

**Bias correction and negative estimates.** The sample \\\eta^2\\ is the
proportion of *sample* variance accounted for and overestimates its
population value. \\\hat{\omega}^2\\ subtracts \\df\_{\text{effect}}
\mathit{MS}\_{\text{error}}\\ from the effect's sum of squares and adds
to the denominator, which gives an estimator of the population
proportion with much smaller bias (Hays, 1994; Olejnik & Algina, 2000;
Kelley, 2007). Its numerator is negative whenever *F* \< 1, and a
negative estimate of either estimand is reported as 0, since \\\omega^2
\ge 0\\ by definition (Olejnik & Algina, 2000, p. 262; Maxwell, Delaney,
& Kelley, 2027, Section 3.9.2).

**What the inputs determine.** A fitted model with one error stratum
determines both estimands of every effect. Sums of squares determine the
total when `SS_total` is given and the partial when `N` is given, with
`df_effect` and the error mean square in either case. An *F* statistic
with its degrees of freedom and *N* determines the partial value only:
the total also depends on the sums of squares of the design's other
effects, and the function cannot tell from *F* whether there are any, so
the *F* interface reports the total as `NA` even for a one-way design,
where it equals the partial value. With `estimand = "all"` an estimand
the inputs do not determine is `NA`, with the reason in the
`estimand_notes` attribute and beneath the printed table; asked for by
name, it is refused with an error of class `dmar_estimand_undetermined`
(see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).
Sums of squares and *F* are two ways of giving the same effect, and a
call that gives both, or a fitted model with summary statistics, or both
`MS_error` and `SS_error`, is refused.

**Designs.** The total \\\omega^2\\ is defined for fixed effects in
designs whose factors are all between subjects, with one error term
(Olejnik & Algina, 2000, p. 266). For an `aovlist` fit, which has an
error term for each stratum, the total is `NA`; designs with random
factors and designs with factors within subjects have \\\omega^2\\
formulas of their own (Olejnik & Algina, 2000, Tables 11 and 12; Olejnik
& Algina, 2003), which this function does not compute. An interaction or
a covariate of an `lm` fit is an effect like any other, with its own sum
of squares and degrees of freedom. The sums of squares of a fitted model
are the sequential (Type I) ones of
[`anova()`](https://rdrr.io/r/stats/anova.html); see the Details of
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)
on the corrected total, on Types II and III, and on the coding of
factors.

**Within-subjects designs.** An effect of an `aovlist` fit is tested in
its own error stratum, and its partial estimate is computed with the
size of that stratum in place of *N*, the size with which
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
and
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
convert the limits of its interval: the stratum's total degrees of
freedom for an effect tested within subjects, and the number of subjects
for an effect in the between-subjects stratum of a mixed design (see the
Details of
[`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)).
Since \\df\_{\text{effect}} (F - 1)\\ estimates the noncentrality
parameter \\\lambda\\ of the effect's test,
\\\hat{\omega}^2\_{\text{partial}}\\ then estimates \\\lambda/(\lambda +
M)\\, with *M* the stratum's size, which is the value the sample partial
\\\eta^2\\ converges to as subjects are added and the quantity whose
interval `ci_eta_squared(estimand = "partial")` reports for the row.
Computed with the number of observations *N*, the estimate of a
within-subjects effect estimates \\\lambda/(\lambda + N)\\, a smaller
quantity, which the sample partial \\\eta^2\\ does not estimate. That
estimate is the partial \\\omega^2\\ that Olejnik and Algina (2003, p.
443) compute for a repeated measures factor, whose *N* is the total
number of observations, and the *F* interface given the number of rows
as `N` reproduces it. For the three IQ occasions of
[`pygmalion`](https://yelleknek.github.io/DMAR/reference/pygmalion.md),
*F*(2, 618) with 930 rows gives .150, below the partial \\\eta^2\\ of
.212 and its interval \[.157, .265\]; with the stratum's 2 + 618 = 620,
the estimate is .209.

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

Kim, H.-Y. (2016). Statistical notes for clinical researchers: Sample
size calculation 3. Comparison of several means using one-way ANOVA.
*Restorative Dentistry & Endodontics, 41*(3), 231–234.
[doi:10.5395/rde.2016.41.3.231](https://doi.org/10.5395/rde.2016.41.3.231)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\omega^2\\ and Chapter 7 on total
and partial \\\omega^2\\ in factorial designs.)

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

[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
for the confidence interval,
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`var_omega_squared`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md)

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
omega_squared(SS_effect = c(LIGHT = 34.716, MATERIAL = 1996.998,
                            "LIGHT:MATERIAL" = 1570.389),
              df_effect = c(1, 3, 3), SS_error = 2261.612,
              df_error = 72, SS_total = 5863.715, N = 80)
#>  effect         omega_squared_total omega_squared_partial F_value df_effect
#>  LIGHT          0.001               0.001                 1.105   1        
#>  MATERIAL       0.323               0.431                 21.192  3        
#>  LIGHT:MATERIAL 0.250               0.370                 16.665  3        
#>  df_error N 
#>  72       80
#>  72       80
#>  72       80

# 2. An F statistic determines the partial value only. Bargmann's
#    (1970) example, used in turn by Venables (1975), Fleishman (1980),
#    and Steiger (2004): a one-way ANOVA with five groups of 11 and an
#    observed F of 11.221 on 4 and 50 df. In a one-way design the
#    partial value is the total.
omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
              estimand = "partial")
#>  effect  omega_squared_partial F_value df_effect df_error N 
#>  overall 0.426                 11.221  4         50       55

# 3. One-way ANOVA from a fitted model (depression_bdi: three
#    treatment arms, 10 per arm, N = 30); the two estimands agree.
fit_one <- aov(bdi_post ~ condition, data = depression_bdi)
omega_squared(fit_one)
#>  effect    omega_squared_total omega_squared_partial F_value df_effect df_error
#>  condition 0.119               0.119                 3.035   2         27      
#>  N 
#>  30

# 4. Two-factor ANOVA for the manipulated expectancy treatment and the
#    measured grade classification (pygmalion data, unequal cell
#    sizes, N = 310). The treatment by grade interaction is weak here
#    (F = 1.19), so the additive model is used.
fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
omega_squared(fit_additive)
#>  effect        omega_squared_total omega_squared_partial F_value df_effect
#>  treatment     0.017               0.017                 6.519   1        
#>  factor(grade) 0.028               0.028                 2.787   5        
#>  df_error N  
#>  303      310
#>  303      310

# 5. The same partial values with their confidence intervals.
ci_omega_squared(fit_additive, estimand = "partial")
#>  effect        omega_squared_partial omega_squared_partial_ci_lower
#>  treatment     0.017                 0.001                         
#>  factor(grade) 0.028                 0.001                         
#>  omega_squared_partial_ci_upper F_value df_effect df_error N  
#>  0.062                          6.519   1         303      310
#>  0.082                          2.787   5         303      310
#> 
#> Confidence level: 95%
```
