# Confidence Interval for Eta Squared (Effect Size for ANOVA)

Computes the point estimate of eta squared (\\\eta^2\\) of each effect
of an analysis of variance and an exact, noncentrality-based confidence
interval for its population value, for the *partial* and the *total*
estimand (see
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)).
With `estimand = "all"`, the default, the function reports both; the
interval from an effect's *F* statistic is an interval for the partial
value, and the total has one only where it equals the partial value, in
a one-way design (see Details). `"total"` or `"partial"` reports one
estimand.

## Usage

``` r
ci_eta_squared(
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

  Optional. A fitted model object of class
  [`aov`](https://rdrr.io/r/stats/aov.html),
  [`lm`](https://rdrr.io/r/stats/lm.html), or `aovlist` (multi-stratum
  aov fit, e.g.\\ `aov(y ~ A + Error(subject/A), data = d)`). For
  multi-stratum fits the function walks every error stratum and returns
  one row per non-`Residuals` effect, identifying the stratum used. Give
  either `object` or the *F* interface, not both.

- estimand:

  Which estimand to report: `"all"` (the default, both, with `NA` where
  the inputs determine no estimate or no interval), `"total"`, or
  `"partial"`. An estimand asked for by name for which the inputs
  determine no interval is refused (see Details).

- F_value:

  Observed *F*-value, or a vector of the *F* statistics of several
  effects tested against the same error term.

- df_effect:

  Numerator degrees of freedom for the effect (one value for each
  effect, or one for all).

- df_error:

  Error (residual) degrees of freedom.

- N:

  Total sample size, the total number of observations (derived
  automatically from a fitted model, via
  [`nobs`](https://rdrr.io/r/stats/nobs.html)`(object)` for
  single-stratum fits, or, for `aovlist` fits, recovered as one more
  than the sum of the effect and residual degrees of freedom across
  every stratum, which is the total number of observations). The
  confidence limits for the noncentrality parameter are converted into
  limits for \\\eta^2\\ with `N`. For an *F* from a design with an
  `Error()` stratum, give the size of the effect's stratum instead (see
  Details): the number of subjects for an effect tested between
  subjects, and for an effect tested within subjects the stratum's total
  degrees of freedom, which is `df_effect + df_error` when the effect is
  the only one tested in its stratum, not the number of rows of the
  data. `N` must be at least `df_effect + df_error`.

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
reported, its estimate and limits, `eta_squared_total`,
`eta_squared_total_ci_lower`, `eta_squared_total_ci_upper`,
`eta_squared_partial`, `eta_squared_partial_ci_lower`, and
`eta_squared_partial_ci_upper`; `stratum` (for an `aovlist` fit); and
`F_value`, `df_effect`, `df_error`, and `N`. With the *F* interface
`effect` is `"overall"` (the names of `F_value`, when it has them, or
`effect_1`, `effect_2`, and so on, for several). `N` is the number of
observations (for the *F* interface, `N` as supplied). The limits of the
*F* interface and of a single-stratum fit are converted with that `N`;
the limits of an effect of an `aovlist` fit are converted with the size
of its own stratum (see Details), so that its row's *F*, degrees of
freedom, and that size, passed to the *F* interface, reproduce them. The
size is not a column of the table: it is the sum of the `Df` column of
the stratum's table in `summary(object)`, plus one for the stratum that
lies between subjects (`id` in `Error(id / occasion)`). The table
records the confidence level of its intervals in its `conf_level`
attribute (`conf_level`, or `1 - alpha_lower - alpha_upper` when the
tail areas are supplied), which the console prints beneath it and
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
reports; the `estimand` attribute records the estimand asked for, and
the `estimand_notes` attribute names each estimate or interval that is
`NA` with the reason, which the console also prints.
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) gives one row
per effect and estimand, with the limits in `ci_lower` and `ci_upper`
(see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

## Details

The function takes a fitted model
([`aov`](https://rdrr.io/r/stats/aov.html),
[`lm`](https://rdrr.io/r/stats/lm.html), or an `aovlist` from
[`aov()`](https://rdrr.io/r/stats/aov.html) with an `Error()` term) or
an *F* statistic with its degrees of freedom and *N*. For factorial and
within-subjects designs it returns one row per effect, the partial value
computed against that effect's own error stratum.

**Point estimates.** The partial value is
\\\hat{\eta}^2\_{\text{partial}} = df\_{\text{effect}} F /
(df\_{\text{effect}} F + df\_{\text{error}})\\, which equals
\\\mathit{SS}\_{\text{effect}}/(\mathit{SS}\_{\text{effect}} +
\mathit{SS}\_{\text{error}})\\; the total value of a fitted model is
\\\mathit{SS}\_{\text{effect}}/\mathit{SS}\_{\text{total}}\\, with the
corrected total sum of squares. They are the estimates of
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
whose Details describe the two estimands and the sums of squares of a
fitted model. An *F* statistic does not determine the total value, which
the *F* interface reports as `NA`.

**Confidence interval.** The interval rests on the two principles
Steiger (2004, pp. 167–168) sets out. A confidence interval for the
noncentrality parameter \\\lambda\\ of the effect's *F* test is found by
inverting the noncentral *F* distribution (via
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md)), and
because \\\eta^2\\ is an increasing function of \\\lambda\\, the limits
for \\\lambda\\ map into limits for \\\eta^2\\:
\$\$\eta^2\_{\text{bound}} =
\frac{\lambda\_{\text{bound}}}{\lambda\_{\text{bound}} + N},\$\$ where
*N* is the number of observations for an effect tested between subjects
and the size of the effect's error stratum for an effect of an `aovlist`
fit (see below). When the observed *F* is below the critical value that
sets the lower limit, the lower limit for \\\lambda\\, and so for
\\\eta^2\\, is 0. In a one-way design the interval is the one
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md) and
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
give, the population \\\eta^2\\ and \\\omega^2\\ being the same
quantity; the functions differ there only in their point estimates.

**What the interval estimates: the partial value.** The population value
of an effect's partial \\\eta^2\\ is the value its sample counterpart
converges to as the sample grows, which is how Steiger (2004, Eq. 35, p.
171) defines the population \\\rho^2\\ of a one-way design. For an
effect tested between subjects, in a one-way or a factorial design,
\\\lambda = N \sigma^2\_{\text{effect}} / \sigma^2_e\\, with
\\\sigma^2\_{\text{effect}}\\ the variance of the effect over all *N*
observations and \\\sigma^2_e\\ the within-cell variance, so the
population partial \\\eta^2\\, \\\sigma^2\_{\text{effect}} /
(\sigma^2\_{\text{effect}} + \sigma^2_e)\\, equals \\\lambda/(\lambda +
N)\\ at every sample size (Steiger, 2004, Eqs. 13, 15, 16, and 19, pp.
169–170, where it is written \\\omega^2\\ and described as the
proportion of the variance that remains after the other effects have
been partialed out; Kelley, 2007, Eqs. 44, 45, and 50, p. 12). The
interval converted with *N* is therefore exact for it, and each effect
of a fitted `aov` or `lm` model is converted with the number of
observations.

**The total value.** The population total \\\eta^2\\ of an effect,
\\\sigma^2\_{\text{effect}}\\ as a proportion of the variance due to
every effect of the design and the error, depends on the noncentrality
parameters of the other effects as well as on the effect's own, and the
distribution of one *F* statistic carries no information about them. The
intervals of Steiger (2004, pp. 169–171), Smithson (2003, pp. 43–44,
where the correlation ratio of a factorial effect is “a squared partial
correlation”), and Fidler and Thompson (2001, Table 8) for the effects
of a fixed-effects factorial design are intervals for the partial value.
In a one-way design the effect and the error are the whole design, the
total and the partial value are the same quantity, and the interval of a
fitted one-way model is reported for both. For an effect of a design
with other effects, and for every effect of an `aovlist` fit (whose
total includes the differences among subjects), no interval is formed
for the total: with `estimand = "all"` its estimate is reported with
`NA` limits, and `estimand = "total"` is refused with an error of class
`dmar_estimand_undetermined` (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).
The *F* interface determines neither the total value nor its interval,
since *F* does not say whether the design has other effects; for a
one-way design give the fitted model, or use `estimand = "partial"`,
whose interval is then the total's too.

**A conversion with the degrees of freedom.** Smithson (2001, Eqs. 8 and
9, p. 618; 2003, Formulas 5.3 and 5.4, p. 43) defines the population
\\\eta^2\\ of a factorial effect in the same way, through the effect's
variance over all *N* cases, but converts the limits with
\\\lambda/(\lambda + u + v + 1)\\, *u* and *v* being the effect and
error degrees of freedom (2001, Eq. 11; 2003, Formula 5.6). That form
comes from regression, where \\u + v + 1\\ is the sample size (2001, Eq.
4, p. 612), and it agrees with his definition only when \\u + v + 1 =
N\\, as in a one-way design. For a factorial effect it is a convention
for the conversion rather than a definition of the parameter: its limits
are exact for \\\lambda/(\lambda + u + v + 1)\\, a larger quantity that
depends on the degrees of freedom of the design's other effects. For the
2 by 7 design with four observations per cell that both Steiger (2004,
p. 171) and Smithson (2001, p. 619) analyze, the 90% limits for the
interaction are \[.008, .316\] with *N* = 56 and \[.009, .346\] with 6 +
42 + 1 = 49, and those for the row effect are \[.009, .234\] with *N* =
56 and \[.012, .280\] with 1 + 42 + 1 = 44. Converted with 44, a 95%
interval for that row effect covers the population value 94.5% of the
time when the value is .02 and 91.9% of the time when it is .40. Fidler
and Thompson (2001, Table 8, p. 596) computed their fixed-effects
factorial intervals with Smithson's syntax, from the *F* values and
their degrees of freedom alone. Supplying \\u + v + 1\\ as `N` with
`estimand = "partial"` reproduces such published limits.

**Within-subjects designs.** An effect of an `aovlist` fit is tested in
its own error stratum, and the number of rows of the data is not a
sample size for that test. Under compound symmetry, Steiger (2004, Eq.
58, p. 177) gives the noncentrality parameter of a one-way
within-subjects design with *n* subjects and *p* occasions as \\\lambda
= n \sum_j (\mu_j - \mu)^2 / \[\sigma^2 (1 - \rho)\]\\. The sums of
squares for the occasions and for their error grow as \\n \sum_j
(\mu_j - \mu)^2\\ and \\n (p - 1) \sigma^2 (1 - \rho)\\, so the value
the sample partial \\\eta^2\\ converges to is \\\lambda/\[\lambda +
n(p - 1)\]\\, where \\n(p - 1)\\ is the stratum's total degrees of
freedom, \\df\_{\text{effect}} + df\_{\text{error}}\\. Under sphericity
the same holds for every effect of a within-subjects stratum, with the
stratum's total degrees of freedom (those of its effects and its error)
in place of \\n(p - 1)\\, and each such effect is converted with that
total. An effect in the between-subjects stratum of a mixed design is
tested among the subjects' means, as an effect of a between-subjects
design is tested among its observations, and is converted with the
number of subjects, one more than that stratum's total degrees of
freedom. Converted with the number of rows, \\np\\, the interval of a
within-subjects effect is for \\\lambda/(\lambda + np)\\, a smaller
quantity, which the sample partial \\\eta^2\\ does not estimate. In the
three-occasion analysis of
[`pygmalion`](https://yelleknek.github.io/DMAR/reference/pygmalion.md)
in
[`vignette("DMAR", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/DMAR.md),
*F*(2, 618) converted with the 930 rows gives the 95% interval \[.111,
.194\], which lies below the estimate of .212; converted with 2 + 618 =
620, the interval is \[.157, .265\].

The literature cited here does not settle the within-subjects case.
Steiger (2004, p. 176) restricts his procedures to completely
between-subjects designs, noting that the noncentrality parameter of a
repeated measures design confounds the effects of the treatments with
the correlation among the observations (pp. 176–177), and the factorial
examples of Smithson (2001, 2003) and Fidler and Thompson (2001) are
between-subjects designs. The conversion of a within-subjects stratum
given here follows from Steiger's Eq. 58 and rests, as the stratum's *F*
test does, on sphericity. A within-subjects partial \\\eta^2\\ also
leaves the differences among subjects out of its denominator, so it is
not comparable with a partial \\\eta^2\\ from a between-subjects design
(Bakeman, 2005, p. 380).

**The *F* interface.** The limits are converted with `N` as supplied.
The function cannot tell a between-subjects factorial *F*, for which the
total sample size properly exceeds \\df\_{\text{effect}} +
df\_{\text{error}} + 1\\, from a within-subjects *F* given with the
number of rows, so it does not warn; for an *F* tested within subjects,
supply the size of the effect's stratum as `N`.

**Designs supported.**

- *Between-subjects ANOVA*: fitted `aov`/`lm`, or *F*, df, and *N*. Each
  effect's limits are converted with the number of observations.

- *Within-subjects or mixed ANOVA*: fitted `aovlist`. Each effect's CI
  is built from its own stratum's *F* and residual *df* and converted
  with the size of that stratum: the number of subjects for the
  between-subjects stratum (the stratum of the outermost term of
  `Error()`, such as `subject` in `Error(subject/time)`) and the total
  degrees of freedom of a within-subjects stratum. The `N` column
  reports the total number of observations across all strata, and the
  `stratum` column identifies which error term each row used.

**Sums of squares in factorial designs.** When a fitted model is
supplied, *F*-values are read from
[`anova()`](https://rdrr.io/r/stats/anova.html) (single-stratum) or
[`summary()`](https://rdrr.io/r/base/summary.html) (multi-stratum), both
of which use Type I (sequential) sums of squares in base R. For balanced
designs Types I, II, and III agree; for unbalanced designs they differ.
If Type II or III *F*-values are required, compute them with e.g.\\
`car::Anova(object, type = 3)` and pass the relevant *F*, degrees of
freedom, and *N* into the *F* interface. See the Details of
[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md)
on the coding of factors that Type III tests depend on.
[`factorial_anova`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md)
applies sum-to-zero contrasts itself and reports partial \\\eta^2\\ and
partial \\\omega^2\\ with their confidence intervals for
between-subjects factorial designs.

## References

Bakeman, R. (2005). Recommended effect size statistics for repeated
measures designs. *Behavior Research Methods, 37*(3), 379–384.
[doi:10.3758/BF03192707](https://doi.org/10.3758/BF03192707)

Bargmann, R. E. (1970). Interpretation and use of a generalized
discriminant function. In R. C. Bose, I. M. Chakravarti, P. C.
Mahalanobis, C. R. Rao, & K. J. C. Smith (Eds.), *Essays in probability
and statistics* (pp. 35–60). University of North Carolina Press.

Fidler, F., & Thompson, B. (2001). Computing correct confidence
intervals for ANOVA fixed- and random-effects effect sizes. *Educational
and Psychological Measurement, 61*(4), 575–604.
[doi:10.1177/0013164401614003](https://doi.org/10.1177/0013164401614003)

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

Smithson, M. (2003). *Confidence intervals*. Thousand Oaks, CA: Sage
Publications.
[doi:10.4135/9781412983761](https://doi.org/10.4135/9781412983761)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

Venables, W. (1975). Calculation of confidence intervals for
noncentrality parameters. *Journal of the Royal Statistical Society,
Series B, 37*(3), 406–412.
[doi:10.1111/j.2517-6161.1975.tb01554.x](https://doi.org/10.1111/j.2517-6161.1975.tb01554.x)

## See also

[`eta_squared`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_nc_F`](https://yelleknek.github.io/DMAR/reference/ci_nc_F.md),
[`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
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
#    50 df, at 90%. F determines the partial value, which in a one-way
#    design is the total.
ci_eta_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55,
               estimand = "partial", conf_level = 0.90)
#>  effect  eta_squared_partial eta_squared_partial_ci_lower
#>  overall 0.473               0.261                       
#>  eta_squared_partial_ci_upper F_value df_effect df_error N 
#>  0.565                        11.221  4         50       55
#> 
#> Confidence level: 90%

# 2. One-way ANOVA from a fitted model: mean IQ gain differs across
#    the six grades of the pygmalion data (N = 310). The total and the
#    partial value coincide, and so do their intervals.
fit_one <- aov(iq_gain ~ factor(grade), data = pygmalion)
ci_eta_squared(fit_one)
#>  effect        eta_squared_total eta_squared_total_ci_lower
#>  factor(grade) 0.118             0.048                     
#>  eta_squared_total_ci_upper eta_squared_partial eta_squared_partial_ci_lower
#>  0.176                      0.118               0.048                       
#>  eta_squared_partial_ci_upper F_value df_effect df_error N  
#>  0.176                        8.159   5         304      310
#> 
#> Confidence level: 95%

# 3. Two-factor ANOVA for the manipulated expectancy treatment and the
#    measured grade classification (pygmalion data, N = 310). The
#    treatment by grade interaction is weak here (F = 1.19), so the
#    additive model is used. The totals are reported without limits.
fit_additive <- aov(iq_8 ~ treatment + factor(grade), data = pygmalion)
ci_eta_squared(fit_additive)
#>  effect        eta_squared_total eta_squared_total_ci_lower
#>  treatment     0.020             <NA>                      
#>  factor(grade) 0.043             <NA>                      
#>  eta_squared_total_ci_upper eta_squared_partial eta_squared_partial_ci_lower
#>  <NA>                       0.021               0.001                       
#>  <NA>                       0.044               0.001                       
#>  eta_squared_partial_ci_upper F_value df_effect df_error N  
#>  0.062                        6.519   1         303      310
#>  0.082                        2.787   5         303      310
#> 
#> Confidence level: 95%
#> 
#> NA eta_squared_total limits: No confidence interval is formed for the total eta squared of an effect in a design with other effects or with an Error() term: the noncentral F interval from the effect's F is an interval for the partial value, and the total depends as well on the other effects of the design (see the Details of ?ci_eta_squared).

# 4. Within-subjects ANOVA. The partial value's interval is computed
#    against the within-subjects error stratum, named in 'stratum'.
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
ci_eta_squared(fit_rm, estimand = "partial")
#>  effect eta_squared_partial eta_squared_partial_ci_lower
#>  time   0.231               0.023                       
#>  eta_squared_partial_ci_upper stratum      F_value df_effect df_error N 
#>  0.417                        subject:time 5.701   2         38       60
#> 
#> Confidence level: 95%

# 5. The size with which the limits are converted. The interaction of
#    the 2 by 7 design with four observations per cell analyzed by
#    Steiger (2004, p. 171), F = 2.50 on 6 and 42 df, at 90%: the
#    total sample size of 56 gives the interval for the population
#    partial eta squared, and 6 + 42 + 1 = 49 gives the limits that
#    Smithson (2001) reports for the same F.
ci_eta_squared(F_value = 2.50, df_effect = 6, df_error = 42, N = 56,
               estimand = "partial", conf_level = 0.90)
#>  effect  eta_squared_partial eta_squared_partial_ci_lower
#>  overall 0.263               0.008                       
#>  eta_squared_partial_ci_upper F_value df_effect df_error N 
#>  0.316                        2.500   6         42       56
#> 
#> Confidence level: 90%
ci_eta_squared(F_value = 2.50, df_effect = 6, df_error = 42, N = 49,
               estimand = "partial", conf_level = 0.90)
#>  effect  eta_squared_partial eta_squared_partial_ci_lower
#>  overall 0.263               0.009                       
#>  eta_squared_partial_ci_upper F_value df_effect df_error N 
#>  0.346                        2.500   6         42       49
#> 
#> Confidence level: 90%

# 6. A within-subjects F from a report: the three IQ occasions of the
#    pygmalion data, F(2, 618) = 83.11 for 310 children. N is the
#    620 degrees of freedom of the occasion stratum, not 930 rows.
ci_eta_squared(F_value = 83.11, df_effect = 2, df_error = 618, N = 620,
               estimand = "partial")
#>  effect  eta_squared_partial eta_squared_partial_ci_lower
#>  overall 0.212               0.157                       
#>  eta_squared_partial_ci_upper F_value df_effect df_error N  
#>  0.265                        83.110  2         618      620
#> 
#> Confidence level: 95%
```
