# Moving From MBESS to DMAR: A Migration Guide and a Tour of the New Methods

``` r

library(DMAR)
```

## Why a Reimagining of MBESS

DMAR is a more general and greatly expanded reimagining of the MBESS
package (*Methods for the Behavioral, Educational, and Social Sciences*;
Kelley, 2007a, *Journal of Statistical Software*; 2007b, *Behavior
Research Methods*). MBESS shipped in 2006 and has been on CRAN ever
since. Two features of MBESS led to a new package. First, the methods it
implements apply well beyond the behavioral, educational, and social
sciences that its name lists, for example in clinical and translational
research, biostatistics, information systems, marketing, organizational
science, and sociology. Second, the package’s conventions (dotted
argument names, mixed-style returns, a `verbose =` flag controlling
print) reflect an earlier era of R style. DMAR ships under a new name to
make the broader scope explicit and to introduce a uniform interface
without changing the MBESS interface.

MBESS remains on CRAN, and scripts written for it continue to run under
it. The purpose of this vignette is to help users who move to DMAR, for
new work or to migrate existing scripts, and to show the methods that
DMAR adds beyond the scope of MBESS.

## Migration Table: The Renames That Matter

The naming convention in DMAR is `snake_case` throughout. Argument names
use underscores rather than dots (`conf_level`, not `conf.level`;
`alpha_level`, not `alpha.level`), and the canonical return is a tidy
`data.frame` with `term` and `value` columns. The table below maps the
most-used MBESS calls to their DMAR equivalents.

| MBESS call | DMAR call | Notes |
|----|----|----|
| `MBESS::ci.smd(ncp, n.1, n.2, conf.level)` | `DMAR::ci_smd(ncp, n_1, n_2, conf_level)` | Same noncentral *t* inversion; tidy `data.frame` return. |
| `MBESS::ci.smd.c(...)` | `DMAR::ci_smd_c(...)` | Glass’s $`g`$ with control-group SD. |
| `MBESS::ci.R2(R2, N, K, conf.level)` | `DMAR::ci_R2(R2, N, p, conf_level)` | `K -> p`; same fixed-vs-random predictors switch. |
| `MBESS::ci.reg.coef(...)` | `DMAR::ci_reg_coef(...)` | Same noncentral / central paths. |
| `MBESS::ci.rc(...)`, `MBESS::ci.src(...)` | `DMAR::ci_rc(...)`, `DMAR::ci_src(...)` | Per-coefficient CIs. |
| `MBESS::ci.cv(cv, n, ...)` | `DMAR::ci_cv(cv, N, ...)` | `n -> N`. Noncentral *t* CI for the coefficient of variation. The limits agree with MBESS’s until the noncentrality parameter at a limit passes 37.62, where R’s noncentral *t* is approximate and DMAR computes exactly (a difference of 0.1 percent for a coefficient of .10 with `N = 200`). |
| `MBESS::ci.pvaf(...)` | `DMAR::ci_pvaf(...)` | Proportion of variance accounted for. The limits agree with MBESS’s. When the lower limit is 0, the `actual_coverage` row is the stated level (`Actual.Coverage` in MBESS is 1 - `alpha.upper` there); [`?ci_pvaf`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md) describes the tail columns of the `ci_lower` row in that case. |
| `MBESS::ci.snr(...)` | `DMAR::ci_snr(...)` | Signal-to-noise CI. |
| `MBESS::ci.srsnr(...)` | `DMAR::ci_srsnr(...)` | Square root of signal-to-noise CI. |
| `MBESS::ci.sm(...)` | `DMAR::ci_sm(...)` | Standardized-mean CI; capitalized `Mean`/`SD` are now `mean`/`sd`. |
| `MBESS::ci.sc(...)`, `MBESS::ci.sc.ancova(...)` | `DMAR::ci_sc(...)`, `DMAR::ci_sc_ancova(...)` | Standardized contrast CIs. |
| `MBESS::ci.c(...)`, `MBESS::ci.c.ancova(...)` | `DMAR::ci_c(...)`, `DMAR::ci_c_ancova(...)` | Unstandardized contrast CIs. |
| `MBESS::ci.rmsea(...)` | `DMAR::ci_rmsea(...)` | Noncentral $`\chi^2`$ inversion for RMSEA; a 90% interval by default (MBESS: 95%; `conf_level = 0.95` matches it). |
| `MBESS::ci.cc(...)` | `DMAR::ci_r(...)` | Correlation CI; `n -> N`. |
| `MBESS::conf.limits.nct(...)` | `DMAR::ci_nc_t(...)` | Noncentral *t*; the observed *t* is `ncp` (no `t.value` alias). |
| `MBESS::conf.limits.ncf(...)` | `DMAR::ci_nc_F(...)` | Noncentral $`F`$. |
| `MBESS::conf.limits.nc.chisq(...)` | `DMAR::ci_nc_chisq(...)` | Noncentral $`\chi^2`$. |
| `MBESS::ss.aipe.smd(delta, conf.level, width, ...)` | `DMAR::ss_aipe_smd(true_smd, conf_level, width, ...)` | Same AIPE planner; tidy return. `degree.of.certainty -> assurance`, with which DMAR returns the smallest per-group *n* at which the exact probability that the interval is no wider than `width` reaches it, a case or so from MBESS’s plan in either direction. |
| `MBESS::ss.aipe.cv(C.of.V, width, mu, sigma, ...)` | `DMAR::ss_aipe_cv(true_cv, width, true_mean, true_sd, ...)` | `C.of.V -> true_cv`, the name of the data generating value in [`ss_aipe_cv_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md); `mu -> true_mean`, `sigma -> true_sd`. DMAR plans at the expected value of the sample coefficient with the exact interval (see [`?ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)), so without `assurance` its sample sizes equal MBESS’s or are 1 to 4 smaller, and occasionally 1 larger (for a coefficient of .25, a width of .10, and 99% confidence: 100 in MBESS, 101 in DMAR). With `assurance` the differences are of that size in nearly all conditions (142 in both for that example with an assurance of .99). |
| `MBESS::ss.aipe.cv.sensitivity(True.C.of.V, Estimated.C.of.V, ...)` | `DMAR::ss_aipe_cv_sensitivity(true_cv, supposed_cv, ...)` | `True.C.of.V -> true_cv`, `Estimated.C.of.V -> supposed_cv`, `mean -> true_mean`, `Specified.N -> specified_N`; `print.iter` (default `TRUE`) is `print_iter` (default `FALSE`). |
| `MBESS::ss.aipe.R2(Population.R2, conf.level, width, ...)` | `DMAR::ss_aipe_R2(true_R2, conf_level, width, ...)` | AIPE planner for $`R^2`$; `Population.R2 -> true_R2`; `K -> p`; `Random.Predictors -> random_predictors` (MBESS’s alias `Random.Regressors` is not carried over); `degree.of.certainty -> assurance`. With `which_width = "Lower"` or `"Upper"` DMAR returns the smallest sample size such that the distance from the expected sample $`R^2`$ to that limit is no larger than `width` at that sample size and at every larger one, where MBESS returns the first sample size of an upward search at which it is no larger, so a plan for `"Upper"` can be far larger than MBESS’s (for a population value of .70, a width of .05, and 10 predictors at 95% confidence: 23 in MBESS, 202 in DMAR). DMAR accepts `assurance` for the full width only: with `"Lower"` or `"Upper"` it is refused (see [`?ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)), where MBESS returns a sample size for `"Lower"`. |
| `MBESS::ss.aipe.reg.coef(...)` | `DMAR::ss_aipe_reg_coef(...)` | AIPE planner for a regression coefficient; `degree.of.certainty -> assurance`. With `assurance`, the central interval of the default `noncentral = FALSE` gets the smallest *N* at which the interval is no wider than `width` with that probability, computed exactly, where MBESS used a chi square approximation that DMAR does not offer for that interval. With `noncentral = TRUE` the default, `assurance_method = "analytic"`, is MBESS’s approximation, and `assurance_method = "simulation"` gives the smallest *N* at which a simulated estimate of the probability reaches `assurance`. The plans are usually a case or more below MBESS’s: the expected squared multiple correlation of the targeted predictor with the others is that of a regression on *p* - 1 predictors, where MBESS took *p*, and the plan is the smallest *N* that meets the width, where MBESS’s search could start above it; the central plans of the default `noncentral = FALSE` with unequal `alpha_lower` and `alpha_upper` take the *t* quantile of each tail, the interval [`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md) reports, where MBESS took both at the average tail. |
| `MBESS::ss.aipe.rmsea(...)` | `DMAR::ss_aipe_rmsea(...)` | AIPE for RMSEA; plans a 90% interval by default (MBESS: 95%; `conf_level = 0.95` matches it). |
| `MBESS::ss.power.R2(...)` | `DMAR::ss_power_R2(...)` | Power-based planner; `alpha.level -> alpha_level`. `Null.R2` is not carried over; the test of an increment is [`ss_power_R2_change()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md). |
| `MBESS::ss.power.reg.coef(...)` | `DMAR::ss_power_reg_coef(...)` | Power for a regression coefficient. |
| `MBESS::cv(C.of.V, mean, sd)` | `DMAR::cv(cv, mean, sd)` | Coefficient of variation; `C.of.V -> cv`; tidy return. |
| `MBESS::mr.cv(data, A, structural.cost, epsilon, sampling.cost, ...)` | `DMAR::mr_cv(data, A, structural_cost, epsilon, sampling_cost, ...)` | Minimum risk estimation of the coefficient of variation. The decision in MBESS’s `Is.Satisfied?` column is printed beneath the table and held in `attr(x, "is_satisfied")`. DMAR estimates the variance of the sample coefficient differently (see [`?mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.md)), so its risk and stopping point differ from MBESS’s for nearly any data. |
| `MBESS::mr.smd(A, structural.cost, epsilon, d, n, sampling.cost, ...)` | `DMAR::mr_smd(A, structural_cost, epsilon, d, n_per_group, sampling_cost, ...)` | Minimum risk estimation of the standardized mean difference; `n -> n_per_group`. The decision in MBESS’s `Is.Satisfied?` column is printed beneath the table and held in `attr(x, "is_satisfied")`. |
| `MBESS::sd.unbiased(...)` | `DMAR::sd_unbiased(...)` | Holtzman-corrected SD. |
| `MBESS::signal.to.noise.R2(R.Square, ...)` | `DMAR::signal_to_noise_R2(R2, ...)` | `R.Square -> R2` per the meaningful-capital rule. |
| `MBESS::smd(...)`, `MBESS::smd.c(...)` | `DMAR::smd(...)`, `DMAR::smd_c(...)` | Standardized mean difference point estimates. |
| `MBESS::HS`, `MBESS::Prime.Time` | [`DMAR::holzinger_swineford`](https://yelleknek.github.io/DMAR/reference/holzinger_swineford.md), [`DMAR::prime_time_achievement`](https://yelleknek.github.io/DMAR/reference/prime_time_achievement.md) | Same data under the documented snake_case names. |

The convention for the rename is: 1. Replace `.` with `_` in function
and argument names (`ci.smd -> ci_smd`, `n.1 -> n_1`,
`conf.level -> conf_level`). 2. Lowercase abbreviations that are not
statistical notation (`R.Square -> R2`, `Mean -> mean`, `SD -> sd`);
keep meaningful capitals (`R2`, `N`, `S`, `F_value`) and the capital
labels of variables (`X`, `Y`, `Z`, as in `true_cor_YX`). 3. Use `p` for
the number of predictors (MBESS sometimes used `K`). 4. Use
`random_predictors` (not `Random.Predictors` or its alias
`Random.Regressors`). 5. Name a value for its role and its quantity in
English: `true_` for a population value given to a planner
(`delta -> true_smd` in
[`ss_aipe_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md),
`C.of.V -> true_cv` in
[`ss_aipe_cv()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)),
and `supposed_` for the planning value a sensitivity function chooses
its design under.

Most DMAR functions return a tidy `data.frame` with `term` and `value`
columns
([`vignette("dmar_output")`](https://yelleknek.github.io/DMAR/articles/dmar_output.md)
describes the exceptions). Scripts that consumed MBESS’s named-list
returns generally need only adjust the extraction (e.g.,
`out$Lower.Conf.Limit.smd` becomes `out$value[out$term == "ci_lower"]`).

## What DMAR Adds Beyond the MBESS Scope

DMAR also implements methods outside the scope of MBESS. The additions
cluster in five areas.

### 1. Maximum Likelihood Multiple Regression With FIML (`mlmr()` / `mlmr_mv()`)

[`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) is an
[`lm()`](https://rdrr.io/r/stats/lm.html)-like front end to full
information maximum likelihood regression. The formula interface and S3
methods mirror [`lm()`](https://rdrr.io/r/stats/lm.html) (`coef`,
`vcov`, `confint`, `summary`, `anova`, `predict`, `update`); the default
confidence intervals are profile likelihood intervals with Wald and
bootstrap as alternatives; and the missing data handling is
`missing = "fiml"` by default. The multivariate sibling
[`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md)
takes `cbind(y1, y2) ~ ...` and models the joint distribution of
correlated outcomes, so that the observed outcomes of a case inform the
estimates for an outcome it is missing. The companion vignette
[`vignette("mlmr", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/mlmr.md)
walks through the missingness scenarios in which FIML actually matters.

### 2. broom-Style Integration (`generics::tidy()` and `generics::glance()`)

DMAR outputs dispatch through the broom-ecosystem generics in the
`generics` package, so `purrr::map_dfr(fits, generics::tidy)` works
across DMAR fits the same way it works across
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html), and other broom-supported
models. The families covered are `mlmr`, `mlmr_mv`, `cfa_1`, the
reliability family, the long-format CI family, the ANOVA effect size CI
family, and the power planner family.

### 3. AIPE Sensitivity Analysis

Every closed-form AIPE planner has a Monte Carlo sensitivity companion
(`ss_aipe_*_sensitivity()`) that simulates from a true population value
to quantify the realized CI width and empirical coverage when the
planning value is wrong. The companion is the recommended workflow when
the planning value comes from a small pilot or from a literature with
publication bias. A 10,000-replication sweep across the planner family,
reporting realized interval width and empirical coverage for each, is
maintained separately from the package.

### 4. ANOVA and ANCOVA Wrappers

DMAR adds tidy entry points for several ANOVA designs that required
manual model fitting in MBESS:
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md) for
the classical ANCOVA with adjusted means and the omnibus
$`\hat\omega^2_{\text{partial}}`$ CI;
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md)
for the fixed-random F-ratio bookkeeping in a two-way crossed design;
[`anova_within_two_way()`](https://yelleknek.github.io/DMAR/reference/anova_within_two_way.md)
for the two-factor within-subjects ANOVA with sphericity corrections per
effect;
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md)
for the mixed-design multivariate ANOVA. The functions return tidy
`data.frame`s suitable for direct piping into reporting tables.

### 5. The Reliability Family

The reliability family (`reliability_alpha`, `reliability_omega` (with a
model implied or observed total-variance denominator; the latter is
`MBESS::ci.reliability`’s “hierarchical” type),
`reliability_omega_categorical`, `reliability_kr20`, `reliability_H`,
and the dispatch wrapper
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md))
returns the point estimate with its standard error and confidence
interval where the chosen `ci_method` gives them (the menu includes the
Feldt, Fisher, Bonett, maximum likelihood, robust maximum likelihood,
and bootstrap intervals;
[`vignette("reliability")`](https://yelleknek.github.io/DMAR/articles/reliability.md)
gives the defaults), the sample size, and the number of items, in a tidy
schema that matches the rest of the package.
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)
provides the single-factor CFA fit on which several of those estimators
depend.

### Other Additions

- [`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md):
  tidy two-sample *t* test with separate variances.
- [`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md):
  exact paired sign-flip randomization test with Monte Carlo fallback.
- [`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md):
  power of Fisher’s exact 2×2 test.
- [`is_orthogonal_set()`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md):
  check orthogonality of an entire contrast matrix.
- [`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
  [`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
  [`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
  [`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md):
  critical values and CIs for the classical multiple-comparison
  procedures, returning tidy `data.frame`s.
- [`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md):
  tidy ICC + CI from a fitted `lmerMod`.
- [`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
  [`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
  [`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md):
  critical values for the classical procedures, with explicit handling
  of the studentized maximum modulus distribution.
- A modernized, ggplot2-based set of plotting functions:
  [`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
  [`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
  [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
  [`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
  and
  [`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md).
- A package-wide `term`/`value` `data.frame` schema; broom-style S3
  methods for many families; `set.seed(113)` as the reproducibility seed
  of the help page examples; `seed = NULL` as the default for every
  bootstrap and Monte Carlo function, with the caller’s generator state
  restored when a function given a seed returns.

## A Short Worked Migration

A small MBESS script computing a noncentral *t* CI on the standardized
mean difference (Cohen’s $`d`$), its companion AIPE sample size plan,
and the realized CI width under the plan looks like this in MBESS:

``` r

# MBESS, classic
library(MBESS)
out_ci  <- ci.smd(ncp = 4, n.1 = 30, n.2 = 30, conf.level = 0.95)
out_ss  <- ss.aipe.smd(delta = 0.5, conf.level = 0.95, width = 0.40)
out_ci$Lower.Conf.Limit.smd
out_ci$Upper.Conf.Limit.smd
```

The DMAR equivalent, with tidy returns and the
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
route into a unified table:

``` r

ci  <- ci_smd(ncp = 4, n_1 = 30, n_2 = 30, conf_level = 0.95)
ss  <- ss_aipe_smd(true_smd = 0.5, conf_level = 0.95, width = 0.40)
ci
```

| term     | value |
|:---------|:------|
| ci_lower | 0.489 |
| smd      | 1.033 |
| ci_upper | 1.569 |

Confidence level: 95%

``` r

ss
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 199   |
| true_smd              | 0.500 |
| width                 | 0.400 |

Confidence level: 95%

``` r

generics::tidy(ci)
#>   term estimate  ci_lower ci_upper conf_level
#> 1  smd 1.032796 0.4891759 1.568559       0.95
```

The two return values compose with `dplyr` summaries and `ggplot2` plots
without further wrapping.

## See Also

- [`vignette("DMAR", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/DMAR.md)
  for the introductory tour.
- [`vignette("mlmr", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/mlmr.md)
  for the FIML regression family.
- [`vignette("reliability", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/reliability.md)
  for the reliability family.
- `NEWS.md` for the per-release changelog.

## References

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*(11), 1547–1562. <https://doi.org/10.1177/0956797617723724>

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. <https://doi.org/10.18637/jss.v020.i08>

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
<https://doi.org/10.3758/BF03192993>

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
<https://doi.org/10.1146/annurev.psych.59.103006.093735>

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
