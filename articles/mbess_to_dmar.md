# Moving From MBESS to DMAR: A Migration Guide and a Tour of the New Methods

``` r
library(DMAR)
```

## Why a Reimplementation and Expansion

The MBESS package (*Methods for the Behavioral, Educational, and Social
Sciences*; Kelley, 2007a, *Journal of Statistical Software*; 2007b,
*Behavior Research Methods*) shipped in 2006 and has been on CRAN ever
since. Its uptake outgrew the original framing in two ways. First, the
methods it implements are used well beyond the behavioral, educational,
and social sciences, including in clinical and translational research,
biostatistics, information systems, marketing, organizational science,
sociology, education, and the methodological literature itself. Second,
the package’s API conventions (dotted argument names, mixed-style
returns, a `verbose =` flag controlling print) reflect an earlier era of
R style. DMAR ships under a new name to make the scope expansion
explicit and to introduce a uniform, modern API without breaking
compatibility with the long-stable MBESS interface.

MBESS itself remains stable on CRAN; researchers with running scripts
that depend on MBESS can continue to use it indefinitely. The purpose of
this vignette is to help users who want to move forward to DMAR, either
for new work or to gradually migrate existing scripts, and to make
visible the methods that DMAR adds beyond MBESS’s scope.

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
| `MBESS::ci.cv(...)` | `DMAR::ci_cv(...)` | CV with McKay/Vangel CIs. |
| `MBESS::ci.pvaf(...)` | `DMAR::ci_pvaf(...)` | Proportion of variance accounted for. |
| `MBESS::ci.snr(...)` | `DMAR::ci_snr(...)` | Signal-to-noise CI. |
| `MBESS::ci.srsnr(...)` | `DMAR::ci_srsnr(...)` | Square root of signal-to-noise CI. |
| `MBESS::ci.sm(...)` | `DMAR::ci_sm(...)` | Standardized-mean CI; capitalized `Mean`/`SD` are now `mean`/`sd`. |
| `MBESS::ci.sc(...)`, `MBESS::ci.sc.ancova(...)` | `DMAR::ci_sc(...)`, `DMAR::ci_sc_ancova(...)` | Standardized contrast CIs. |
| `MBESS::ci.c(...)`, `MBESS::ci.c.ancova(...)` | `DMAR::ci_c(...)`, `DMAR::ci_c_ancova(...)` | Unstandardized contrast CIs. |
| `MBESS::ci.rmsea(...)` | `DMAR::ci_rmsea(...)` | Noncentral-$`\chi^2`$ inversion for RMSEA. |
| `MBESS::ci.cc(...)` | `DMAR::ci_cc(...)` | Correlation CI; `r` and `n` arguments. |
| `MBESS::conf.limits.nct(...)` | `DMAR::conf_limits_nct(...)` | Noncentral *t*; `t.value -> t_value`. |
| `MBESS::conf.limits.ncf(...)` | `DMAR::conf_limits_ncf(...)` | Noncentral $`F`$. |
| `MBESS::conf.limits.nc.chisq(...)` | `DMAR::conf_limits_nc_chisq(...)` | Noncentral $`\chi^2`$. |
| `MBESS::ss.aipe.smd(delta, conf.level, width, ...)` | `DMAR::ss_aipe_smd(delta, conf_level, width, ...)` | Same AIPE planner; tidy return. |
| `MBESS::ss.aipe.R2(...)` | `DMAR::ss_aipe_R2(...)` | AIPE planner for $`R^2`$; `K -> p`; `random.regressors` argument retired in favor of `random_predictors`. |
| `MBESS::ss.aipe.reg.coef(...)` | `DMAR::ss_aipe_reg_coef(...)` | AIPE planner for a regression coefficient. |
| `MBESS::ss.aipe.rmsea(...)` | `DMAR::ss_aipe_rmsea(...)` | AIPE for RMSEA. |
| `MBESS::ss.power.R2(...)` | `DMAR::ss_power_R2(...)` | Power-based planner; `alpha.level -> alpha_level`. |
| `MBESS::ss.power.reg.coef(...)` | `DMAR::ss_power_reg_coef(...)` | Power for a regression coefficient. |
| `MBESS::cv(mean, sd)` | `DMAR::cv(mean, sd)` | Coefficient of variation; tidy return. |
| `MBESS::sd.unbiased(...)` | `DMAR::sd_unbiased(...)` | Holtzman-corrected SD. |
| `MBESS::signal.to.noise.R2(R.Square, ...)` | `DMAR::signal_to_noise_R2(R2, ...)` | `R.Square -> R2` per the meaningful-capital rule. |
| `MBESS::smd(...)`, `MBESS::smd.c(...)` | `DMAR::smd(...)`, `DMAR::smd_c(...)` | Standardized mean difference point estimates. |
| [`MBESS::HS`](https://rdrr.io/pkg/MBESS/man/HS.html), `MBESS::Prime.Time` | [`DMAR::holzinger_swineford`](https://yelleknek.github.io/DMAR/reference/holzinger_swineford.md) (`HS_Data` alias), [`DMAR::prime_time_achievement`](https://yelleknek.github.io/DMAR/reference/prime_time_achievement.md) (`Prime_Time` alias) | Same data; documented canonical names; short aliases. |

The convention for the rename is: 1. Replace `.` with `_` in function
and argument names (`ci.smd -> ci_smd`, `n.1 -> n_1`,
`conf.level -> conf_level`). 2. Lowercase abbreviations that are not
statistical notation (`R.Square -> R2`, `Mean -> mean`, `SD -> sd`);
keep meaningful capitals (`R2`, `N`, `S`, `Lambda`, `F_value`). 3. Use
`p` for the number of predictors (MBESS sometimes used `K`). 4. Use
`random_predictors` (not `random.regressors`).

The DMAR return is a tidy `data.frame` with stable column schemas across
the package. Scripts that consumed MBESS’s named-list returns generally
need only adjust the extraction (e.g., `out$Lower.Conf.Limit.smd`
becomes `out$value[out$term == "lower_limit"]`).

## What DMAR Adds Beyond the MBESS Scope

The migration story is only half the story. The reason to move forward
is what DMAR does that MBESS does not. The additions cluster in five
areas.

### 1. Maximum Likelihood Multiple Regression with FIML (`mlmr()` / `mlmr_mv()`)

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
correlated outcomes, which is the case where the FIML advantage over
listwise deletion is largest. The companion vignette
[`vignette("mlmr", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/mlmr.md)
walks through the missingness scenarios in which FIML actually matters.

### 2. broom-Style Integration (`generics::tidy()` and `generics::glance()`)

DMAR outputs dispatch through the broom-ecosystem generics in the
`generics` package, so `purrr::map_dfr(fits, generics::tidy)` works
across DMAR fits the same way it works across
[`lm()`](https://rdrr.io/r/stats/lm.html),
[`glm()`](https://rdrr.io/r/stats/glm.html), and other broom-supported
models. The families covered as of this release: `mlmr`, `mlmr_mv`,
`cfa_1`, the reliability family, the long-format CI family, the ANOVA
effect size CI family, and the power planner family.

### 3. AIPE Sensitivity Analysis

Every closed-form AIPE planner has a Monte Carlo sensitivity companion
(`ss_aipe_*_sensitivity()`) that simulates from a true population value
to quantify the realized CI width and empirical coverage when the
planning value is wrong. The companion is the recommended workflow when
the planning value comes from a small pilot or from a literature with
publication bias. The
[`vignette("aipe_simulation_study", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/aipe_simulation_study.md)
reports a 10,000-replication sweep across the planner family.

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

### 5. Reliability with Proper CIs

The reliability family (`reliability_alpha`, `reliability_omega` (with a
model implied or observed total-variance denominator; the latter is
[`MBESS::ci.reliability`](https://rdrr.io/pkg/MBESS/man/ci.reliability.html)’s
“hierarchical” type), `reliability_omega_categorical`,
`reliability_kr20`, `reliability_H`, and the dispatch wrapper
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md))
returns the point estimate alongside the Feldt/Bonett/Fisher CI, the
delta method SE when applicable, the sample size, and the number of
items, in a tidy schema that matches the rest of the package.
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
  [`dunnett_ci()`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md),
  [`tukey_kramer_ci()`](https://yelleknek.github.io/DMAR/reference/tukey_kramer_ci.md),
  [`scheffe_ci()`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md):
  critical values and CIs for the classical multiple-comparison
  procedures, returning tidy `data.frame`s.
- [`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md):
  tidy ICC + CI from a fitted `lmerMod`.
- [`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
  [`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
  [`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md):
  critical values for the classical procedures, with explicit handling
  of the Studentized maximum modulus distribution.
- A modernized, ggplot2-based set of plotting functions:
  [`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
  [`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
  [`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
  [`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
  and
  [`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md).
- A package-wide `term`/`value` `data.frame` schema; broom-style S3
  methods for many families; `set.seed(113)` as the package-wide
  reproducibility seed for examples and tests; `seed = NULL` as the
  default for every bootstrap and Monte Carlo function, with explicit
  save/restore of `.Random.seed` when the user supplies a seed.

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
ss  <- ss_aipe_smd(delta = 0.5, conf_level = 0.95, width = 0.40)
ci
```

| term        | value |
|:------------|:------|
| lower_limit | 0.489 |
| smd         | 1.03  |
| upper_limit | 1.57  |

Confidence level: 95%

``` r
ss
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 199   |
| supposed_smd          | 0.5   |
| width                 | 0.4   |

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
