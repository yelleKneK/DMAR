# Marginal and Conditional \\R^2\\ for a Mixed-Effects Model

Computes the Nakagawa and Schielzeth (2013) marginal and conditional
coefficients of determination from a fitted mixed-effects model,
returned in tidy long form. The marginal \\R^2\\ is the proportion of
total variance explained by the fixed effects alone; the conditional
\\R^2\\ is the proportion explained by the fixed and random effects
together. The two quantities are the mixed-effects companions to the
intraclass correlation returned by
[`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md):
where the ICC isolates the share of variance attributable to a single
grouping factor, the marginal and conditional \\R^2\\ summarize how much
of the outcome variance the fixed and random parts of the model account
for.

## Usage

``` r
R2_mixed_effects(
  model,
  conf_level = 0.95,
  ci_method = c("none", "boot"),
  B = 1000,
  seed = NULL
)
```

## Arguments

- model:

  A fitted mixed-effects model. A
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) fit (class `lmerMod`)
  is the primary target; an
  [`lme`](https://rdrr.io/pkg/nlme/man/lme.html) fit is also accepted.

- conf_level:

  Confidence level. Default `0.95`. Used only when a bootstrap interval
  is requested (see `ci_method`).

- ci_method:

  Interval method for the two \\R^2\\ values. The default `"none"`
  returns point estimates only. `"boot"` adds a parametric bootstrap
  percentile interval (`lmerMod` models only, via
  [`bootMer`](https://rdrr.io/pkg/lme4/man/bootMer.html)); see
  *Details*.

- B:

  Number of bootstrap replications when `ci_method = "boot"`, a single
  whole number. Default `1000`, because every replication simulates a
  response and refits the model. The fewest accepted is 400 at the
  default `conf_level` (200 at 0.90, and 2000 at 0.99, which the default
  does not meet); a smaller `B` is refused before anything is computed
  (class `dmar_count_invalid`). See *Bootstrap replications* under
  *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- seed:

  Optional integer seed for the bootstrap. Default `NULL`, under which
  the bootstrap draws from, and advances, the caller's current random
  number stream. When supplied, the seed is set inside the function and
  the caller's RNG state is restored on exit.

## Value

A `data.frame` with rows `"R2_marginal"` and `"R2_conditional"` in the
`value` column. When `ci_method = "boot"`, lower- and upper-limit rows
for each quantity are appended (`R2_marginal_ci_lower`,
`R2_marginal_ci_upper`, `R2_conditional_ci_lower`,
`R2_conditional_ci_upper`), the confidence level is carried on the
object, and the attributes `B`, `B_used`, `n_dropped`, `B_required`, and
`bootstrap_status` record the replications (see *Bootstrap replications*
under Details).

## Details

**Variance decomposition.** Writing \\\sigma^2_f\\ for the variance of
the fixed-effect linear predictor, \\\sigma^2_r\\ for the variance
attributable to the random effects, and \\\sigma^2\_\varepsilon\\ for
the residual variance, \$\$R^2\_{\mathrm{marginal}} \\=\\
\frac{\sigma^2_f}{\sigma^2_f + \sigma^2_r + \sigma^2\_\varepsilon},
\qquad R^2\_{\mathrm{conditional}} \\=\\ \frac{\sigma^2_f +
\sigma^2_r}{\sigma^2_f + \sigma^2_r + \sigma^2\_\varepsilon}.\$\$ The
fixed-effect variance is \\\sigma^2_f =
\mathrm{var}(\mathbf{X}\boldsymbol{\beta})\\, the variance of the fitted
fixed-effect linear predictor across the observations. The residual
variance is \\\sigma^2\_\varepsilon =
\mathrm{sigma}(\mathrm{model})^2\\.

**Random-effect variance.** For a random-intercept model the
random-effect variance is the sum of the variance components read off
[`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html). For a model with
random slopes the variance contributed by a random-effects term depends
on the values of the associated covariates, so the sum of the diagonal
variance components is not correct on its own. This function uses the
Johnson (2014) extension: for each random-effects term with design
matrix \\\mathbf{Z}\\ and estimated covariance matrix
\\\boldsymbol{\Sigma}\\, its contribution is the mean over the
observations of the quadratic form \\\mathbf{z}\_i^\top
\boldsymbol{\Sigma}\\ \mathbf{z}\_i\\, that is,
\\\tfrac{1}{n}\\\mathrm{tr}(\mathbf{Z}\boldsymbol{\Sigma}\mathbf{Z}^\top)\\,
and \\\sigma^2_r\\ is the sum of these contributions across all
random-effects terms. For a random-intercept term this reduces to the
intercept variance component, so the two paths agree.

**Scope.** The decomposition here is the one appropriate for a Gaussian
(identity-link) linear mixed model, which is what `lmer` and `lme` fit.
Generalized linear mixed models introduce a distribution-specific
variance term and are not handled by this function.

**The bootstrap interval.** The default `ci_method = "none"` reports the
two point estimates alone, so the bootstrap is what to ask for when the
marginal and conditional \\R^2\\ are to be reported with an interval and
the refits it costs are affordable. With `ci_method = "boot"` the
interval comes from a parametric bootstrap
([`bootMer`](https://rdrr.io/pkg/lme4/man/bootMer.html)): each of the
`B` replicates (1000 by default) simulates a new response vector from
the fitted model, refits the model, and recomputes the two \\R^2\\
values. The unit of resampling is therefore a whole simulated data set
drawn from the estimated model, not a resampled set of cases. Only the
percentile interval is offered: the limits are the empirical quantiles
of the `B` bootstrap values (Efron & Tibshirani, 1993); there is no BCa
or bootstrap standard error variant. Raising `B` tightens the Monte
Carlo error of the reported limits. Bootstrap results vary from run to
run; supply `seed` for reproducibility.

**Bootstrap replications.** The limits of each interval are quantiles of
the `B` replicated values, so the interval is run only when `B` puts at
least 10 replicates beyond each limit (400 at the 95 percent level). A
replication in which either \\R^2\\ is undefined (a refit that fails, or
one whose variance components cannot be computed) is dropped and
counted, never replaced, and the interval is computed from the rest. The
result records `B` (requested), `B_used` (the number the interval rests
on), `n_dropped`, and `B_required` (the fewest replications a request
may ask for; that count is held to `B`, not to `B_used`) as attributes,
with `bootstrap_status` (`"complete"`, `"dropped"`, `"exploratory"`,
`"insufficient_count"`, or `"insufficient_share"`). One warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the 2.5 percent of `B` that lies beyond one limit at the 95
percent level, because they could then decide by themselves where that
limit falls. When more than half were dropped, the four limit rows are
`NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`. Such a replicate is dropped for both
quantities, so the two intervals rest on the same replicates.
`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration; every such call gives a
warning of class `dmar_bootstrap_exploratory` and the result is not one
to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

## References

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Johnson, P. C. D. (2014). Extension of Nakagawa & Schielzeth's
\\R^2\_{GLMM}\\ to random slopes models. *Methods in Ecology and
Evolution, 5*(9), 944–946.
[doi:10.1111/2041-210X.12225](https://doi.org/10.1111/2041-210X.12225)

Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
obtaining \\R^2\\ from generalized linear mixed-effects models. *Methods
in Ecology and Evolution, 4*(2), 133–142.
[doi:10.1111/j.2041-210x.2012.00261.x](https://doi.org/10.1111/j.2041-210x.2012.00261.x)

## See also

[`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`ss_power_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html),
[`VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html)

Other agreement and measurement:
[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
[`limits_of_agreement()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md),
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)

Other mixed models:
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                  data = lme4::sleepstudy)

# Marginal R2 is the proportion of variance the fixed effects account for;
# conditional R2 adds what the random effects account for, so the gap
# between the two is what the subject-level terms buy.
R2_mixed_effects(fit)
#>  term           value
#>  R2_marginal    0.279
#>  R2_conditional 0.799
#> 
#> Confidence level: 95%

# A parametric bootstrap percentile interval for both quantities. Each
# replication refits the model. B = 20 keeps the example quick and is
# below the 400 a reported interval needs; the option lets it run and
# the warning says so. A reported interval deserves the default
# B = 1000. The seed makes the limits reproducible and leaves the
# caller's generator state as it was.
old <- options(dmar.bootstrap_exploratory = TRUE)
R2_mixed_effects(fit, ci_method = "boot", B = 20, seed = 113)
#> Warning: R2_mixed_effects(): B = 20 is below the 400 replications a 95% interval is reported from (10 beyond each limit). The limits are shown because options(dmar.bootstrap_exploratory = TRUE) is set. With 20 replications each limit rests on about 0.5 of them and will move from seed to seed; use the result to explore or to demonstrate, and report an interval only from B of at least 400. B, B_used, and B_required are recorded in the result.
#>  term                    value
#>  R2_marginal             0.279
#>  R2_marginal_ci_lower    0.151
#>  R2_marginal_ci_upper    0.470
#>  R2_conditional          0.799
#>  R2_conditional_ci_lower 0.710
#>  R2_conditional_ci_upper 0.855
#> 
#> Confidence level: 95%
#> 
#> Bootstrap: 20 replications, 20 used (0 dropped; at least 400 must be requested). Exploratory: below the number required for a reported interval.
options(old)
```
