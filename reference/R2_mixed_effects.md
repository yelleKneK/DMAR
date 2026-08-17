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
  seed = NULL,
  ...
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

  Number of bootstrap replications when `ci_method = "boot"`. Default
  `1000`.

- seed:

  Optional integer seed for the bootstrap. Default `NULL`, meaning the
  caller's current RNG state is used and left unchanged. When supplied,
  the seed is set inside the function and the caller's RNG state is
  restored on exit.

- ...:

  Currently unused.

## Value

A `data.frame` with rows `"R2_marginal"` and `"R2_conditional"` in the
`value` column. When `ci_method = "boot"`, lower- and upper-limit rows
for each quantity are appended and the confidence level is carried on
the object.

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
or bootstrap standard error variant. Replicates whose refit fails are
dropped, and the interval is computed from the replications that return
a value. The default `B = 1000` is adequate for the central quantiles a
percentile interval uses; raising it tightens the Monte Carlo error of
the reported limits. Bootstrap results vary from run to run; supply
`seed` for reproducibility.

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

# A bootstrap interval is available through ci_method = "boot". It refits
# the model once per replication, so it is not run here; the call is
#   R2_mixed_effects(fit, ci_method = "boot", B = 1000, seed = 113)
# where B = 1000 is the default and seed is supplied because the limits
# otherwise move from run to run. Raise B when the Monte Carlo error of
# the reported limits needs to be smaller.
```
