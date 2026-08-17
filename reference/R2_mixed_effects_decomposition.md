# R-Squared Measures for Mixed-Effects Models

Computes the Rights and Sterba (2019) integrative framework of R-squared
measures for a fitted two-level linear mixed-effects (multilevel) model.
The model implied outcome variance is fully decomposed into five
sources: variance due to level-1 predictors via fixed slopes (\\f_1\\),
level-2 predictors via fixed slopes (\\f_2\\), predictors via random
slope (co)variation (\\v\\), cluster-specific outcome means via random
intercept variation (\\m\\), and level-1 residuals (\\\sigma^2\\).
Proportions of the total, within-cluster, and between-cluster outcome
variance attributable to combinations of these sources give the family
of R-squared measures.

## Usage

``` r
R2_mixed_effects_decomposition(model)
```

## Arguments

- model:

  A fitted two-level model of class `merMod` (from lme4, e.g.\\
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html)) or `lme` (from
  nlme). The model must use numeric predictors (only the cluster
  variable may be a factor) and must not contain
  [`I()`](https://rdrr.io/r/base/AsIs.html) terms; create any
  transformed predictors as their own columns first.

## Value

A `data.frame` (`dmar_tbl`) with columns `term` and `value`. When the
level-1 predictors are cluster-mean-centered, the full set of 12
measures is returned, named `total_f1`, `total_f2`, `total_v`,
`total_m`, `total_f`, `total_fv`, `total_fvm`, `within_f1`, `within_v`,
`within_fv`, `between_f2`, and `between_m`; otherwise the five
total-variance measures `total_f`, `total_v`, `total_m`, `total_fv`, and
`total_fvm` are returned (the within/between split requires
cluster-mean-centering). The returned object carries the
source-by-target variance decomposition in `attr(x, "decomposition")`.

## Details

The companion
[`R2_mixed_effects`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md)
returns the Nakagawa and Schielzeth marginal and conditional R-squared;
this function contains those two as special cases (`total_f` and
`total_fvm`) within the fuller source decomposition.

The measure superscripts index the variance sources in the numerator and
the subscripts index the outcome variance in the denominator: `total_*`
use the total outcome variance, `within_*` the within-cluster variance
(\\f_1 + v + \sigma^2\\), and `between_*` the between-cluster variance
(\\f_2 + m\\). `total_fvm` is the omnibus measure (all explained sources
over the total variance) and, for a random-intercept model, coincides
with the Nakagawa and Schielzeth conditional R-squared computed by
[`R2_mixed_effects`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md);
`total_f` coincides with their marginal R-squared. See Rights and Sterba
(2019, Table 1) for the definitions.

The measures were derived under the assumption that the fitted model
uses cluster-mean-centering of the level-1 predictors (with the cluster
means entered as level-2 predictors). When that centering is not
detected, only the total-variance measures are returned, matching the
reference implementation.

Fitting the model requires lme4 (for a `merMod` fit) or nlme (for an
`lme` fit) to be installed.

## References

Rights, J. D., & Sterba, S. K. (2019). Quantifying explained variance in
multilevel models: An integrative framework for defining R-squared
measures. *Psychological Methods, 24*(3), 309–338.
[doi:10.1037/met0000184](https://doi.org/10.1037/met0000184)

Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
obtaining \\R^2\\ from generalized linear mixed-effects models. *Methods
in Ecology and Evolution, 4*(2), 133–142.
[doi:10.1111/j.2041-210x.2012.00261.x](https://doi.org/10.1111/j.2041-210x.2012.00261.x)

## See also

[`R2_mixed_effects`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
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
R2_mixed_effects_decomposition(fit)
#>  term      value 
#>  total_f   0.279 
#>  total_v   0.0892
#>  total_m   0.432 
#>  total_fv  0.368 
#>  total_fvm 0.799 
```
