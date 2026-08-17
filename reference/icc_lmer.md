# Intraclass Correlation From a Fitted `lme4` Mixed-Effects Model

Reads the variance-component decomposition off a fitted
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) model and returns the
implied intraclass correlation (ICC) for the specified grouping factor
along with a Bonett (2002) Fisher-\\L\\-transform confidence interval,
in tidy long form. The bridge function between DMAR's classical- ANOVA
path
([`variance_components_mls`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md)) and
the modern mixed-effects path.

## Usage

``` r
icc_lmer(fit, group = NULL, conf_level = 0.95)
```

## Arguments

- fit:

  A fitted [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) model (class
  `lmerMod`) with at least one random-effects grouping factor.

- group:

  Character name of the grouping factor whose ICC is wanted. If `NULL`
  (default), the first random-effects grouping factor is used. For
  three-level models, supply the target level explicitly.

- conf_level:

  Confidence level. Default `0.95`.

## Value

A `data.frame` with rows for the point estimate of the ICC, the variance
components (between-group and residual), the implied total variance, the
Bonett (2002) CI lower/upper limits, and the cluster-level effective
sample size used in the CI.

## Details

**Definition.** For a two-level model with random intercept by group
\\j\\ and within-group residual, \$\$\rho \\=\\
\frac{\sigma^2_b}{\sigma^2_b + \sigma^2_w},\$\$ read directly from
`VarCorr(fit)`. For three-level models, the user specifies which level
(`group`) provides the variance contribution; the denominator is the
total variance summed across all variance components.

**Bonett (2002) CI.** The CI is built on the Fisher-style
\\L\\-transformation \$\$L \\=\\ \tfrac{1}{2} \log\left(\frac{1 +
(k - 1) \rho}{1 - \rho}\right),\$\$ with variance \\k / (2 (k - 1) (n -
2))\\, where \\k\\ is the average cluster size and \\n\\ is the number
of clusters. Back- transformation keeps the bounds in \\\[0, 1\]\\.

**Limitations.** The Bonett CI is built on a balanced / approximately
balanced approximation; for severely unbalanced designs the profile
likelihood CI from `confint(fit, ...)` is preferable.

## References

Bonett, D. G. (2002). Sample size requirements for estimating intraclass
correlations with desired precision. *Statistics in Medicine, 21*(9),
1331–1335. [doi:10.1002/sim.1108](https://doi.org/10.1002/sim.1108)

Donner, A. (1986). A review of inference procedures for the intraclass
correlation coefficient in the one-way random effects model.
*International Statistical Review, 54*(1), 67–82.

Snijders, T. A. B., & Bosker, R. J. (2012). *Multilevel analysis: An
introduction to basic and advanced multilevel modeling* (2nd ed.). Sage.

## See also

[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md),
[`variance_components_mls`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md),
[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html)

Other agreement and measurement:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
[`limits_of_agreement()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md),
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
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
# Twenty groups of six, with a group-level standard deviation of 0.7 on
# top of within-group noise with a standard deviation of 1, so the data
# generating ICC is 0.7^2 / (0.7^2 + 1) = 0.329. The estimate below comes
# with a wide interval: 20 clusters is not many, and the interval is what
# keeps that fact visible.
set.seed(113)
n_grp <- 20; n_per <- 6
grp <- factor(rep(1:n_grp, each = n_per))
y   <- rnorm(n_grp * n_per, 0, 1) + rep(rnorm(n_grp, 0, 0.7), each = n_per)
d <- data.frame(y, grp)
fit <- lme4::lmer(y ~ 1 + (1 | grp), data = d)
icc_lmer(fit)
#>  term                 value 
#>  icc                  0.174 
#>  sigma2_between       0.248 
#>  sigma2_within        1.18  
#>  total_variance       1.42  
#>  lower_limit          0.0174
#>  upper_limit          0.377 
#>  n_clusters           20    
#>  average_cluster_size 6     
#> 
#> Confidence level: 95%
```
