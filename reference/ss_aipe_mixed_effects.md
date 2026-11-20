# AIPE Sample Size Planning for a Fixed Effect in a Two-Level Mixed-Effects Model

Computes the smallest number of clusters (i.e., Level 2 units) for which
the confidence interval for the fixed-effect slope of one predictor has
an expected full width no larger than `width`, in a balanced two-level
model with a random intercept. The predictor may be measured on the
individuals (i.e., Level 1) or on the clusters (i.e., Level 2), and may
be centered at its cluster mean, centered at its grand mean, or left as
it is.
[`ss_aipe_mixed_effects_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md)
simulates the same design, analysis, and interval, so the two can be
used together: plan here, then check the plan there, including under
values other than the ones planned with.

## Usage

``` r
ss_aipe_mixed_effects(
  true_total_var,
  true_var_X,
  true_icc,
  width,
  cluster_size = 20L,
  conf_level = 0.95,
  predictor_level = 1,
  centering = c("cluster_mean", "grand_mean", "none"),
  true_icc_X = 0
)
```

## Arguments

- true_total_var:

  Population variance of the outcome given the predictor (a planning
  value): the variance of the random intercept (between clusters) plus
  the residual variance (within clusters). It is the name the total
  variance has in the cluster randomized planners
  ([`ss_aipe_crd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crd.md)).

- true_var_X:

  Population variance of the predictor (a planning value). For an
  individual-level (i.e., Level 1) predictor it is the total variance
  over individuals, between and within clusters together (see
  `true_icc_X`); for a cluster-level (i.e., Level 2) predictor it is the
  variance over clusters.

- true_icc:

  Population intraclass correlation of the outcome given the predictor
  (a planning value): the share of `true_total_var` that lies between
  clusters. In \[0, 1).

- width:

  Desired full width of the confidence interval for the slope.

- cluster_size:

  Number of individuals (i.e., Level 1 units) in each cluster (i.e.,
  Level 2 unit), the same in every cluster. Default `20L`.

- conf_level:

  Confidence level. Default `0.95`.

- predictor_level:

  `1` (the default) for an individual-level predictor, which varies
  within a cluster, or `2` for a cluster-level predictor, which is the
  same for every individual of a cluster.

- centering:

  How the predictor enters the model. `"cluster_mean"` (the default)
  subtracts each cluster's mean from its individuals' values, so the
  slope is the within-cluster slope; it applies to an individual-level
  predictor only. `"grand_mean"` subtracts the overall mean, and
  `"none"` leaves the predictor as it is. These two give the same slope
  and the same interval (only the intercept differs), and both are
  offered so that the call can say what the analysis will do.

- true_icc_X:

  Population intraclass correlation of an individual-level predictor (a
  planning value): the share of `true_var_X` that lies between clusters.
  Default `0`, a predictor whose cluster means differ by sampling alone.
  It is not used with `predictor_level = 2`.

## Value

A `data.frame` with rows for the necessary number of clusters
`necessary_n_clusters`, the implied total sample size `total_N`
(`necessary_n_clusters * cluster_size`), the desired `width`, the
`expected_width` at the necessary number of clusters, the degrees of
freedom `df` of the interval's *t* quantile, and the inputs `true_icc`,
`true_icc_X`, `cluster_size`, and `predictor_level`. The centering is
recorded in `attr(x, "centering")`.

## Details

**The model.** For individual \\i\\ of cluster \\k\\, \$\$y\_{ki} =
\gamma_0 + \beta x\_{ki} + u_k + e\_{ki},\$\$ with \\u_k \sim N(0,
\tau^2)\\ the random intercept of the cluster and \\e\_{ki} \sim N(0,
\sigma^2)\\ the residual of the individual, in \\K\\ clusters of \\m\\
individuals each. With \\V\\ for `true_total_var` and \\\rho\\ for
`true_icc`, \\\tau^2 = \rho V\\ and \\\sigma^2 = (1 - \rho) V\\. The
slope has the same value within and between clusters. Writing
\\\sigma^2_X\\ for `true_var_X` and \\\rho_X\\ for `true_icc_X`, an
individual-level predictor has variance \\\rho_X \sigma^2_X\\ between
clusters and \\\sigma^2_W = (1 - \rho_X)\sigma^2_X\\ within them.

**The interval planned for.** The interval is \\\hat\beta \pm t\_{1 -
\alpha/2,\\ \nu}\\ \widehat{SE}(\hat\beta)\\, with \\\nu = K(m - 1) -
1\\ for an individual-level predictor and \\\nu = K - 2\\ for a
cluster-level predictor (the regression of \\K\\ cluster means). Its
expected width is \\2\\ t\_{1 - \alpha/2,\\ \nu}\\
E\[\widehat{SE}(\hat\beta)\]\\, and the function returns the smallest
\\K\\ at which that is no larger than `width`. No value of the slope is
needed: \\\beta\\ does not appear in any of the expressions below.

**The expected standard error.** Each expression follows from the sums
of squares of a balanced design.

- *Individual-level predictor centered at its cluster mean.* The
  estimate is the within-cluster slope, and for a normally distributed
  predictor the expected standard error is exactly
  \$\$E\[\widehat{SE}(\hat\beta)\] = \frac{\sigma}{\sigma_W \sqrt{K
  (m - 1) - 1}}.\$\$ Clustering of the outcome enters only through
  \\\sigma\\: the design effect \\1 + (m - 1)\rho\\ that inflates the
  variance of a cluster-level quantity does not appear, because
  centering removes the differences among clusters from both variables.

- *Cluster-level predictor.* The estimate is the slope of the regression
  of the cluster means of the outcome on the predictor, and for a
  normally distributed predictor the expected standard error is exactly
  \$\$E\[\widehat{SE}(\hat\beta)\] = \frac{\sqrt{\tau^2 + \sigma^2 /
  m}}{\sigma_X \sqrt{K - 2}}.\$\$ Here the number of clusters is what
  buys precision: as \\m\\ grows the numerator approaches \\\tau\\ and
  not 0.

- *Individual-level predictor centered at its grand mean, or not
  centered.* The estimate pools the within-cluster and the
  between-cluster information about the one slope, each weighted by its
  precision, and to first order \$\$E\[\widehat{SE}(\hat\beta)\]^2
  \approx \left\[ \frac{\\K (m - 1) - 1\\\\ \sigma^2_W}{\sigma^2} +
  \frac{(K - 2)\\(\rho_X \sigma^2_X + \sigma^2_W / m)} {\tau^2 +
  \sigma^2 / m} \right\]^{-1},\$\$ which reduces to the two expressions
  above when either part is absent. It is never larger than the value
  for the within-cluster slope, and the gain is what the assumption of
  one slope within and between clusters buys; when that assumption is in
  doubt, plan for the cluster-mean-centered analysis.

In 4,000 simulated data sets for each of six designs (clusters of 5 to
20, 9 to 44 clusters, intraclass correlations of .05 to .30), the mean
width of the interval was within 0.3 percent of the expected width in
the four designs with an exact expression, and 0.5 and 1.3 percent
narrower than it in the two with a pooled slope, so the first-order
expression errs toward a slightly larger number of clusters.
[`ss_aipe_mixed_effects_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md)
reports the mean width and the coverage for any design.

**Scope.** One predictor, a random intercept, a fixed slope, and
clusters of equal size. For random slopes, cross-level interactions, or
unequal cluster sizes the variance changes and a Monte Carlo planner
should be used (Schoemann, Boulton, & Short, 2017). The two-level
random-intercept model is treated at length by Snijders and Bosker
(2012).

## References

Schoemann, A. M., Boulton, A. J., & Short, S. D. (2017). Determining
power and sample size for simple and complex mediation models. *Social
Psychological and Personality Science, 8*(4), 379–386.
[doi:10.1177/1948550617715068](https://doi.org/10.1177/1948550617715068)

Snijders, T. A. B., & Bosker, R. J. (2012). *Multilevel analysis: An
introduction to basic and advanced multilevel modeling* (2nd ed.). Sage.

## See also

[`ss_aipe_mixed_effects_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md)
to check a plan by simulation;
[`ss_power_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md),
[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_change()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md),
[`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects_sensitivity.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. An individual-level predictor centered at its cluster mean (the
#    default): clusters of 20, intraclass correlation .10, unit
#    variances, and a desired width of .20.
ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
                      width = 0.20, cluster_size = 20)
#>  term                 value
#>  necessary_n_clusters 19   
#>  total_N              380  
#>  width                0.200
#>  expected_width       0.197
#>  df                   360  
#>  true_icc             0.100
#>  true_icc_X           0    
#>  cluster_size         20   
#>  predictor_level      1    
#> 
#> Confidence level: 95%

# 2. The same predictor not centered, when 30 percent of its variance
#    lies between clusters.
ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
                      width = 0.20, cluster_size = 20,
                      centering = "none", true_icc_X = 0.30)
#>  term                 value
#>  necessary_n_clusters 23   
#>  total_N              460  
#>  width                0.200
#>  expected_width       0.200
#>  df                   436  
#>  true_icc             0.100
#>  true_icc_X           0.300
#>  cluster_size         20   
#>  predictor_level      1    
#> 
#> Confidence level: 95%

# 3. A cluster-level predictor: far more clusters are needed for the
#    same width, because each cluster supplies one value of it.
ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
                      width = 0.20, cluster_size = 20,
                      predictor_level = 2, centering = "grand_mean")
#>  term                 value
#>  necessary_n_clusters 61   
#>  total_N              1220 
#>  width                0.200
#>  expected_width       0.198
#>  df                   59   
#>  true_icc             0.100
#>  true_icc_X           0    
#>  cluster_size         20   
#>  predictor_level      2    
#> 
#> Confidence level: 95%
```
