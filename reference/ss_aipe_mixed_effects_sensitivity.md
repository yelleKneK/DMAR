# Sensitivity Analysis for Sample Size Planning From the AIPE Perspective for a Mixed-Effects Fixed Effect

Simulates the design that
[`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md)
plans and reports what its confidence interval for the fixed-effect
slope delivers: the mean width, the share of intervals no wider than
desired, and the coverage. The two functions take the same description
of the predictor (`predictor_level`, `centering`, the intraclass
correlation of the predictor), use the same expressions, and form the
same interval, so a plan made by the first can be checked by the second,
at the planning values or at population values that differ from them.

On each replication the function draws `n_clusters` clusters (i.e.,
Level 2 units) of `cluster_size` individuals (i.e., Level 1 units) from
\$\$y\_{ki} = \beta x\_{ki} + u_k + e\_{ki},\$\$ with \\u_k \sim N(0,
\rho V)\\ and \\e\_{ki} \sim N(0, (1 - \rho) V)\\, where \\\beta\\ is
`true_reg_coef`, \\V\\ is `true_total_var`, and \\\rho\\ is `true_icc`.
An individual-level (i.e., Level 1) predictor is \\x\_{ki} = b_k +
w\_{ki}\\, with variance \\\rho_X \sigma^2_X\\ between clusters and
\\(1 - \rho_X)\sigma^2_X\\ within them (`true_var_X`, `true_icc_X`); a
cluster-level (i.e., Level 2) predictor is \\x_k \sim N(0,
\sigma^2_X)\\, the same for every individual of a cluster. The predictor
is then centered as `centering` says, the model is fit, and the interval
\\\hat\beta \pm t\_{1 - \alpha/2,\\ \nu}\\\widehat{SE}(\hat\beta)\\ is
recorded, with the degrees of freedom \\\nu\\ of
[`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md).

## Usage

``` r
ss_aipe_mixed_effects_sensitivity(
  true_total_var = NULL,
  true_var_X = NULL,
  true_icc = NULL,
  true_reg_coef = 0,
  supposed_total_var = NULL,
  supposed_var_X = NULL,
  supposed_icc = NULL,
  width,
  cluster_size = 20L,
  specified_n_clusters = NULL,
  conf_level = 0.95,
  predictor_level = 1,
  centering = c("cluster_mean", "grand_mean", "none"),
  true_icc_X = 0,
  supposed_icc_X = true_icc_X,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_total_var:

  Population variance of the outcome given the predictor, the variance
  of the random intercept plus the residual variance: the value the data
  are generated from.

- true_var_X:

  Population variance of the predictor the data are generated from: over
  individuals for an individual-level predictor (between and within
  clusters together), over clusters for a cluster-level predictor.

- true_icc:

  Population intraclass correlation of the outcome given the predictor
  (the share of `true_total_var` between clusters).

- true_reg_coef:

  Population value of the fixed-effect slope, the regression coefficient
  of the outcome on the predictor (\\\beta\\ in the model above; default
  `0`). The width of the interval does not depend on it; the coverage is
  evaluated against it.

- supposed_total_var, supposed_var_X, supposed_icc:

  Supposed values, the planning values passed to
  [`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md)
  (not estimates computed from data), which may differ from the
  population values above; supply all three, or `specified_n_clusters`,
  but not both.

- width:

  Desired full width of the confidence interval for the slope.

- cluster_size:

  Number of individuals in each cluster (the same in every cluster).

- specified_n_clusters:

  Number of clusters to evaluate, in place of the number the planning
  values give.

- conf_level:

  Confidence level (default `0.95`).

- predictor_level, centering:

  The level of the predictor (`1`, the default, or `2`) and how it is
  centered (`"cluster_mean"`, `"grand_mean"`, or `"none"`), as in
  [`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md).
  The same values are used to plan and to analyze each simulated data
  set.

- true_icc_X:

  Population intraclass correlation of an individual-level predictor:
  the share of `true_var_X` between clusters in the data generated.
  Default `0`. Not used with `predictor_level = 2`.

- supposed_icc_X:

  The supposed value of `true_icc_X`, the planning value passed to
  [`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md).
  Default: `true_icc_X`.

- G:

  Number of Monte Carlo replications.

- print_iter:

  Logical.

- filename:

  Optional path for a comma separated file recording every replication
  (the fixed effect estimate, the two confidence limits, the interval
  width, and two indicators of whether the interval missed
  `true_reg_coef` below or above): nothing is written when `filename` is
  `NULL` (the default), a new file with a header row is created
  otherwise, an existing file at that path is appended to, and a
  throwaway run should point it at `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with rows for the mean, median, and standard deviation of
the estimated slope (`mean_reg_coef`, `median_reg_coef`, `sd_reg_coef`)
and of the interval width, the share of intervals no wider than `width`
(`pct_ci_less_width`), the shares that missed `true_reg_coef` below and
above and their sum, the design simulated (`n_clusters`, `cluster_size`,
`total_N`, `predictor_level`), the `expected_width` that
[`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md)
gives for that design at the population values, and the inputs. The
centering is recorded in `attr(x, "centering")`.

## Details

**How each case is fit.** Two of the three analyses have closed forms in
a balanced design and are computed directly. The slope of an
individual-level predictor centered at its cluster mean is the
within-cluster regression slope, with the residual mean square on
\\K(m - 1) - 1\\ degrees of freedom. The slope of a cluster-level
predictor is the slope of the regression of the \\K\\ cluster means of
the outcome on the predictor, on \\K - 2\\ degrees of freedom. In both,
restricted maximum likelihood for the random-intercept model gives the
same estimate and standard error when its estimate of the intercept
variance is positive. An individual-level predictor that is centered at
its grand mean or not centered has no closed form, because its weights
depend on the estimated variance components; it is fit by restricted
maximum likelihood with lme4, which must then be installed, and its
interval uses \\K(m - 1) - 1\\ degrees of freedom.

**What to look at.** `mean_ci_width` against `width` says whether the
plan delivers the width on average, and `pct_ci_less_width` how often.
`total_type_I_error` against \\1 - \\ `conf_level` says whether the
interval covers as it should; with few clusters and a cluster-level
predictor this is where a problem would show. Replications in which no
interval could be formed are counted in the record of the simulation.

## The record of the simulation

Every sensitivity function of the package summarizes its replications
over one denominator: the replications that produced a confidence
interval (for a power sensitivity function, a test). The proportion of
intervals that met the width target (`pct_ci_less_width`), the
proportions of intervals that missed the true value, and the summaries
of the widths are all taken over those replications. A replication
without an interval is in none of these rows. It is counted, and the
returned table carries the counts in its attribute `simulation`, a list
read with `attr(x, "simulation")`:

- `attempted`:

  the number of replications asked for (`G`).

- `returned`:

  the number that produced an interval (a test).

- `failed`:

  the number that did not (`attempted - returned`).

- `denominator`:

  the number of replications the rows named above are taken over, which
  is `returned`.

- `width_target_probability`:

  the proportion of the replications with an interval whose interval met
  the width target, the value of the `pct_ci_less_width` row; `NA` for a
  function or a call without a width target, and when no replication
  produced an interval.

- `width_target_mc_ci`:

  the exact binomial (Clopper-Pearson) limits for that proportion, at
  the confidence level of the table. They describe the error of the
  simulation, which shrinks as `G` grows, and nothing else.

A function that records more about its run adds elements to the list,
and its help page names them. The line "Replications" beneath the
printed table gives `G`, and the three counts when a replication was
without an interval. With many replications without an interval the rows
describe only the studies in which the interval could be formed, so
`failed` is to be reported beside them.

## References

McNeish, D., & Kelley, K. (2019). Fixed effects versus mixed effects
models for clustered data: Reviewing the approaches, disentangling the
differences, and making recommendations. *Psychological Methods, 24*,
20–35. [doi:10.1037/met0000182](https://doi.org/10.1037/met0000182)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 15 and 16 on mixed-effects models.)

## See also

[`ss_aipe_mixed_effects`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`icc_lmer`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other AIPE sample size planning:
[`ss_aipe_c_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_sensitivity.md),
[`ss_aipe_cliff_delta()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta.md),
[`ss_aipe_cliff_delta_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cliff_delta_sensitivity.md),
[`ss_aipe_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md),
[`ss_aipe_crossed()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed.md),
[`ss_aipe_crossed_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_crossed_sensitivity.md),
[`ss_aipe_equivalence_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r.md),
[`ss_aipe_equivalence_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_r_sensitivity.md),
[`ss_aipe_equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd.md),
[`ss_aipe_equivalence_smd_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_equivalence_smd_sensitivity.md),
[`ss_aipe_icc()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md),
[`ss_aipe_icc_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc_sensitivity.md),
[`ss_aipe_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md),
[`ss_aipe_indirect_effect_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect_sensitivity.md),
[`ss_aipe_mean()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean.md),
[`ss_aipe_mean_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mean_sensitivity.md),
[`ss_aipe_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md),
[`ss_aipe_omega_squared_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared_sensitivity.md),
[`ss_aipe_partial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r.md),
[`ss_aipe_partial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_partial_r_sensitivity.md),
[`ss_aipe_pcm_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_pcm_sensitivity.md),
[`ss_aipe_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r.md),
[`ss_aipe_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_r_sensitivity.md),
[`ss_aipe_reliability_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md),
[`ss_aipe_semipartial_r()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_aipe_semipartial_r_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md),
[`ss_aipe_two_phase()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase.md),
[`ss_aipe_two_phase_cost()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_cost.md),
[`ss_aipe_two_phase_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_two_phase_sensitivity.md)

Other mixed models:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`R2_mixed_effects_decomposition()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects_decomposition.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Plan, then check the plan (few replications, for a fast example).
ss_aipe_mixed_effects(true_total_var = 1, true_var_X = 1, true_icc = 0.10,
                      width = 0.40, cluster_size = 10)
#>  term                 value
#>  necessary_n_clusters 10   
#>  total_N              100  
#>  width                0.400
#>  expected_width       0.400
#>  df                   89   
#>  true_icc             0.100
#>  true_icc_X           0    
#>  cluster_size         10   
#>  predictor_level      1    
#> 
#> Confidence level: 95%
ss_aipe_mixed_effects_sensitivity(
  true_total_var = 1, true_var_X = 1, true_icc = 0.10,
  true_reg_coef = 0.30,
  supposed_total_var = 1, supposed_var_X = 1, supposed_icc = 0.10,
  width = 0.40, cluster_size = 10,
  G = 50, print_iter = FALSE, seed = 113
)
#>  term               value
#>  mean_reg_coef      0.311
#>  median_reg_coef    0.313
#>  sd_reg_coef        0.112
#>  mean_ci_width      0.398
#>  median_ci_width    0.393
#>  sd_ci_width        0.039
#>  pct_ci_less_width  0.560
#>  pct_ci_miss_low    0.020
#>  pct_ci_miss_high   0.020
#>  total_type_I_error 0.040
#>  n_clusters         10   
#>  cluster_size       10   
#>  total_N            100  
#>  predictor_level    1    
#>  expected_width     0.400
#>  true_total_var     1    
#>  true_var_X         1    
#>  true_icc           0.100
#>  true_icc_X         0    
#>  true_reg_coef      0.300
#>  supposed_total_var 1    
#>  supposed_var_X     1    
#>  supposed_icc       0.100
#>  supposed_icc_X     0    
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%

# A cluster-level predictor, at a number of clusters chosen by hand.
ss_aipe_mixed_effects_sensitivity(
  true_total_var = 1, true_var_X = 1, true_icc = 0.10,
  true_reg_coef = 0.30,
  specified_n_clusters = 25, cluster_size = 10, width = 0.40,
  predictor_level = 2, centering = "grand_mean",
  G = 50, print_iter = FALSE, seed = 113
)
#>  term               value
#>  mean_reg_coef      0.287
#>  median_reg_coef    0.304
#>  sd_reg_coef        0.090
#>  mean_ci_width      0.380
#>  median_ci_width    0.370
#>  sd_ci_width        0.084
#>  pct_ci_less_width  0.680
#>  pct_ci_miss_low    0.020
#>  pct_ci_miss_high   0.020
#>  total_type_I_error 0.040
#>  n_clusters         25   
#>  cluster_size       10   
#>  total_N            250  
#>  predictor_level    2    
#>  expected_width     0.376
#>  true_total_var     1    
#>  true_var_X         1    
#>  true_icc           0.100
#>  true_icc_X         0    
#>  true_reg_coef      0.300
#>  supposed_total_var <NA> 
#>  supposed_var_X     <NA> 
#>  supposed_icc       <NA> 
#>  supposed_icc_X     <NA> 
#>  width              0.400
#>  conf_level         0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
