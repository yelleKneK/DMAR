# Sample Sizes for a Cluster Randomized Design: Accuracy in Parameter Estimation for the Difference Between Condition Means

Finds the number of clusters, the number of individuals per cluster, or
both, for a two-condition cluster randomized design (CRD), so that the
confidence interval for the difference between the condition means has a
desired width or is as narrow as a budget allows. Four of the functions
fix one of the two sizes and take a desired width or a budget; the two
`both` functions take a budget and return the sizes that give the
narrowest interval, or take a desired width and return the least costly
sizes.

## Usage

``` r
ss_aipe_crd_n_clusters_fixed_width(
  width,
  cluster_size,
  pr_treat,
  true_var_between = NULL,
  true_var_within = NULL,
  true_total_var = NULL,
  true_icc = NULL,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  cluster_cost = NULL,
  individual_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_n_individuals_fixed_width(
  width,
  n_clusters,
  pr_treat,
  true_var_between = NULL,
  true_var_within = NULL,
  true_total_var = NULL,
  true_icc = NULL,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  cluster_cost = NULL,
  individual_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_n_clusters_fixed_budget(
  budget,
  cluster_size,
  cluster_cost = 0,
  individual_cost = 1,
  pr_treat = NULL,
  true_var_between = NULL,
  true_var_within = NULL,
  true_total_var = NULL,
  true_icc = NULL,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_n_individuals_fixed_budget(
  budget,
  n_clusters,
  cluster_cost = 0,
  individual_cost = 1,
  pr_treat = NULL,
  true_var_between = NULL,
  true_var_within = NULL,
  true_total_var = NULL,
  true_icc = NULL,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_both_fixed_budget(
  budget,
  cluster_cost = 0,
  individual_cost = 1,
  pr_treat,
  true_var_between = NULL,
  true_var_within = NULL,
  true_total_var = NULL,
  true_icc = NULL,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_both_fixed_width(
  width,
  cluster_cost = 0,
  individual_cost = 1,
  pr_treat,
  true_var_between = NULL,
  true_var_within = NULL,
  true_total_var = NULL,
  true_icc = NULL,
  true_R2_between = 0,
  true_R2_within = 0,
  n_covariates = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)
```

## Arguments

- width:

  The desired width of the confidence interval of the unstandardized
  means difference

- cluster_size:

  The number of individuals in each cluster: a single whole number of at
  least 2. The searches over the cluster size start at 2, and with one
  individual in a cluster the variance between clusters cannot be told
  from the variance within them

- pr_treat:

  The proportion of treatment clusters: a single number strictly between
  0 and 1

- true_var_between:

  The population residual variance in the between level before
  accounting for the covariate (a planning value), 0 or more

- true_var_within:

  The population residual variance in the within level before accounting
  for the covariate (a planning value), 0 or more

- true_total_var:

  The population total residual variance before accounting for the
  covariate (a planning value), positive. Any two of `true_var_between`,
  `true_var_within`, `true_total_var`, and `true_icc` determine the
  other two; a third or fourth must agree with them (see the section on
  refused calls)

- true_icc:

  The population intraclass correlation of the dependent variable (a
  planning value), the share of `true_total_var` that lies between
  clusters, in \[0, 1\]

- true_R2_between:

  The population proportion of the between-level variance that the
  cluster-level covariates explain (a planning value in \[0, 1); default
  0). It must be 0 when `n_covariates = 0`, since there is then no
  covariate at the cluster level to explain it

- true_R2_within:

  The population proportion of the within-level variance that the
  covariates explain (a planning value in \[0, 1); default 0, no
  covariate)

- n_covariates:

  The number of covariates at the cluster (between) level, not counting
  the treatment indicator: a single whole number of 0 or more. The
  degrees of freedom of the interval are
  `n_clusters - 2 - n_covariates`. The default, 0, is a model with the
  treatment indicator alone. It is not `p`, the total number of
  predictors, which here is `p = n_covariates + 1`.

- assurance:

  The degree of assurance, which is the value with which confidence can
  be placed that describes the likelihood of obtaining a confidence
  interval less than the value specified (e.g., .80, .90, .95)

- conf_level:

  The desired level of confidence for the confidence interval

- cluster_cost:

  The cost of adding one cluster, whatever the number of individuals in
  it: a single finite number of 0 or more. Where the costs choose the
  design (the two budget functions that fix a size and the two `both`
  functions), `cluster_cost` and `individual_cost` cannot both be 0

- individual_cost:

  The cost of adding one individual to a cluster: a single finite number
  of 0 or more, and positive where the costs choose the cluster size
  (`ss_aipe_crd_n_individuals_fixed_budget` and the two `both`
  functions), since otherwise clusters of every size cost the same

- diff_size:

  Difference cluster size specification. The differences in cluster
  sizes can be specified in two ways, and the specified vector is
  recycled across the clusters. First, users may specify differences as
  integers, which can be negative or positive; the resulting cluster
  sizes add the specified values to the estimated cluster size. For
  example, if the cluster size is 25, the number of clusters is 10, and
  `diff_size = c(-1, 0, 1)`, the cluster sizes will be 24, 25, 26, 24,
  25, 26, 24, 25, 26, and 24. Second, users may specify multipliers of
  the cluster size as positive decimals; at least one value must be
  non-integer, which is what selects the multiplicative form. The
  resulting cluster sizes multiply the estimated cluster size by the
  specified values and round to the nearest integer. For example, if the
  cluster size is 25, the number of clusters is 10, and
  `diff_size = c(0.8, 1, 1.2)`, the cluster sizes will be 20, 25, 30,
  20, 25, 30, 20, 25, 30, and 20. In either form a resulting cluster
  size below 1 is set to 1. If `NULL`, the cluster size is equal across
  clusters. Anything other than `NULL` or a vector of one or more finite
  numbers is refused

- n_clusters:

  The number of clusters: a single whole number of at least
  `3 + n_covariates`, since the interval has
  `n_clusters - 2 - n_covariates` degrees of freedom. With `pr_treat` it
  must put a cluster in each condition (`round(n_clusters * pr_treat)`
  of them are treated)

- budget:

  The amount available for the study, in the units of `cluster_cost` and
  `individual_cost`: a single finite number of 0 or more. A budget that
  does not buy the smallest design the function considers is refused
  (see the section on refused calls)

## Value

The `ss_aipe_crd_n_clusters_fixed_width` and
`ss_aipe_crd_n_clusters_fixed_budget` functions provide the number of
clusters. The `ss_aipe_crd_n_individuals_fixed_width` and
`ss_aipe_crd_n_individuals_fixed_budget` functions provide the cluster
size. The `ss_aipe_crd_both_fixed_budget` and
`ss_aipe_crd_both_fixed_width` provide the number of clusters and the
cluster size, respectively.

Each returns a table with the columns `term` and `value`. The rows,
where they apply, are `necessary_n_clusters` (the number of clusters,
where it is solved for), `necessary_cluster_size` (the number of
individuals in each cluster, where it is solved for; with `diff_size`
the cluster sizes and their frequencies are printed ahead of the table
instead), `width` (the desired width, in the functions that plan from
one), the width of the interval the design gives, and `budget` (the cost
of the design, when the two costs are supplied). Without `assurance` the
width row is `expected_width`, the expected width of the interval. With
`assurance` it is `assured_width`, the width that the interval does not
exceed with the requested probability, followed by the row `assurance`
that echoes that probability.

A size that was given is an argument of the call and is not repeated in
the table: only a size the function solved for has a row, and its name
begins with `necessary_`.
[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) report the
planned size as `estimate` (the number of clusters where it was solved
for, otherwise the cluster size) and the desired width as `width`, which
is `NA` for a plan made from a budget; `glance()` carries
`expected_width` or `assured_width` as a column of that name.

## Details

Here are the functions' descriptions:

- `ss_aipe_crd_n_clusters_fixed_width`:

  Find the number of clusters given a specified width of the confidence
  interval and the cluster size

- `ss_aipe_crd_n_individuals_fixed_width`:

  Find the cluster size given a specified width of the confidence
  interval and the number of clusters

- `ss_aipe_crd_n_clusters_fixed_budget`:

  Find the number of clusters given a budget and the cluster size

- `ss_aipe_crd_n_individuals_fixed_budget`:

  Find the cluster size given a budget and the number of clusters

- `ss_aipe_crd_both_fixed_budget`:

  Find the sample size combinations (the number of clusters and that
  cluster size) providing the narrowest confidence interval given the
  fixed budget

- `ss_aipe_crd_both_fixed_width`:

  Find the sample size combinations (the number of clusters and that
  cluster size) providing the lowest cost given the specified width of
  the confidence interval

The searches over the number of clusters consider at least
`ceiling(2 / min(pr_treat, 1 - pr_treat)) + n_covariates` clusters,
which puts at least two clusters in each condition after rounding (the
rule of MBESS; Kelley, 2007a, 2007b) and gives the interval at least two
degrees of freedom: 4 clusters at `pr_treat = 0.5`, 20 at 0.1, and 200
at 0.01. From there `ss_aipe_crd_n_clusters_fixed_width` returns the
smallest number of clusters whose width is no larger than `width` (a
width equal to the target meets it, so asking for the width a plan
reports returns that plan), and the two `both` functions compare numbers
of clusters from there with cluster sizes of 2 or more. The searches
consider at most 1,000,000 clusters, the bound the power planners use
for a sample size: the width falls toward 0 as clusters are added, so
every positive target is met by some number, but a target too narrow for
any study would otherwise be searched for without end. The widths at
100, 200, ..., 1000 clusters bracket the target, and the bracket is
searched one cluster at a time, which gives the smallest number when the
width falls as clusters are added, as it does with equal cluster sizes
and no `assurance`. At an `assurance` below .5 the width need not fall
as the number of clusters grows, so the search may return more than the
smallest number that meets the target; the same can happen with unequal
cluster sizes (`diff_size`), whose width depends on how the sizes fall
between the conditions. The cluster size searches consider 2 to 100,000
individuals per cluster.

## Refused calls

Each refusal below is an error of class `dmar_argument_invalid` that
names the argument or arguments as `arg` and, where the cause needs a
name, gives it as `reason`.

- `n_covariates` must be a whole number of 0 or more, `cluster_size` a
  whole number of at least 2, and `n_clusters` a whole number of at
  least `3 + n_covariates`; `pr_treat` must be strictly between 0 and 1,
  and `true_R2_between` and `true_R2_within` in \[0, 1). A nonzero
  `true_R2_between` with `n_covariates = 0` is refused with `reason`
  `"R2_between_without_covariates"`; `true_R2_within` may be nonzero
  without them. A number of clusters that, at `pr_treat`, leaves a
  condition without a cluster is refused with `reason`
  `"arm_without_clusters"`. `diff_size` must be `NULL` or a vector of
  one or more finite numbers.

- Each cost, and the budget, must be a single finite number of 0 or
  more. Costs that are both 0 are refused where the costs choose the
  design (`reason` `"costs_all_zero"`), and `individual_cost = 0` where
  they choose the cluster size (`reason` `"individual_cost_zero"`): a
  search over designs that all cost the same has no end. In the two
  functions that plan from a width and a fixed size the costs only price
  the design, and costs of 0 give a budget of 0.

- A budget that does not buy the smallest design the function can return
  is refused with `reason` `"budget_infeasible"` and the cost of that
  design as `min_cost`: in `ss_aipe_crd_n_clusters_fixed_budget` 4
  clusters of `cluster_size`, or, with `pr_treat`, the fewest clusters
  that put a cluster in each condition and leave a degree of freedom (6
  at `pr_treat = 0.1`); in `ss_aipe_crd_n_individuals_fixed_budget`
  `n_clusters` clusters of 2; and in `ss_aipe_crd_both_fixed_budget` the
  fewest clusters its search considers, given above, each of 2. When
  `ss_aipe_crd_n_clusters_fixed_budget` is given `pr_treat` and the
  budget buys 4 clusters or more but a design with no finite width, it
  is refused with `reason` `"arm_without_clusters"` (a condition without
  a cluster) or `"too_few_clusters"` (no degrees of freedom), with the
  budget that buys the fewest clusters that would do.

- `ss_aipe_crd_n_individuals_fixed_width` refuses a width that no
  cluster size up to 100,000 reaches with `n_clusters` clusters, with
  `reason` `"width_not_reachable"`; it gives as `width_limit` the width
  the interval falls toward as the clusters grow, which the variance
  between clusters sets, and says whether a larger cluster size would
  reach the target. `ss_aipe_crd_n_clusters_fixed_width` and
  `ss_aipe_crd_both_fixed_width` refuse a width that no design of up to
  1,000,000 clusters meets (with clusters of up to 100,000 individuals
  in the latter), with `reason` `"n_clusters_beyond_reach"` and the
  width there as `width_at_max`.

- The variances are needed where a width is computed: two of
  `true_var_between`, `true_var_within`, `true_total_var`, and
  `true_icc`, any two (`reason` `"underidentified_variances"` otherwise,
  as for a between-cluster variance of 0 with an intraclass correlation
  of 0). Each must be `NULL` or a single finite number (`reason`
  `"variance_input_invalid"`). A variance must be 0 or more, the total
  variance positive (whether it is given or is the sum of the two
  variances, so two variances of 0 are refused), and the intraclass
  correlation in \[0, 1\] (`reason` `"variance_out_of_range"`). A third
  or fourth value must agree with the two that determine the model to
  within about 1e-14 of its scale (`reason` `"inconsistent_variances"`,
  with the value they imply in the message), so a rounded value is
  refused; give two and leave the others `NULL`. With `diff_size` the
  variance within clusters must be positive (`reason`
  `"within_variance_zero"`). These refusals carry the four values as
  given in `supplied`.

## References

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
[doi:10.3758/BF03192993](https://doi.org/10.3758/BF03192993)

Pornprasertmanit, S., & Schneider, W. J. (2014). Accuracy in parameter
estimation in cluster randomized designs. *Psychological Methods,
19*(3), 356–379.
[doi:10.1037/a0037036](https://doi.org/10.1037/a0037036)

## See also

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Examples for each function
ss_aipe_crd_n_clusters_fixed_width(width = 0.3, cluster_size = 30,
  pr_treat = 0.5, true_var_between = 0.25, true_var_within = 0.75)
#>  term                 value
#>  necessary_n_clusters 191  
#>  width                0.300
#>  expected_width       0.299
#> 
#> Confidence level: 95%

ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 250,
  pr_treat = 0.5, true_var_between = 0.25, true_var_within = 0.75)
#>  term                   value
#>  necessary_cluster_size 7    
#>  width                  0.300
#>  expected_width         0.298
#> 
#> Confidence level: 95%

ss_aipe_crd_n_clusters_fixed_budget(budget = 10000, cluster_size = 20,
  cluster_cost = 20, individual_cost = 1)
#>  term                 value
#>  necessary_n_clusters 250  
#>  budget               10000
#> 
#> Confidence level: 95%

ss_aipe_crd_n_individuals_fixed_budget(budget = 10000, n_clusters = 30,
  cluster_cost = 20, individual_cost = 1,
  pr_treat = 0.5, true_var_between = 0.05, true_var_within = 0.95, assurance = 0.8)
#>  term                   value
#>  necessary_cluster_size 313  
#>  assured_width          0.380
#>  assurance              0.800
#>  budget                 9990 
#> 
#> Confidence level: 95%

ss_aipe_crd_both_fixed_budget(budget = 10000, cluster_cost = 30, individual_cost = 1,
  pr_treat = 0.5, true_var_between = 0.25, true_var_within = 0.75)
#>  term                   value
#>  necessary_n_clusters   250  
#>  necessary_cluster_size 10   
#>  expected_width         0.284
#>  budget                 10000
#> 
#> Confidence level: 95%

ss_aipe_crd_both_fixed_width(width = 0.3, cluster_cost = 0, individual_cost = 1,
  pr_treat = 0.5, true_var_between = 0.25, true_var_within = 0.75)
#>  term                   value
#>  necessary_n_clusters   430  
#>  necessary_cluster_size 2    
#>  width                  0.300
#>  expected_width         0.300
#>  budget                 860  
#> 
#> Confidence level: 95%

# Examples for different cluster size
set.seed(113)
ss_aipe_crd_n_clusters_fixed_width(width = 0.3, cluster_size = 30,
  pr_treat = 0.5, true_var_between = 0.25, true_var_within = 0.75,
  diff_size = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))
#>   cluster_size freq
#> 1           27   21
#> 2           28   22
#> 3           29   21
#> 4           30   63
#> 5           31   22
#> 6           32   21
#> 7           33   21
#>  term                 value
#>  necessary_n_clusters 191  
#>  width                0.300
#>  expected_width       0.299
#> 
#> Confidence level: 95%

# Examples for different number of clusters
ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 250,
  pr_treat = 0.5, true_var_between = 0.25, true_var_within = 0.75,
  diff_size = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
#>   cluster_size freq
#> 1            4   32
#> 2            6   62
#> 3            7   62
#> 4            8   63
#> 5           10   31
#>  term           value
#>  width          0.300
#>  expected_width 0.300
#> 
#> Confidence level: 95%
```
