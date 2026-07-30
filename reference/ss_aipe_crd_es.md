# Find Target Sample Sizes for the Accuracy in Standardized Conditions Means Estimation in CRD

Find target sample sizes (the number of clusters, cluster size, or both)
for the accuracy in standardized conditions means estimation in CRD. If
users wish to seek for both types of sample sizes simultaneously, an
additional constraint is required, such as a desired width or a desired
budget. This function uses the likelihood-based confidence interval
(Cheung, 2009) by the `OpenMx` package (Boker et al., 2011). See further
details at Pornprasertmanit and Schneider (2014).

## Usage

``` r
ss_aipe_crd_es_n_clusters_fixed_width(
  width,
  n_individuals,
  es,
  es_type = 1,
  icc_Y,
  pr_treat,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  nrep = 1000,
  icc_Z = NULL,
  seed = NULL,
  multicore = FALSE,
  num_proc = NULL,
  clus_cost = NULL,
  indiv_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_es_n_individuals_fixed_width(
  width,
  n_clusters,
  es,
  es_type = 1,
  icc_Y,
  pr_treat,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  nrep = 1000,
  icc_Z = NULL,
  seed = NULL,
  multicore = FALSE,
  num_proc = NULL,
  clus_cost = NULL,
  indiv_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_es_n_clusters_fixed_budget(
  budget,
  n_individuals,
  clus_cost,
  indiv_cost,
  nrep = NULL,
  pr_treat = NULL,
  icc_Y = NULL,
  es = NULL,
  es_type = 1,
  num_predictors = 0,
  icc_Z = NULL,
  R2_within = NULL,
  R2_between = NULL,
  assurance = NULL,
  seed = NULL,
  multicore = FALSE,
  num_proc = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_es_n_individuals_fixed_budget(
  budget,
  n_clusters,
  clus_cost,
  indiv_cost,
  nrep = NULL,
  pr_treat = NULL,
  icc_Y = NULL,
  es = NULL,
  es_type = 1,
  num_predictors = 0,
  icc_Z = NULL,
  R2_within = NULL,
  R2_between = NULL,
  assurance = NULL,
  seed = NULL,
  multicore = FALSE,
  num_proc = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_es_both_fixed_budget(
  budget,
  clus_cost = 0,
  indiv_cost = 1,
  es,
  es_type = 1,
  icc_Y,
  pr_treat,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  nrep = 1000,
  icc_Z = NULL,
  seed = NULL,
  multicore = FALSE,
  num_proc = NULL,
  diff_size = NULL
)

ss_aipe_crd_es_both_fixed_width(
  width,
  clus_cost = 0,
  indiv_cost = 1,
  es,
  es_type = 1,
  icc_Y,
  pr_treat,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  nrep = 1000,
  icc_Z = NULL,
  seed = NULL,
  multicore = FALSE,
  num_proc = NULL,
  diff_size = NULL
)
```

## Arguments

- width:

  The desired width of the confidence interval of the unstandardized
  means difference

- n_individuals:

  The number of individuals in each cluster (cluster size)

- es:

  The amount of effect size

- es_type:

  The type of effect size. There are only three possible options: 0 =
  the effect size using total standard deviation, 1 = the effect size
  using the individual-level standard deviation (level 1), 2 = the
  effect size using the cluster-level standard deviation (level 2)

- icc_Y:

  The intraclass correlation of the dependent variable

- pr_treat:

  The proportion of treatment clusters

- R2_between:

  The proportion of variance explained in the between level (used when
  `covariate = TRUE`)

- R2_within:

  The proportion of variance explained in the within level (used when
  `covariate = TRUE`)

- num_predictors:

  The number of predictors used in the between level

- assurance:

  The degree of assurance, which is the value with which confidence can
  be placed that describes the likelihood of obtaining a confidence
  interval less than the value specified (e.g., .80, .90, .95)

- conf_level:

  The desired level of confidence for the confidence interval

- nrep:

  The number of replications used in a priori Monte Carlo simulation

- icc_Z:

  The intraclass correlation of the covariate (used when
  `covariate = TRUE`). If `icc_Z = 0`, the within-level covariate will
  be only used. If `icc_Z = 1`, the between-level covariate will be only
  used

- seed:

  An optional integer seed for the a priori Monte Carlo simulation. The
  default `NULL` uses the current state of the random number generator
  and leaves it unchanged, so repeated calls reflect the genuine
  sampling variability of the simulation. Supply an integer for
  reproducible results, in which case the generator state is restored on
  exit

- multicore:

  Use multiple processors within a computer. Specify as `TRUE` to use it

- num_proc:

  The number of processors to be used when `multicore = TRUE`. If it is
  not specified, the package will use the maximum number of processors
  in a machine

- clus_cost:

  The cost of collecting a new cluster regardless of the number of
  individuals collected in each cluster

- indiv_cost:

  The cost of collecting a new individual

- diff_size:

  Difference cluster size specification. The difference in cluster sizes
  can be specified in two ways. First, users may specify cluster size as
  integers, which can be negative or positive. The resulting cluster
  sizes will be based on the estimated cluster size adding by the
  specified vectors. For example, if the cluster size is 25, the number
  of clusters is 10, and the specified different cluster size is
  `c(-1, 0, 1)`, the cluster sizes will be 24, 25, 26, 24, 25, 26, 24,
  25, 26, and 24. Second, users may specify cluster size as positive
  decimals. The resulting cluster size will be based on the estimated
  cluster size multiplied by the specified vectors. For example, if the
  cluster size is 25, the number of clusters is 10, and the specified
  different cluster size is `c(-1, 0, 1)`, the cluster sizes will be 24,
  25, 26, 24, 25, 26, 24, 25, 26, and 24. If `NULL`, the cluster size is
  equal across clusters

- n_clusters:

  The desired number of clusters

- budget:

  The desired amount of budget

## Value

The `ss_aipe_crd_es_n_clusters_fixed_width` and
`ss_aipe_crd_es_n_clusters_fixed_budget` functions provide the number of
clusters. The `ss_aipe_crd_es_n_individuals_fixed_width` and
`ss_aipe_crd_es_n_individuals_fixed_budget` functions provide the
cluster size. The `ss_aipe_crd_es_both_fixed_budget` and
`ss_aipe_crd_es_both_fixed_width` provide the number of clusters and the
cluster size, respectively.

## Details

Here are the functions' descriptions:

- `ss_aipe_crd_es_n_clusters_fixed_width`:

  Find the number of clusters given a specified width of the confidence
  interval and the cluster size

- `ss_aipe_crd_es_n_individuals_fixed_width`:

  Find the cluster size given a specified width of the confidence
  interval and the number of clusters

- `ss_aipe_crd_es_n_clusters_fixed_budget`:

  Find the number of clusters given a budget and the cluster size

- `ss_aipe_crd_es_n_individuals_fixed_budget`:

  Find the cluster size given a budget and the number of clusters

- `ss_aipe_crd_es_both_fixed_budget`:

  Find the sample size combinations (the number of clusters and that
  cluster size) providing the narrowest confidence interval given the
  fixed budget

- `ss_aipe_crd_es_both_fixed_width`:

  Find the sample size combinations (the number of clusters and that
  cluster size) providing the lowest cost given the specified width of
  the confidence interval

## References

Boker, S. M., Neale, M. C., Maes, H. H., Wilde, M., Spiegel, M., Brick,
T. R., ... Fox, J. (2011). OpenMx: An open source extended structural
equation modeling framework. *Psychometrika, 76*(2), 306–317.
[doi:10.1007/s11336-010-9200-6](https://doi.org/10.1007/s11336-010-9200-6)

Cheung, M. W.-L. (2009). Constructing approximate confidence intervals
for parameters with structural equation models. *Structural Equation
Modeling, 16*(2), 267–294.
[doi:10.1080/10705510902751291](https://doi.org/10.1080/10705510902751291)

Pornprasertmanit, S., & Schneider, W. J. (2010). *Efficient sample size
for power and desired accuracy in Cohen's d estimation in two-group
cluster randomized design* (Master Thesis). Illinois State University,
Normal, IL.

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
# Every planner here runs an OpenMx likelihood-based Monte Carlo at each
# step of a sample size search. How long that takes depends on what is
# being searched over. Searching over cluster size, or over a budget, is
# cheap and those calls run under automated checking. Searching over the
# NUMBER OF CLUSTERS to hit a target confidence interval width evaluates
# many candidate designs and takes minutes to tens of minutes even at a
# small nrep, so those two calls are shown but not run; the package's tests
# exercise them.
# \donttest{
# Cluster size needed for a target width, given the number of clusters.
ss_aipe_crd_es_n_individuals_fixed_width(width = 0.3, 250, es = 0.5,
  es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 20)
#>  term                                      value
#>  cluster_size                              2    
#>  exp_width_of_individual-level_effect_size 0.236
#> 
#> Confidence level: 95%

# The same design questions under a budget rather than a target width.
ss_aipe_crd_es_n_clusters_fixed_budget(budget = 1000, n_individuals = 20,
  clus_cost = 0, indiv_cost = 1, nrep = 20, pr_treat = 0.5, icc_Y = 0.25, es = 0.5)
#>  term                                      value
#>  necessary_n_clusters                      50   
#>  exp_width_of_individual-level_effect_size 0.338
#>  budget                                    1000 
#> 
#> Confidence level: 95%

ss_aipe_crd_es_n_individuals_fixed_budget(budget = 1000, n_clusters = 200,
  clus_cost = 0, indiv_cost = 1, nrep = 20, pr_treat = 0.5, icc_Y = 0.25, es = 0.5)
#>  term                                      value
#>  cluster_size                              5    
#>  exp_width_of_individual-level_effect_size 0.203
#>  budget                                    1000 
#> 
#> Confidence level: 95%

# Both the number of clusters and the cluster size, under a budget.
ss_aipe_crd_es_both_fixed_budget(budget = 1000, clus_cost = 5, indiv_cost = 1, es = 0.5,
  es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 20)
#>  term                                      value
#>  necessary_n_clusters                      111  
#>  cluster_size                              4    
#>  exp_width_of_individual-level_effect_size 0.274
#>  budget                                    999  
#> 
#> Confidence level: 95%
# }

if (FALSE) { # \dontrun{
# Number of clusters needed for a target width, given the cluster size.
ss_aipe_crd_es_n_clusters_fixed_width(width = 0.3, n_individuals = 20, es = 0.5,
  es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 20)

# Both quantities under a target width.
ss_aipe_crd_es_both_fixed_width(width = 0.5, clus_cost = 5, indiv_cost = 1, es = 0.5,
  es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 20)

# Unequal cluster sizes: diff_size gives each cluster's deviation from
# n_individuals (additive) or its multiplicative factor.
ss_aipe_crd_es_n_clusters_fixed_width(width = 0.3, n_individuals = 20, es = 0.5,
  es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 20,
  diff_size = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))

ss_aipe_crd_es_n_clusters_fixed_width(width = 0.3, n_individuals = 20, es = 0.5,
  es_type = 1, icc_Y = 0.25, pr_treat = 0.5, nrep = 20,
  diff_size = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
} # }
```
