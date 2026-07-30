# Find Target Sample Sizes for the Accuracy in Unstandardized Conditions Means Estimation in CRD

Find target sample sizes (the number of clusters, cluster size, or both)
for the accuracy in unstandardized conditions means estimation in CRD.
If users wish to seek for both types of sample sizes simultaneously, an
additional constraint is required, such as a desired width or a desired
budget.

## Usage

``` r
ss_aipe_crd_n_clusters_fixed_width(
  width,
  n_individuals,
  pr_treat,
  tau_Y = NULL,
  sigma2_Y = NULL,
  total_var = NULL,
  icc_Y = NULL,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  clus_cost = NULL,
  indiv_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_n_individuals_fixed_width(
  width,
  n_clusters,
  pr_treat,
  tau_Y = NULL,
  sigma2_Y = NULL,
  total_var = NULL,
  icc_Y = NULL,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  clus_cost = NULL,
  indiv_cost = NULL,
  diff_size = NULL
)

ss_aipe_crd_n_clusters_fixed_budget(
  budget,
  n_individuals,
  clus_cost = 0,
  indiv_cost = 1,
  pr_treat = NULL,
  tau_Y = NULL,
  sigma2_Y = NULL,
  total_var = NULL,
  icc_Y = NULL,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_n_individuals_fixed_budget(
  budget,
  n_clusters,
  clus_cost = 0,
  indiv_cost = 1,
  pr_treat = NULL,
  tau_Y = NULL,
  sigma2_Y = NULL,
  total_var = NULL,
  icc_Y = NULL,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_both_fixed_budget(
  budget,
  clus_cost = 0,
  indiv_cost = 1,
  pr_treat,
  tau_Y = NULL,
  sigma2_Y = NULL,
  total_var = NULL,
  icc_Y = NULL,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)

ss_aipe_crd_both_fixed_width(
  width,
  clus_cost = 0,
  indiv_cost = 1,
  pr_treat,
  tau_Y = NULL,
  sigma2_Y = NULL,
  total_var = NULL,
  icc_Y = NULL,
  R2_between = 0,
  R2_within = 0,
  num_predictors = 0,
  assurance = NULL,
  conf_level = 0.95,
  diff_size = NULL
)
```

## Arguments

- width:

  The desired width of the confidence interval of the unstandardized
  means difference

- n_individuals:

  The number of individuals in each cluster (cluster size)

- pr_treat:

  The proportion of treatment clusters

- tau_Y:

  The residual variance in the between level before accounting for the
  covariate

- sigma2_Y:

  The residual variance in the within level before accounting for the
  covariate

- total_var:

  The total residual variance before accounting for the covariate

- icc_Y:

  The intraclass correlation of the dependent variable

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

The `ss_aipe_crd_n_clusters_fixed_width` and
`ss_aipe_crd_n_clusters_fixed_budget` functions provide the number of
clusters. The `ss_aipe_crd_n_individuals_fixed_width` and
`ss_aipe_crd_n_individuals_fixed_budget` functions provide the cluster
size. The `ss_aipe_crd_both_fixed_budget` and
`ss_aipe_crd_both_fixed_width` provide the number of clusters and the
cluster size, respectively.

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

## References

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
ss_aipe_crd_n_clusters_fixed_width(width = 0.3, n_individuals = 30,
  pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#>  term                               value
#>  necessary_n_clusters               191  
#>  exp_width_of_unstd_conditions_diff 0.299
#> 
#> Confidence level: 95%

ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 250,
  pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#>  term                               value
#>  cluster_size                       7    
#>  exp_width_of_unstd_conditions_diff 0.298
#> 
#> Confidence level: 95%

ss_aipe_crd_n_clusters_fixed_budget(budget = 10000, n_individuals = 20,
  clus_cost = 20, indiv_cost = 1)
#>  term                 value
#>  necessary_n_clusters 250  
#>  budget               10000
#> 
#> Confidence level: 95%

ss_aipe_crd_n_individuals_fixed_budget(budget = 10000, n_clusters = 30,
  clus_cost = 20, indiv_cost = 1,
  pr_treat = 0.5, tau_Y = 0.05, sigma2_Y = 0.95, assurance = 0.8)
#>  term                                              value
#>  cluster_size                                      313  
#>  width_of_unstd_conditions_diff_with_0.8_assurance 0.38 
#>  budget                                            9990 
#> 
#> Confidence level: 95%

ss_aipe_crd_both_fixed_budget(budget = 10000, clus_cost = 30, indiv_cost = 1,
  pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#>  term                               value
#>  necessary_n_clusters               250  
#>  cluster_size                       10   
#>  exp_width_of_unstd_conditions_diff 0.284
#>  budget                             10000
#> 
#> Confidence level: 95%

ss_aipe_crd_both_fixed_width(width = 0.3, clus_cost = 0, indiv_cost = 1,
  pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75)
#>  term                               value
#>  necessary_n_clusters               430  
#>  cluster_size                       2    
#>  exp_width_of_unstd_conditions_diff 0.3  
#>  budget                             860  
#> 
#> Confidence level: 95%

# Examples for different cluster size
set.seed(113)
ss_aipe_crd_n_clusters_fixed_width(width = 0.3, n_individuals = 30,
  pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75,
  diff_size = c(-2, 1, 0, 2, -1, 3, -3, 0, 0))
#>   cluster_size freq
#> 1           27   21
#> 2           28   22
#> 3           29   21
#> 4           30   63
#> 5           31   22
#> 6           32   21
#> 7           33   21
#>  term                               value
#>  necessary_n_clusters               191  
#>  exp_width_of_unstd_conditions_diff 0.299
#> 
#> Confidence level: 95%

# Examples for different number of clusters
ss_aipe_crd_n_individuals_fixed_width(width = 0.3, n_clusters = 250,
  pr_treat = 0.5, tau_Y = 0.25, sigma2_Y = 0.75,
  diff_size = c(0.6, 1.2, 0.8, 1.4, 1, 1, 1.1, 0.9))
#>   cluster_size freq
#> 1            4   32
#> 2            6   62
#> 3            7   62
#> 4            8   63
#> 5           10   31
#>  term                               value
#>  exp_width_of_unstd_conditions_diff 0.3  
#> 
#> Confidence level: 95%
```
