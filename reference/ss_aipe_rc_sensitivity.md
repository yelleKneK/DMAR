# Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the Unstandardized Regression Coefficient

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation Perspective for the unstandardized
regression coefficient. It is
[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
under the name that pairs with
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md);
[`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md)
is the counterpart for the standardized coefficient.

## Usage

``` r
ss_aipe_rc_sensitivity(
  true_var_Y = NULL,
  true_cov_YX = NULL,
  true_cov_XX = NULL,
  supposed_var_Y = NULL,
  supposed_cov_YX = NULL,
  supposed_cov_XX = NULL,
  specified_N = NULL,
  which_predictor = 1,
  width = NULL,
  noncentral = FALSE,
  conf_level = 0.95,
  assurance = NULL,
  G = 1000,
  print_iter = TRUE,
  filename = NULL,
  assurance_method = c("analytic", "simulation"),
  seed = NULL
)
```

## Arguments

- true_var_Y:

  Population variance of the dependent variable (*Y*)

- true_cov_YX:

  Population covariances vector between the *p* predictor variables and
  the dependent variable (*Y*)

- true_cov_XX:

  Population covariance matrix of the *p* predictor variables

- supposed_var_Y:

  Supposed variance of the dependent variable (*Y*), the planning value
  from which the sample size is planned

- supposed_cov_YX:

  Supposed covariances vector (planning values) between the *p*
  predictor variables and the dependent variable (*Y*)

- supposed_cov_XX:

  Supposed covariance matrix (planning values) of the *p* predictor
  variables

- specified_N:

  Directly specified sample size (instead of planning one from the
  supposed covariance structure)

- which_predictor:

  identifies which of the *p* predictors is of interest

- width:

  Desired confidence interval width for the regression coefficient of
  interest, in the units of the data

- noncentral:

  specify with a `TRUE` or `FALSE` statement whether or not the
  noncentral approach to sample size planning should be used; the
  noncentral plan is for a standardized coefficient, so a plan from the
  supposed values with `noncentral = TRUE` needs supposed variances of 1
  for *Y* and the targeted predictor and is refused otherwise (see
  [`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md))

- conf_level:

  desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate)

- assurance:

  degree of certainty that the obtained confidence interval will be
  sufficiently narrow (i.e., the probability that the observed interval
  will be no larger than desired)

- G:

  the number of generations (i.e., replications) of the simulation
  within the function

- print_iter:

  specify with a `TRUE/FALSE` statement if the iteration number should
  be printed as the simulation within the function runs

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the unstandardized coefficient estimate `reg_coef_j`, its confidence
  limits `ci_lower` and `ci_upper`, the observed \\R^2\\ `R2`, the
  standard error `se_reg_coef_j`, and the *t* statistic `t_value`),
  appended when the file already exists and created otherwise; the
  default `NULL` writes nothing, and a throwaway run that wants the file
  should point it at `tempfile(fileext = ".csv")`.

- assurance_method:

  Passed to
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
  which plans the sample size from the supposed values with `assurance`:
  `"analytic"` (the default) for the smallest sample size at which the
  exact probability of an interval no wider than `width` is at least
  `assurance` with the central interval, and for the published
  approximation of MBESS, which does not compute that probability, with
  `noncentral = TRUE`; `"simulation"`, offered with `noncentral = TRUE`
  only, for the smallest sample size at which a simulated estimate of
  the probability (made by the planner from its default `G` of 20000
  samples with the random number stream that `seed` sets) is at least
  `assurance`. With `noncentral = FALSE` `"simulation"` is refused with
  a condition of class `dmar_argument_invalid`. Not used with
  `specified_N` or without `assurance`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis, the table of
[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md):
mean / median / SD summaries of the realized unstandardized regression
coefficient (`mean_reg_coef_j`, `median_reg_coef_j`, `sd_reg_coef_j`),
the realized interval widths, and the realized squared multiple
correlation coefficient; the proportion of intervals at or below the
planning target (`pct_ci_less_width`); the tail-specific and overall
empirical non-coverage rates (`pct_ci_miss_low`, `pct_ci_miss_high`,
`total_type_I_error`), all proportions on the 0 to 1 scale; and the
input echoes (`total_N`, `p`, `which_predictor`, `true_reg_coef_j`,
`supposed_reg_coef_j`, `width`, `conf_level`, and, when one was
supplied, `assurance`). See
[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md)
for the full row list and the scale of each value.

## Details

Direct specification of `true_cov_YX` and `true_cov_XX` is necessary,
even if one is interested in a single regression coefficient, so that
the covariance/correlation structure can be specified when the
simulation within the function runs.

## Note

Note that when the true and supposed covariance structures agree
(`true_cov_YX` equals `supposed_cov_YX` and `true_cov_XX` equals
`supposed_cov_XX`), the results are not literally from a sensitivity
analysis, rather the function performs a standard simulation study. A
simulation study can be helpful in order to determine if the sample size
procedure under or overestimates necessary sample size. See
`ss_aipe_reg_coef_sensitivity` in DMAR for more details.

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

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

## See also

[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md),
[`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md),
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Sensitivity analysis for an unstandardized regression coefficient
# with two correlated predictors. G is kept small here so the example
# runs quickly; raise G for a stable Monte Carlo summary.
set.seed(113)
cov_XX <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
cor_YX <- c(0.4, 0.3)
cov_YX <- cor_YX
ss_aipe_rc_sensitivity(
  true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = cov_XX,
  supposed_var_Y = 1, supposed_cov_YX = cov_YX, supposed_cov_XX = cov_XX,
  which_predictor = 1, width = 0.20, conf_level = 0.95,
  G = 20, print_iter = FALSE
)
#>  term                value
#>  mean_reg_coef_j     0.351
#>  median_reg_coef_j   0.349
#>  sd_reg_coef_j       0.043
#>  mean_ci_width       0.199
#>  median_ci_width     0.200
#>  sd_ci_width         0.009
#>  pct_ci_less_width   0.500
#>  pct_ci_miss_low     0.000
#>  pct_ci_miss_high    0.000
#>  total_type_I_error  0.000
#>  mean_R2             0.211
#>  median_R2           0.213
#>  sd_R2               0.040
#>  total_N             345  
#>  p                   2    
#>  which_predictor     1    
#>  true_reg_coef_j     0.341
#>  supposed_reg_coef_j 0.341
#>  width               0.200
#>  conf_level          0.950
#> 
#> Replications: 20
#> 
#> Confidence level: 95%
```
