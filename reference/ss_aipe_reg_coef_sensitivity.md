# Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the Unstandardized Regression Coefficient

This function performs a sensitivity analysis when planning sample size
from the Accuracy in Parameter Estimation Perspective for the
unstandardized regression coefficient.
[`ss_aipe_rc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc_sensitivity.md)
is the same analysis under the name that pairs with
[`ss_aipe_rc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rc.md),
and
[`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md)
is its counterpart for the standardized regression coefficient.

## Usage

``` r
ss_aipe_reg_coef_sensitivity(
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

  Population covariances vector between the `p` predictor variables and
  the dependent variable (*Y*)

- true_cov_XX:

  Population covariance matrix of the `p` predictor variables

- supposed_var_Y:

  Supposed variance of the dependent variable (*Y*), the planning value
  from which the sample size is planned

- supposed_cov_YX:

  Supposed covariances vector (planning values) between the `p`
  predictor variables and the dependent variable (`Y`)

- supposed_cov_XX:

  Supposed covariance matrix (planning values) of the `p` predictor
  variables

- specified_N:

  Directly specified sample size (instead of planning one from the
  supposed covariance structure)

- which_predictor:

  Identifies which of the *p* predictors is of interest

- width:

  Desired confidence interval width for the regression coefficient of
  interest, in the units of the data

- noncentral:

  Specify with a `TRUE` or `FALSE` statement whether or not the
  noncentral approach to sample size planning should be used. The
  noncentral plan is for a standardized coefficient, so when the sample
  size is planned from the supposed values, `noncentral = TRUE` needs
  supposed variances of 1 for *Y* and the targeted predictor, where the
  two scales coincide; otherwise the call is refused (with a condition
  of class `dmar_argument_invalid`), and
  [`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md)
  is the function for the standardized coefficient.

- conf_level:

  Desired level of confidence for the computed interval (i.e., 1 - the
  Type I error rate)

- assurance:

  Degree of certainty that the obtained confidence interval will be
  sufficiently narrow

- G:

  The number of generations (i.e., replications) of the simulation
  within the function

- print_iter:

  Specify with a `TRUE`/`FALSE` statement if the iteration number should
  be printed as the simulation within the function runs

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the coefficient estimate `reg_coef_j`, its confidence limits
  `ci_lower` and `ci_upper`, the observed \\R^2\\ `R2`, the standard
  error `se_reg_coef_j`, and the *t* statistic `t_value`), appended when
  the file already exists and created otherwise; the default `NULL`
  writes nothing, and a throwaway run that wants the file should point
  it at `tempfile(fileext = ".csv")`.

- assurance_method:

  Passed to
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
  when the sample size is planned from the supposed values with
  `assurance`. `"analytic"` (the default) plans, for the central
  interval, the smallest sample size at which the exact probability of
  an interval no wider than `width` is at least `assurance`, and, for
  the noncentral interval, the approximation that MBESS uses, which does
  not compute that probability. `"simulation"`, offered with
  `noncentral = TRUE` only, plans the smallest sample size at which a
  simulated estimate of the probability (made by the planner from its
  default `G` of 20000 samples with the random number stream that `seed`
  sets) is at least `assurance`; with `noncentral = FALSE` it is refused
  with a condition of class `dmar_argument_invalid`, as
  [`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
  refuses it. Not used with `specified_N` or without `assurance`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. The `term` entries
are: `mean_reg_coef_j`, `median_reg_coef_j`, `sd_reg_coef_j` (summaries
of the realized regression coefficient point estimates);
`mean_ci_width`, `median_ci_width`, `sd_ci_width` (summaries of the
realized interval widths); `pct_ci_less_width` (proportion of intervals
at or below the planning target `width`); `pct_ci_miss_low` and
`pct_ci_miss_high` (tail-specific empirical non-coverage of the
population coefficient); `total_type_I_error` (overall empirical
non-coverage, the sum of the two tails); `mean_R2`, `median_R2`, `sd_R2`
(summaries of the realized squared multiple correlation coefficient);
and the input echoes `total_N` (the sample size evaluated), `p`,
`which_predictor`, `true_reg_coef_j` and `supposed_reg_coef_j` (the
population and planning values of the targeted coefficient implied by
the supplied covariance structures), `width`, `conf_level`, and
`assurance` (present only when an assurance was supplied). The
proportion and Type I error rows are proportions on the 0 to 1 scale,
not percentages, so `total_type_I_error` is the sum of `pct_ci_miss_low`
and `pct_ci_miss_high`.

## Details

Direct specification of `true_cov_YX` and `true_cov_XX` is necessary,
even if one is interested in a single regression coefficient, so that
the covariance/correlation structure can be specified when the
simulation within the function runs.

**The coefficient.** The coefficient is the unstandardized one, the
*j*th element of \\\Sigma\_{XX}^{-1}\sigma\_{XY}\\, where
\\\Sigma\_{XX}\\ is `true_cov_XX` and \\\sigma\_{XY}\\ is `true_cov_YX`;
call it \\b_j\\. It is the true value the intervals are judged against,
and the supposed coefficient echoed beside it is the same element
computed from `supposed_cov_XX` and `supposed_cov_YX`. Every replication
is fit in the units of the data, and a sample size planned from the
supposed values is planned by
[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md)
with the supposed standard deviations of *Y* and of the targeted
predictor (the square roots of `supposed_var_Y` and of the *j*th
diagonal element of `supposed_cov_XX`), so that `width` is a width in
the units of the data. The standardized coefficient, \\\beta_j = b_j
\sigma\_{X_j} / \sigma_Y\\ (Kelley & Maxwell, 2003, p. 316), is the
subject of
[`ss_aipe_src_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_src_sensitivity.md).

The interval of each replication is that of
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)
with the standard error of the fitted regression, central or, with
`noncentral = TRUE`, noncentral.

## Note

Note that when the true and supposed covariance structures agree
(`true_cov_YX` equals `supposed_cov_YX` and `true_cov_XX` equals
`supposed_cov_XX`), the results are not literally from a sensitivity
analysis, rather the function performs a standard simulation study. A
simulation study can be helpful in order to determine if the sample size
procedure under or overestimates necessary sample size.

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
# with two predictors at a modest R squared. The Monte Carlo loop is
# run with a small number of generations (G) here so the example is
# fast; use a larger G (for example G = 1000) in real applications.
set.seed(113)
cov_XX <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
cov_YX <- c(0.4, 0.3)
ss_aipe_reg_coef_sensitivity(
  true_var_Y = 1, true_cov_YX = cov_YX, true_cov_XX = cov_XX,
  supposed_var_Y = 1, supposed_cov_YX = cov_YX, supposed_cov_XX = cov_XX,
  which_predictor = 1, width = 0.20, conf_level = 0.95,
  G = 100, print_iter = FALSE
)
#>  term                value
#>  mean_reg_coef_j     0.348
#>  median_reg_coef_j   0.349
#>  sd_reg_coef_j       0.049
#>  mean_ci_width       0.199
#>  median_ci_width     0.200
#>  sd_ci_width         0.010
#>  pct_ci_less_width   0.530
#>  pct_ci_miss_low     0.030
#>  pct_ci_miss_high    0.010
#>  total_type_I_error  0.040
#>  mean_R2             0.204
#>  median_R2           0.197
#>  sd_R2               0.042
#>  total_N             345  
#>  p                   2    
#>  which_predictor     1    
#>  true_reg_coef_j     0.341
#>  supposed_reg_coef_j 0.341
#>  width               0.200
#>  conf_level          0.950
#> 
#> Replications: 100
#> 
#> Confidence level: 95%

# The outcome in units whose variance is 4. The unstandardized
# coefficient of the first predictor is about .68 units of Y per unit of
# the predictor, and the 0.55 width is in those units. For the same
# structure ss_aipe_src_sensitivity() reports the standardized
# coefficient, .68 times 1 / 2, about .34, with a width in standardized
# units.
ss_aipe_reg_coef_sensitivity(
  true_var_Y = 4, true_cov_YX = c(0.8, 0.6), true_cov_XX = cov_XX,
  specified_N = 200, which_predictor = 1, width = 0.55,
  G = 100, print_iter = FALSE, seed = 113
)
#>  term                value
#>  mean_reg_coef_j     0.691
#>  median_reg_coef_j   0.690
#>  sd_reg_coef_j       0.131
#>  mean_ci_width       0.525
#>  median_ci_width     0.522
#>  sd_ci_width         0.041
#>  pct_ci_less_width   0.730
#>  pct_ci_miss_low     0.040
#>  pct_ci_miss_high    0.030
#>  total_type_I_error  0.070
#>  mean_R2             0.203
#>  median_R2           0.204
#>  sd_R2               0.052
#>  total_N             200  
#>  p                   2    
#>  which_predictor     1    
#>  true_reg_coef_j     0.681
#>  supposed_reg_coef_j 0.681
#>  width               0.550
#>  conf_level          0.950
#> 
#> Replications: 100
#> 
#> Confidence level: 95%
ss_aipe_src_sensitivity(
  true_var_Y = 4, true_cov_YX = c(0.8, 0.6), true_cov_XX = cov_XX,
  specified_N = 200, which_predictor = 1, width = 0.28,
  G = 100, print_iter = FALSE, seed = 113
)
#>  term                    value
#>  mean_std_reg_coef_j     0.346
#>  median_std_reg_coef_j   0.348
#>  sd_std_reg_coef_j       0.065
#>  mean_ci_width           0.270
#>  median_ci_width         0.271
#>  sd_ci_width             0.007
#>  pct_ci_less_width       0.920
#>  pct_ci_miss_low         0.040
#>  pct_ci_miss_high        0.020
#>  total_type_I_error      0.060
#>  mean_R2                 0.203
#>  median_R2               0.204
#>  sd_R2                   0.052
#>  total_N                 200  
#>  p                       2    
#>  which_predictor         1    
#>  true_std_reg_coef_j     0.341
#>  supposed_std_reg_coef_j 0.341
#>  width                   0.280
#>  conf_level              0.950
#> 
#> Replications: 100
#> 
#> Confidence level: 95%
```
