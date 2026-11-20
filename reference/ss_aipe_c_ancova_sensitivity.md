# Sensitivity Analysis for Sample Size Planning for the (Unstandardized) Contrast in Randomized ANCOVA From the Accuracy in Parameter Estimation (AIPE) Perspective

Performs a sensitivity analysis when planning sample size from the
Accuracy in Parameter Estimation (AIPE) Perspective for the
(unstandardized) contrast in randomized ANCOVA design.

## Usage

``` r
ss_aipe_c_ancova_sensitivity(
  true_var_error_ancova = NULL,
  supposed_var_error_ancova = NULL,
  true_var_error_anova = NULL,
  supposed_var_error_anova = NULL,
  true_cor,
  supposed_cor = NULL,
  G = 10000,
  true_mean_Y,
  true_sd_Y,
  true_mean_X,
  true_sd_X,
  c_weights,
  width,
  conf_level = 0.95,
  assurance = NULL,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_var_error_ancova:

  population error variance of the ANCOVA model

- supposed_var_error_ancova:

  the supposed error variance of the ANCOVA model, the planning value
  given to
  [`ss_aipe_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova.md)
  to choose the sample size

- true_var_error_anova:

  population error variance of the ANOVA model (i.e., excluding the
  covariate)

- supposed_var_error_anova:

  the supposed error variance of the ANOVA model (i.e., excluding the
  covariate), a planning value used with `supposed_cor`

- true_cor:

  population correlation coefficient of the response and the covariate

- supposed_cor:

  the supposed correlation coefficient of the response and the
  covariate, a planning value used with `supposed_var_error_anova`

- G:

  number of generations (i.e., replications) of the simulation

- true_mean_Y:

  vector that contains the response's population mean of each group

- true_sd_Y:

  the population standard deviation of the response

- true_mean_X:

  the population mean of the covariate

- true_sd_X:

  the population standard deviation of the covariate

- c_weights:

  the contrast weights

- width:

  the desired full width of the obtained confidence interval

- conf_level:

  the desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- assurance:

  parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be NULL or between zero and unity)

- filename:

  an optional path for a comma separated file recording every
  replication (the realized contrast, its full and covariate-ignoring
  standard errors and their ratio, the interval width, and the three
  non-coverage indicators): nothing is written when `filename` is `NULL`
  (the default), a new file with a header row is created otherwise, an
  existing file at that path is appended to, and a throwaway run should
  point it at `tempfile(fileext = ".csv")`

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. The `term` entries
are: `mean_contrast`, `median_contrast`, `sd_contrast` (summaries of the
realized unstandardized contrast); `mean_ci_width`, `median_ci_width`,
`sd_ci_width` (summaries of the realized interval widths);
`pct_ci_less_width` (proportion of intervals narrower than the planning
target `width`); `pct_ci_miss_low` and `pct_ci_miss_high` (tail-specific
empirical non-coverage of the population contrast); `total_type_I_error`
(overall empirical non-coverage, the sum of the two tails);
`mean_se_ratio` (mean ratio of the contrast standard error that ignores
the covariate-imbalance term to the full ANCOVA standard error); and the
input echoes `n_per_group`, `total_N`, `true_contrast` (the population
contrast implied by `true_mean_Y` and `c_weights`),
`supposed_var_error_ancova` (as supplied or as resolved from
`supposed_var_error_anova` and `supposed_cor`), `true_cor`, `width`,
`conf_level`, and `assurance` (present only when an assurance was
supplied). The proportion rows are on the 0 to 1 scale, not percentages.
The per-replication vectors (`contrast`, `se_contrast`,
`se_contrast_restricted`, `width_obs`) are not returned; they are
written to the comma separated file named by `filename` when one is
supplied.

## Details

The arguments `true_mean_Y`, `true_mean_X`, `true_sd_Y`, and `true_sd_X`
are used to generate random data in the simulations for the sensitivity
analysis. The value of `true_sd_Y` should be the same as the square root
of `true_var_error_anova`.

So far this function is based on one-covariate randomized ANCOVA design
only. The argument `true_mean_X` should be a single number, because it
is assumed that the population mean of the covariate is equal across
groups in randomized ANCOVA.

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

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Monte Carlo sensitivity sweep; G is small here so the example runs quickly.
# Raise G (e.g., G = 1000 or more) for a stable sensitivity analysis.
set.seed(113)
ss_aipe_c_ancova_sensitivity(true_var_error_ancova = 30,
                             supposed_var_error_ancova = 30,
                             true_cor = .2, true_mean_Y = c(10, 12, 15, 13),
                             true_mean_X = 2, true_sd_X = 1.3, true_sd_Y = 2,
                             c_weights = c(1, 0, -1, 0), width = 3, G = 25)
#>  term                      value 
#>  mean_contrast             -4.985
#>  median_contrast           -4.987
#>  sd_contrast               0.287 
#>  mean_ci_width             1.069 
#>  median_ci_width           1.061 
#>  sd_ci_width               0.038 
#>  pct_ci_less_width         1.000 
#>  pct_ci_miss_low           0.040 
#>  pct_ci_miss_high          0.040 
#>  total_type_I_error        0.080 
#>  mean_se_ratio             0.999 
#>  n_per_group               104   
#>  total_N                   416   
#>  true_contrast             -5    
#>  supposed_var_error_ancova 30    
#>  true_cor                  0.200 
#>  width                     3     
#>  conf_level                0.950 
#> 
#> Replications: 25
#> 
#> Confidence level: 95%

ss_aipe_c_ancova_sensitivity(true_var_error_anova = 36,
                             supposed_var_error_anova = 36,
                             true_cor = .2, supposed_cor = .2,
                             true_mean_Y = c(10, 12, 15, 13), true_mean_X = 2,
                             true_sd_X = 1.3, true_sd_Y = 6,
                             c_weights = c(1, 0, -1, 0), width = 3,
                             assurance = NULL, G = 25)
#>  term                      value 
#>  mean_contrast             -4.933
#>  median_contrast           -4.936
#>  sd_contrast               0.797 
#>  mean_ci_width             2.989 
#>  median_ci_width           2.991 
#>  sd_ci_width               0.065 
#>  pct_ci_less_width         0.600 
#>  pct_ci_miss_low           0.040 
#>  pct_ci_miss_high          0.000 
#>  total_type_I_error        0.040 
#>  mean_se_ratio             0.999 
#>  n_per_group               119   
#>  total_N                   476   
#>  true_contrast             -5    
#>  supposed_var_error_ancova 34.560
#>  true_cor                  0.200 
#>  width                     3     
#>  conf_level                0.950 
#> 
#> Replications: 25
#> 
#> Confidence level: 95%
```
