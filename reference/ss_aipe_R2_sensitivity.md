# Sensitivity Analysis for Sample Size Planning With the Goal of Accuracy in Parameter Estimation (I.e., a Narrow Observed Confidence Interval)

Given `supposed_R2` and `true_R2`, one can perform a sensitivity
analysis to determine the effect of a misspecified population squared
multiple correlation coefficient using the Accuracy in Parameter
Estimation (AIPE) approach to sample size planning. The function
evaluates the effect of a misspecified `true_R2` on the width of
obtained confidence intervals.

## Usage

``` r
ss_aipe_R2_sensitivity(
  true_R2 = NULL,
  supposed_R2 = NULL,
  width = NULL,
  p = NULL,
  random_predictors = TRUE,
  specified_N = NULL,
  assurance = NULL,
  conf_level = 0.95,
  generate_random_predictors = TRUE,
  true_cor_YX = 0.3,
  true_cor_XX = 0.3,
  G = 10000,
  print_iter = TRUE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_R2:

  Value of the population squared multiple correlation coefficient

- supposed_R2:

  Supposed value of the squared multiple correlation coefficient, the
  planning value from which the sample size is planned

- width:

  Full confidence interval width of interest

- p:

  Number of predictors

- random_predictors:

  Whether or not the sample size procedure and the simulation itself
  should be based on random (set to `TRUE`) or fixed predictors (set to
  `FALSE`)

- specified_N:

  Selected sample size to use in order to determine distributional
  properties at a given value of sample size. Give `supposed_R2` or
  `specified_N`, not both; a call that gives both is refused with a
  condition of class `dmar_argument_invalid`.

- assurance:

  Parameter to ensure confidence interval width with a specified degree
  of certainty

- conf_level:

  Confidence interval coverage (symmetric coverage)

- generate_random_predictors:

  Specify whether the simulation should be based on random (default) or
  fixed regressors.

- true_cor_YX:

  Value of the correlation between *y* (dependent variable) and each of
  the *x* variables (independent variables)

- true_cor_XX:

  Value of the correlation among the *x* variables (independent
  variables)

- G:

  Number of generations (i.e., replications) of the simulation

- print_iter:

  Should the iteration number (between 1 and `G`) during the run of the
  function

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the confidence limits `ci_lower` and `ci_upper`, the observed \\R^2\\
  `R2`, the one-sided widths `ci_width_lower` and `ci_width_upper`,
  measured from the observed \\R^2\\ to each limit, and the full width
  `ci_width`, the quantities the `*_ci_width_lower`, `*_ci_width_upper`,
  and `*_ci_width` rows summarize), overwriting any file already at that
  path; the default `NULL` writes nothing, and a throwaway run that
  wants the file should point it at `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the run reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. The `term` entries
are: `mean_ci_lower`, `median_ci_lower`, `sd_ci_lower`, `mean_ci_upper`,
`median_ci_upper`, `sd_ci_upper` (summaries of the realized confidence
limits); `mean_R2`, `median_R2`, `sd_R2` (summaries of the observed
\\R^2\\); `mean_ci_width_lower`, `median_ci_width_lower`,
`sd_ci_width_lower`, `mean_ci_width_upper`, `median_ci_width_upper`,
`sd_ci_width_upper` (summaries of the one-sided widths, measured from
the observed \\R^2\\ to each limit); `mean_ci_width`, `median_ci_width`,
`sd_ci_width` (summaries of the full interval widths);
`pct_ci_less_width` (proportion of intervals with width at or below the
planning target `width`); `pct_ci_miss_low` and `pct_ci_miss_high`
(tail-specific empirical non-coverage of `true_R2`);
`total_type_I_error` (overall empirical non-coverage, the sum of the two
tails); and the input echoes `total_N` (the sample size evaluated), `p`,
`true_R2`, `supposed_R2` (NA when `specified_N` was supplied instead),
`width`, `conf_level`, and `assurance` (present only when an assurance
was supplied). The proportion rows are on the 0 to 1 scale, not
percentages. The number of replications on which a confidence interval
could not be obtained is `attr(x, "simulation")$failed` (see the section
on the record of the simulation).

## Details

When `supposed_R2`=`true_R2`, the results are that of a simulation study
when all assumptions are satisfied. Rather than specifying
`supposed_R2`, one can specify `specified_N` to determine the results of
a particular sample size (when doing this `supposed_R2` cannot be
specified).

The sample size estimation procedure technically assumes multivariate
normal variables (`p`+1) with fixed predictors (`x`/independent
variables), yet the function assumes random multivariate normal
predictors (having a `p`+1 multivariate distribution). As Gatsonis and
Sampson (1989) note in the context of statistical power analysis (recall
this function is used in the context of precision), there is little
difference in the outcome.

In the behavioral, educational, and social sciences, predictor variables
are almost always random, and thus `random_predictors` should generally
be used. `random_predictors=TRUE` specifies how both the sample size
planning procedure and the confidence intervals are calculated based on
the random predictors/regressors. The internal simulation generates
random or fixed predictors/regressors based on whether variables
predictor variables are random or fixed. However, when
`random_predictors=FALSE`, only the sample size planning procedure and
the confidence intervals are calculated based on the parameter. The
parameter `generate_random_predictors` (where the default is `TRUE` so
that random predictors/regressors are generated) allows random or fixed
predictor variables to be generated. Because the sample size planning
procedure and the internal simulation are both specified, for purposes
of sensitivity analysis random/fixed can be crossed to examine the
effects of specifying sample size based on one but using it on data
based on the other.

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

Algina, J. & Olejnik, S. (2000). Determining sample size for accurate
estimation of the squared multiple correlation coefficient.
*Multivariate Behavioral Research, 35*, 119–137.
[doi:10.1207/s15327906mbr3501_5](https://doi.org/10.1207/s15327906mbr3501_5)

Gatsonis, C. & Sampson, A. R. (1989). Multiple Correlation: Exact power
and sample size calculations. *Psychological Bulletin, 106*(3), 516–524.

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals, *Multivariate Behavioral Research, 43*(4),
524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Steiger, J. H., & Fouladi, R. T. (1992). R2: A computer program for
interval estimation, power calculations, sample size estimation, and
hypothesis testing in multiple regression. *Behavior Research Methods,
Instruments, & Computers, 24*(4), 581–582.
[doi:10.3758/BF03203611](https://doi.org/10.3758/BF03203611)

## See also

[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The planner used .4 for a population squared multiple correlation
# coefficient that is really .5. G = 25 keeps the example quick; a
# reported sensitivity analysis deserves the default G of 10000.
ss_aipe_R2_sensitivity(true_R2 = .5, supposed_R2 = .4, width = .10, p = 5,
                       conf_level = 0.95, G = 25, print_iter = FALSE,
                       seed = 113)
#>  term                  value
#>  mean_ci_lower         0.455
#>  median_ci_lower       0.452
#>  sd_ci_lower           0.020
#>  mean_ci_upper         0.548
#>  median_ci_upper       0.545
#>  sd_ci_upper           0.019
#>  mean_R2               0.505
#>  median_R2             0.502
#>  sd_R2                 0.019
#>  mean_ci_width_lower   0.050
#>  median_ci_width_lower 0.050
#>  sd_ci_width_lower     0.001
#>  mean_ci_width_upper   0.043
#>  median_ci_width_upper 0.043
#>  sd_ci_width_upper     0.001
#>  mean_ci_width         0.093
#>  median_ci_width       0.093
#>  sd_ci_width           0.002
#>  pct_ci_less_width     1.000
#>  pct_ci_miss_low       0.000
#>  pct_ci_miss_high      0.000
#>  total_type_I_error    0.000
#>  total_N               889  
#>  p                     5    
#>  true_R2               0.500
#>  supposed_R2           0.400
#>  width                 0.100
#>  conf_level            0.950
#> 
#> Replications: 25
#> 
#> Confidence level: 95%
```
