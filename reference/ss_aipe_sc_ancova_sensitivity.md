# Sensitivity Analysis for the Sample Size Planning Method for Standardized ANCOVA Contrast

Sensitivity analysis for the sample size planning method with the goal
to obtain sufficiently narrow confidence intervals for standardized
ANCOVA complex contrasts.

## Usage

``` r
ss_aipe_sc_ancova_sensitivity(
  true_std_contrast = NULL,
  supposed_std_contrast = NULL,
  c_weights,
  width = NULL,
  specified_n_per_group = NULL,
  true_mean_X = 0,
  true_sd_X = 1,
  true_cor,
  divisor = "s_ancova",
  assurance = NULL,
  conf_level = 0.95,
  G = 10000,
  print_iter = TRUE,
  filename = NULL,
  assurance_method = c("analytic", "simulation"),
  seed = NULL
)
```

## Arguments

- true_std_contrast:

  the population standardized ANCOVA contrast

- supposed_std_contrast:

  the supposed standardized ANCOVA contrast, the planning value given to
  [`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md)
  to choose the sample size

- c_weights:

  the contrast weights

- width:

  the desired full width of the obtained confidence interval. Required
  when `supposed_std_contrast` is given, because the sample size is
  planned from it; a call without it is refused with a condition of
  class `dmar_argument_invalid`.

- specified_n_per_group:

  selected sample size per group to use in order to determine
  distributional properties of a given value of sample size

- true_mean_X:

  the population mean for the covariate

- true_sd_X:

  the population standard deviation of the covariate

- true_cor:

  the population correlation coefficient between the response and the
  covariate

- divisor:

  which error standard deviation to be used in standardizing the
  contrast; the value can be either `"s_ancova"` or `"s_anova"`

- assurance:

  parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty.
  When the sample size is planned from `supposed_std_contrast`, it is
  passed to
  [`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md)
  and must be `NULL` or greater than .5 and less than 1; with
  `specified_n_per_group` it is only reported in the result.

- conf_level:

  the desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- G:

  number of generations (i.e., replications) of the simulation

- print_iter:

  to print the current value of the iterations

- filename:

  Optional path of a CSV file to receive the per-replication results
  (the observed standardized contrast, `std_contrast`, the full and
  one-sided interval widths, the tail and overall misses, and the
  confidence limits), appended when the file already exists and created
  otherwise; the default `NULL` writes nothing, and a throwaway run that
  wants the file should point it at `tempfile(fileext = ".csv")`.

- assurance_method:

  Passed to
  [`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md),
  which plans the sample size from `supposed_std_contrast` with
  `assurance`: `"analytic"` (the default) for the published plan of Lai
  and Kelley (2012), which does not compute the probability it attains
  and can attain less than `assurance`, and `"simulation"` for the
  smallest sample size per group at which the probability of an interval
  no wider than `width`, estimated by the planner from its default `G`
  of 10000 simulated samples with the random number stream that `seed`
  sets, is at least `assurance`. Not used with `specified_n_per_group`
  or without `assurance`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. The `term` entries
are: `mean_std_contrast`, `median_std_contrast`, `sd_std_contrast`
(summaries of the realized standardized ANCOVA contrast);
`mean_ci_width`, `median_ci_width`, `sd_ci_width` (summaries of the full
interval widths); `mean_ci_width_lower` and `mean_ci_width_upper` (mean
one-sided widths, measured from the observed contrast to each limit);
`pct_ci_less_width` (proportion of intervals at or below the target
width); `pct_ci_miss_low` and `pct_ci_miss_high` (tail-specific
empirical non-coverage of `true_std_contrast`); `total_type_I_error`
(overall empirical non-coverage, the sum of the two tails); and the
input echoes `n_per_group` (the size simulated at, planned or
specified), `total_N`, `true_std_contrast`, `supposed_std_contrast` (NA
when `specified_n_per_group` was supplied instead), `true_cor`, `width`,
`conf_level`, and `assurance` (present only when an assurance was
supplied). The proportion and Type I error rows are proportions on the 0
to 1 scale, not percentages.

## Details

The sample size planning method this function is based on is developed
in the context of simple (i.e., one-response-one-covariate) ANCOVA model
and randomized design (i.e., same population covariate mean across
groups).

An ANCOVA contrast can be standardized in at least two ways: (a) divided
by the error standard deviation of the ANOVA model, (b) divided by the
error standard deviation of the ANCOVA model. This function can be used
to analyze both types of standardized ANCOVA contrasts.

The population mean and standard deviation of the covariate does not
affect the sample size planning procedure; they can be specified as any
values that are considered as reasonable by the user.

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

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical models. In L. L. Harlow, S.
A. Mulaik, & J. H. Steiger (Eds.), *What if there were no significance
tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ss_aipe_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova.md),
[`ss_aipe_sc_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Sensitivity analysis for a standardized ANCOVA contrast across
# three groups, contrast (-1, 0, 1), a covariate-outcome correlation
# of 0.4, and a planning target width of 0.5. G = 50 keeps the
# example quick; a reported sensitivity analysis deserves the default of
# G = 10000 replications.
set.seed(113)
ss_aipe_sc_ancova_sensitivity(
  true_std_contrast = 0.5, supposed_std_contrast = 0.5,
  c_weights = c(-1, 0, 1),
  width = 0.5, true_cor = 0.4,
  conf_level = 0.95, G = 50, print_iter = FALSE
)
#>  term                  value
#>  mean_std_contrast     0.481
#>  median_std_contrast   0.467
#>  sd_std_contrast       0.119
#>  mean_ci_width         0.500
#>  median_ci_width       0.499
#>  sd_ci_width           0.002
#>  mean_ci_width_lower   0.250
#>  mean_ci_width_upper   0.249
#>  pct_ci_less_width     0.580
#>  pct_ci_miss_low       0.000
#>  pct_ci_miss_high      0.020
#>  total_type_I_error    0.020
#>  n_per_group           126  
#>  total_N               378  
#>  true_std_contrast     0.500
#>  supposed_std_contrast 0.500
#>  true_cor              0.400
#>  width                 0.500
#>  conf_level            0.950
#> 
#> Replications: 50
#> 
#> Confidence level: 95%
```
