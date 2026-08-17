# Sensitivity Analysis for the Sample Size Planning Method for Standardized ANCOVA Contrast

Sensitivity analysis for the sample size planning method with the goal
to obtain sufficiently narrow confidence intervals for standardized
ANCOVA complex contrasts.

## Usage

``` r
ss_aipe_sc_ancova_sensitivity(
  true_psi = NULL,
  estimated_psi = NULL,
  c_weights,
  desired_width = NULL,
  n_per_group = NULL,
  mu_x = 0,
  sigma_x = 1,
  rho,
  divisor = "s_ancova",
  assurance = NULL,
  conf_level = 0.95,
  G = 10000,
  print_iter = TRUE,
  save = FALSE,
  filename = "ss_aipe_sc_ancova_sensitivity_result.csv",
  ...
)
```

## Arguments

- true_psi:

  the population standardized ANCOVA contrast

- estimated_psi:

  the estimated standardized ANCOVA contrast

- c_weights:

  the contrast weights

- desired_width:

  the desired full width of the obtained confidence interval

- n_per_group:

  selected sample size to use in order to determine distributional
  properties of a given value of sample size

- mu_x:

  the population mean for the covariate

- sigma_x:

  the population standard deviation of the covariate

- rho:

  the population correlation coefficient between the response and the
  covariate

- divisor:

  which error standard deviation to be used in standardizing the
  contrast; the value can be either `"s_ancova"` or `"s_anova"`

- assurance:

  parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be `NULL` or between zero and unity)

- conf_level:

  the desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- G:

  number of generations (i.e., replications) of the simulation

- print_iter:

  to print the current value of the iterations

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

- ...:

  allows one to potentially include parameter values for inner functions

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo sensitivity analysis across `G` replications. The `term` entries
are: `mean_psi`, `median_psi`, `sd_psi` (summaries of the realized
standardized ANCOVA contrast); `mean_ci_width`, `median_ci_width`,
`sd_ci_width` (summaries of the full interval widths);
`mean_ci_width_lower` and `mean_ci_width_upper` (mean one-sided widths,
measured from the observed contrast to each limit); `pct_ci_less_w`
(proportion of intervals at or below the target width);
`pct_ci_miss_low` and `pct_ci_miss_high` (tail-specific empirical
non-coverage of `true_psi`); `total_type_I_error` (overall empirical
non-coverage, the sum of the two tails); and the input echoes
`n_per_group`, `total_N`, `true_psi`, `estimated_psi` (NA when
`n_per_group` was supplied instead), `rho`, `width`, `conf_level`, and
`assurance` (present only when an assurance was supplied). The
proportion and Type I error rows are proportions on the 0 to 1 scale,
not percentages.

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
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

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
# of 0.4, and a planning target width of 0.5. Sizes are kept small
# here so the Monte Carlo sweep runs quickly; raise G for a stable
# estimate in practice.
set.seed(113)
ss_aipe_sc_ancova_sensitivity(
  true_psi = 0.5, estimated_psi = 0.5,
  c_weights = c(-1, 0, 1),
  desired_width = 0.5, rho = 0.4,
  conf_level = 0.95, G = 50, print_iter = FALSE
)
#>  term                value  
#>  mean_psi            0.481  
#>  median_psi          0.467  
#>  sd_psi              0.119  
#>  mean_ci_width       0.5    
#>  median_ci_width     0.499  
#>  sd_ci_width         0.00249
#>  mean_ci_width_lower 0.25   
#>  mean_ci_width_upper 0.249  
#>  pct_ci_less_w       0.58   
#>  pct_ci_miss_low     0      
#>  pct_ci_miss_high    0.02   
#>  total_type_I_error  0.02   
#>  n_per_group         126    
#>  total_N             378    
#>  true_psi            0.5    
#>  estimated_psi       0.5    
#>  rho                 0.4    
#>  width               0.5    
#>  conf_level          0.95   
#> 
#> Confidence level: 95%
```
