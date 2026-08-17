# Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the Coefficient of Variation

Quantifies how much misspecification of the population coefficient of
variation can distort an AIPE-based sample size plan. Given a true
(population) value of the coefficient of variation and the value that
was used in planning, the function simulates draws of size *N* from a
normal population, computes the confidence interval for the coefficient
of variation on each replication, and summarizes how often the realized
interval width is below the desired target and how often the interval
covers the population value. This is the standard sensitivity-analysis
workflow described in Kelley (2007) and Maxwell, Delaney, and Kelley
(2027, Section 3.11 on sample size planning).

## Usage

``` r
ss_aipe_cv_sensitivity(
  true_cv = NULL,
  estimated_cv = NULL,
  width = NULL,
  assurance = NULL,
  mean = 100,
  specified_N = NULL,
  conf_level = 0.95,
  G = 1000,
  print_iter = FALSE,
  save = FALSE,
  filename = "ss_aipe_cv_sensitivity_result.csv"
)
```

## Arguments

- true_cv:

  Population coefficient of variation (the data generating value)

- estimated_cv:

  Coefficient of variation used to plan the study (the value the
  researcher guessed when invoking
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md));
  must be positive. Supply this or `specified_N` but not both.

- width:

  Desired (full) width of the two-sided confidence interval for the
  population coefficient of variation

- assurance:

  Probability with which the realized interval should be no wider than
  `width` (must be `NULL` or strictly between 0 and 1; `NULL` means plan
  to the *expected* width without an assurance constraint)

- mean:

  Population mean used by the simulator to generate data (the standard
  deviation is determined by `mean` and `true_cv`, since CV = sigma/mu).
  Default 100.

- specified_N:

  Pre-specified sample size to evaluate (use this when you want the
  sensitivity results at a fixed *N* rather than at the *N* that
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
  would recommend); incompatible with `estimated_cv`.

- conf_level:

  Desired confidence level (i.e., 1 - Type I error rate); default 0.95

- G:

  Number of Monte Carlo replications; defaults to 1000. Increase (e.g.,
  5000 or 10000) for stable Type I error estimates.

- print_iter:

  Logical. If `TRUE` the simulation prints the iteration index after
  each replication (helpful for long runs); default `FALSE`.

- save:

  Logical. If `TRUE` the per-replication results are appended to a CSV
  file at `filename`; default `FALSE`.

- filename:

  Path used when `save = TRUE`; default
  `"ss_aipe_cv_sensitivity_result.csv"` in the current working
  directory.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo results across the `G` replications. The `term` entries are:
`"mean_cv"`, `"median_cv"`, `"sd_cv"` (mean / median / SD of the *G*
observed sample coefficients of variation); `"mean_ci_width"`,
`"median_ci_width"`, `"sd_ci_width"` (corresponding summaries of the
realized interval widths); `"pct_ci_less_w"` (proportion of intervals at
or below the planning width `width`); `"pct_ci_miss_low"` and
`"pct_ci_miss_high"` (tail-specific non-coverage);
`"total_type_I_error"` (overall empirical non-coverage of `true_cv`);
plus the input echoes `"total_N"` (the sample size evaluated),
`"true_cv"`, `"estimated_cv"` (NA when `specified_N` was supplied
instead), `"width"`, `"conf_level"`, and `"assurance"` (present only
when an assurance was supplied). The proportion rows are on the 0 to 1
scale, not percentages, so `total_type_I_error` is the sum of
`pct_ci_miss_low` and `pct_ci_miss_high`.

## Details

Sample size planning for the coefficient of variation under the Accuracy
in Parameter Estimation framework chooses *N* so that the expected (or,
with assurance, the high-probability) confidence interval width is no
larger than `width` (Kelley, 2007). Because the procedure assumes the
planning value `estimated_cv` matches the population value `true_cv`, in
practice the realized width will deviate from the planned width whenever
the planning value is wrong. This sensitivity analysis quantifies the
deviation by Monte Carlo simulation: the planned *N* is obtained from
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
with `estimated_cv`, then samples are drawn from the *true* population
(with coefficient of variation `true_cv`) and the realized confidence
interval widths are summarized.

For a discussion of AIPE-based sample size planning more generally and
how sensitivity analyses guard against misspecification, see Maxwell,
Delaney, & Kelley (2027, Section 3.5).

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3.)

## See also

[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md),
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Textbook scenario (Kelley, 2007). A researcher plans a study to estimate the
# coefficient of variation of reaction times with a 95% confidence interval no
# wider than .10. They guess from prior work that the population coefficient of
# variation is around .25 and apply ss_aipe_cv() to obtain a planned N.

# Question 1: how does the realized interval width behave when the planning
# value is correct (well-specified case)?
set.seed(113)
ss_aipe_cv_sensitivity(
  true_cv      = .25,
  estimated_cv = .25,
  width            = .10,
  assurance        = NULL,
  conf_level       = .95,
  G                = 200,
  print_iter       = FALSE
)
#> Warning: The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.
#> Warning: The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.
#> Warning: The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.
#> Warning: The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.
#> Warning: The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.
#> Warning: The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.
#> Warning: The observed noncentrality parameter exceeds 37.62 in magnitude, which is the limit at which R's pt()/qt() can return accurate noncentral t probabilities. Results may be inaccurate; use caution.
#>  term               value 
#>  mean_cv            0.247 
#>  median_cv          0.247 
#>  sd_cv              0.0239
#>  mean_ci_width      0.0976
#>  median_ci_width    0.0979
#>  sd_ci_width        0.0108
#>  pct_ci_less_w      0.59  
#>  pct_ci_miss_low    0.015 
#>  pct_ci_miss_high   0.03  
#>  total_type_I_error 0.045 
#>  total_N            60    
#>  true_cv            0.25  
#>  estimated_cv       0.25  
#>  width              0.1   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%

# Question 2: what happens if the planning value is materially smaller than
# the true coefficient of variation (a common direction of misspecification,
# since planning values are often optimistic)? The intervals will be wider on
# average than the target and the pct_ci_less_w will fall.
set.seed(113)
ss_aipe_cv_sensitivity(
  true_cv      = .35,
  estimated_cv = .25,
  width            = .10,
  assurance        = NULL,
  conf_level       = .95,
  G                = 200,
  print_iter       = FALSE
)
#>  term               value 
#>  mean_cv            0.345 
#>  median_cv          0.346 
#>  sd_cv              0.0347
#>  mean_ci_width      0.144 
#>  median_ci_width    0.145 
#>  sd_ci_width        0.0176
#>  pct_ci_less_w      0.005 
#>  pct_ci_miss_low    0.015 
#>  pct_ci_miss_high   0.03  
#>  total_type_I_error 0.045 
#>  total_N            60    
#>  true_cv            0.35  
#>  estimated_cv       0.25  
#>  width              0.1   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%
```
