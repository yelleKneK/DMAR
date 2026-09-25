# A Priori Monte Carlo Simulation for Sample Size Planning for SEM Targeted Effects

Conduct a priori Monte Carlo simulation to empirically study the effects
of (mis)specifications of input information on the calculated sample
size. Random data are generated from the true covariance matrix but fit
to the proposed model, whereas sample size is calculated based on the
input covariance matrix and proposed model.

## Usage

``` r
ss_aipe_sem_path_sensitivity(
  model,
  est_Sigma,
  true_Sigma = est_Sigma,
  which_path,
  desired_width,
  N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 100,
  filename = NULL,
  ...
)
```

## Arguments

- model:

  A single character string giving the free analysis model in lavaan
  model syntax (see
  [`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)),
  the model that would be fit to the data. The target path must carry a
  parameter label so it can be referred to by name, for example
  `"f2 ~ b*f1"` labels the structural path `b`. The model may or may not
  be the true data generating model

- est_Sigma:

  the covariance matrix used to calculate sample size, may or may not be
  the true covariance matrix. The row names and column names of
  `est_Sigma` should be the same as the observed variables in `model`

- true_Sigma:

  the true population covariance matrix, which will be used to generate
  random data for the simulation study. The row names and column names
  of `true_Sigma` should be the same as the observed variables in
  `model`

- which_path:

  the parameter label of the targeted path, given as a character string,
  for example `"b"` for the path labeled `f2 ~ b*f1` in `model`

- desired_width:

  desired confidence interval width for the model parameter of interest

- N:

  the sample size of random data. If it is `NULL`, it will be determined
  by the sample size planning method

- conf_level:

  confidence level (i.e., 1- Type I error rate)

- assurance:

  the assurance that the confidence interval obtained in a particular
  study will be no wider than desired (must be `NULL` or a value between
  0.50 and 1)

- G:

  number of replications in the Monte Carlo simulation

- filename:

  an optional path for a comma separated file recording every
  replication (the estimate of the targeted path, its standard error,
  the two confidence limits, and the interval width): nothing is written
  when `filename` is `NULL` (the default), a new file with a header row
  is created otherwise, an existing file at that path is appended to,
  and a throwaway run should point it at `tempfile(fileext = ".csv")`

- ...:

  allows one to potentially include parameter values for inner functions

## Value

A `data.frame` with columns `term` and `value` summarizing the a priori
Monte Carlo study. The `term` entries are: `"mean_path"`,
`"median_path"`, `"sd_path"` (summaries of the realized estimates of the
targeted path across the converged replications); `"mean_ci_width"`,
`"median_ci_width"`, `"sd_ci_width"` (summaries of the realized interval
widths); `"pct_ci_less_w"` (proportion of realized widths at or below
`desired_width`); `"pct_ci_miss_low"` and `"pct_ci_miss_high"`
(tail-specific empirical non-coverage of the population path);
`"total_type_I_error"` (overall empirical non-coverage, the sum of the
two tails); and the echoes `"suc_rep"` (number of converged
replications), `"total_N"` (the *N* evaluated), `"true_path"` (the
population value of the targeted path under `true_Sigma`), `"width"`,
`"conf_level"`, and `"assurance"` (present only when an assurance was
supplied). The proportion rows are on the 0 to 1 scale, not percentages.

## Details

This function implements the sample size planning methods proposed in
Lai and Kelley (2011). It calls
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)
to plan the sample size and to identify the targeted path, then fits the
analysis model to data simulated from `true_Sigma` with
[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html). The analysis model is
written in lavaan model syntax with the targeted path given a parameter
label; see
[`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html) for
the syntax and [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) for the
fitting machinery. The population covariance matrices are most naturally
produced by
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md) from
a fully fixed population model. This function requires lavaan and MASS
to be installed.

## Note

Occasionally a replication fails to converge when the analysis model is
fit to a simulated data set. Such replications are not counted toward
the `G` converged replications; the simulation draws fresh data and
continues. A safety cap stops the loop after `20 * G` attempts, and a
single warning is issued if fewer than `G` replications converged.

## References

Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
targeted effects in structural equation modeling: Sample size planning
for narrow confidence intervals. *Psychological Methods, 16*(2),
127–148. [doi:10.1037/a0021764](https://doi.org/10.1037/a0021764)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

## See also

[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md),
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The planning values a researcher would bring to ss_aipe_sem_path(): two
# factors, each measured by three indicators, joined by a structural path
# of 0.5, with every residual variance at 0.5.
planning_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f2 ~ 0.5*f1
  f1 ~~ 1*f1
  f2 ~~ 0.75*f2
  y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
  y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
"

# The population the study will actually sample from: the same structural
# path, but noisier indicators than the planning values assumed, with
# every residual variance at 0.8.
true_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f2 ~ 0.5*f1
  f1 ~~ 1*f1
  f2 ~~ 0.75*f2
  y1 ~~ 0.8*y1; y2 ~~ 0.8*y2; y3 ~~ 0.8*y3
  y4 ~~ 0.8*y4; y5 ~~ 0.8*y5; y6 ~~ 0.8*y6
"

# The free analysis model, with the targeted path labeled "b".
analysis_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f2 ~ b*f1
"

est_Sigma <- cov_sem(planning_model)$sigma_theta
true_Sigma <- cov_sem(true_model)$sigma_theta

# The sample size planned from the optimistic measurement quality is
# evaluated against the population that actually holds. Notice that the
# mean realized width (mean_ci_width) exceeds the desired 0.30, and that
# pct_ci_less_w, the proportion of intervals no wider than desired, falls
# well below the one half that the planned sample size would deliver if
# the planning values were right. G = 20 keeps the example quick; a
# reported sensitivity study deserves the default G = 100 or more.
set.seed(113)
ss_aipe_sem_path_sensitivity(model = analysis_model, est_Sigma = est_Sigma,
                             true_Sigma = true_Sigma, which_path = "b",
                             desired_width = 0.30, G = 20)
#>  term               value 
#>  mean_path          0.488 
#>  median_path        0.495 
#>  sd_path            0.0713
#>  mean_ci_width      0.338 
#>  median_ci_width    0.335 
#>  sd_ci_width        0.0372
#>  pct_ci_less_w      0.2   
#>  pct_ci_miss_low    0     
#>  pct_ci_miss_high   0.05  
#>  total_type_I_error 0.05  
#>  suc_rep            20    
#>  total_N            264   
#>  true_path          0.5   
#>  width              0.3   
#>  conf_level         0.95  
#> 
#> Confidence level: 95%
```
