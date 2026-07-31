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
  save = FALSE,
  filename = "ss_aipe_sem_path_sensitivity_result.csv",
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

- save:

  option to save simulation results. It can be saved with `save = TRUE`
  outside of the printed results

- filename:

  the name of the file that simulation results will be saved to

- ...:

  allows one to potentially include parameter values for inner functions

## Value

A `data.frame` with columns `term` and `value` summarizing the a priori
Monte Carlo study. The `term` entries are `"total_N"` (the planned *N*),
`"desired_width"` (target CI width), `"mean_width"` and `"median_width"`
(realized CI width distribution across the `G` replications),
`"width_less_than_desired"` (proportion of realized widths at or below
`desired_width`), `"type_I_err_upper"` and `"type_I_err_lower"`
(tail-specific empirical Type I error rates), `"type_I_err"` (overall
empirical Type I error rate), `"conf_level"`, and `"suc_rep"` (number of
converged replications).

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
# \donttest{
set.seed(113)
pop_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f2 ~ 0.5*f1
  f1 ~~ 1*f1
  f2 ~~ 0.75*f2
  y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
  y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
"
Sigma <- cov_sem(pop_model)$sigma_theta
analysis_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f2 ~ b*f1
"
ss_aipe_sem_path_sensitivity(model = analysis_model, est_Sigma = Sigma,
                             true_Sigma = Sigma, which_path = "b",
                             desired_width = 0.30, N = 150, G = 25)
#>  term                    value
#>  total_N                 150  
#>  desired_width           0.3  
#>  mean_width              0.412
#>  median_width            0.41 
#>  width_less_than_desired 0    
#>  type_I_err_upper        0    
#>  type_I_err_lower        0    
#>  type_I_err              0    
#>  conf_level              0.95 
#>  suc_rep                 25   
#> 
#> Confidence level: 95%
# }
```
