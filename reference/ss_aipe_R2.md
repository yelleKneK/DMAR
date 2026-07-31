# Sample Size Planning for Accuracy in Parameter Estimation for the Multiple Correlation Coefficient

Determines necessary sample size for the multiple correlation
coefficient so that the confidence interval for the population multiple
correlation coefficient is sufficiently narrow. Optionally, there is a
certainty parameter that allows one to be a specified percent certain
that the observed interval will be no wider than desired.

## Usage

``` r
ss_aipe_R2(
  population_R2 = NULL,
  conf_level = 0.95,
  width = NULL,
  random_predictors = TRUE,
  which_width = "Full",
  p = NULL,
  assurance = NULL,
  verify_ss = FALSE,
  tol = 1e-09,
  ...
)
```

## Arguments

- population_R2:

  Value of the population multiple correlation coefficient

- conf_level:

  Confidence interval level (e.g., .95, .99, .90); 1-Type I error rate

- width:

  Width of the confidence interval (see `which_width`)

- random_predictors:

  Whether or not the predictor variables are random (set to `TRUE`) or
  are fixed (set to `FALSE`)

- which_width:

  Defines the width that `width` refers to

- p:

  The number of predictor variables

- assurance:

  Value with which confidence can be placed that describes the
  likelihood of obtaining a confidence interval less than the value
  specified (e.g, .80, .90, .95)

- verify_ss:

  Evaluates numerically via an internal Monte Carlo simulation the exact
  sample size given the specifications

- tol:

  The tolerance of the iterative function `conf_limits_nct` for
  convergence

- ...:

  For modifying the parameters of functions this function calls upon

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` value
is `"necessary_N"` and `value` is the necessary total sample size *N*
given the input specifications.

## Details

This function determines a necessary sample size so that the expected
confidence interval width for the squared multiple correlation
coefficient is sufficiently narrow (when `assurance=NULL`) so that the
obtained confidence interval is no larger than the value specified with
some desired degree of certainty (i.e., a probability that the obtained
width is less than the specified width). The method depends on whether
or not the regressors are regarded as fixed or random. This is the case
because the distribution theory for the two cases is different and thus
the confidence interval procedure is conditional on the type of
regressors. The default methods are approximate but can be made exact
with the specification of `verify_ss=TRUE`, which performs an a priori
Monte Carlo simulation study. Kelley (2008) and Kelley & Maxwell (2008)
detail the methods used in the function, with the former focusing on
random regressors and the latter on fixed regressors.

It is recommended that the option `verify_ss` should always be used!
Doing so uses the method implied sample size as an estimate and then
evaluates with an internal Monte Carlo simulation (i.e., via
"brute-force" methods) the exact sample size given the goals specified.
When `verify_ss=TRUE`, the default number of iterations is 10,000 but
this can be changed by specifying G=5000 (or some other value; 10000 is
the recommended) When `verify_ss=TRUE` is specified, an internal
function `verify_ss_aipe_r2` calls upon the `ss_aipe_R2_sensitivity`
function for purposes of the internal Monte Carlo simulation study. See
the `verify_ss_aipe_r2` function for arguments that can be passed from
`ss_aipe_R2` to `verify_ss_aipe_r2`.

## Note

With `verify_ss = TRUE` the function can take some time to converge
(e.g., several minutes to a quarter hour) because the closed form
approximation is followed by an a priori Monte Carlo simulation. The
default `verify_ss = FALSE` returns the closed form approximation only
and is essentially instantaneous.

## References

Algina, J. & Olejnik, S. (2000). Determining sample size for accurate
estimation of the squared multiple correlation coefficient.
*Multivariate Behavioral Research, 35*, 119–137.
[doi:10.1207/s15327906mbr3501_5](https://doi.org/10.1207/s15327906mbr3501_5)

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
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
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md),
[`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Closed-form planner under random predictors (the typical case).
#    Sample size sufficient for the expected CI width on rho^2 to be .10.
ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, random_predictors = TRUE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_ncf for the meaning of the clamp.
#>  term        value
#>  necessary_N 773  
#> 
#> Confidence level: 95%

# 2. Same target under fixed predictors (planned dosing levels,
#    factorial covariates, etc.). Required N is typically smaller.
ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, random_predictors = FALSE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_ncf for the meaning of the clamp.
#>  term        value
#>  necessary_N 584  
#> 
#> Confidence level: 95%

# 3. Adding assurance (.85): the realized CI width will be no larger
#    than the target width in 85 percent of replications, not just on
#    average. Required N grows accordingly.
ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, assurance = .85,
           random_predictors = TRUE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_ncf for the meaning of the clamp.
#>  term        value
#>  necessary_N 815  
#> 
#> Confidence level: 95%

# 4. verify_ss = TRUE runs a Monte Carlo verification of the
#    closed-form approximation. Slow (minutes), so wrapped in
#    \donttest{}; uncomment to run.
# \donttest{
ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, random_predictors = TRUE,
           verify_ss = TRUE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_ncf for the meaning of the clamp.
#>  term        value
#>  necessary_N 771  
#> 
#> Confidence level: 95%
ss_aipe_R2(population_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, assurance = .85,
           random_predictors = FALSE, verify_ss = TRUE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in conf_limits_ncf() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?conf_limits_ncf for the meaning of the clamp.
#>  term        value
#>  necessary_N 622  
#> 
#> Confidence level: 95%
# }
```
