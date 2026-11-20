# Sample Size Planning for Accuracy in Parameter Estimation for the Multiple Correlation Coefficient

Determines necessary sample size for the multiple correlation
coefficient so that the confidence interval for the population multiple
correlation coefficient is sufficiently narrow. Optionally, there is a
certainty parameter that allows one to be a specified percent certain
that the observed interval will be no wider than desired.

## Usage

``` r
ss_aipe_R2(
  true_R2 = NULL,
  conf_level = 0.95,
  width = NULL,
  random_predictors = TRUE,
  which_width = "Full",
  p = NULL,
  assurance = NULL,
  verify_ss = FALSE,
  tol = 1e-09,
  G_initial = 500,
  G = 10000
)
```

## Arguments

- true_R2:

  Population squared multiple correlation coefficient the sample size is
  planned for. The planner treats this value as the truth;
  [`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md)
  shows what the plan delivers when the truth (`true_R2` there) differs
  from the value planned with (`supposed_R2` there).

- conf_level:

  Confidence interval level (e.g., .95, .99, .90); 1-Type I error rate

- width:

  Width of the confidence interval (see `which_width`)

- random_predictors:

  Whether or not the predictor variables are random (set to `TRUE`) or
  are fixed (set to `FALSE`)

- which_width:

  Which width `width` refers to: `"Full"` (the default), the full width
  of the confidence interval; `"Lower"`, the distance from the expected
  sample \\R^2\\ to the lower limit; or `"Upper"`, the distance from the
  upper limit to the expected sample \\R^2\\. An abbreviation such as
  `"L"` stands for the value it begins. With `"Lower"` or `"Upper"` the
  plan is for that distance at the expected sample \\R^2\\: the function
  returns the smallest *N* such that the distance is no larger than
  `width` at that *N* and at every larger sample size. (The search
  begins at \\N = 2p + 3\\, so no such plan is smaller.) The distance to
  the upper limit does not fall steadily as *N* grows. With several
  predictors and a small *N* the expected sample \\R^2\\ is well above
  the population value and the upper limit is close to that expected
  value or below it, so the distance is small or negative; it rises to a
  maximum as *N* grows and then falls. With `true_R2 = .70`,
  `width = .05`, `p = 10`, and `which_width = "Upper"` the distance is
  .043 at *N* = 23, is larger than .05 at every *N* from 27 to 201, and
  the function returns 202.
  [`expected_R2`](https://yelleknek.github.io/DMAR/reference/expected_R2.md)
  and [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md)
  give the distance at any *N*. A very small `width` (.0025 or less in
  the calls examined) can be narrower than can be planned with R's
  distribution functions, and such a half-width call is refused with an
  error that names `width`. To confirm that the distance stays within
  `width` the search needs the limits at sample sizes several times the
  plan, and there the noncentrality parameter is in the millions, beyond
  the range in which R's noncentral *F* and beta distribution functions
  ([`pf`](https://rdrr.io/r/stats/Fdist.html) and
  [`pbeta`](https://rdrr.io/r/stats/Beta.html)) are accurate (their help
  says the code is meant for moderate values of the noncentrality
  parameter), or the sample sizes pass 10,000,000, the largest the
  search tries. The refusal does not mean that no sample size meets
  `width`. `assurance` can be used only with `"Full"`: given with
  `"Lower"` or `"Upper"` it is refused with an error, because the
  assurance step does not deliver the requested probability for a
  half-width. `verify_ss = TRUE` is refused with a half-width as well.

- p:

  The number of predictor variables

- assurance:

  Value with which confidence can be placed that describes the
  likelihood of obtaining a confidence interval less than the value
  specified (e.g, .80, .90, .95). It can be used only with
  `which_width = "Full"`.

- verify_ss:

  Evaluates numerically via an internal Monte Carlo simulation the exact
  sample size given the specifications

- tol:

  The tolerance of the iterative function `ci_nc_t` for convergence

- G_initial:

  Used only with `verify_ss = TRUE`: the number of Monte Carlo
  replications for each candidate *N* in the coarse search that brackets
  the sample size (default 500).

- G:

  Used only with `verify_ss = TRUE`: the number of Monte Carlo
  replications in the final pass that confirms the sample size near that
  bracket (default 10000).

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

**What the assurance plan delivers.** The assurance procedure is not
exact. In simulation with a population \\R^2\\ of .10, 4 predictors, and
width .20 the realized probability was .86 for a requested .90 and .73
for a requested .80; with a population \\R^2\\ of .40 and 2 predictors
it was 1.00 for a requested .90. The largest shortfall found was at a
population \\R^2\\ of .10 with 5 predictors and width .30, where the
realized probability was .67 for a requested .80. (Each value is the
proportion of 2,000 simulated samples, at the planned sample sizes of
172, 156, 227, and 75, whose interval was no wider than the target; its
Monte Carlo standard error is .008 to .011 for the three that fell
short.) The shortfall appears where the population \\R^2\\ is small and
there are several predictors. These are the plans that were simulated,
not a limit on the shortfall, so a plan that was not simulated can fall
further short. Check a plan with
[`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md),
whose `pct_ci_less_width` row is the realized probability at the planned
sample size, or use `verify_ss = TRUE`.

Using `verify_ss = TRUE` is recommended. It takes the method implied
sample size as a starting value and then finds, by Monte Carlo
simulation, the sample size that meets the goals specified. When
`verify_ss=TRUE`, the default number of iterations is 10,000 but this
can be changed by specifying G=5000 (or some other value; 10000 is the
recommended). When `verify_ss=TRUE` is specified, an internal function
`verify_ss_aipe_r2` calls upon the `ss_aipe_R2_sensitivity` function for
purposes of the internal Monte Carlo simulation study, with `G_initial`
replications for each candidate *N* in the coarse search that brackets
the answer and `G` replications in the final pass that confirms the
sample size near that bracket.

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
[`ci_nc_t`](https://yelleknek.github.io/DMAR/reference/ci_nc_t.md),
[`ss_aipe_R2_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Closed form planner under random predictors (the typical case).
#    Sample size sufficient for the expected CI width on rho^2 to be .10.
ss_aipe_R2(true_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, random_predictors = TRUE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in ci_nc_F() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?ci_nc_F for the meaning of the clamp.
#>  term        value
#>  necessary_N 773  
#> 
#> Confidence level: 95%

# 2. The same target under fixed predictors (planned dosing levels,
#    factorial covariates, and the like) needs a smaller N, since fixed
#    predictors contribute no sampling variability of their own.
ss_aipe_R2(true_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, random_predictors = FALSE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in ci_nc_F() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?ci_nc_F for the meaning of the clamp.
#>  term        value
#>  necessary_N 584  
#> 
#> Confidence level: 95%

# 3. An assurance of .85, so that the realized width is no larger than the
#    target in 85 percent of replications rather than only on average,
#    needs a larger N than the expected width plan in (1).
ss_aipe_R2(true_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, assurance = .85,
           random_predictors = TRUE)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in ci_nc_F() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?ci_nc_F for the meaning of the clamp.
#>  term        value
#>  necessary_N 815  
#> 
#> Confidence level: 95%

# 4. verify_ss = TRUE follows the closed form approximation with an a
#    priori Monte Carlo simulation of the realized width, starting from
#    the closed form answer and returning the sample size the simulation
#    settles on, which is what a plan meant to be defended deserves. The
#    coarse search runs G_initial replications per candidate N and the final pass
#    runs G; the small counts here keep the example quick, and a reported
#    plan deserves the defaults of G_initial = 500 and G = 10000.
set.seed(113)
ss_aipe_R2(true_R2 = .50, conf_level = .95, width = .10,
           which_width = "Full", p = 5, random_predictors = TRUE,
           verify_ss = TRUE, G_initial = 10, G = 30)
#> Warning: During the iterative sample size search, the noncentral F lower-limit clamp in ci_nc_F() fired in 6 intermediate evaluations. The returned sample size accounts for this; see ?ci_nc_F for the meaning of the clamp.
#>  term        value
#>  necessary_N 766  
#> 
#> Confidence level: 95%
```
