# Confidence Interval for a Contrast in a Fixed Effects ANOVA

Computes the confidence interval for an unstandardized contrast of means
in a fixed effects analysis of variance, so a focused comparison among
groups (a pairwise difference or any weighted combination of the means)
is reported with its precision and in the units of the response.
Homogeneity of variance is assumed, as in the ANOVA on which `s_anova`
is based.

## Usage

``` r
ci_c(
  means = NULL,
  s_anova = NULL,
  c_weights = NULL,
  n = NULL,
  N = NULL,
  psi = NULL,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  df_error = NULL,
  ...
)
```

## Arguments

- means:

  A vector of the group means or the means of the particular level of
  the effect (for fixed effect designs)

- s_anova:

  The standard deviation of the errors from the ANOVA model (i.e., the
  square root of the mean square error)

- c_weights:

  The contrast weights (choose weights so that the positive *c*-weights
  sum to 1 and the negative *c*-weights sum to -1; i.e., use fractional
  values not integers)

- n:

  Sample sizes *per group* or level of the particular factor (if length
  1 it is assumed that the per group/level sample sizes are equal)

- N:

  Total sample size

- psi:

  The contrast effect, obtained by multiplying the *j*th mean by the
  *j*th contrast weight.

- conf_level:

  Confidence interval coverage (i.e., 1- Type I error rate); default is
  .95

- alpha_lower:

  Type I error for the lower confidence limit

- alpha_upper:

  Type I error for the upper confidence limit

- df_error:

  The degrees of freedom for the error. In one-way designs, this is
  simply *N*-length (means) and need not be specified; it must be
  specified if the design has multiple factors.

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 3-row `data.frame` with columns `term` and `value`. The `term` values
are `"lower_limit"` (the lower confidence limit on the population
contrast), `"contrast"` (the estimated unstandardized contrast), and
`"upper_limit"` (the upper limit).

## Note

Be sure to use the standard deviation and not the error variance for
`s_anova`, not the square of this value (the error variance) which would
come from the source table (i.e., use the root mean square error, not
the mean square error).

Be sure to use fractional *c*-weights when doing complex contrasts (not
integers) to specify `c_weights`. For example, in an ANCOVA of four
groups, if the user wants to compare the mean of group 1 and 2 with the
mean of group 3 and 4, `c_weights` should be specified as c(0.5, 0.5,
-0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the
contrast weights is zero.

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means.)

Steiger, J. H. (2004). Beyond the *F* Test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`ci_sc`](https://yelleknek.github.io/DMAR/reference/ci_sc.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Here is a four group example. Suppose that the means of groups 1--4 are 2, 4, 9,
# and 13, respectively. Further, let the error variance be .64 and thus the standard
# deviation would be .80 (note we use the standard deviation in the function, not the
# variance). The contrast of interest here is the average of groups 1 and 4 versus the
# average of groups 2 and 3.
ci_c(means = c(2, 4, 9, 13), s_anova = .80, c_weights = c(.5, -.5, -.5, .5),
n = c(3, 3, 3, 3), N = 12, conf_level = .95)
#>  term        value  
#>  lower_limit -0.0651
#>  contrast    1      
#>  upper_limit 2.07   
#> 
#> Confidence level: 95%

# Here is an example with two groups.
ci_c(means = c(1.6, 0), s_anova = .80, c_weights = c(1, -1),
n = c(10, 10), N = 20, conf_level = .95)
#>  term        value
#>  lower_limit 0.848
#>  contrast    1.6  
#>  upper_limit 2.35 
#> 
#> Confidence level: 95%

# An example given by Maxwell, Delaney, & Kelley (2027) :
# 24 subjects of mild hypertensives are assigned to one of four treatments: drug
# therapy, biofeedback, dietary modification, and a treatment combining all the
# three previous treatments. Subjects' blood pressure is measured two weeks
# after the termination of treatment. Now we want to form a 95% level
# confidence interval for the difference in blood pressure between subjects
# who received drug treatment and those who received biofeedback treatment

## Drug group's mean = 94; group size=4
## Biofeedback group's mean = 91; group size=6
## Diet group's mean = 92; group size=5
## Combination group's mean = 83; group size=5
## Mean Square Within (i.e., 'error_variance') = 67.375

ci_c(means = c(94, 91, 92, 83), s_anova = sqrt(67.375), c_weights = c(1, -1, 0, 0),
n = c(4, 6, 5, 5), N = 20, conf_level = .95)
#>  term        value
#>  lower_limit -8.23
#>  contrast    3    
#>  upper_limit 14.2 
#> 
#> Confidence level: 95%
```
