# Sample Size Planning for Accuracy in Parameter Estimation for Reliability Coefficients

Computes the necessary sample size for the confidence interval on a
population reliability coefficient (coefficient alpha or coefficient
omega, depending on the assumed measurement model) to have expected
width no larger than `width`, or (when `assurance` is supplied) to be no
wider than `width` with the specified probability. This is the accuracy
in parameter estimation (AIPE) counterpart to power based planning for
reliability and is the companion of
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md).
The closed form planner uses the Terry and Kelley (2012) formulas under
the user-selected measurement model and confidence interval type; when
`assurance` is supplied, the function follows the closed form with an
internal Monte Carlo simulation to find the smallest *N* delivering the
requested assurance.

Coefficient alpha (Guttman, 1945; subsequently popularized by Cronbach,
1951) is returned for the parallel and tau equivalent (*i.e.*, the so
called True Score) measurement models; coefficient omega (McDonald,
1999) is returned for the congeneric model.

## Usage

``` r
ss_aipe_reliability(
  model = NULL,
  type = NULL,
  width = NULL,
  S = NULL,
  conf_level = 0.95,
  assurance = NULL,
  data = NULL,
  i = NULL,
  cor_est = NULL,
  lambda = NULL,
  psi_square = NULL,
  initial_iter = 500,
  final_iter = 5000,
  start_ss = NULL,
  verbose = FALSE
)
```

## Arguments

- model:

  The measurement model assumed for the population. Accepts
  (case-sensitive aliases shown in parentheses): `"Parallel"`
  (`"parallel"`, `"SB"`, `"Spearman Brown"`, `"Spearman-Brown"`, `"sb"`)
  for the strictly parallel items model; `"True Score"`
  (`"True Score Equivalent"`, `"True-Score Equivalent"`, `"Equivalent"`,
  `"Tau Equivalent"`, `"tau-equivalent"`, `"Tau-Equivalent"`,
  `"True-Score"`, `"true-score"`, `"true score"`, `"Cronbach"`,
  `"cronbach"`, `"Chronbach"`, `"alpha"`) for the tau equivalent model
  (in which case the function plans for coefficient alpha); or
  `"Congeneric"` (`"congeneric"`, `"omega"`, `"Omega"`) for the
  congeneric model (in which case the function plans for coefficient
  omega).

- type:

  The method used to construct the confidence interval on the
  reliability coefficient: either `"Factor Analytic"` (McDonald, 1999),
  available for all three measurement models, or `"Normal Theory"` (van
  Zyl, Neudecker, & Nel, 2000), available for the parallel and tau
  equivalent models only.

- width:

  The desired full width of the two-sided confidence interval.

- S:

  A symmetric population covariance (or correlation) matrix among the
  items, used to imply the population reliability and its sampling
  distribution.

- conf_level:

  Confidence level (i.e., \\1 - \alpha\\, where \\\alpha\\ is the Type I
  error rate). Default `0.95`.

- assurance:

  Optional probability with which the realized interval is to be no
  wider than `width`. When `NULL` (the default), the planner targets the
  *expected* width; when supplied (e.g., 0.80, 0.85, 0.95), the function
  follows the closed form with a Monte Carlo search.

- data:

  A data set from which the population covariance matrix should be
  inferred.

- i:

  Number of items.

- cor_est:

  The presumed inter-item correlation. One value for the parallel and
  tau equivalent models.

- lambda:

  Vector of population factor loadings.

- psi_square:

  Vector of population unique (error) variances.

- initial_iter:

  Number of Monte Carlo iterations used in the initial assurance search.

- final_iter:

  Number of Monte Carlo iterations used in the final assurance
  verification.

- start_ss:

  Optional starting sample size for the iterative assurance search.

- verbose:

  If `TRUE`, prints the current sample size and empirical assurance at
  each step of the Monte Carlo search.

## Value

A `data.frame` with columns `term` and `value`. Without `assurance` the
data frame has a single row, `"necessary_N"`, giving the necessary *N*.
With `assurance` supplied, the data frame has five rows: `"necessary_N"`
(necessary *N*), `"width"` (echo of the target width),
`"specified_assurance"` (echo of the requested probability),
`"empirical_assurance"` (the assurance achieved at the returned *N* in
the Monte Carlo verification), and `"final_iter"` (number of Monte Carlo
iterations used).

## Details

The Monte Carlo assurance search simulates covariance matrices from the
population implied by the inputs and, at each candidate sample size,
computes the realized confidence interval width with the same machinery
the estimation side of the package uses. For `type = "Factor Analytic"`
the single-factor model is fit by maximum likelihood (equal loadings for
the parallel and tau equivalent models, free loadings for the congeneric
model) and the interval is the delta method Wald interval on the model
implied reliability, the interval of
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)`(denominator = "model_implied", ci_method = "ml")`
and of
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)`(estimator = "model_implied", ci_method = "ml")`.
For `type = "Normal Theory"` the interval uses the van Zyl, Neudecker,
and Nel (2000) closed form standard error for the tau equivalent model
and its compound symmetry simplification for the parallel model, the
same closed form behind `reliability_alpha(ci_method = "ml")`. The
congeneric model has no normal theory form, so `type = "Normal Theory"`
is an error there.

## Note

Not all of the items can be entered into the function to represent the
population values. For example, either 'data' can be used, or `S`, or
`i`, `cor_est`, and `psi_square`, or `i`, `lambda`, and `psi_square`.
With a large number of iterations (`final_iter`) this function may take
considerable time.

## Warning

In some conditions the factor analytic model fit by
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) (via
lavaan) may fail to converge, and you may see a non-convergence message
from lavaan. The Monte Carlo assurance search treats a non-converged
iteration as missing and continues, so a few such messages do not
invalidate the result. Frequent non-convergence usually means the model
is poorly determined by the data, for example because of a small sample
size, a low number of iterations, or a poorly behaved covariance matrix.

## References

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

McDonald, R. P. (1999). *Test theory: A unified treatment*. Mahwah, NJ:
Lawrence Erlbaum Associates.

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000). On the distribution
of the maximum likelihood estimator of Cronbach's alpha. *Psychometrika,
65*(3), 271–280.
[doi:10.1007/BF02296146](https://doi.org/10.1007/BF02296146)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Expected confidence interval width (closed form, no Monte Carlo search).
ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
  i = 6, cor_est = .3, psi_square = .2, conf_level = .95, assurance = NULL)
#>  term        value
#>  necessary_N 38   
#> 
#> Confidence level: 95%

# \donttest{
# The assurance cases run a Monte Carlo search; the iteration counts below are
# reduced so the example runs quickly. Raise initial_iter and final_iter for a
# production plan, and set.seed() for a reproducible search.
set.seed(113)

# Same population, now targeting an assurance.
ss_aipe_reliability(model = "Parallel", type = "Normal Theory", width = .1,
  i = 6, cor_est = .3, psi_square = .2, conf_level = .95, assurance = .85,
  initial_iter = 50, final_iter = 200)
#>  term                value
#>  necessary_N         57   
#>  width               0.1  
#>  specified_assurance 0.85 
#>  empirical_assurance 0.875
#>  final_iter          200  
#> 
#> Confidence level: 95%

# True score (tau equivalent) model. Here psi_square is a vector of length i
# (number of items), while cor_est stays a single value.
ss_aipe_reliability(model = "True Score", type = "Normal Theory", width = .1,
  i = 5, cor_est = .3, psi_square = c(.2, .3, .3, .2, .3), conf_level = .95,
  assurance = .85, initial_iter = 50, final_iter = 200)
#>  term                value
#>  necessary_N         115  
#>  width               0.1  
#>  specified_assurance 0.85 
#>  empirical_assurance 0.85 
#>  final_iter          200  
#> 
#> Confidence level: 95%

# Congeneric model with the factor analytic approach (coefficient omega).
# Each Monte Carlo iteration fits a one factor model, so this is the slow case.
ss_aipe_reliability(model = "Congeneric", type = "Factor Analytic", width = .15,
  i = 4, lambda = c(.8, .7, .7, .8), psi_square = c(.4, .5, .5, .4),
  conf_level = .95, assurance = .80, initial_iter = 10, final_iter = 20)
#>  term                value
#>  necessary_N         60   
#>  width               0.15 
#>  specified_assurance 0.8  
#>  empirical_assurance 0.8  
#>  final_iter          20   
#> 
#> Confidence level: 95%

# Planning from a presumed population correlation matrix among the items.
pop_mat <- rbind(
  c(1.0000000, 0.3813850, 0.4216370, 0.3651484, 0.4472136),
  c(0.3813850, 1.0000000, 0.4020151, 0.3481553, 0.4264014),
  c(0.4216370, 0.4020151, 1.0000000, 0.3849002, 0.4714045),
  c(0.3651484, 0.3481553, 0.3849002, 1.0000000, 0.4082483),
  c(0.4472136, 0.4264014, 0.4714045, 0.4082483, 1.0000000))
ss_aipe_reliability(model = "True Score", type = "Normal Theory", width = .15,
  S = pop_mat, conf_level = .95, assurance = .85, initial_iter = 50,
  final_iter = 200)
#>  term                value
#>  necessary_N         117  
#>  width               0.15 
#>  specified_assurance 0.85 
#>  empirical_assurance 0.855
#>  final_iter          200  
#> 
#> Confidence level: 95%
# }
```
