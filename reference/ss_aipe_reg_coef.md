# Sample Size Planning for a Single Regression Coefficient (AIPE)

Computes the necessary sample size for the confidence interval on a
targeted regression coefficient \\\beta_j\\ (or its standardized
counterpart) in a multiple regression with *p* predictors to be no wider
than a user-specified value. This is the accuracy in parameter
estimation (AIPE) framework of Kelley and Maxwell (2003), targeted at a
specific coefficient rather than at the omnibus \\R^2\\. The
`noncentral = TRUE` variant inverts the noncentral *t* distribution of
the standardized \\b_j\\ under joint multivariate normality of the
predictors; the default central-*t* variant uses a closed form Wald
style approximation. Optionally, supplying `assurance` returns the
larger *N* that guarantees the realized width with the specified
probability rather than just on average.

## Usage

``` r
ss_aipe_reg_coef(
  rho2_Y_X = NULL,
  rho2_j_X_without_j = NULL,
  p = NULL,
  b_j = NULL,
  width,
  which_width = "Full",
  sigma_Y = 1,
  sigma_X = 1,
  rho_XX = NULL,
  rho_YX = NULL,
  which_predictor = NULL,
  noncentral = FALSE,
  alpha_lower = NULL,
  alpha_upper = NULL,
  conf_level = 0.95,
  assurance = NULL
)
```

## Arguments

- rho2_Y_X:

  Population value of \\\rho^2\_{Y \cdot X_1, \ldots, X_p}\\, the
  squared multiple correlation of the outcome *Y* with all *p*
  predictors.

- rho2_j_X_without_j:

  Population value of \\\rho^2\_{X_j \cdot X\_{-j}}\\, the squared
  multiple correlation when the *j*th predictor is regressed on the
  remaining \\p - 1\\ predictors. Quantifies the multicollinearity faced
  by the targeted coefficient.

- p:

  The number of predictor variables.

- b_j:

  The (unstandardized) regression coefficient for the *j*th predictor,
  the predictor of interest.

- width:

  Desired (full) width of the two-sided confidence interval on
  \\\beta_j\\.

- which_width:

  Which portion of the confidence interval `width` refers to. Only
  `"Full"` is currently implemented.

- sigma_Y:

  Population standard deviation of *Y*.

- sigma_X:

  Population standard deviation of the *j*th predictor.

- rho_XX:

  Population correlation matrix for the *p* predictor variables. If
  supplied with `rho_YX`, `rho2_Y_X` and `rho2_j_X_without_j` are
  derived from the covariance structure.

- rho_YX:

  Length-*p* vector of population correlations between *Y* and the *p*
  predictors.

- which_predictor:

  Which of the *p* predictors is the targeted coefficient.

- noncentral:

  If `TRUE`, plans using the exact noncentral *t* sampling distribution
  of the standardized \\b_j\\ (Kelley, 2007). If `FALSE` (the default),
  uses the central *t* approximation. The noncentral path requires
  `sigma_Y = sigma_X = 1` (*i.e.*, a standardized solution).

- alpha_lower:

  Type I error rate for the lower confidence limit.

- alpha_upper:

  Type I error rate for the upper confidence limit.

- conf_level:

  Confidence level (i.e., \\1 - \alpha\\, where \\\alpha\\ is the Type I
  error rate). Default `0.95`. Mutually exclusive with `alpha_lower` and
  `alpha_upper`.

- assurance:

  Optional probability with which the realized confidence interval is to
  be no wider than `width`. When `NULL` (the default), the planning
  targets the *expected* width.

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` value
is `"necessary_N"` and `value` is the necessary total sample size *N*
given the input specifications.

## Details

**Calling conventions.** The function offers several mutually exclusive
ways to supply the population information needed to plan \\N\\; the user
picks the one most aligned with their available planning values. Specify
exactly one of:

- *Covariance structure path.* Supply `rho_XX` and `rho_YX` along with
  `which_predictor`. The function derives \\\rho^2\_{Y \cdot X}\\ and
  \\\rho^2\_{X_j \cdot X\_{-j}}\\ from the covariance structure and also
  computes the population \\b_j\\ for consistency checks.

- *Squared multiple correlations path.* Supply `rho2_Y_X`,
  `rho2_j_X_without_j`, `p`, and `b_j` directly when the user already
  has these planning values from prior research and does not need to
  specify the full covariance structure.

- *Standardized solution.* For the `noncentral = TRUE` path, set
  `sigma_Y = sigma_X = 1`; the result returns the standardized
  regression coefficient sample size.

**Noncentral vs.\\ central planning.** The central *t* closed form is
fast and adequate at moderate to large *N*; the `noncentral = TRUE` path
additionally accounts for the noncentral *t* sampling distribution of
the standardized \\b_j\\ and is preferred when planning at small to
moderate *N* or when reporting planning that will be matched against the
noncentral CI from
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md).

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

## See also

[`ss_aipe_reg_coef_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef_sensitivity.md),
[`conf_limits_nct`](https://yelleknek.github.io/DMAR/reference/conf_limits_nct.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Covariance structure path: supply the population correlation
#    matrix and the population YX cross-correlations. Five predictors
#    in an exchangeable structure (all pairwise correlations 0.5,
#    all Y-X correlations 0.3).
rho_YX <- c(.3, .3, .3, .3, .3)
rho_XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5),
                c(.5, .5, 1, .5, .5), c(.5, .5, .5, 1, .5),
                c(.5, .5, .5, .5, 1))

# Closed-form (central t) planning, targeting the first predictor's
# standardized coefficient.
ss_aipe_reg_coef(width = .10, which_width = "Full",
                 sigma_Y = 1, sigma_X = 1,
                 rho_XX = rho_XX, rho_YX = rho_YX,
                 which_predictor = 1, noncentral = FALSE,
                 conf_level = .95)
#>  term        value
#>  necessary_N 2185 
#> 
#> Confidence level: 95%

# Adding assurance (.85): the realized CI width will be no larger than
# 0.10 in 85 percent of replications. Required N grows accordingly.
ss_aipe_reg_coef(width = .10, which_width = "Full",
                 sigma_Y = 1, sigma_X = 1,
                 rho_XX = rho_XX, rho_YX = rho_YX,
                 which_predictor = 1, noncentral = FALSE,
                 conf_level = .95, assurance = .85)
#>  term        value
#>  necessary_N 2259 
#> 
#> Confidence level: 95%

# Exact noncentral t planning. Required N differs at small to
# moderate samples.
ss_aipe_reg_coef(width = .10, which_width = "Full",
                 sigma_Y = 1, sigma_X = 1,
                 rho_XX = rho_XX, rho_YX = rho_YX,
                 which_predictor = 1, noncentral = TRUE,
                 conf_level = .95)
#>  term        value
#>  necessary_N 2191 
#> 
#> Confidence level: 95%

# 2. Squared multiple correlations path: when the user has planning
#    values for rho^2_Y.X and rho^2_j.X_-j directly (e.g., from a
#    prior power analysis), without specifying the full covariance
#    structure. b_j is required on this path.
ss_aipe_reg_coef(rho2_Y_X = 0.30, rho2_j_X_without_j = 0.20,
                 p = 5, b_j = 0.25,
                 width = .15, which_width = "Full",
                 sigma_Y = 1, sigma_X = 1,
                 noncentral = FALSE, conf_level = .95)
#>  term        value
#>  necessary_N 607  
#> 
#> Confidence level: 95%
```
