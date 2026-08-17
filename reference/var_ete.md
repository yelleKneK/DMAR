# Variance of the Estimated Treatment Effect in Two-Group ANCOVA With Heterogeneous Slopes

Computes the variance of the estimated treatment effect (ETE) at a
chosen covariate value in a two-group analysis of covariance with
heterogeneity of regression and a random covariate, following Li,
McLouth, and Delaney (2020). When the two groups' slopes differ, the
treatment effect is a function of the covariate, and its sampling
variance at the sample grand mean, one standard deviation from the mean,
or a fixed covariate value must account for the covariate being a random
variable rather than a set of fixed constants; the fixed-constant
formulas understate or misstate that variability. This is the
reimplementation of `var.ete()` from MBESS, contributed there by Li Li.

## Usage

``` r
var_ete(
  sigma2,
  sigma2_Z,
  n_1,
  n_2,
  beta_1,
  beta_2,
  mu_Z = 0,
  fixed_value = 0,
  type = c("sample", "population"),
  covariate_value = c("sample_mean", "sd", "fixed")
)
```

## Arguments

- sigma2:

  Residual error variance: the population value when
  `type = "population"`, the sample estimate when `type = "sample"`.

- sigma2_Z:

  Variance of the random covariate: population value or sample estimate,
  matching `type`.

- n_1, n_2:

  Sample sizes of the two groups (each must exceed 3; the formulas
  involve \\n - 3\\ in denominators).

- beta_1, beta_2:

  Slopes of the covariate in group 1 and group 2: population values or
  sample estimates, matching `type`.

- mu_Z:

  Mean of the covariate (population value or sample mean, matching
  `type`). Defaults to 0. Used when `covariate_value = "fixed"`.

- fixed_value:

  The fixed covariate value at which the treatment effect is assessed
  when `covariate_value = "fixed"`. Defaults to 0.

- type:

  `"sample"` (default) for the unbiased estimate of the variance from
  sample slopes and variances, or `"population"` for the variance from
  population values.

- covariate_value:

  Where the treatment effect is assessed: `"sample_mean"` (default) at
  the sample grand mean of the covariate, `"sd"` at the grand mean plus
  or minus one sample standard deviation, or `"fixed"` at `fixed_value`.

## Value

A `data.frame` (class `dmar_tbl`) with one row, `term = "var_ete"`,
whose `value` is the variance of the estimated treatment effect at the
chosen covariate value. The `type` and `covariate_value` choices are
recorded as attributes of the same names.

## Details

Randomized experiments with a covariate commonly probe the simple
treatment effect at the grand mean and one standard deviation either
side of it when the slopes differ across groups. The variance
expressions here treat the covariate as normally distributed rather than
fixed, which Li, McLouth, and Delaney (2020) show can change the
estimated standard error substantially when heterogeneity of regression
is strong. The square root of the returned value is the standard error
used for a confidence interval or test of the treatment effect at the
chosen covariate value.

At the sample grand mean of the covariate, writing \\N = n_1 + n_2\\,
the population variance (their Equation 10) is \$\$\mathrm{Var} =
\sigma^2 C_0 + \frac{(\beta_1 - \beta_2)^2 \sigma^2_Z}{N}, \qquad C_0 =
\frac{1}{n_1} + \frac{1}{n_2} + \frac{n_2}{N n_1 (n_1 - 3)} +
\frac{n_1}{N n_2 (n_2 - 3)},\$\$ and with `type = "sample"` the returned
value is their unbiased estimator (Equation C.7), which subtracts
\\\sigma^2 \\(N-3)/(n_1-3) + (N-3)/(n_2-3)\\ / \\N (N-1)\\\\ so that
plugging in sample estimates does not overstate the variance. The
`covariate_value = "sd"` expressions are their Equations 12 and C.9,
which add the variance contribution of estimating the covariate's
standard deviation, and the `"fixed"` expressions are their Equations 14
and C.10.

The two `"fixed"` estimands differ in where the deviation of
`fixed_value` is measured from. With `type = "sample"` the deviation is
taken from the sample grand mean (Equation C.10), so `mu_Z` should be
the sample mean of the covariate. With `type = "population"` the
deviation is taken from a known population mean (Equation 14);
evaluating the treatment effect at that known mean itself, as in the
paper's worked example, sets `fixed_value = mu_Z`, which zeroes the
deviation term.

## References

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
[doi:10.3758/BF03192993](https://doi.org/10.3758/BF03192993)

Li, L., McLouth, C. J., & Delaney, H. D. (2020). Analysis of covariance
in randomized experiments with heterogeneity of regression and a random
covariate: The variance of the estimated treatment effect at selected
covariate values. *Multivariate Behavioral Research, 55*(6), 926–940.
[doi:10.1080/00273171.2019.1693953](https://doi.org/10.1080/00273171.2019.1693953)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 on heterogeneity of regression.)

## See also

[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md) for the
model whose treatment effect this variance describes;
[`regions_of_significance`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md)
for the companion question of where a moderated effect is
distinguishable from zero.

Other variance utilities:
[`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
[`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md),
[`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
[`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
[`var_smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Pygmalion data (Maxwell, Delaney, & Kelley, 2027): the treatment
# effect of the "Bloomer" expectation at the covariate grand mean,
# with heterogeneous pre-IQ slopes.
data(pygmalion)
fit <- lm(iq_8 ~ iq_pre * treatment, data = pygmalion)
s2  <- sum(residuals(fit)^2) / fit$df.residual
var_ete(sigma2 = s2, sigma2_Z = var(pygmalion$iq_pre),
        n_1 = sum(pygmalion$treatment == "Bloomer"),
        n_2 = sum(pygmalion$treatment == "Control"),
        beta_1 = coef(fit)["iq_pre"] + coef(fit)["iq_pre:treatmentBloomer"],
        beta_2 = coef(fit)["iq_pre"])
#>  term    value
#>  var_ete 3.52 
```
