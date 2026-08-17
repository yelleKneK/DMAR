# Variance of the Mediated (Indirect) Effect \\ab\\

Computes the asymptotic variance of the product of two regression
coefficients \\\hat a \hat b\\ (the mediated/indirect effect in a simple
three-variable mediator model: \\X \to M \to Y\\) under four competing
formulas: Sobel (1982) first-order, Aroian (1947) / Goodman (1960)
second-order, and the full second-order delta method with optional
cross-product covariance. All four are reported in a single output so
the user can see the relative contributions of the higher-order terms.

## Usage

``` r
var_indirect_effect(a, b, var_a, var_b, cov_ab = 0)
```

## Arguments

- a, b:

  Anticipated population (or estimated) regression coefficients for \\X
  \to M\\ and \\M \to Y\\ (typically on the standardized scale).

- var_a, var_b:

  Variances (squared standard errors) of \\\hat a\\ and \\\hat b\\. For
  standardized regression with no covariates, \\\mathrm{Var}(\hat a)
  \approx (1 - a^2)/(n - 2)\\ and \\\mathrm{Var}(\hat b) \approx (1 -
  b^2)/\\(n - 3)(1 - a^2)\\\\, the standardized-model variance
  accounting for the correlation the \\a\\ path induces between the
  predictors of \\Y\\.

- cov_ab:

  Optional covariance between \\\hat a\\ and \\\hat b\\. Defaults to 0
  (the assumption underlying the standard Sobel formula). In practice
  the two estimators are nearly uncorrelated when the controls in the
  \\Y\\-on-\\M\\ regression are uncorrelated with the \\M\\-on-\\X\\
  regression's predictors.

## Value

A `data.frame` with rows for the four variance formulas; columns are
`term` and `value`.

## Details

**Sobel (1982) first-order.** The delta method variance of \\\hat a \hat
b\\ under independent \\\hat a\\, \\\hat b\\ is
\$\$\mathrm{Var}\_{\mathrm{Sobel}}(\hat a \hat b) \\=\\ a^2
\mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a).\$\$ This is the most
cited form and is the variance used by the standard Sobel *z*-test
(Sobel, 1982).

**Aroian (1947).** Aroian retains the second-order term:
\$\$\mathrm{Var}\_{\mathrm{Aroian}}(\hat a \hat b) \\=\\ a^2
\mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a) + \mathrm{Var}(\hat
a)\\\mathrm{Var}(\hat b).\$\$ Aroian shows this is exact under joint
normality of the two independent estimators.

**Goodman (1960).** Goodman's "unbiased" variance subtracts the
second-order term instead of adding it:
\$\$\mathrm{Var}\_{\mathrm{Goodman}}(\hat a \hat b) \\=\\ a^2
\mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a) - \mathrm{Var}(\hat
a)\\\mathrm{Var}(\hat b).\$\$ For small variances the three forms agree
to leading order; they diverge for noisy \\\hat a\\, \\\hat b\\.

**Second-order delta method (with covariance).** When \\\hat a\\ and
\\\hat b\\ share variability (e.g., they are both estimated from the
same regression of \\Y\\ on \\X\\ and \\M\\), the cross-product
covariance term enters: \$\$\mathrm{Var}(\hat a \hat b) \\\approx\\ a^2
\mathrm{Var}(\hat b) + b^2 \mathrm{Var}(\hat a) + 2 a b \cdot
\mathrm{Cov}(\hat a, \hat b).\$\$ MacKinnon et al.\\ (2002) show this
matters in models with covariates that simultaneously load on \\M\\ and
\\Y\\.

**Connection to
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md).**
The Sobel (delta method) variance is what
[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)
builds on under `method = "closed_form"` for AIPE planning; this
function makes the alternative formulas available for explicit
comparison.

## References

Aroian, L. A. (1947). The probability function of the product of two
normally distributed variables. *The Annals of Mathematical Statistics,
18*(2), 265–271.

Goodman, L. A. (1960). On the exact variance of products. *Journal of
the American Statistical Association, 55*(292), 708–713.

Lachowicz, M. J., Preacher, K. J., & Kelley, K. (2018). A novel measure
of effect size for mediation analysis. *Psychological Methods, 23*,
244–261. [doi:10.1037/met0000165](https://doi.org/10.1037/met0000165)

MacKinnon, D. P., Lockwood, C. M., Hoffman, J. M., West, S. G., &
Sheets, V. (2002). A comparison of methods to test mediation and other
intervening variable effects. *Psychological Methods, 7*(1), 83–104.
[doi:10.1037/1082-989X.7.1.83](https://doi.org/10.1037/1082-989X.7.1.83)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Preacher, K. J., & Kelley, K. (2011). Effect size measures for mediation
models: Quantitative strategies for communicating indirect effects.
*Psychological Methods, 16*(2), 93–115.
[doi:10.1037/a0022658](https://doi.org/10.1037/a0022658)

Sobel, M. E. (1982). Asymptotic confidence intervals for indirect
effects in structural equation models. *Sociological Methodology, 13*,
290–312.

Tofighi, D., & Kelley, K. (2020). Improved inference in mediation
analysis: Introducing the model-based constrained optimization
procedure. *Psychological Methods, 25*, 496–515.
[doi:10.1037/met0000259](https://doi.org/10.1037/met0000259)

## See also

[`ss_aipe_indirect_effect`](https://yelleknek.github.io/DMAR/reference/ss_aipe_indirect_effect.md)

Other variance utilities:
[`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
[`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`var_ete()`](https://yelleknek.github.io/DMAR/reference/var_ete.md),
[`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
[`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
[`var_smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. a = 0.40, b = 0.40, var_a = 0.02, var_b = 0.02, no covariance:
var_indirect_effect(a = 0.40, b = 0.40, var_a = 0.02, var_b = 0.02)
#>  term                   value 
#>  var_sobel              0.0064
#>  var_aroian             0.0068
#>  var_goodman            0.006 
#>  var_delta_second_order 0.0064

# 2. With a positive covariance between a-hat and b-hat:
var_indirect_effect(a = 0.40, b = 0.40,
                     var_a = 0.02, var_b = 0.02, cov_ab = 0.005)
#>  term                   value 
#>  var_sobel              0.0064
#>  var_aroian             0.0068
#>  var_goodman            0.006 
#>  var_delta_second_order 0.008 
```
