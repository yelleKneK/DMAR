# Asymptotic Variance of the Semipartial (Part) Correlation Coefficient

Computes the large-sample variance of the sample semipartial (also known
as the part) correlation coefficient \\r\_{Y(X \cdot Z_1 \cdots Z_J)}\\
under multivariate normality. The function provides two variances: the
asymptotic variance under the alternative (the partial-correlation form
applied by analogy) and, when the full-model coefficient of multiple
determination \\R^2\_{Y \cdot X Z_1 \cdots Z_J}\\ is supplied, the exact
null-hypothesis variance derived from the multiple-regression *F*-test
for the unique contribution of \\X\\ (Cohen, Cohen, West, & Aiken,
2003).

## Usage

``` r
var_semipartial_r(r_sp, n, J = 1, R2_full = NULL)
```

## Arguments

- r_sp:

  Sample semipartial correlation coefficient \\r\_{Y(X \cdot Z_1 \cdots
  Z_J)}\\, with the controlled variables partialled out of \\X\\ only
  (not \\Y\\). Must be in \\\[-1, 1\]\\.

- n:

  Total sample size.

- J:

  Number of variables partialled out of \\X\\ (i.e., the count of \\Z_1,
  \ldots, Z_J\\); must be at least 1. Defaults to 1.

- R2_full:

  Optional coefficient of multiple determination \\R^2\_{Y \cdot X Z_1
  \cdots Z_J}\\ from the full model that includes \\X\\ and all
  controls. When supplied, the null-hypothesis variance \\(1 - R^2)/(n -
  J - 2)\\ is returned instead of the alternative-side asymptotic
  variance. Must be in \\\[0, 1\]\\.

## Value

A one-row `data.frame` with columns `term` (either `"var_semipartial_r"`
or `"var_semipartial_r_under_null"`) and `value` (the requested
variance).

## Details

**Background.** The squared semipartial \\r^2\_{Y(X \cdot Z)}\\ equals
the increase in \\R^2\\ when \\X\\ is added to a model already
containing the controls \\Z_1, \ldots, Z_J\\, i.e., the unique variance
in \\Y\\ attributable to \\X\\. Unlike the partial, the semipartial is
on the original scale of \\Y\\ rather than on the partialled scale,
which makes it the natural effect size companion to standardized
regression coefficients in multiple-regression reports (Cohen et al.,
2003).

**Asymptotic variance (default).** Under multivariate normality the
semipartial admits the same large-sample form as the partial (Fisher,
1924, applied by analogy): \$\$\mathrm{Var}(\hat r\_{Y(X \cdot Z)})
\\\approx\\ \frac{(1 - \rho^2\_{Y(X \cdot Z)})^2}{n - J - 1}.\$\$ The
function evaluates this with \\\hat r\_{sp}\\ substituted for
\\\rho\_{sp}\\. This is the appropriate quantity for Wald-style
inference and for AIPE-style precision planning analogous to that of the
partial correlation. Aloe and Becker (2012) develop the asymptotic
variance of the semipartial as a function of the full population
correlation structure, and Yuan and Chan (2011) give exact higher-order
results for the closely related standardized regression coefficients;
the present approximation matches the leading \\1/n\\ behavior.

**Null-hypothesis variance (when `R2_full` is supplied).** In multiple
regression the unique contribution of \\X\\ is tested with \$\$F \\=\\
\frac{r^2\_{Y(X \cdot Z)}\\(n - J - 2)}{1 - R^2\_{Y \cdot X Z}} \\\sim\\
F(1,\\ n - J - 2)\$\$ under \\H_0\\: \rho\_{Y(X \cdot Z)} = 0\\ (Cohen
et al., 2003, equation 3.7.3). Equivalently \\t = \hat
r\_{sp}\\\sqrt{(n - J - 2)/(1 - R^2\_{Y \cdot X Z})}\\ is a
*t*-statistic on \\n - J - 2\\ degrees of freedom, so the
*under-the-null* variance of \\\hat r\_{sp}\\ is
\$\$\mathrm{Var}\_0(\hat r\_{sp}) \\=\\ \frac{1 - R^2\_{Y \cdot X
Z}}{n - J - 2}.\$\$ Supplying `R2_full` returns this null variance,
which is the standard ingredient for testing the significance of \\X\\'s
unique contribution.

## References

Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
predictors in meta-analysis. *Journal of Educational and Behavioral
Statistics, 37*(2), 278–297.
[doi:10.3102/1076998610396901](https://doi.org/10.3102/1076998610396901)

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum.

Fisher, R. A. (1924). The distribution of the partial correlation
coefficient. *Metron, 3*, 329–332.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on contrasts, Chapter 5 on multiple
comparisons, and Chapter 9 on ANCOVA.)

Yuan, K.-H., & Chan, W. (2011). Biases and standard errors of
standardized regression coefficients. *Psychometrika, 76*(4), 670–690.
[doi:10.1007/s11336-011-9224-6](https://doi.org/10.1007/s11336-011-9224-6)

## See also

[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
[`var_R2`](https://yelleknek.github.io/DMAR/reference/var_R2.md),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Olkin-Finn-style asymptotic variance of a semipartial r = .25 with
#     n = 100 and J = 3 controls in X.
var_semipartial_r(r_sp = 0.25, n = 100, J = 3)
#>  term              value  
#>  var_semipartial_r 0.00916

# With R^2 of the full model also supplied, the function returns the
#     null-hypothesis variance used in the F-test for X's unique
#     contribution.
var_semipartial_r(r_sp = 0.25, n = 100, J = 3, R2_full = 0.42)
#>  term                         value  
#>  var_semipartial_r_under_null 0.00611
```
