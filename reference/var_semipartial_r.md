# Asymptotic Variance of the Semipartial (Part) Correlation Coefficient

Computes a variance for the sample semipartial (also known as the part)
correlation coefficient \\r\_{Y(X \cdot Z_1 \cdots Z\_{p-1})}\\ under
multivariate normality. The function provides two quantities, both of
which need the full-model coefficient of multiple determination
\\R^2\_{Y \cdot X Z_1 \cdots Z\_{p-1}}\\: the squared standard error in
the *t* test of \\X\\'s coefficient (the default when `R2_full` is
supplied; Cohen, Cohen, West, & Aiken, 2003), and the first-order
normal-theory variance of the sample semipartial correlation
(`under_null = FALSE`), which is the quantity for a Wald confidence
interval and for sample size planning. Without `R2_full` the function
returns the first-order variance for the special case in which the
variables partialed out explain none of the outcome, and warns.

## Usage

``` r
var_semipartial_r(r_sp, N, p = 2, R2_full = NULL, under_null = TRUE)
```

## Arguments

- r_sp:

  Sample semipartial correlation coefficient \\r\_{Y(X \cdot Z_1 \cdots
  Z\_{p-1})}\\, with the other variables partialed out of \\X\\ only
  (not \\Y\\). Must be in \\\[-1, 1\]\\.

- N:

  Total sample size.

- p:

  Number of predictors in the regression of *Y* on *X* and the variables
  partialed out, so `p - 1` variables are partialed out and `p` is at
  least 2. A categorical variable with *a* levels counts as *a* - 1
  predictors. It is the `p` of
  [`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
  [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md), and
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md):
  the correlation concerns one predictor of a regression with `p`
  predictors. In the notation of a model comparison, \\q\\ is the number
  of predictors removed from the full model to give the restricted model
  (\\q = p\_{full} - p\_{restricted}\\, the predictors being tested),
  and the number of variables partialed out is \\p - q\\; this
  correlation concerns one predictor, so \\q = 1\\ and \\p - 1\\
  variables are partialed out. Defaults to 2 (one variable partialed
  out).

- R2_full:

  Coefficient of multiple determination \\R^2\_{Y \cdot X Z_1 \cdots
  Z\_{p-1}}\\ from the full model that includes \\X\\ and every variable
  partialed out (`R2_Y_X` in
  [`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md);
  the squared semipartial correlation is `R2_full` minus the squared
  multiple correlation of the model without *X*). Must be in \\\[0,
  1\]\\ and cannot be smaller than `r_sp` squared, because the squared
  semipartial correlation is the part of the full-model \\R^2\\ that
  \\X\\ adds; a smaller value is refused with an error of class
  `dmar_semipartial_r_incoherent`, unless the two values could have come
  from rounding to two decimals (see Details). When `NULL` (the
  default), the value returned is the first-order variance that holds
  when the variables partialed out explain none of the outcome, with a
  warning of class `dmar_semipartial_r_no_R2_full`.

- under_null:

  Logical, used only when `R2_full` is supplied. If `TRUE` (the
  default), the function returns \\(1 - R^2)/(N - p - 1)\\, the squared
  standard error in the *t* test of \\X\\'s coefficient. If `FALSE`, it
  returns the first-order normal-theory variance of the sample
  semipartial correlation, evaluated at `r_sp` and `R2_full`.

## Value

A one-row `data.frame` with columns `term` (either `"var_semipartial_r"`
or `"var_semipartial_r_under_null"`) and `value` (the requested
variance).

## Details

**Background.** The squared semipartial \\r^2\_{Y(X \cdot Z)}\\ equals
the increase in \\R^2\\ when \\X\\ is added to a model already
containing the variables partialed out, \\Z_1, \ldots, Z\_{p-1}\\, that
is, the unique variance in \\Y\\ attributable to \\X\\. Unlike the
partial, the semipartial is on the original scale of \\Y\\ rather than
on the partialed scale, which makes it the natural effect size companion
to standardized regression coefficients in multiple regression reports
(Cohen et al., 2003). Aloe and Becker (2012) propose it as an effect
size for regression predictors in meta-analysis.
[`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md)
plans the power of the test of that change in \\R^2\\.

**First-order variance (`under_null = FALSE`).** Write \\\rho\_{sp}\\
for the population semipartial correlation, \\P^2\\ for the population
full-model squared multiple correlation, and \\P_r^2 = P^2 -
\rho\_{sp}^2\\, the squared multiple correlation of \\Y\\ with \\Z_1,
\ldots, Z\_{p-1}\\ alone. Under multivariate normality, to first order
in \\1/N\\, \$\$N\\\mathrm{Var}(\hat r\_{sp}) \\\to\\ (1 -
\rho\_{sp}^2)^2 \\-\\ P_r^2\\(1 - 2 \rho\_{sp}^2).\$\$ The variance
depends on the population through \\\rho\_{sp}\\ and \\P_r^2\\ only.
Three cases show how it behaves. When the variables partialed out
explain none of the outcome (\\P_r^2 = 0\\) it is \\(1 -
\rho\_{sp}^2)^2\\, the form of a simple or partial correlation. When
\\\rho\_{sp} = 0\\ it is \\1 - P^2\\. For \\P_r^2 \> 0\\, it is smaller
than \\(1 - \rho\_{sp}^2)^2\\ when \\\rho\_{sp}^2 \< 1/2\\ and larger
when \\\rho\_{sp}^2 \> 1/2\\. Because the variance depends on \\P_r^2\\,
no function of \\r\_{sp}\\ alone stabilizes it, and Fisher's *Z*
transformation does not apply to a semipartial correlation as it does to
a partial correlation. The expression is the first-order variance of
Aloe and Becker (2012, Equation 4), which they obtain from Alf and Graf
(1999), with \\N - p\\ in place of their \\n\\: their numerator, \\P^4 -
2 P^2 + P\_{(f)}^2 + 1 - P\_{(f)}^4\\ with \\P\_{(f)}^2\\ the squared
multiple correlation of the model without \\X\\, equals the expression
above with \\P_r^2 = P\_{(f)}^2\\. The result follows from the delta
method variance of the difference between two nested squared multiple
correlations (Model A of Olkin & Finn, 1995, for one variable partialed
out; Alf & Graf, 1999, for any number), because \\r\_{sp}^2\\ is that
difference and so \\\mathrm{Var}(\hat r\_{sp}) \approx
\mathrm{Var}(\widehat{\Delta R^2}) / (4 \rho\_{sp}^2)\\.

The function evaluates the expression at `r_sp` and `R2_full` and
divides by \\N - p\\. The number of variables partialed out does not
enter the first-order variance. The divisor \\N - p\\ is a convention,
kept for agreement with
[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md)
and because, with estimated values, a Wald interval formed with it has
coverage closer to the nominal level than one formed with \\N\\; it is
not the most accurate estimate of the variance itself (in simulation
neither \\N - 1\\ nor \\N - p\\ was closer to the observed variance in
every condition).

The approximation has limits that simulation shows and that no theory
given here corrects: it understates the variance when \\\|\rho\_{sp}\|\\
is large and \\N - p\\ is small (in one simulated condition,
\\\rho\_{sp} = .95\\, \\P_r^2 = 0\\, \\p = 2\\, and \\N = 40\\, it was
about 39 percent of the variance observed), and the sample semipartial
correlation is biased toward zero by roughly \\\rho\_{sp} (p - 1) /
(2N)\\. A Wald interval formed from it also misses the population value
too often when the full-model \\R^2\\ is near 1 and \\N\\ is small, for
any \\\rho\_{sp}\\ (in one simulated condition, \\\rho\_{sp} = .1\\, a
full-model \\R^2\\ of .99, \\p = 4\\, and \\N = 12\\, it missed in .13
of samples). When \\\|\rho\_{sp}\|\\ is above about .9 and the variables
partialed out explain some of the outcome, the interval misses more
often on one side than the other even with several hundred observations
(in one simulated condition, \\\rho\_{sp} = .95\\, a full-model \\R^2\\
of .97, \\p = 2\\, and \\N = 395\\, it missed in .064 of samples, .051
of them on one side).

**Values rounded to two decimals.** Reported values are often rounded,
and a rounded pair can have `R2_full` slightly below `r_sp` squared: a
semipartial correlation reported as .46 has a square of .2116, and the
full-model \\R^2\\ beside it may be reported as .21. When
`R2_full + .005` is at least \\(\|r\_{sp}\| - .005)^2\\ (with the
difference in parentheses not below 0), the pair is consistent with
rounding. The function then gives a warning of class
`dmar_semipartial_r_rounded`, which carries the two values, and
proceeds: the squared semipartial correlation cannot exceed the
full-model \\R^2\\, so the variables partialed out are taken to explain
none of the outcome (\\P_r^2 = 0\\) in the first-order variance, and the
squared standard error of the *t* test uses `R2_full` as given. A pair
further apart than rounding allows is refused.
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md)
and
[`ss_aipe_semipartial_r_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r_sensitivity.md)
apply the same rule to their values.

**Without `R2_full`.** The function returns \\(1 - r\_{sp}^2)^2 / (N -
p)\\, the first-order variance at \\P_r^2 = 0\\. When the variables
partialed out explain some of the outcome, this value is too large if
\\\rho\_{sp}^2 \< 1/2\\ and too small if \\\rho\_{sp}^2 \> 1/2\\, which
is why the call warns. At \\r\_{sp} = 0\\ it is \\1/(N - p)\\ whatever
the full-model \\R^2\\, where the first-order variance is \\(1 -
P^2)/N\\.

**Squared standard error of the *t* test (`R2_full` supplied,
`under_null = TRUE`).** In multiple regression the unique contribution
of \\X\\ is tested with \$\$F \\=\\ \frac{r^2\_{Y(X \cdot Z)}\\(N - p -
1)}{1 - R^2\_{Y \cdot X Z}} \\\sim\\ F(1,\\ N - p - 1)\$\$ under
\\H_0\\: \rho\_{Y(X \cdot Z)} = 0\\ (Cohen et al., 2003). Equivalently
\\t = \hat r\_{sp}\\\sqrt{(N - p - 1)/(1 - R^2\_{Y \cdot X Z})}\\ is a
*t*-statistic on \\N - p - 1\\ degrees of freedom, and the function
returns the square of its denominator, \$\$\frac{1 - R^2\_{Y \cdot X
Z}}{N - p - 1},\$\$ the squared standard error in the *t* test of
\\X\\'s coefficient. It uses the sample \\R^2\\ and equals the
first-order variance at \\\rho\_{sp} = 0\\; it is the ingredient for
testing the significance of \\X\\'s unique contribution, not a variance
for an interval about a nonzero \\r\_{sp}\\.

With `under_null = TRUE` the value is \\(1 - R^2\_{full})/(N - p - 1)\\,
which equals \\r\_{sp}^2/t^2\\ for the *t* statistic of the predictor's
coefficient, by the identity \\r\_{sp} = t\sqrt{(1 - R^2\_{full})/(N -
p - 1)}\\ (Aloe & Becker, 2012, Equation 2).

## References

Alf, E. F., Jr., & Graf, R. G. (1999). Asymptotic confidence limits for
the difference between two squared multiple correlations: A simplified
approach. *Psychological Methods, 4*(1), 70–75.
[doi:10.1037/1082-989X.4.1.70](https://doi.org/10.1037/1082-989X.4.1.70)

Aloe, A. M., & Becker, B. J. (2012). An effect size for regression
predictors in meta-analysis. *Journal of Educational and Behavioral
Statistics, 37*(2), 278–297.
[doi:10.3102/1076998610396901](https://doi.org/10.3102/1076998610396901)

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on contrasts, Chapter 5 on multiple
comparisons, and Chapter 9 on ANCOVA.)

Olkin, I., & Finn, J. D. (1995). Correlations redux. *Psychological
Bulletin, 118*(1), 155–164.
[doi:10.1037/0033-2909.118.1.155](https://doi.org/10.1037/0033-2909.118.1.155)

## See also

[`var_partial_r`](https://yelleknek.github.io/DMAR/reference/var_partial_r.md),
[`var_R2`](https://yelleknek.github.io/DMAR/reference/var_R2.md),
[`ci_reg_coef`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ss_aipe_semipartial_r`](https://yelleknek.github.io/DMAR/reference/ss_aipe_semipartial_r.md),
[`ss_power_R2_change`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_change.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# First-order variance of a semipartial r = .25 with N = 100, three
#     variables partialed out of X (p = 4 predictors in the full model),
#     and a full-model R^2 of .42.
var_semipartial_r(r_sp = 0.25, N = 100, p = 4, R2_full = 0.42,
                  under_null = FALSE)
#>  term              value
#>  var_semipartial_r 0.006

# By default, with R^2 of the full model supplied, the function returns
#     the squared standard error in the t test of X's unique
#     contribution.
var_semipartial_r(r_sp = 0.25, N = 100, p = 4, R2_full = 0.42)
#>  term                         value
#>  var_semipartial_r_under_null 0.006

# A reported pair, r_sp = .46 with a full-model R^2 of .21, in which
#     rounding has left R^2 below .46^2 = .2116: the value is returned
#     with a warning that says how the pair was read. The handler keeps
#     the message and lets the call finish, so both are shown.
noted <- new.env()
rounded <- withCallingHandlers(
  var_semipartial_r(r_sp = 0.46, N = 100, p = 2, R2_full = 0.21),
  dmar_semipartial_r_rounded = function(w) {
    noted$message <- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
rounded
#>  term                         value
#>  var_semipartial_r_under_null 0.008
noted$message
#> [1] "'R2_full' (0.21) is below 'r_sp' squared (0.2116) by no more than rounding to two decimals allows. The squared semipartial correlation cannot exceed the full-model R2, so one of the two was rounded. The value returned does not depend on r_sp and uses R2_full as given."
```
