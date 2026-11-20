# Asymptotic Variance of Coefficient Alpha (Cronbach, Guttman)

Computes two large-sample approximations to the variance of the sample
coefficient alpha (Guttman, 1945; Cronbach, 1951) under multivariate
normality of the *J* item scores and compound symmetry of their
covariance matrix (equal item variances and equal covariances, the
parallel-items model): the asymptotic variance of van Zyl, Neudecker,
and Nel (2000), \\2J(1 - \alpha)^2 / ((J - 1)(N - 1))\\, and the
variance that follows from Bonett (2002), \\2J(1 - \alpha)^2 / ((J -
1)(N - 2))\\.

## Usage

``` r
var_alpha(alpha, N, n_items)
```

## Arguments

- alpha:

  Population coefficient alpha. Numeric scalar in \\\[0, 1)\\.

- N:

  Total sample size (number of *respondents*).

- n_items:

  Number of items contributing to the composite alpha coefficient (\\J\\
  in the formulas); at least 2.

## Value

A `data.frame` with columns `term` and `value` and two rows:
`var_alpha_van_zyl`, the asymptotic variance of van Zyl, Neudecker, and
Nel (2000), and `var_alpha_bonett`, the variance that follows from
Bonett (2002).

## Details

For an interval estimate of alpha from data see
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
and
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md);
for sample size planning see
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md).

**van Zyl, Neudecker, and Nel (2000).** Under multivariate normality and
compound symmetry, van Zyl et al. show that
\\n\\\mathrm{Var}(\hat\alpha)\\ tends to \\2J(1 - \alpha)^2 / (J - 1)\\
as \\n\\ grows, where \\n = N - 1\\ is the degrees of freedom of the
item covariance matrix (their Equation 13). The asymptotic variance is
therefore \$\$\mathrm{Var}(\hat\alpha) \\\approx\\ \frac{2 J (1 -
\alpha)^2}{(J - 1)(N - 1)}.\$\$ Without compound symmetry the asymptotic
variance is a function of the whole item covariance matrix (their
Equation 21);
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
uses that general form for its `"ml"` interval, and at a compound
symmetric covariance matrix it equals the form above.

**Bonett (2002).** Bonett (2003, Equations 2 and 4), citing Bonett
(2002), uses \\2J / ((J - 1)(N - 2))\\ as the large-sample variance of
\\\ln(1 - \hat\alpha)\\ under the same assumptions. Carried to the alpha
scale by the delta method it gives \$\$\mathrm{Var}(\hat\alpha)
\\\approx\\ \frac{2 J (1 - \alpha)^2}{(J - 1)(N - 2)}.\$\$

**Comparing the two.** The two forms differ only in the denominator: the
form with \\N - 2\\ is the form with \\N - 1\\ times \\(N - 1)/(N -
2)\\, and a sample size chosen to reach a target variance is one larger
with the \\N - 2\\ denominator than with \\N - 1\\. Both are
large-sample approximations. Under multivariate normality and compound
symmetry the exact variance of \\\hat\alpha\\ for \\N \> 5\\ is
\$\$\frac{2 (1 - \alpha)^2 n (nJ - 2)}{(J - 1)(n - 2)^2 (n - 4)}, \quad
n = N - 1\$\$ (van Zyl, Neudecker, & Nel, 2000, Equation 12). It is
larger than both forms, and the form with \\N - 2\\ is the closer of the
two: with \\J = 5\\ the exact variance is 2.53 times that form at \\N =
10\\, 1.27 times at \\N = 30\\, and 1.07 times at \\N = 100\\.

## References

Bonett, D. G. (2002). Sample size requirements for testing and
estimating coefficient alpha. *Journal of Educational and Behavioral
Statistics, 27*(4), 335–340.
[doi:10.3102/10769986027004335](https://doi.org/10.3102/10769986027004335)

Bonett, D. G. (2003). Sample size requirements for comparing two alpha
coefficients. *Applied Psychological Measurement, 27*(1), 72–74.
[doi:10.1177/0146621602239477](https://doi.org/10.1177/0146621602239477)

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika, 16*(3), 297–334.

Guttman, L. (1945). A basis for analyzing test-retest reliability.
*Psychometrika, 10*(4), 255–282.

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

McDonald, R. P. (1999). *Test theory: A unified treatment*. Lawrence
Erlbaum.

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

[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)

Other variance utilities:
[`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`var_ete()`](https://yelleknek.github.io/DMAR/reference/var_ete.md),
[`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md),
[`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
[`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
[`var_smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Variance of alpha = 0.80 from a 10-item test with N = 100.
var_alpha(alpha = 0.80, N = 100, n_items = 10)
#>  term              value
#>  var_alpha_van_zyl 0.001
#>  var_alpha_bonett  0.001

# 2. Variance shrinks with N and grows as alpha moves away from 1:
var_alpha(alpha = 0.90, N = 50,  n_items = 5)
#>  term              value
#>  var_alpha_van_zyl 0.001
#>  var_alpha_bonett  0.001
var_alpha(alpha = 0.90, N = 500, n_items = 5)
#>  term              value
#>  var_alpha_van_zyl 0.000
#>  var_alpha_bonett  0.000
```
