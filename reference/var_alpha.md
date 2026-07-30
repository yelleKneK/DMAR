# Asymptotic Variance of Coefficient Alpha (Cronbach, Guttman)

Computes the asymptotic variance of the sample coefficient alpha
(Guttman, 1945; Cronbach, 1951) under multivariate normality of the *p*
item scores, using the closed form derived by van Zyl, Neudecker, & Nel
(2000) under the assumption that the population covariance matrix has
equal off-diagonals (the parallel-items model), along with the simpler
Bonett (2002) approximation \\2p/(p-1) \cdot (1 - \alpha)^2 / (n - 1)\\
that is widely used in planning.

## Usage

``` r
var_alpha(alpha, n, p_items)
```

## Arguments

- alpha:

  Population coefficient alpha. Numeric scalar in \\\[0, 1)\\.

- n:

  Sample size (number of *respondents*).

- p_items:

  Number of items contributing to the composite alpha coefficient.

## Value

A `data.frame` with rows for the van Zyl, Neudecker, & Nel (2000)
exact-under-parallel-items variance and the Bonett (2002) simpler
approximation; columns are `term` and `value`.

## Details

Companion to the existing reliability-coefficient infrastructure
([`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md))
and a building block for AIPE planning around alpha.

**van Zyl-Neudecker-Nel (2000) variance.** Under multivariate normality
and the parallel-items model (all items have equal variances and equal
pairwise covariances), the asymptotic variance of the maximum likelihood
estimator of alpha is \$\$\mathrm{Var}(\hat\alpha) \\=\\ \frac{2 p (1 -
\alpha)^2}{(p - 1)(n - 2)}.\$\$ This is one of two closed forms van Zyl
et al. derive; the more general (non-parallel) form involves matrix
expressions and is implemented separately by the existing reliability
infrastructure.

**Bonett (2002) approximation.** Bonett (2002) gives the easy planning
form \$\$\mathrm{Var}(\hat\alpha) \\\approx\\ \frac{2 p}{(p - 1)} \cdot
\frac{(1 - \alpha)^2}{n - 1}.\$\$ This differs from the van Zyl form
only in the denominator (\\n-1\\ vs.\\ \\n-2\\) and converges to the
same value for moderate \\n\\. Bonett's version is what most sample size
tables use.

**When to use which.** For inference (a CI on \\\alpha\\), the van Zyl
form is preferable, especially at small \\n\\; for sample size
*planning* the difference is immaterial and the Bonett form is widely
cited and easier to invert.

## References

Bonett, D. G. (2002). Sample size requirements for testing and
estimating coefficient alpha. *Journal of Educational and Behavioral
Statistics, 27*(4), 335–340.
[doi:10.3102/10769986027004335](https://doi.org/10.3102/10769986027004335)

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
# 1. Variance of alpha = 0.80 from a 10-item test with n = 100.
var_alpha(alpha = 0.80, n = 100, p_items = 10)
#>  term              value   
#>  var_alpha_van_zyl 0.000907
#>  var_alpha_bonett  0.000898

# 2. Variance shrinks with n and grows as alpha moves away from 1:
var_alpha(alpha = 0.90, n = 50,  p_items = 5)
#>  term              value   
#>  var_alpha_van_zyl 0.000521
#>  var_alpha_bonett  0.00051 
var_alpha(alpha = 0.90, n = 500, p_items = 5)
#>  term              value   
#>  var_alpha_van_zyl 5.02e-05
#>  var_alpha_bonett  5.01e-05
```
