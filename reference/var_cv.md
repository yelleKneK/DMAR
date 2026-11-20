# Asymptotic Variance of the Coefficient of Variation

Computes the large-sample variance of the sample coefficient of
variation \\k = s / M\\, the sample standard deviation divided by the
sample mean, for a random sample from a normal population: \\\kappa^2
(1/2 + \kappa^2) / (N - 1)\\. It takes the population coefficient of
variation and the sample size and returns the variance. Companion to
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) and
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md).

## Usage

``` r
var_cv(true_cv, N)
```

## Arguments

- true_cv:

  Population coefficient of variation \\\kappa = \sigma / \mu\\. A
  single finite number greater than 0. For \\\kappa\\ above about 0.4
  the value returned is too small, and more so the larger \\\kappa\\ is
  (see Details); for an interval use
  [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md), which
  does not use this variance.

- N:

  Total sample size. A single finite number of at least 3.

## Value

A `data.frame` with the columns `term` and `value` and one row,
`var_cv_normal`: the first-order variance of the sample coefficient of
variation under normality (see Details).

## Details

**Estimand.** The variance, over repeated samples of \\N\\ independent
observations from a normal population with mean \\\mu \> 0\\ and
standard deviation \\\sigma\\, of the sample coefficient of variation
\\k = s / M\\, evaluated at the population value \\\kappa = \sigma /
\mu\\ given as `true_cv`. Strictly, the sample coefficient has no finite
variance, because the sample mean can be arbitrarily close to zero; the
expression below describes its spread when that event has negligible
probability, that is, when \\\sqrt{N} / \kappa\\ is large. Substituting
a sample coefficient for `true_cv` gives a plug-in estimate of the
variance, which carries the sampling error of the coefficient it was
computed from.

**The variance.** The first-order (delta method) variance is
\$\$\mathrm{Var}(k) \\\approx\\ \frac{\kappa^2}{N - 1}
\left(\frac{1}{2} + \kappa^2\right).\$\$ It follows from
\\\mathrm{Var}(s) \approx \sigma^2 / \\2 (N - 1)\\\\, \\\mathrm{Var}(M)
= \sigma^2 / N\\, and the independence of \\s\\ and \\M\\ under
normality, which give \\\mathrm{Var}(k) \approx \kappa^2 / \\2 (N -
1)\\ + \kappa^4 / N\\. That expression is Expression 5 of Sokal and
Braumann (1980, p. 52) without its third-order term. The function uses
\\N - 1\\ under both terms, which agrees with it to first order in
\\1/N\\ and is larger by \\\kappa^4 / \\N (N - 1)\\\\. The form usually
given, their Expression 7, has \\N\\ under both terms: \\\kappa^2 (1/2 +
\kappa^2) / N\\.

**Accuracy under normality.** The comparisons in this paragraph are with
the variance of \\k\\ when the sample mean is restricted to lie within 8
standard errors of \\\mu\\ and to be at least \\\mu / 50\\. The
restriction excludes a probability below \\10^{-9}\\ in every case
reported, and the figures given are the same for any lower bound from
\\\mu / 1000\\ to \\\mu / 5\\. For \\\kappa \le 0.3\\ and \\N \ge 10\\
the value is too large. The excess is largest for a small \\\kappa\\,
where it is that of \\\sigma^2 / \\2 (N - 1)\\\\ over the exact variance
of \\s\\: to one decimal place, 3.0 percent at \\N = 10\\, 1.4 percent
at \\N = 20\\, 0.9 percent at \\N = 30\\, and 0.5 percent at \\N = 50\\.
For \\\kappa\\ above about 0.4 the value is too small, because of the
terms of higher order that a first-order variance omits: by 5.7 percent
at \\\kappa = 0.5\\ and \\N = 10\\, and by 1.7 percent at \\\kappa =
0.5\\ and \\N = 30\\.

**Normality.** The value is for a normal population. For any population
with a finite fourth moment, the large-sample variance of \\k\\ is
\\\xi^2 / N\\, with \$\$\xi^2 \\=\\ \frac{\kappa^2}{2} + \kappa^4 -
\kappa^3 \gamma_1 + \frac{\kappa^2 \gamma_2}{4},\$\$ where \\\gamma_1\\
is the skewness and \\\gamma_2\\ the excess kurtosis of the population
(Kelley, Darku, & Chattopadhyay, 2018, Equation 30, which is written in
central moments). Under normality both are 0 and \\\xi^2 = \kappa^2
(1/2 + \kappa^2)\\. For other populations the normal-theory value can be
far off in either direction: \\\kappa^2 (1/2 + \kappa^2)\\ is 39 percent
larger than \\\xi^2\\ for a gamma population with \\\kappa = 0.8\\, and
56 percent smaller for a Laplace (double exponential) population with
\\\kappa = 0.3\\. These compare the limiting values, so a larger sample
does not remove the difference.
[`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.md) estimates
\\\xi^2\\ from the scores, without assuming a normal population.

**Use.** Neither
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) nor
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
uses this variance: the interval inverts the noncentral *t* distribution
of \\\sqrt{N} / k\\, and the planner (Kelley, 2007) chooses the sample
size for that interval. The value serves as a quick description of the
sampling variability of \\k\\ and for comparison with a simulated
variance. McKay (1932) and Vangel (1996) approximate the distribution of
a function of \\k\\ by a chi square distribution, from which confidence
limits follow; the variance here is not taken from those approximations.

## References

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Kelley, K., Darku, F. B., & Chattopadhyay, B. (2018). Accuracy in
parameter estimation for a general class of effect sizes: A sequential
approach. *Psychological Methods, 23*(2), 226–243.
[doi:10.1037/met0000127](https://doi.org/10.1037/met0000127)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

McKay, A. T. (1932). Distribution of the coefficient of variation and
the extended *t* distribution. *Journal of the Royal Statistical
Society, 95*(4), 695–698.

Sokal, R. R., & Braumann, C. A. (1980). Significance tests for
coefficients of variation and variability profiles. *Systematic Zoology,
29*(1), 50–66.
[doi:10.1093/sysbio/29.1.50](https://doi.org/10.1093/sysbio/29.1.50)

Vangel, M. G. (1996). Confidence intervals for a normal coefficient of
variation. *The American Statistician, 50*(1), 21–26.
[doi:10.1080/00031305.1996.10473537](https://doi.org/10.1080/00031305.1996.10473537)

## See also

[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md),
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md),
[`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md),
[`mr_cv`](https://yelleknek.github.io/DMAR/reference/mr_cv.md)

Other variance utilities:
[`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
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
# The variance of the sample coefficient of variation in samples of 30
# from a normal population whose coefficient is .20. The variance is
# small, so it is printed to six decimal places.
v <- var_cv(true_cv = 0.20, N = 30)
print(v, digits = 6)
#>  term          value   
#>  var_cv_normal 0.000745

# Its square root is the standard error of the sample coefficient, shown
# here beside the standard deviation of 2000 simulated coefficients.
set.seed(113)
k <- replicate(2000, {
  y <- rnorm(30, mean = 10, sd = 2)
  sd(y) / mean(y)
})
c(normal_theory = sqrt(v$value), simulated = sd(k))
#> normal_theory     simulated 
#>    0.02729153    0.02689039 

# The same comparison for a skewed population, a gamma distribution
# with a coefficient of .80 and samples of 200: the normal-theory
# standard error is too large.
k_gamma <- replicate(2000, {
  y <- rgamma(200, shape = 1 / 0.80^2)
  sd(y) / mean(y)
})
c(normal_theory = sqrt(var_cv(true_cv = 0.80, N = 200)$value),
  simulated = sd(k_gamma))
#> normal_theory     simulated 
#>    0.06055024    0.05082959 
```
