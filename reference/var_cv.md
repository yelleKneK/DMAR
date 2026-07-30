# Asymptotic Variance of the Coefficient of Variation

Computes the asymptotic variance of the sample coefficient of variation
\\\hat\kappa = s / \bar Y\\ under normality, using McKay's (1932)
original noncentral *t*-based approximation and Vangel's (1996)
refinement. Companion to
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md) and
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md).

## Usage

``` r
var_cv(cv, n)
```

## Arguments

- cv:

  Population coefficient of variation \\\kappa = \sigma / \mu\\. Numeric
  scalar in \\(0, \infty)\\. When \\\kappa\\ is very large (\> 0.5),
  both McKay's and Vangel's approximations degrade and exact methods
  ([`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md))
  should be preferred.

- n:

  Sample size.

## Value

A `data.frame` with rows for the McKay (1932) and Vangel (1996)
approximations; columns are `term` and `value`.

## Details

**McKay (1932).** The classical large-sample variance of the sample CV
under normality is \$\$\mathrm{Var}(\hat\kappa) \\\approx\\
\frac{\kappa^2}{n - 1} \cdot \left(\frac{1}{2} + \kappa^2\right).\$\$
This is exact up to \\O(1/n)\\ and is what most planning tables use. It
begins to drift when \\\kappa \> 0.3\\ or so.

**Vangel (1996).** Vangel showed that a small-sample correction that
adjusts the McKay form for the noncentral *t* mean factor gives
substantially better coverage of CIs derived from the variance:
\$\$\mathrm{Var}\_{\mathrm{Vangel}}(\hat\kappa) \\\approx\\
\frac{\kappa^2}{n - 1} \cdot \left(\frac{1}{2} + \kappa^2 \cdot
\frac{n + 1}{n - 1}\right).\$\$ The two forms coincide in the
large-\\n\\ limit. We report both so the user can see the magnitude of
the small-sample correction.

## References

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3.)

McKay, A. T. (1932). Distribution of the coefficient of variation and
the extended *t* distribution. *Journal of the Royal Statistical
Society, 95*(4), 695–698.

Vangel, M. G. (1996). Confidence intervals for a normal coefficient of
variation. *The American Statistician, 50*(1), 21–26.
[doi:10.1080/00031305.1996.10473537](https://doi.org/10.1080/00031305.1996.10473537)

## See also

[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md),
[`ss_aipe_cv_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv_sensitivity.md)

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
# 1. CV = 0.20 in a sample of 30:
var_cv(cv = 0.20, n = 30)
#>  term          value   
#>  var_cv_mckay  0.000745
#>  var_cv_vangel 0.000749

# 2. The Vangel correction grows with cv (becomes non-trivial
#        for kappa > 0.3):
var_cv(cv = 0.50, n = 30)
#>  term          value  
#>  var_cv_mckay  0.00647
#>  var_cv_vangel 0.00661
```
