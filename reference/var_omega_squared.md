# Asymptotic Variance of Omega Squared (ANOVA Effect Size)

Computes the large-sample (delta method) variance of the sample
\\\hat\omega^2\\ (Hays' 1994 bias-corrected estimator) in a
fixed-effects ANOVA. Fleishman (1980, Eq. 22, p. 669) gives the exact
variance of the unbiased estimator of the signal-to-noise ratio \\f^2 =
\sigma^2_a / \sigma^2_e\\ under the noncentral *F* sampling distribution
of the observed *F* statistic; because \\\omega^2 = f^2 / (1 + f^2)\\
(his Eq. 8), the delta method carries that variance to the \\\omega^2\\
scale with the Jacobian \\\mathrm{d}\omega^2/\mathrm{d}f^2 = (1 -
\omega^2)^2\\. Fleishman gives no variance on the \\\omega^2\\ scale
himself, so the transfer is this package's step rather than his. The
result is the natural companion to
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
(CI) and
[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md)
(point estimate).

## Usage

``` r
var_omega_squared(
  population_omega_squared = NULL,
  df_effect = NULL,
  df_error = NULL,
  N = NULL,
  object = NULL
)
```

## Arguments

- population_omega_squared:

  Population \\\omega^2\\. Numeric scalar in \\\[0, 1)\\. Ignored when
  `object` is supplied.

- df_effect:

  Numerator degrees of freedom for the effect. Ignored when `object` is
  supplied.

- df_error:

  Error (residual) degrees of freedom. Ignored when `object` is
  supplied.

- N:

  Total sample size. Ignored when `object` is supplied.

- object:

  Optional fitted [`aov`](https://rdrr.io/r/stats/aov.html) or
  [`lm`](https://rdrr.io/r/stats/lm.html) object. When supplied, the
  function loops over the non-`Residuals` rows of
  [`anova`](https://rdrr.io/r/stats/anova.html)`(object)` and returns
  one row per effect, plugging in the sample \\\hat\omega^2_p\\ for each
  as the working population value.

## Value

A 1-row `data.frame` with columns `term` (`"var_omega_squared"`) and
`value` (the asymptotic variance).

## Details

**Derivation.** In a fixed-effects ANOVA with numerator df \\df_1\\ and
denominator df \\df_2\\, the observed *F* statistic follows a noncentral
*F* distribution with noncentrality \\\lambda = df_1 (F - 1)\\ when
\\\hat\omega^2 = df_1 (F - 1) / \[df_1 (F - 1) + N\]\\ is the population
value (Hays, 1994). The asymptotic variance of \\\hat\omega^2\\ is
obtained by the delta method on this relationship (Fleishman, 1980),
yielding: \$\$\mathrm{Var}(\hat\omega^2) \\\approx\\ \frac{2 \cdot df_1
\cdot (df_2 - 2) (1 - \omega^2)^2 (1 + \lambda^\*/df_1)^2} {N^2 (df_2 -
4)},\$\$ with \\\lambda^\* = \omega^2 N / (1 - \omega^2)\\ the
noncentrality implied by the population value. This is the form used by
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
when constructing a Wald-style interval; the noncentral *F* CI is
generally preferred.

**Caveats.** The variance is a delta method approximation; it becomes
inaccurate when \\df_2\\ is small (\\\< 10\\), when \\\omega^2\\ is near
the boundaries 0 or 1, or when the residual distribution is
heavy-tailed. For small-sample inference, the noncentral *F* CI
([`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md))
is preferred over a Wald-style interval built on this variance.

## References

Fleishman, A. I. (1980). Confidence intervals for correlation ratios.
*Educational and Psychological Measurement, 40*(3), 659–670.

Hays, W. L. (1994). *Statistics* (5th ed.). Fort Worth, TX: Harcourt
Brace College Publishers.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*, 137–152.
[doi:10.1037/a0028086](https://doi.org/10.1037/a0028086)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\\eta^2\\, Chapter 7 on factorial
designs, and Chapter 11 on generalized \\\eta^2\\ for within-subjects
designs.)

Olejnik, S., & Algina, J. (2003). Generalized eta and omega squared
statistics: Measures of effect size for some common research designs.
*Psychological Methods, 8*(4), 434–447.
[doi:10.1037/1082-989X.8.4.434](https://doi.org/10.1037/1082-989X.8.4.434)

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`omega_squared`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ss_aipe_omega_squared`](https://yelleknek.github.io/DMAR/reference/ss_aipe_omega_squared.md)

Other variance utilities:
[`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
[`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`var_ete()`](https://yelleknek.github.io/DMAR/reference/var_ete.md),
[`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md),
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
[`var_smd()`](https://yelleknek.github.io/DMAR/reference/var_smd.md),
[`var_smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. One way ANOVA: 3 groups (df_effect = 2), 60 total (df_error = 57),
#        population omega^2 = 0.10.
var_omega_squared(population_omega_squared = 0.10,
                   df_effect = 2,
                   df_error = 57,
                   N = 60)
#>  term              value  
#>  var_omega_squared 0.00632

# 2. Per-effect variance from a fitted lm() / aov() (pygmalion data:
#        expectancy treatment x grade, N = 310):
fit_factorial <- aov(iq_8 ~ treatment * factor(grade), data = pygmalion)
var_omega_squared(object = fit_factorial)
#>  effect                  omega_squared_partial var_omega_squared df_effect
#>  treatment               0.0176                0.000239          1        
#>  factor(grade)           0.0281                0.000441          5        
#>  treatment:factor(grade) 0.00307               0.000145          5        
#>  df_error N  
#>  298      310
#>  298      310
#>  298      310
```
