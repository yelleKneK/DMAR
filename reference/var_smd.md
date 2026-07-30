# Variance of Cohen's *d* and Hedges' *g*

Computes the variance of the sample standardized mean difference
(Cohen's *d*) under bivariate normality and homogeneous variances, using
the *exact* noncentral *t* sampling distribution (Hedges, 1981) as the
default and reporting the large-sample Hedges-Olkin (1985) approximation
alongside for comparison. Optionally returns the variance of Hedges'
*g*, the bias-corrected counterpart of *d*.

## Usage

``` r
var_smd(delta, n_1, n_2 = NULL, unbiased = FALSE)
```

## Arguments

- delta:

  Population standardized mean difference. Numeric scalar or vector.

- n_1:

  Sample size in group 1. Scalar or vector.

- n_2:

  Sample size in group 2. Scalar or vector. Defaults to `n_1` (balanced
  design).

- unbiased:

  Logical. If `TRUE`, returns the variance of Hedges' *g* (the
  bias-corrected estimator); if `FALSE` (the default), returns the
  variance of Cohen's *d*.

## Value

A `data.frame` with rows for the exact (noncentral- *t*) variance and
the Hedges-Olkin large-sample approximation. Columns are `term`
(`"var_smd_exact"` or `"var_smd_approx"`) and `value`.

## Details

`var_smd()` is a stand-alone variance utility: most R packages return
only the Hedges-Olkin large-sample approximation and conflate the
variances of *d* and *g* (Goulet-Pelletier & Cousineau, 2018). The drift
between the exact and approximate forms becomes non-trivial below about
\\n = 30\\ per group and matters whenever `var_smd()` feeds into
meta-analytic weighting, AIPE planning, or a Wald-style standard-error
report.

**Exact noncentral *t* form.** For \\\hat d = (\bar Y_1 - \bar
Y_2)/s_p\\ with pooled \\s_p\\, the rescaled statistic \\t = \hat d
\sqrt{n_1 n_2 / (n_1 + n_2)}\\ follows a noncentral *t* with
\\\mathit{df} = n_1 + n_2 - 2\\ degrees of freedom and noncentrality
parameter \\\lambda = \delta \sqrt{n_1 n_2 / (n_1 + n_2)}\\. The
variance of a noncentral *t* is (Johnson, Kotz, & Balakrishnan, 1995,
Sec.\\ 31.3) \$\$\mathrm{Var}(t) \\=\\ \frac{\mathit{df}\\(1 +
\lambda^2)}{\mathit{df} - 2} \\-\\ \lambda^2 \\ c(\mathit{df})^{2},\$\$
where \\c(\mathit{df}) = \sqrt{\mathit{df}/2}\\
\Gamma((\mathit{df}-1)/2)\\/\\\Gamma(\mathit{df}/2)\\; dividing by the
design factor \\n_1 n_2 / (n_1 + n_2)\\ returns \\\mathrm{Var}(\hat
d)\\. For Hedges' *g*, multiply the result by \\J(\mathit{df})^2\\ where
\\J(\mathit{df}) = 1/c(\mathit{df})\\ is the Hedges-Olkin (1985)
bias-correction factor (see
[`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md)).

**Hedges-Olkin large-sample approximation.** The frequently quoted
approximation (Hedges & Olkin, 1985, equation 8) is
\$\$\mathrm{Var}(\hat d) \\\approx\\ \frac{n_1 + n_2}{n_1 n_2} \\+\\
\frac{\delta^2}{2(n_1 + n_2 - 2)}.\$\$ This approaches the exact form
only as the degrees of freedom grow: even at \\\delta = 0\\ it returns
\\1/(n_1 n_2 / (n_1 + n_2))\\ while the exact noncentral *t* variance is
\\\[\mathit{df}/(\mathit{df} - 2)\]/(n_1 n_2 / (n_1 + n_2))\\, so the
approximation is biased downward by a factor of \\(\mathit{df} -
2)/\mathit{df}\\, and the downward bias grows with \\\delta\\ and small
\\n\\. Goulet-Pelletier & Cousineau (2018) document the drift and
recommend the exact form for \\n \< 30\\ per group.

**Companions.** `var_smd()` is the variance partner of
[`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md)
(mean) and
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) (CI).
For design-stage AIPE planning that solves for \\n\\ given a target CI
width on *d*, see
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md).

## References

Goulet-Pelletier, J.-C., & Cousineau, D. (2018). A review of effect
sizes and their confidence intervals, Part I: The Cohen's *d* family.
*The Quantitative Methods for Psychology, 14*(4), 242–265.
[doi:10.20982/tqmp.14.4.p242](https://doi.org/10.20982/tqmp.14.4.p242)

Hedges, L. V. (1981). Distribution theory for Glass's estimator of
effect size and related estimators. *Journal of Educational Statistics,
6*(2), 107–128.

Hedges, L. V., & Olkin, I. (1985). *Statistical methods for
meta-analysis*. Academic Press.

Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). *Continuous
univariate distributions, volume 2* (2nd ed.). Wiley.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons and Chapter 3
on one-way ANOVA.)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`expected_smd`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`ss_aipe_smd`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)

Other variance utilities:
[`var_alpha()`](https://yelleknek.github.io/DMAR/reference/var_alpha.md),
[`var_cv()`](https://yelleknek.github.io/DMAR/reference/var_cv.md),
[`var_ete()`](https://yelleknek.github.io/DMAR/reference/var_ete.md),
[`var_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/var_indirect_effect.md),
[`var_omega_squared()`](https://yelleknek.github.io/DMAR/reference/var_omega_squared.md),
[`var_r()`](https://yelleknek.github.io/DMAR/reference/var_r.md),
[`var_smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/var_smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Balanced design, delta = 0.5, n = 20 per group.
var_smd(delta = 0.5, n_1 = 20)
#>  term           value
#>  var_smd_exact  0.109
#>  var_smd_approx 0.103

# 2. Hedges-Olkin approximation drifts from exact when n is small.
var_smd(delta = 0.5, n_1 = 5)
#>  term           value
#>  var_smd_exact  0.56 
#>  var_smd_approx 0.416
var_smd(delta = 0.5, n_1 = 50)
#>  term           value 
#>  var_smd_exact  0.0422
#>  var_smd_approx 0.0413

# 3. Variance of Hedges' g (bias-corrected).
var_smd(delta = 0.5, n_1 = 20, unbiased = TRUE)
#>  term           value 
#>  var_smd_exact  0.105 
#>  var_smd_approx 0.0992
```
