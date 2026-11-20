# Exact Expected Value of the Sample Pearson Correlation Given \\\rho\\ and \\N\\

Computes \\\mathrm{E}\[r \mid \rho, N\]\\, the exact expected value of
the sample Pearson product-moment correlation coefficient under
bivariate normality, using the Olkin-Pratt (1958) closed form built on
Hotelling's (1953) exact density. The sample \\r\\ is downward-biased as
an estimator of the population \\\rho\\, a fact known since Soper (1913)
and Fisher (1915). The exact expected value requires summing a Gauss
hypergeometric series (see Details), which this function does, so that
users can see the exact bias and decide whether to apply the Olkin-Pratt
(1958) unbiased estimator of \\\rho\\ (see Details).

## Usage

``` r
expected_r(true_cor, N)
```

## Arguments

- true_cor:

  Population correlation coefficient \\\rho\\. A numeric scalar or
  vector in the interval \\\[-1, 1\]\\.

- N:

  Sample size on which the Pearson \\r\\ would be computed. A scalar or
  vector of integers with \\N \ge 4\\. (At \\N = 3\\ the bias is defined
  but the exact formula is numerically delicate; we require \\N \ge 4\\
  for well-conditioned computation.) When `true_cor` and `N` are both
  vectors they must have the same length or recycle cleanly.

## Value

A `data.frame` with one row per (`true_cor`, `N`) input and the columns

- `true_cor`, the input population correlation,

- `N`, the input sample size,

- `expected_r`, \\\mathrm{E}\[r \mid \rho, N\]\\,

- `bias`, \\\rho - \mathrm{E}\[r \mid \rho, N\]\\ (the amount by which
  the sample \\r\\ underestimates \\\rho\\ on average),

- `relative_bias`, `bias / true_cor` when \\\rho \ne 0\\, and `NA` when
  \\\rho = 0\\.

## Details

`expected_r()` is especially useful at the *design stage* of a study,
where sample size planning typically proceeds from an assumed value of
the population correlation \\\rho\\. Because the value of \\r\\ a
researcher should expect to observe on average is *smaller in absolute
value* than the assumed \\\rho\\, plugging \\\rho\\ directly into
sampling-distribution machinery (precision of \\r\\, width of a
confidence interval, power of a test of \\H_0\\: \rho = 0\\)
systematically over-promises on the realized precision or power.
Substituting `expected_r(true_cor, N)` for the bare \\\rho\\ in such
planning calculations corrects the leading-order over-promise. See
**Examples** for a design-stage walk-through.

**The exact formula.** Under bivariate normality, if \\r\\ is the sample
Pearson correlation in a sample of size \\N\\, \$\$\mathrm{E}\[r \mid
\rho, N\] \\=\\ \rho \\\cdot\\ {}\_2F_1\\\left(\tfrac{1}{2},\\
\tfrac{1}{2};\\ \tfrac{N+1}{2};\\ \rho^2 \right) \\\cdot\\
\frac{\Gamma\\\left(\tfrac{N}{2}\right)^2}
{\Gamma\\\left(\tfrac{N-1}{2}\right)\\\Gamma\\\left(\tfrac{N+1}{2}\right)}.\$\$
Here \\{}\_2F_1(a, b; c; z) = \sum\_{k=0}^{\infty} \frac{(a)\_k
(b)\_k}{(c)\_k\\ k!}\\ z^k\\ is the Gauss hypergeometric function with
Pochhammer symbols \\(x)\_k = x(x+1)\cdots(x+k-1)\\ (Hotelling, 1953,
Section 7 for the moments of \\r\\, Section 3, equation 25, for the
exact density; Olkin & Pratt, 1958, equation 3.2, for this closed form).
For \\\rho^2 \< 1\\ the series converges absolutely; this implementation
sums by a numerically stable forward recurrence and stops when the
relative contribution of the next term falls below a tolerance. The
`bias` column is computed directly, as \\\rho\\ times one minus the
product of the series and the gamma ratio, formed from their logarithms,
and not as the difference \\\rho - \mathrm{E}\[r\]\\. At a large \\N\\
the two numbers share their leading digits (about eight of them at \\N =
10^8\\), and the subtraction would leave the bias only the digits after
those.

**Tuning the series convergence (rarely needed).** Two internal tuning
constants control the \\{}\_2F_1\\ series summation:

- `dmar.expected_r.tol`, relative tolerance for stopping the series
  (default `1e-15`). The summation stops once a term is below this
  fraction of the sum of the terms after the leading 1, the only part of
  the series on which the bias depends.

- `dmar.expected_r.max_iter`, maximum number of series terms (default
  `5000L`). If the last term summed is still not below
  `dmar.expected_r.tol` times the whole series, a non-convergence
  warning is issued.

Users who need different values (e.g., for very near-boundary
\\\|\rho\|\\ where the series converges slowly) can set them via
[`options()`](https://rdrr.io/r/base/options.html), for example
`options(dmar.expected_r.tol = 1e-12)`. They are deliberately hidden
from the function signature so as not to clutter the everyday user's
view of the call.

**Why the bias is present even though \\s^2\\ is unbiased for
\\\sigma^2\\.** The downward bias arises because \\r\\ is a nonlinear
function of unbiased sample moments. By Jensen's inequality and the
concavity of the square root in the denominator of \\r\\, the
expectation of the ratio is not the ratio of the expectations. The bias
is largest when \\N\\ is small or \\\|\rho\|\\ is moderate; as \\N \to
\infty\\, \\\mathrm{E}\[r\] \to \rho\\.

**Sign and magnitude.** The bias \\\rho - \mathrm{E}\[r\]\\ has the same
sign as \\\rho\\ and is approximately \\\rho(1 - \rho^2)/\[2(N - 1)\]\\
to leading order (Fisher, 1915; Hotelling, 1953; Ghosh, 1966); the exact
formula above incorporates all higher-order corrections. For \\\rho =
0.5\\, \\N = 10\\, the bias is about \\+0.021\\; for \\\rho = 0.5\\, \\N
= 30\\, about \\+0.0065\\.

**Olkin-Pratt (1958) unbiased estimator of \\\rho\\.** The companion to
this expected-value calculation is the Olkin-Pratt unbiased estimator of
\\\rho\\ given an observed \\r\\: \$\$\tilde{\rho}\_{\text{OP}}(r, N)
\\=\\ r \\\cdot\\ {}\_2F_1\\\left(\tfrac{1}{2},\\ \tfrac{1}{2};\\
\tfrac{N-2}{2};\\ 1 - r^2 \right).\$\$ Olkin & Pratt (1958) prove
\\\mathrm{E}\[\tilde{\rho}\_{\text{OP}}(r, N) \mid \rho, N\] = \rho\\
exactly, for every \\\rho \in (-1, 1)\\ and \\N \ge 4\\. The commonly
quoted first-order approximation \\\tilde{\rho} \approx r\\\[1 + (1 -
r^2)/(2(N - 3))\]\\ (Olkin, 1967) is the truncation of the OP series at
\\k = 1\\; the full series is what makes the estimator exactly unbiased.

**Size of the bias relative to sampling variability.** The standard
deviation of \\r\\ is approximately \\(1 - \rho^2)/\sqrt{N - 1}\\, so
the bias of a single \\r\\ is about \\\rho/(2\sqrt{N - 1})\\ of its
standard deviation, whether the coefficient is tested, given a
confidence interval, or used as a reliability, validity, or planning
value; at \\\rho = 0.8\\ and \\N = 20\\ the exact ratio is 0.09. When
\\K\\ correlations from samples of the same \\N\\ and the same \\\rho\\
are averaged, as in a meta-analysis, the bias of the average equals the
bias of a single \\r\\, while the standard error of the average falls as
\\1/\sqrt{K}\\; the bias is then about \\\rho\sqrt{K}/(2\sqrt{N - 1})\\
standard errors of the average, and at \\\rho = 0.8\\, \\N = 20\\, and
\\K = 50\\ the exact ratio is 0.63.

**Connection to Fisher's Z transform.** The variance-stabilizing
transform \\z = \tanh^{-1}(r)\\, proposed in passing in Fisher (1915, p.
521) and developed in Fisher (1921), has approximate variance
\\1/(N-3)\\ regardless of \\\rho\\, but \\\mathrm{E}\[z\]\\ also carries
a small-sample bias of order \\1/N\\ (Hotelling, 1953, Section 8).
Hotelling (1953, Sections 9–10) gives bias-adjusted and
variance-stabilized refinements of \\Z\\, e.g. \\z - (3z + r)/(4n)\\.

## References

Anderson, T. W. (2003). *An introduction to multivariate statistical
analysis* (3rd ed.), Section 4.2. Wiley.

Fisher, R. A. (1915). Frequency distribution of the values of the
correlation coefficient in samples from an indefinitely large
population. *Biometrika, 10*(4), 507–521.

Fisher, R. A. (1921). On the "probable error" of a coefficient of
correlation deduced from a small sample. *Metron, 1*, 3–32.

Ghosh, B. K. (1966). Asymptotic expansions for the moments of the
distribution of correlation coefficient. *Biometrika, 53*(1/2), 258–262.

Hotelling, H. (1953). New light on the correlation coefficient and its
transforms. *Journal of the Royal Statistical Society, Series B, 15*(2),
193–232. (Discussion, pp.\\ 225–232.)

Olkin, I. (1967). Correlations revisited. In J. C. Stanley (Ed.),
*Improving experimental design and statistical analysis* (pp.\\
102–128). Rand McNally.

Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain
correlation coefficients. *The Annals of Mathematical Statistics,
29*(1), 201–211.

Soper, H. E. (1913). On the probable error of the correlation
coefficient to a second approximation. *Biometrika, 9*(1/2), 91–115.

Soper, H. E., Young, A. W., Cave, B. M., Lee, A., & Pearson, K. (1917).
On the distribution of the correlation coefficient in small samples.
Appendix II to the papers of "Student" and R. A. Fisher: A cooperative
study. *Biometrika, 11*(4), 328–413.

Stuart, A., & Ord, J. K. (1994). *Kendall's advanced theory of
statistics, Vol.\\ 1: Distribution theory* (6th ed.), Section 16.32.
Edward Arnold.

## See also

[`ci_r`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`cor`](https://rdrr.io/r/stats/cor.html),
[`cor.test`](https://rdrr.io/r/stats/cor.test.html)

Other effect size estimates:
[`cles()`](https://yelleknek.github.io/DMAR/reference/cles.md),
[`cliff_delta()`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
[`correction_for_attenuation()`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
[`eta_squared()`](https://yelleknek.github.io/DMAR/reference/eta_squared.md),
[`eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/eta_squared_generalized.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. A single value: true_cor = 0.5, N = 10. The sample r is downwardly
#        biased by about 0.021, roughly 4% of true_cor.
expected_r(true_cor = 0.5, N = 10)
#>  true_cor N  expected_r bias  relative_bias
#>  0.500    10 0.479      0.021 0.043        

# 2. Bias as a function of N for fixed true_cor. The bias is roughly
#        true_cor * (1 - true_cor^2) / (2(N - 1)) to leading order; as N
#        grows it shrinks toward zero.
expected_r(true_cor = 0.5, N = c(5, 10, 20, 50, 100, 500))
#>  true_cor N   expected_r bias  relative_bias
#>  0.500    5   0.452      0.048 0.097        
#>  0.500    10  0.479      0.021 0.043        
#>  0.500    20  0.490      0.010 0.020        
#>  0.500    50  0.496      0.004 0.008        
#>  0.500    100 0.498      0.002 0.004        
#>  0.500    500 0.500      0.000 0.001        

# 3. Bias as a function of true_cor for fixed N. The bias is zero at
#        true_cor = 0 and true_cor = +/- 1, and largest near
#        true_cor = +/- 0.6.
expected_r(true_cor = seq(0, 0.95, by = 0.05), N = 10)
#>  true_cor N  expected_r bias  relative_bias
#>  0.000    10 0.000      0.000 <NA>         
#>  0.050    10 0.047      0.003 0.054        
#>  0.100    10 0.095      0.005 0.054        
#>  0.150    10 0.142      0.008 0.053        
#>  0.200    10 0.190      0.010 0.052        
#>  0.250    10 0.237      0.013 0.051        
#>  0.300    10 0.285      0.015 0.050        
#>  0.350    10 0.333      0.017 0.049        
#>  0.400    10 0.381      0.019 0.047        
#>  0.450    10 0.430      0.020 0.045        
#>  0.500    10 0.479      0.021 0.043        
#>  0.550    10 0.528      0.022 0.040        
#>  0.600    10 0.578      0.022 0.037        
#>  0.650    10 0.628      0.022 0.034        
#>  0.700    10 0.678      0.022 0.031        
#>  0.750    10 0.730      0.020 0.027        
#>  0.800    10 0.782      0.018 0.023        
#>  0.850    10 0.835      0.015 0.018        
#>  0.900    10 0.889      0.011 0.013        
#>  0.950    10 0.944      0.006 0.007        

# 4. The Olkin-Pratt unbiased estimator: invert the bias for an
#        observed sample r.
set.seed(113)
x <- rnorm(20); y <- 0.4 * x + sqrt(1 - 0.4^2) * rnorm(20)
r_obs <- cor(x, y)
r_obs
#> [1] 0.0439591

# Olkin-Pratt unbiased estimator of rho:
op_unbiased <- function(r, N, tol = 1e-15, max_iter = 5000) {
  z <- 1 - r^2; c_par <- (N - 2) / 2
  s <- 1; term <- 1
  for (k in seq_len(max_iter)) {
    term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
    s <- s + term
    if (abs(term) < tol * abs(s)) break
  }
  r * s
}
op_unbiased(r_obs, N = 20)
#> [1] 0.04535046

# 5. Design-stage use: a study planned around true_cor = 0.4 with N = 30.
#        Naive plug-in says we expect to observe r = 0.40 on average,
#        but the realized expected r is smaller, and that gap matters
#        for any precision- or power-based sample size calculation that
#        plugs in true_cor as if it were the expected sample r.
cor_planned <- 0.4
n_planned   <- 30
expected_r(true_cor = cor_planned, N = n_planned)
#>  true_cor N  expected_r bias  relative_bias
#>  0.400    30 0.394      0.006 0.015        

# 6. Tuning constants are hidden from the signature but tunable
#        through options() for the rare cases that need them (e.g.,
#        very near-boundary |true_cor| where the 2F1 series converges
#        slowly). The defaults rarely need to be changed.
options(dmar.expected_r.tol = 1e-12,
        dmar.expected_r.max_iter = 20000L)
expected_r(true_cor = 0.999, N = 5)
#>  true_cor N expected_r bias  relative_bias
#>  0.999    5 0.999      0.000 0.000        
options(dmar.expected_r.tol = NULL,
        dmar.expected_r.max_iter = NULL)   # restore defaults
```
