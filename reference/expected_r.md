# Exact Expected Value of the Sample Pearson Correlation Given \\\rho\\ and \\n\\

Computes \\\mathrm{E}\[r \mid \rho, n\]\\, the exact expected value of
the sample Pearson product-moment correlation coefficient under
bivariate normality, using the Olkin-Pratt (1958) closed form built on
Hotelling's (1953) exact density. The sample \\r\\ is downward-biased as
an estimator of the population \\\rho\\, a fact known since Soper (1913)
and Fisher (1915); although this exact-bias equation has been available
for more than half a century, it has historically not been easy to
implement in general-purpose statistical software and is correspondingly
rarely used in applied work. This function makes the exact bias visible
so that users can decide whether to apply the Olkin-Pratt (1958)
unbiased estimator of \\\rho\\ (see Details).

## Usage

``` r
expected_r(rho, n)
```

## Arguments

- rho:

  Population correlation coefficient. A numeric scalar or vector in the
  interval \\\[-1, 1\]\\.

- n:

  Sample size on which the Pearson \\r\\ would be computed. A scalar or
  vector of integers with \\n \ge 4\\. (At \\n = 3\\ the bias is defined
  but the exact formula is numerically delicate; we require \\n \ge 4\\
  for well-conditioned computation.) When `rho` and `n` are both vectors
  they must have the same length or recycle cleanly.

## Value

A `data.frame` with one row per (`rho`, `n`) input and the columns

- `rho`, the input population correlation,

- `n`, the input sample size,

- `expected_r`, \\\mathrm{E}\[r \mid \rho, n\]\\,

- `bias`, \\\rho - \mathrm{E}\[r \mid \rho, n\]\\ (the amount by which
  the sample \\r\\ underestimates \\\rho\\ on average),

- `relative_bias`, `bias / rho` when \\\rho \ne 0\\, and `NA` when
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
Substituting `expected_r(rho, n)` for the bare \\\rho\\ in such planning
calculations corrects the leading-order over-promise. See **Examples**
for a design-stage walk-through.

**The exact formula.** Under bivariate normality, if \\r\\ is the sample
Pearson correlation in a sample of size \\n\\, \$\$\mathrm{E}\[r \mid
\rho, n\] \\=\\ \rho \\\cdot\\ {}\_2F_1\\\left(\tfrac{1}{2},\\
\tfrac{1}{2};\\ \tfrac{n+1}{2};\\ \rho^2 \right) \\\cdot\\
\frac{\Gamma\\\left(\tfrac{n}{2}\right)^2}
{\Gamma\\\left(\tfrac{n-1}{2}\right)\\\Gamma\\\left(\tfrac{n+1}{2}\right)}.\$\$
Here \\{}\_2F_1(a, b; c; z) = \sum\_{k=0}^{\infty} \frac{(a)\_k
(b)\_k}{(c)\_k\\ k!}\\ z^k\\ is the Gauss hypergeometric function with
Pochhammer symbols \\(x)\_k = x(x+1)\cdots(x+k-1)\\ (Hotelling, 1953,
Section 7 for the moments of \\r\\, Section 3, equation 25, for the
exact density; Olkin & Pratt, 1958, equation 3.2, for this closed form).
For \\\rho^2 \< 1\\ the series converges absolutely; this implementation
sums by a numerically stable forward recurrence and stops when the
relative contribution of the next term falls below a tolerance.

**Tuning the series convergence (rarely needed).** Two internal tuning
constants control the \\{}\_2F_1\\ series summation:

- `DMAR.expected_r.tol`, relative tolerance for stopping the series
  (default `1e-15`).

- `DMAR.expected_r.max_iter`, maximum number of series terms before
  issuing a non-convergence warning (default `5000L`).

Both have sensible defaults; advanced users who need different values
(e.g., for very near-boundary \\\|\rho\|\\ where the series converges
slowly) can set them via
[`options()`](https://rdrr.io/r/base/options.html), for example
`options(DMAR.expected_r.tol = 1e-12)`. They are deliberately hidden
from the function signature so as not to clutter the everyday user's
view of the call.

**Why the bias is present even though \\s^2\\ is unbiased for
\\\sigma^2\\.** The downward bias arises because \\r\\ is a nonlinear
function of unbiased sample moments. By Jensen's inequality and the
concavity of the square root in the denominator of \\r\\, the
expectation of the ratio is not the ratio of the expectations. The bias
is largest when \\n\\ is small or \\\|\rho\|\\ is moderate; as \\n \to
\infty\\, \\\mathrm{E}\[r\] \to \rho\\.

**Sign and magnitude.** The bias \\\rho - \mathrm{E}\[r\]\\ has the same
sign as \\\rho\\ and is approximately \\\rho(1 - \rho^2)/\[2(n - 1)\]\\
to leading order (Fisher, 1915; Hotelling, 1953; Ghosh, 1966); the exact
formula above incorporates all higher-order corrections. For \\\rho =
0.5\\, \\n = 10\\, the bias is about \\+0.021\\; for \\\rho = 0.5\\, \\n
= 30\\, about \\+0.0065\\.

**Olkin-Pratt (1958) unbiased estimator of \\\rho\\.** The companion to
this expected-value calculation is the Olkin-Pratt unbiased estimator of
\\\rho\\ given an observed \\r\\: \$\$\tilde{\rho}\_{\text{OP}}(r, n)
\\=\\ r \\\cdot\\ {}\_2F_1\\\left(\tfrac{1}{2},\\ \tfrac{1}{2};\\
\tfrac{n-2}{2};\\ 1 - r^2 \right).\$\$ Olkin & Pratt (1958) prove
\\\mathrm{E}\[\tilde{\rho}\_{\text{OP}}(r, n) \mid \rho, n\] = \rho\\
exactly, for every \\\rho \in (-1, 1)\\ and \\n \ge 4\\. The commonly
quoted first-order approximation \\\tilde{\rho} \approx r\\\[1 + (1 -
r^2)/(2(n - 3))\]\\ (Olkin, 1967) is the truncation of the OP series at
\\k = 1\\; the full series is what makes the estimator exactly unbiased.

Although the unbiased estimator is straightforward to apply, it is
rarely used because in most downstream uses (significance testing,
Fisher's \\Z\\ confidence intervals, structural-equation models) the
bias is small relative to other sources of uncertainty. Where it does
matter, meta-analyses with many small samples, reliability / validity
coefficients estimated from short calibration samples, design-stage
estimates feeding into AIPE sample size machinery the correction is well
worth applying.

**Connection to Fisher's Z transform.** The variance-stabilizing
transform \\z = \tanh^{-1}(r)\\, proposed in passing in Fisher (1915, p.
521) and developed in Fisher (1921), has approximate variance
\\1/(n-3)\\ regardless of \\\rho\\, but \\\mathrm{E}\[z\]\\ also carries
a small-sample bias of order \\1/n\\ (Hotelling, 1953, Section 8).
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
[`eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/eta_squared_partial.md),
[`expected_partial_r()`](https://yelleknek.github.io/DMAR/reference/expected_partial_r.md),
[`expected_smd()`](https://yelleknek.github.io/DMAR/reference/expected_smd.md),
[`nnt_from_smd()`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md),
[`omega_squared()`](https://yelleknek.github.io/DMAR/reference/omega_squared.md),
[`omega_squared_partial()`](https://yelleknek.github.io/DMAR/reference/omega_squared_partial.md),
[`probability_of_superiority_paired()`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md),
[`proportion_of_superiority()`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
[`responder_analysis()`](https://yelleknek.github.io/DMAR/reference/responder_analysis.md),
[`smd_trimmed()`](https://yelleknek.github.io/DMAR/reference/smd_trimmed.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. A single value: rho = 0.5, n = 10. The sample r is downwardly
#        biased by about 0.021, roughly 4% of rho.
expected_r(rho = 0.5, n = 10)
#>  rho n  expected_r bias   relative_bias
#>  0.5 10 0.479      0.0213 0.0427       

# 2. Bias as a function of n for fixed rho. The bias is roughly
#        rho * (1 - rho^2) / (2(n - 1)) to leading order; as n grows it
#        shrinks toward zero.
expected_r(rho = 0.5, n = c(5, 10, 20, 50, 100, 500))
#>  rho n   expected_r bias     relative_bias
#>  0.5 5   0.452      0.0483   0.0966       
#>  0.5 10  0.479      0.0213   0.0427       
#>  0.5 20  0.49       0.01     0.02         
#>  0.5 50  0.496      0.00385  0.0077       
#>  0.5 100 0.498      0.0019   0.0038       
#>  0.5 500 0.5        0.000376 0.000752     

# 3. Bias as a function of rho for fixed n. The bias is zero at
#        rho = 0 and rho = +/- 1, and largest near rho = +/- 0.6.
expected_r(rho = seq(0, 0.95, by = 0.05), n = 10)
#>  rho  n  expected_r bias    relative_bias
#>  0    10 0          0       <NA>         
#>  0.05 10 0.0473     0.00269 0.0538       
#>  0.1  10 0.0946     0.00535 0.0535       
#>  0.15 10 0.142      0.00794 0.053        
#>  0.2  10 0.19       0.0104  0.0522       
#>  0.25 10 0.237      0.0128  0.0512       
#>  0.3  10 0.285      0.015   0.05         
#>  0.35 10 0.333      0.017   0.0486       
#>  0.4  10 0.381      0.0187  0.0469       
#>  0.45 10 0.43       0.0202  0.0449       
#>  0.5  10 0.479      0.0213  0.0427       
#>  0.55 10 0.528      0.0221  0.0402       
#>  0.6  10 0.578      0.0224  0.0374       
#>  0.65 10 0.628      0.0223  0.0343       
#>  0.7  10 0.678      0.0215  0.0308       
#>  0.75 10 0.73       0.0202  0.0269       
#>  0.8  10 0.782      0.0181  0.0226       
#>  0.85 10 0.835      0.0152  0.0179       
#>  0.9  10 0.889      0.0113  0.0126       
#>  0.95 10 0.944      0.00636 0.00669      

# 4. The Olkin-Pratt unbiased estimator: invert the bias for an
#        observed sample r.
set.seed(113)
x <- rnorm(20); y <- 0.4 * x + sqrt(1 - 0.4^2) * rnorm(20)
r_obs <- cor(x, y)
r_obs
#> [1] 0.0439591

# Olkin-Pratt unbiased estimator of rho:
op_unbiased <- function(r, n, tol = 1e-15, max_iter = 5000) {
  z <- 1 - r^2; c_par <- (n - 2) / 2
  s <- 1; term <- 1
  for (k in seq_len(max_iter)) {
    term <- term * ((k - 0.5)^2) / ((c_par + k - 1) * k) * z
    s <- s + term
    if (abs(term) < tol * abs(s)) break
  }
  r * s
}
op_unbiased(r_obs, n = 20)
#> [1] 0.04535046

# 5. Design-stage use: a study planned around rho = 0.4 with n = 30.
#        Naive plug-in says we expect to observe r = 0.40 on average,
#        but the realized expected r is smaller, and that gap matters
#        for any precision- or power-based sample size calculation that
#        plugs in rho as if it were the expected sample r.
rho_planned <- 0.4
n_planned   <- 30
expected_r(rho = rho_planned, n = n_planned)
#>  rho n  expected_r bias    relative_bias
#>  0.4 30 0.394      0.00581 0.0145       

# 6. Tuning constants are hidden from the signature but tunable
#        through options() for the rare cases that need them (e.g.,
#        very near-boundary |rho| where the 2F1 series converges
#        slowly). The defaults rarely need to be changed.
options(DMAR.expected_r.tol = 1e-12,
        DMAR.expected_r.max_iter = 20000L)
expected_r(rho = 0.999, n = 5)
#>  rho   n expected_r bias     relative_bias
#>  0.999 5 0.999      0.000493 0.000494     
options(DMAR.expected_r.tol = NULL,
        DMAR.expected_r.max_iter = NULL)   # restore defaults
```
