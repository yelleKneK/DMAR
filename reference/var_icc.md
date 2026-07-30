# Asymptotic Variance of the Intraclass Correlation Coefficient

Computes the asymptotic (large-sample) variance of the intraclass
correlation coefficient (ICC) for any of the six classical Shrout-Fleiss
(1979) forms, given a value of the population ICC, the number of
subjects \\n\\, and the number of raters \\k\\. The single-rater forms
use Smith's (1956) one-way-ANOVA asymptotic variance, equivalently the
Fisher (1925) information-matrix result, and the average-of-\\k\\ forms
use a delta method transformation through the Spearman-Brown relation.

## Usage

``` r
var_icc(rho, n, k, type = "ICC(1,1)")
```

## Arguments

- rho:

  The intraclass correlation coefficient at which the asymptotic
  variance is evaluated, on the scale matching `type`: for the
  single-rater types (`"ICC(1,1)"`, `"ICC(2,1)"`, `"ICC(3,1)"`) supply
  the single-rater ICC; for the average-of-\\k\\ types (`"ICC(1,k)"`,
  `"ICC(2,k)"`, `"ICC(3,k)"`) supply the average-of-\\k\\ ICC. Must lie
  in \\\[0, 1\]\\.

  **Population value or sample value?** The formula is derived in terms
  of the unknown population value \\\rho\\. In applied work \\\rho\\ is
  never known, so two conventions are common: (i) for *prospective*
  variance calculation (e.g., when planning a study and asking how
  precise an estimator will be at plausible truths), supply your
  *anticipated* population value motivated by prior literature, theory,
  or a pilot; (ii) for *post hoc* variance calculation (e.g.,
  constructing a Wald standard error or a meta-analytic weight from an
  observed sample), supply the *sample* estimate \\\hat\rho\\ as a
  plug-in for the population value. The two postures look identical at
  the call site but have different conceptual content; the plug-in
  case (ii) yields a consistent (not exact) estimator of the variance.

- n:

  Number of subjects (targets) rated.

- k:

  Number of raters (or repeated measurements) per subject.

- type:

  Which Shrout-Fleiss ICC variant the population value represents. One
  of `"ICC(1,1)"`, `"ICC(2,1)"`, `"ICC(3,1)"`, `"ICC(1,k)"`,
  `"ICC(2,k)"`, `"ICC(3,k)"`. The shorthand aliases used in
  [`icc`](https://yelleknek.github.io/DMAR/reference/icc.md) (`"1"`,
  `"2"`, `"3"`, `"1k"`, `"2k"`, `"3k"`) are also accepted.

## Value

A one-row `data.frame` with columns `term` (always `"var_icc"`) and
`value` (the asymptotic variance).

## Details

**Single-rater forms.** For the single-rater intraclass correlation in a
balanced design with \\n\\ subjects and \\k\\ measurements per subject,
the standard large-sample variance derived from the Fisher information
matrix is (Smith, 1956; Donner, 1986; Searle, 1971, ch. 11):
\$\$\mathrm{Var}(\hat\rho) \\\approx\\ \frac{2 (1 - \rho)^2 \bigl(1 +
(k-1)\rho\bigr)^2} {n\\k\\(k - 1)}.\$\$ This expression is exact under
the one-way random-effects model (`ICC(1,1)`) and serves as the standard
large-sample approximation for the two-way single-rater forms
(`ICC(2,1)`, `ICC(3,1)`) as well; differences with the exact two-way
variance vanish at order \\1/n\\ (Bonett, 2002; Burdick & Graybill,
1992). For small \\n\\ or moderate \\k\\ where exact two-way inference
matters, prefer the *F*-distribution-based confidence intervals returned
by [`icc`](https://yelleknek.github.io/DMAR/reference/icc.md), which
follow Shrout and Fleiss (1979) directly.

**Average-of-*k* forms.** Applying the Spearman-Brown transformation
\\\rho_k = k\rho / \[1 + (k - 1)\rho\]\\ together with the delta method
gives the closed-form \$\$\mathrm{Var}(\hat\rho_k) \\\approx\\
\frac{2\\k\\(1 - \rho_k)^2}{n\\(k - 1)},\$\$ expressed directly in the
average-level ICC \\\rho_k\\ (so the user need not invert Spearman-Brown
when working with reliability of composites). The reduction to this form
follows from the substitutions \\1 - \rho = k(1 - \rho_k)/\[k -
(k-1)\rho_k\]\\ and \\1 + (k-1)\rho = k / \[k - (k-1)\rho_k\]\\.

**Use cases.** The asymptotic variance is the natural ingredient for
Wald-style inference, sample size planning for the width of an ICC
confidence interval (compare with Bonett, 2002, which uses a
Fisher-style transformation), and meta-analytic synthesis of ICCs across
studies (the inverse of `value` weights each study). For confidence
intervals themselves, prefer
[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md), which uses
the exact *F*-distribution inversion of Shrout and Fleiss (1979, pp.
425–426).

## References

Bonett, D. G. (2002). Sample size requirements for estimating intraclass
correlations with desired precision. *Statistics in Medicine, 21*(9),
1331–1335. [doi:10.1002/sim.1108](https://doi.org/10.1002/sim.1108)

Burdick, R. K., & Graybill, F. A. (1992). *Confidence Intervals on
Variance Components*. Marcel Dekker.

Donner, A. (1986). A review of inference procedures for the intraclass
correlation coefficient in the one-way random effects model.
*International Statistical Review, 54*(1), 67–82.

Fisher, R. A. (1925). *Statistical Methods for Research Workers*. Oliver
& Boyd.

McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some
intraclass correlation coefficients. *Psychological Methods, 1*(1),
30–46.
[doi:10.1037/1082-989X.1.1.30](https://doi.org/10.1037/1082-989X.1.1.30)

Searle, S. R. (1971). *Linear Models*. Wiley.

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428.

Smith, C. A. B. (1956). On the estimation of intraclass correlation.
*Annals of Human Genetics, 21*(4), 363–373.

## See also

[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md),
[`var_R2`](https://yelleknek.github.io/DMAR/reference/var_R2.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Single-rater one-way ICC at rho = .60 with 30 subjects and 4 raters.
var_icc(rho = 0.60, n = 30, k = 4, type = "ICC(1,1)")
#>  term    value  
#>  var_icc 0.00697

# Same study, but expressed at the average-of-4-rater level.
#     The Spearman-Brown-transformed population value is
#     rho_k = 4*.6 / (1 + 3*.6) = 0.857
var_icc(rho = 0.857, n = 30, k = 4, type = "ICC(1,k)")
#>  term    value  
#>  var_icc 0.00182

# Two-way mixed-model consistency ICC (uses the same one-way
#     asymptotic variance as a large-sample approximation; see Details).
var_icc(rho = 0.60, n = 30, k = 4, type = "ICC(3,1)")
#>  term    value  
#>  var_icc 0.00697
```
