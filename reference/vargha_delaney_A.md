# Vargha and Delaney's *A* (Stochastic-Superiority Effect Size)

Computes Vargha and Delaney's (2000) *A*, the probability that a
randomly drawn observation from the first sample exceeds a randomly
drawn observation from the second (with tied pairs counted as half),
together with its asymptotic standard error from DeLong, DeLong, and
Clarke-Pearson (1988) and a confidence interval on the population *A*.
*A* is equivalent to the receiver-operating-characteristic
area-under-the-curve (AUC) and to the common-language effect size
(McGraw & Wong, 1992) under a continuous-response assumption, and is a
robust, scale-free ordinal effect size that does not require equal
variances or normality.

## Usage

``` r
vargha_delaney_A(
  x,
  y = NULL,
  data = NULL,
  conf_level = 0.95,
  ci_method = c("logit", "wald")
)
```

## Arguments

- x:

  Either a numeric vector of observations from group 1, or a two-sided
  formula of the form `outcome ~ group` (in which case `data` must be
  supplied and the grouping variable must have exactly two levels).

- y:

  Numeric vector of observations from group 2. Ignored when `x` is a
  formula.

- data:

  Optional data frame containing the variables named in the formula `x`.

- conf_level:

  Confidence coverage for a symmetric interval (default `0.95`).

- ci_method:

  Either `"logit"` (default; Wald interval on the logit of *A* with
  back-transformation, recommended for finite samples; Newcombe, 2006b)
  or `"wald"` (untransformed Wald on the original scale, clipped to \[0,
  1\]).

## Value

A one-row `data.frame` with columns `A` (point estimate), `se`
(DeLong-DeLong-Clarke-Pearson standard error), `lower_limit` and
`upper_limit` (confidence limits at `conf_level`), `z_value` and
`p_value` (Wald test of \\H_0\\: A = 0.5\\, i.e., stochastic equality),
`n_1` and `n_2` (group sample sizes), and `ci_method`.

## Details

**Definition.** For independent samples \\X_1, \ldots, X\_{n_1}\\ and
\\Y_1, \ldots, Y\_{n_2}\\, \$\$A = \Pr(X \> Y) + \tfrac{1}{2}\\\Pr(X =
Y).\$\$ Values of \\A = 0.5\\ indicate stochastic equality; \\A \> 0.5\\
indicates that group 1 tends to score higher. Qualitative magnitude
labels for *A* are not reported here, in keeping with the DMAR
convention of reporting effect sizes as numbers with confidence
intervals.

**Sample estimate.** Equivalent rank-based computation: \$\$\hat A =
\frac{\bar R_X - (n_1 + 1)/2}{n_2},\$\$ where \\\bar R_X\\ is the mean
rank of the first sample in the pooled ranking with mid-ranks for ties
(Vargha & Delaney, 2000, p. 109). Equivalently, \\\hat A = U / (n_1
n_2)\\, with \\U\\ the Mann-Whitney *U*-statistic counting \\X_i \>
Y_j\\ (tied pairs at \\1/2\\).

**Standard error.** The function uses the DeLong-DeLong-Clarke-Pearson
(1988) U-statistic variance estimator, which is unbiased under sampling
from any joint distribution (no parametric or homoscedasticity
assumption). Defining the placement components \$\$V\_{10}(X_i) =
\frac{1}{n_2}\sum\_{j} \psi(X_i, Y_j), \qquad V\_{01}(Y_j) =
\frac{1}{n_1}\sum\_{i} \psi(X_i, Y_j),\$\$ with \\\psi(x, y) = 1,
\tfrac{1}{2}, 0\\ as \\x \> y, =, \<\\, the variance estimate is
\$\$\widehat{\mathrm{Var}}(\hat A) = \frac{S^2\_{10}}{n_1} +
\frac{S^2\_{01}}{n_2},\$\$ where \\S^2\_{10}\\ and \\S^2\_{01}\\ are the
sample variances of the \\V\_{10}\\ and \\V\_{01}\\ placement
components. This is identical to the (single-curve) DeLong AUC variance
and is the standard nonparametric variance for the Mann-Whitney
functional (Brunner & Munzel, 2000).

**Confidence interval.** `ci_method = "logit"` (the default) constructs
a Wald interval on \\\mathrm{logit}(A) = \log\\A/(1-A)\\\\ using the
delta method standard error \\\widehat{\mathrm{SE}}(\hat A)/\\\hat A(1 -
\hat A)\\\\ and back- transforms with the inverse logit. Newcombe
(2006a, 2006b) showed in extensive coverage simulations that logit-Wald
has notably better small- sample coverage than untransformed Wald, while
remaining simple and free of iteration. `ci_method = "wald"` returns the
untransformed Wald interval, clipped to \\\[0, 1\]\\.

## References

Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher
problem: Asymptotic theory and a small-sample approximation.
*Biometrical Journal, 42*(1), 17–25.
[doi:10.1002/(SICI)1521-4036(200001)42:1\<17::AID-BIMJ17\>3.0.CO;2-U](https://doi.org/10.1002/%28SICI%291521-4036%28200001%2942%3A1%3C17%3A%3AAID-BIMJ17%3E3.0.CO%3B2-U)

DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988). Comparing
the areas under two or more correlated receiver operating characteristic
curves: A nonparametric approach. *Biometrics, 44*(3), 837–845.

Hanley, J. A., & McNeil, B. J. (1982). The meaning and use of the area
under a receiver operating characteristic (ROC) curve. *Radiology,
143*(1), 29–36.

McGraw, K. O., & Wong, S. P. (1992). A common language effect size
statistic. *Psychological Bulletin, 111*(2), 361–365.
[doi:10.1037/0033-2909.111.2.361](https://doi.org/10.1037/0033-2909.111.2.361)

Newcombe, R. G. (2006a). Confidence intervals for an effect size measure
based on the Mann-Whitney statistic. Part 1: General issues and
tail-area-based methods. *Statistics in Medicine, 25*(4), 543–557.
[doi:10.1002/sim.2323](https://doi.org/10.1002/sim.2323)

Newcombe, R. G. (2006b). Confidence intervals for an effect size measure
based on the Mann-Whitney statistic. Part 2: Asymptotic methods and
evaluation. *Statistics in Medicine, 25*(4), 559–573.
[doi:10.1002/sim.2324](https://doi.org/10.1002/sim.2324)

Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the
CL common language effect size statistics of McGraw and Wong. *Journal
of Educational and Behavioral Statistics, 25*(2), 101–132.
[doi:10.3102/10769986025002101](https://doi.org/10.3102/10769986025002101)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two numeric vectors.
set.seed(113)
x <- rnorm(40, mean = 0.6)
y <- rnorm(40, mean = 0)
vargha_delaney_A(x, y)
#>         A         se lower_limit upper_limit  z_value     p_value n_1 n_2
#> 1 0.66625 0.06130329   0.5376598   0.7741021 2.711926 0.006689345  40  40
#>   ci_method
#> 1     logit

# Formula interface on the pygmalion field experiment. The first factor
# level (Control) forms group 1, so A below 0.5 says a randomly drawn
# control child tends to score below a child from the bloomer group.
vargha_delaney_A(iq_8 ~ treatment, data = pygmalion)
#>           A         se lower_limit upper_limit   z_value    p_value n_1 n_2
#> 1 0.4292429 0.04119075   0.3510403   0.5111451 -1.717791 0.08583467 246  64
#>   ci_method
#> 1     logit

# Wald (untransformed) interval rather than logit.
vargha_delaney_A(x, y, ci_method = "wald")
#>         A         se lower_limit upper_limit  z_value     p_value n_1 n_2
#> 1 0.66625 0.06130329   0.5460978   0.7864022 2.711926 0.006689345  40  40
#>   ci_method
#> 1      wald
```
