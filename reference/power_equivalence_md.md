# Power of the Two One-Sided Tests Procedure (TOST) for Equivalence

Computes the power of the Schuirmann (1987) two one-sided tests
procedure , the probability that a \\(1 - 2\alpha)\\ confidence interval
for the mean difference (or ratio, on the log scale) lies entirely
within the equivalence interval \\\[\theta_1, \theta_2\]\\, by numerical
integration over the chi distribution of the sample standard deviation.

## Usage

``` r
power_equivalence_md(
  alpha_level,
  logscale,
  ltheta1,
  ltheta2,
  ldiff,
  sigma,
  n,
  nu
)
```

## Arguments

- alpha_level:

  Type I error rate for each of the two one-sided tests (typically
  `0.05`). The full equivalence test uses a \\(1 - 2\alpha)\\ confidence
  interval.

- logscale:

  Logical. If `TRUE`, treatment means are compared on the logarithmic
  scale; `ltheta1`, `ltheta2`, and `ldiff` are expected as ratios
  (untransformed) and are log-transformed internally.

- ltheta1:

  Lower limit of the equivalence interval (on the original scale; logged
  internally if `logscale = TRUE`).

- ltheta2:

  Upper limit of the equivalence interval (on the original scale; logged
  internally if `logscale = TRUE`).

- ldiff:

  True difference in treatment means (or ratio on the log scale).

- sigma:

  \\\sqrt{\mathrm{error\\ variance}}\\; root-MSE from an ANOVA. On the
  log scale, this is the coefficient of variation.

- n:

  Number of subjects per treatment (or total subjects in a crossover
  design).

- nu:

  Degrees of freedom associated with `sigma`.

## Value

A one-row `data.frame` with columns `term` (`"power"`) and `value` (the
computed power, in \\\[0, 1\]\\).

## Details

The computation conditions on the error standard deviation the study
will actually observe. Given that value, whether the confidence interval
fits inside the equivalence interval is an ordinary normal probability,
and the power is that probability averaged over the chi distribution the
error standard deviation follows on `nu` degrees of freedom. The
averaging is carried out on a unit-free scale, with the equivalence
limits expressed in standard errors and the error standard deviation as
a multiple of `sigma`, so the power depends on the design rather than on
the units of the response: multiplying `ltheta1`, `ltheta2`, `ldiff`,
and `sigma` by a common factor leaves the answer unchanged.

For Phillips's (1990) original example (regular-scale two-period
crossover with \\\theta_1 = -0.2\\, \\\theta_2 = 0.2\\, CV \\= 0.20\\,
\\\delta = 0.05\\, \\n = 24\\, \\\nu = 22\\), this function reproduces
the published value of \\0.8029678\\ (Phillips, 1990, Table 1, 5th row,
5th column).

## Note

See the legacy `MBESS` package (Kelley, 2007a, 2007b) for additional
details and discussion.

## References

Diletti, E., Hauschke, D., & Steinijans, V. W. (1991). Sample size
determination of bioequivalence assessment by means of confidence
intervals. *International Journal of Clinical Pharmacology, Therapy and
Toxicology, 29*(1), 1–8.

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
[doi:10.3758/BF03192993](https://doi.org/10.3758/BF03192993)

Phillips, K. F. (1990). Power of the two one-sided tests procedure in
bioequivalence. *Journal of Pharmacokinetics and Biopharmaceutics,
18*(2), 139–144.
[doi:10.1007/BF01063556](https://doi.org/10.1007/BF01063556)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`power_equivalence_md_plot`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`power_density_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md)

Other equivalence testing:
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`tost_c()`](https://yelleknek.github.io/DMAR/reference/tost_c.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Phillips (1990) Table 1, 5th row, 5th column. Expected: 0.8029678.
power_equivalence_md(alpha_level = .05, logscale = FALSE,
                     ltheta1 = -.2, ltheta2 = .2, ldiff = .05,
                     sigma = .20, n = 24, nu = 22)
#>  term  value
#>  power 0.803

# Diletti (1991) Table 1, on the log scale (ratio of test to reference).
# Expected: 0.7922796.
power_equivalence_md(alpha_level = .05, logscale = TRUE,
                     ltheta1 = .8, ltheta2 = 1.25, ldiff = 1.05,
                     sigma = .20, n = 18, nu = 16)
#>  term  value
#>  power 0.792
```
