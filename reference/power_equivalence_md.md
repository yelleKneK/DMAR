# Power of the Two One-Sided Tests Procedure (TOST) for Equivalence

Computes the power of the Schuirmann (1987) two one-sided tests
procedure , the probability that a \\(1 - 2\alpha)\\ confidence interval
for the mean difference (or ratio, on the log scale) lies entirely
within the equivalence interval \\\[\theta_1, \theta_2\]\\
(`lower_bound` and `upper_bound`), by numerical integration over the chi
distribution of the sample standard deviation.

## Usage

``` r
power_equivalence_md(
  alpha_level,
  logscale,
  lower_bound,
  upper_bound,
  true_effect,
  true_sd,
  n_per_treatment,
  df_error
)
```

## Arguments

- alpha_level:

  Type I error rate for each of the two one-sided tests (typically
  `0.05`). The full equivalence test uses a \\(1 - 2\alpha)\\ confidence
  interval.

- logscale:

  Logical. If `TRUE`, treatment means are compared on the logarithmic
  scale; `lower_bound`, `upper_bound`, and `true_effect` are then given
  as ratios of test to reference on the original scale (untransformed)
  and are log-transformed internally. If `FALSE`, they are differences
  of the means.

- lower_bound:

  Lower limit of the equivalence interval on the original scale: a ratio
  when `logscale = TRUE` (logged internally), a difference otherwise.

- upper_bound:

  Upper limit of the equivalence interval on the original scale: a ratio
  when `logscale = TRUE` (logged internally), a difference otherwise.

- true_effect:

  The true (population) effect: the ratio of the treatment means on the
  original scale when `logscale = TRUE` (logged internally), the
  difference in treatment means otherwise.

- true_sd:

  \\\sqrt{\mathrm{error\\ variance}}\\; root-MSE from an ANOVA. On the
  log scale, this is the coefficient of variation.

- n_per_treatment:

  Number of subjects who receive each treatment: the size of each group
  in a parallel design, and the number of subjects in a crossover
  design, where every subject receives both treatments. It is the same
  quantity in both designs, the \\n\\ of the standard error \\\sigma
  \sqrt{2 / n}\\ of the estimated difference, so in a parallel design it
  is the size of one group and never the total. Must be greater than 1.

- df_error:

  Error degrees of freedom associated with `true_sd`.

## Value

A one-row `data.frame` with columns `term` (`"power"`) and `value` (the
computed power, in \\\[0, 1\]\\).

## Details

The computation conditions on the error standard deviation the study
will actually observe. Given that value, whether the confidence interval
fits inside the equivalence interval is an ordinary normal probability,
and the power is that probability averaged over the chi distribution the
error standard deviation follows on `df_error` degrees of freedom. The
averaging is carried out on a unit-free scale, with the equivalence
limits expressed in standard errors and the error standard deviation as
a multiple of `true_sd`, so the power depends on the design rather than
on the units of the response: with `logscale = FALSE`, multiplying
`lower_bound`, `upper_bound`, `true_effect`, and `true_sd` by a common
factor leaves the answer unchanged.

**The scale of the bounds and the true effect.** This function and
[`power_equivalence_md_plot`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)
take `lower_bound`, `upper_bound`, and `true_effect` on the original
scale: ratios of test to reference with `logscale = TRUE`, which are
logged before the computation, and differences with `logscale = FALSE`,
which are used as given.
[`power_density_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md)
takes them on the analysis scale and transforms nothing, so
`power_equivalence_md(logscale = TRUE, lower_bound = .8, ...)`
corresponds to
`power_density_equivalence_md(lower_bound = log(.8), ...)`.

For Phillips's (1990) original example (regular-scale two-period
crossover with \\\theta_1 = -0.2\\, \\\theta_2 = 0.2\\, CV \\= 0.20\\,
\\\delta = 0.05\\, \\n = 24\\, \\\nu = 22\\, which are `lower_bound`,
`upper_bound`, `true_sd`, `true_effect`, `n_per_treatment`, and
`df_error`), this function reproduces the published value of
\\0.8029678\\ (Phillips, 1990, Table 1, 5th row, 5th column).

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
18*(2), 137–144.
[doi:10.1007/BF01063556](https://doi.org/10.1007/BF01063556)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`power_equivalence_md_plot`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`power_density_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md)

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Table 1 of Phillips, 1990, fifth row and fifth column, where the
# published power is 0.8029678.
power_equivalence_md(alpha_level = .05, logscale = FALSE,
                     lower_bound = -.2, upper_bound = .2,
                     true_effect = .05, true_sd = .20,
                     n_per_treatment = 24, df_error = 22)
#>  term  value
#>  power 0.803

# Table 1 of Diletti et al., 1991, on the log scale, so the limits and
# the true difference are ratios of test to reference. The published
# power is 0.7922796.
power_equivalence_md(alpha_level = .05, logscale = TRUE,
                     lower_bound = .8, upper_bound = 1.25,
                     true_effect = 1.05, true_sd = .20,
                     n_per_treatment = 18, df_error = 16)
#>  term  value
#>  power 0.792
```
