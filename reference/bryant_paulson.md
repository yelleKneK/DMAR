# The Bryant-Paulson Generalized Studentized Range Distribution

Distribution function (`pbryant_paulson`), quantile/critical-value
function (`qbryant_paulson`), and density (`dbryant_paulson`) for the
Bryant-Paulson generalized studentized range, the sampling distribution
of the studentized range of covariate-*adjusted* means in the analysis
of covariance (ANCOVA) when the covariate(s) are *random*. These are the
analysis-of-covariance analogs of
[`ptukey`](https://rdrr.io/r/stats/Tukey.html) /
[`qtukey`](https://rdrr.io/r/stats/Tukey.html) and supply the critical
values needed for Tukey-Kramer-type simultaneous confidence intervals on
(and tests of) contrasts of adjusted means.

## Usage

``` r
pbryant_paulson(q, n_covariates, n_groups, df, lower_tail = TRUE)

qbryant_paulson(prob, n_covariates, n_groups, df, lower_tail = TRUE)

dbryant_paulson(q, n_covariates, n_groups, df)
```

## Arguments

- q:

  Vector of quantiles (values of the generalized studentized range
  statistic).

- n_covariates:

  The number of random covariates, \\p\\ (\\p \ge 0\\). With
  `n_covariates = 0` the distribution is exactly the ordinary
  studentized range, computed as in
  [`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
  ([`ptukey`](https://rdrr.io/r/stats/Tukey.html) and
  [`qtukey`](https://rdrr.io/r/stats/Tukey.html) where those are
  accurate, the exact distribution where they are not).

- n_groups:

  The number of groups (treatments) being compared, \\k\\ (\\k \ge 2\\);
  this is the “sample size” of the range.

- df:

  The error degrees of freedom of the ANCOVA model, \\\nu\\. In a
  one-way ANCOVA with \\N\\ total observations, \\k\\ groups, and \\p\\
  covariates, \\\nu = N - k - p\\. `Inf` gives the limiting distribution
  (see Details).

- lower_tail:

  Logical; if `TRUE` (default) probabilities are \\P(Q \le q)\\,
  otherwise \\P(Q \> q)\\.

- prob:

  Vector of probabilities. For `qbryant_paulson` this is the cumulative
  probability (e.g., `0.95` returns the upper 5% critical value).

## Value

Numeric vectors. `pbryant_paulson` returns cumulative (or upper-tail)
probabilities, `dbryant_paulson` returns density values, and
`qbryant_paulson` returns critical values (quantiles) of the generalized
studentized range. Results are recycled to the length of the longest of
`q`/`prob` and the parameter arguments.

## Details

**The statistic.** In a balanced ANCOVA with \\k\\ groups and \\p\\
random covariates, let \\\hat\theta_i\\ be the adjusted group means and
\\\hat\sigma\_{y \mid x}\\ the square root of the ANCOVA error mean
square (on \\\nu\\ degrees of freedom). The Bryant-Paulson statistic is
the studentized range of the adjusted means, \$\$Q \\=\\ \frac{\max_i
\hat\theta_i - \min_i \hat\theta_i}{\hat\sigma\_{y\mid x}\sqrt{K_1 -
K_2}},\$\$ where \\K_1 - K_2\\ is the design constant that scales the
variance of a single adjusted mean (for a one-way design with \\n\\ per
group, \\K_1 - K_2 = 1/n\\). The studentizer uses only this
“between-only” standard error: the extra sampling variability induced by
having to *estimate* the covariate adjustment from random covariates is
carried by the distribution of \\Q\\ itself, not by a per-comparison
standard-error correction. This is what distinguishes the procedure from
naively applying Tukey's method to adjusted means.

**The distribution.** Bryant and Paulson (1976) give the exact CDF of
\\Q_p\\ in their Equation (17), a single integral over a variable that
combines the \\\chi^2\_\nu\\ error estimate with a random
covariate-shrinkage factor \\\delta\\ that, by their Equations
(11)–(12), has a \\\mathrm{Beta}((\nu+1)/2,\\ p/2)\\ distribution.
Carrying out the error integral with the studentized-range routine
`ptukey` reduces Equation (17) to the equivalent one-dimensional form
\$\$P(Q_p \le q) \\=\\ \int_0^1 \mathrm{ptukey}\\\left(q\sqrt{\delta};\\
k,\\ \nu\right)\\ f\_{\mathrm{Beta}}\\\left(\delta;\\
\tfrac{\nu+1}{2},\\ \tfrac{p}{2}\right) d\delta,\$\$ which this package
evaluates (the reduction is exact; see the source-code comments in
`R/bryant_paulson.R` for the one-line derivation from Bryant and
Paulson's p. 634 conditioning argument). When \\p = 0\\ the factor
\\\delta\\ degenerates at 1 and \\Q_p\\ is exactly the ordinary
studentized range (Bryant and Paulson, 1976, Sec. 1), so the code
short-circuits to the studentized range distribution function. The
integral is evaluated with
[`integrate`](https://rdrr.io/r/stats/integrate.html) after the change
of variables \\\delta = 1 - u^2\\, which removes the endpoint
singularity of the Beta weight at \\\delta = 1\\ when \\p = 1\\ and
makes the quadrature converge in a few subdivisions; `qbryant_paulson`
inverts it with [`uniroot`](https://rdrr.io/r/stats/uniroot.html).
Bryant and Bruvold (1980) later showed the same distribution and
critical values remain valid when the covariates are *not* identically
distributed across groups (their grouped-covariate model, Eq. 1.3), and
added the Duncan multiple-range extension.

**Accuracy at small df.** The studentized range distribution function
inside the integral is [`ptukey`](https://rdrr.io/r/stats/Tukey.html)
where that is accurate and is evaluated directly elsewhere, without
`ptukey`, by integrating the probability integral of the range against
the \\\chi^2\_\nu\\ error density. `ptukey` loses accuracy at small
\\\nu\\, increasingly so as \\k\\ grows (at \\\nu = 3\\, \\k = 20\\ its
probability error reaches \\\approx 3\times10^{-4}\\, enough to move the
critical value by about 0.2, and it is larger at \\\nu = 2\\), and again
above 10000 df; the direct evaluation is used for fewer than 9 error df
with up to 6 groups, fewer than 11 with up to 10, fewer than 15 with up
to 20, fewer than 21 with up to 50, fewer than 31 with up to 100, fewer
than 41 with up to 250, and, with more groups, below a bound that rises
by steps to 111 with more than 3000 groups; above 10000 and up to
\\10^8\\ df; and, with more than 50 groups, wherever the probability is
below .75 (.80 with more than 1000 groups;
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
gives the whole rule and how it was measured). At probabilities above
.999 or below .001 the quantile also comes from the exact distribution,
with or without covariates (with covariates every node of the integral
over the covariate-shrinkage weight then uses the exact studentized
range: through `ptukey` the quantile at .9999 with one covariate, 10
groups, and 11 df was 12.977072 instead of 12.977142), and with no
covariates also wherever `qtukey` stops short of its root. The direct
path costs a fraction of a second per probability.

**Large df.** As \\\nu\\ grows, \\\delta\\ tends to 1 and \\Q_p\\ to the
studentized range with infinite df, whatever the number of covariates;
with `df = Inf` the functions return that limit, computed as for
`n_covariates = 0`. At large \\\nu\\ the weight of the integral lies
within a few multiples of \\\sqrt{p/\nu}\\ of \\u = 0\\, where \\\delta
= 1 - u^2\\, so the integral over \\u\\ stops at the square root of the
upper \\10^{-15}\\ quantile of the \\\mathrm{Beta}(p/2,\\ (\nu+1)/2)\\
distribution of \\u^2\\ rather than at 1, which leaves out weight
\\10^{-15}\\. Over the whole unit interval the quadrature missed the
weight from about \\4 \times 10^7\\ df with one covariate (\\5.6 \times
10^7\\ with six), and the distribution function came out near 0. The
weight is computed on the log scale, since with many covariates its
factors overflow and underflow separately. At \\10^8\\ df the
distribution function is still below its infinite-df limit by about
\\10^{-8}\\ (\\7 \times 10^{-9}\\ at 3.7 with one covariate and four
groups).

**Validation.** The implementation reproduces Bryant and Paulson's
(1976) Table 1 exactly, to the two decimal places tabled, over the whole
of its range: both tail areas (\\\alpha = .05\\ and \\\alpha = .01\\),
all three covariate counts (\\p = 1, 2, 3\\), every tabled number of
groups (\\k = 2, \ldots, 8, 10, 12, 16, 20\\), and every tabled error
degrees of freedom (\\\nu = 2, \ldots, 8, 10, 12, 14, 16, 18, 20, 24,
30, 40, 60, 120\\), which is 1188 critical values in all. The corners of
the table are included: \\q\_{.01;\\1,2,2} = 19.09\\,
\\q\_{.01;\\3,20,2} = 73.01\\, \\q\_{.01;\\3,20,3} = 33.13\\, and
\\q\_{.05;\\1,6,14} = 4.83\\ (the value used in the Bryant and Bruvold,
1980 worked example). Two of the 1188 entries, \\q\_{.01;\\2,8,3}\\ and
\\q\_{.01;\\2,20,4}\\, have exact values of 23.165013 and 19.745008,
each roughly \\10^{-5}\\ above the 23.165 and 19.745 half-way points, so
they round to 23.17 and 19.75; the 1976 table rounds them down, to 23.16
and 19.74. Both values were confirmed to fourteen significant figures by
two independent high-order quadrature engines that share no code with
the package implementation. A large-scale simulation of the Bryant and
Paulson statistic confirms the computed values independently, and the
Bryant and Bruvold (1980) Table 2 Duncan ranges are reproduced as well.
See the package tests.

## References

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*, 631–638.
[doi:10.1093/biomet/63.3.631](https://doi.org/10.1093/biomet/63.3.631)

Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures
in the analysis of covariance. *Journal of the American Statistical
Association, 75*(372), 874–880.
[doi:10.2307/2287175](https://doi.org/10.2307/2287175)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
for the simultaneous confidence intervals these critical values produce;
[`ptukey`](https://rdrr.io/r/stats/Tukey.html) and
[`qtukey`](https://rdrr.io/r/stats/Tukey.html) for the ordinary
(fixed-covariate or no-covariate) studentized range.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Critical value from the worked example of Bryant and Bruvold (1980):
# k = 6 panels, p = 1 covariate, nu = 14 error df, alpha = .05. The quantile
# is found by inverting the distribution function with uniroot; it is 4.83,
# the entry in Table 1 of Bryant and Paulson (1976).
qbryant_paulson(0.95, n_covariates = 1, n_groups = 6, df = 14)
#> [1] 4.829856

# The ordinary Tukey value (ignoring that the covariate is random and
# estimated) is smaller, 4.64, so it yields intervals that are too narrow:
qtukey(0.95, nmeans = 6, df = 14)
#> [1] 4.638538

# How much too narrow: the Bryant-Paulson area beyond the Tukey value is the
# familywise error rate Tukey's method actually delivers in this design.
# With no covariate the two distributions coincide, so the area is .0499,
# the nominal .05 up to the rounding of 4.64 itself. Each additional random
# covariate carries more estimation uncertainty, stretches the distribution
# to the right, and pushes the rate up, to .063, .076, and .091.
pbryant_paulson(4.64, n_covariates = 0:3, n_groups = 6, df = 14,
                lower_tail = FALSE)
#> [1] 0.04990724 0.06273456 0.07633320 0.09054790

# The p = 0 entry above is exactly the ordinary studentized range:
ptukey(4.64, nmeans = 6, df = 14, lower.tail = FALSE)
#> [1] 0.04990724
```
