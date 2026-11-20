# Density Underlying the TOST Power Calculation

Evaluates the integrand whose integral over \\(0, \mathrm{upper})\\
yields the power of the Schuirmann (1987) two one-sided tests procedure;
see
[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md).
Useful for plotting the integrand and for diagnostic work.

## Usage

``` r
power_density_equivalence_md(
  sd_values,
  alpha_level,
  lower_bound,
  upper_bound,
  true_effect,
  true_sd,
  n_per_treatment,
  df_error
)
```

## Arguments

- sd_values:

  Numeric vector of \\\sigma\\ values at which to evaluate the
  integrand.

- alpha_level:

  Type I error rate for each of the two one-sided tests.

- lower_bound:

  Lower limit of the equivalence interval on the analysis scale: a
  difference, or the logarithm of a ratio.

- upper_bound:

  Upper limit of the equivalence interval on the analysis scale: a
  difference, or the logarithm of a ratio.

- true_effect:

  The true (population) effect on the analysis scale: the difference in
  treatment means, or the logarithm of their ratio.

- true_sd:

  \\\sqrt{\mathrm{error\\ variance}}\\.

- n_per_treatment:

  Number of subjects who receive each treatment (the size of each group
  in a parallel design; every subject in a crossover design), as in
  [`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md).

- df_error:

  Error degrees of freedom for `true_sd`.

## Value

A `data.frame` with one row per supplied `sd_values`, and columns
`sd_values` and `power_density`.

## Details

**The analysis scale.** Unlike
[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
which takes the bounds and the true effect on the original scale and
logs ratios when `logscale = TRUE`, this function has no `logscale`
argument and transforms nothing: `lower_bound`, `upper_bound`, and
`true_effect` are on the scale of the analysis, differences of the
means, or, for means compared on the log scale, logarithms of the ratios
(`log(0.8)` for a ratio of 0.8).

## Note

See the legacy `MBESS` package (Kelley, 2007a, 2007b) for additional
details and discussion.

## References

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

[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Density at a single value of true_sd:
power_density_equivalence_md(sd_values = 0.10, alpha_level = .05,
                             lower_bound = -.2, upper_bound = .2,
                             true_effect = .05, true_sd = .20,
                             n_per_treatment = 24, df_error = 22)
#>   sd_values power_density
#> 1       0.1    0.02297896

# Vectorized over a grid:
grid <- power_density_equivalence_md(
  sd_values = seq(0.01, 0.40, length.out = 50),
  alpha_level = .05, lower_bound = -.2, upper_bound = .2,
  true_effect = .05, true_sd = .20, n_per_treatment = 24, df_error = 22
)
head(grid)
#>    sd_values power_density
#> 1 0.01000000  3.625481e-22
#> 2 0.01795918  7.451634e-17
#> 3 0.02591837  1.498398e-13
#> 4 0.03387755  3.634158e-11
#> 5 0.04183673  2.582758e-09
#> 6 0.04979592  8.170958e-08

# On the log scale the bounds and the true effect are logarithms of the
# ratios: the ratios .8, 1.25, and 1.05 of power_equivalence_md(logscale =
# TRUE) are given here as their logarithms.
power_density_equivalence_md(sd_values = 0.10, alpha_level = .05,
                             lower_bound = log(.8), upper_bound = log(1.25),
                             true_effect = log(1.05), true_sd = .20,
                             n_per_treatment = 18, df_error = 16)
#>   sd_values power_density
#> 1       0.1     0.1317922
```
