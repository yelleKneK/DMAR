# Density Underlying the TOST Power Calculation

Evaluates the integrand whose integral over \\(0, \mathrm{upper})\\
yields the power of the Schuirmann (1987) two one-sided tests procedure;
see
[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md).
Useful for plotting the integrand and for diagnostic work.

## Usage

``` r
power_density_equivalence_md(
  power_sigma,
  alpha_level,
  theta1,
  theta2,
  diff,
  sigma,
  n,
  nu
)
```

## Arguments

- power_sigma:

  Numeric vector of \\\sigma\\ values at which to evaluate the
  integrand.

- alpha_level:

  Type I error rate for each of the two one-sided tests.

- theta1:

  Lower limit of the equivalence interval on the appropriate scale
  (regular or log).

- theta2:

  Upper limit of the equivalence interval on the appropriate scale
  (regular or log).

- diff:

  True difference in treatment means (ratio on the log scale) on the
  appropriate scale.

- sigma:

  \\\sqrt{\mathrm{error\\ variance}}\\.

- n:

  Number of subjects per treatment.

- nu:

  Degrees of freedom for `sigma`.

## Value

A `data.frame` with one row per supplied `power_sigma`, and columns
`power_sigma` and `power_density`.

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
18*(2), 139–144.
[doi:10.1007/BF01063556](https://doi.org/10.1007/BF01063556)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)

Other equivalence testing:
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`tost_c()`](https://yelleknek.github.io/DMAR/reference/tost_c.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Density at a single value of sigma:
power_density_equivalence_md(power_sigma = 0.10, alpha_level = .05,
                             theta1 = -.2, theta2 = .2, diff = .05,
                             sigma = .20, n = 24, nu = 22)
#>   power_sigma power_density
#> 1         0.1    0.02297896

# Vectorized over a grid:
grid <- power_density_equivalence_md(
  power_sigma = seq(0.01, 0.40, length.out = 50),
  alpha_level = .05, theta1 = -.2, theta2 = .2, diff = .05,
  sigma = .20, n = 24, nu = 22
)
head(grid)
#>   power_sigma power_density
#> 1  0.01000000  3.625481e-22
#> 2  0.01795918  7.451634e-17
#> 3  0.02591837  1.498398e-13
#> 4  0.03387755  3.634158e-11
#> 5  0.04183673  2.582758e-09
#> 6  0.04979592  8.170958e-08
```
