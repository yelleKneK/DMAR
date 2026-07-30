# Plot TOST Equivalence-Test Power Curves Over a Range of True Differences

For each sample size in `n`, draws power as a function of the true mean
difference (or ratio, on the log scale), evaluated at 201 equally spaced
points across the equivalence interval. Returns a ggplot2 object; the
underlying numerical grid is attached as `attr(<plot>, "power_grid")`.

## Usage

``` r
power_equivalence_md_plot(
  alpha_level,
  logscale,
  theta1,
  theta2,
  sigma,
  n,
  nu,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- alpha_level:

  Type I error rate for each of the two one-sided tests.

- logscale:

  Logical. If `TRUE`, the means are compared on the logarithmic scale.

- theta1:

  Lower limit of the equivalence interval.

- theta2:

  Upper limit of the equivalence interval.

- sigma:

  \\\sqrt{\mathrm{error\\ variance}}\\.

- n:

  Vector of sample sizes (one curve per element).

- nu:

  Vector of degrees of freedom for `sigma`, the same length as `n`.

- title:

  Optional plot title (default `"Power of TOST"`).

- subtitle:

  Optional subtitle (typically a reference like `"Phillips, Figure 3"`).

## Value

A `ggplot` object. The 201-row power grid (column 1: true difference;
remaining columns: power for each `n`) is attached as
`attr(<plot>, "power_grid")`.

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

[`power_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_density_equivalence_md`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md)

Other equivalence testing:
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`tost_c()`](https://yelleknek.github.io/DMAR/reference/tost_c.md),
[`tost_r()`](https://yelleknek.github.io/DMAR/reference/tost_r.md),
[`tost_smd()`](https://yelleknek.github.io/DMAR/reference/tost_smd.md)

Other plotting:
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
[`plot_cfa_k()`](https://yelleknek.github.io/DMAR/reference/plot_cfa_k.md),
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`plot_irt_information()`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md),
[`plot_mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md),
[`plot_regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
[`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# \donttest{
# Phillips (1990) Figure 3 reproduction.
n  <- c(9, 12, 18, 24, 30, 40, 60)
nu <- c(7, 10, 16, 22, 28, 38, 58)
power_equivalence_md_plot(
  alpha_level = .05, logscale = FALSE,
  theta1 = -.2, theta2 = .2, sigma = .20,
  n = n, nu = nu,
  subtitle = "Phillips Figure 3"
)


# Diletti (1991) Figure 1c on the log scale.
n_d  <- c(8, 12, 18, 24, 30, 40, 60)
nu_d <- c(6, 10, 16, 22, 28, 38, 58)
power_equivalence_md_plot(
  alpha_level = .05, logscale = TRUE,
  theta1 = .8, theta2 = 1.25, sigma = .20,
  n = n_d, nu = nu_d,
  subtitle = "Diletti, Figure 1c"
)

# }
```
