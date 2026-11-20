# Plot Contrasts Against an Equivalence Region

Draws a forest-style plot of one or more contrast estimates with their
100(1 - 2\\\alpha\\)% confidence intervals against the equivalence
region \\(-\delta_L, \delta_U)\\ and the noninferiority bound
\\-\delta_L\\, colored by the six-way verdict of
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md):
an interval entirely inside the region is equivalent; entirely above
\\\delta_U\\, superior; entirely below \\-\delta_L\\, inferior; a lower
limit above \\-\delta_L\\ with an upper limit at or past \\\delta_U\\,
noninferior only; an upper limit below \\\delta_U\\ with a lower limit
at or below \\-\delta_L\\, not relevantly superior only; and an interval
containing both bounds, inconclusive. The geometry *is* the decision
rule, so the verdict for each contrast can be read from where its
interval falls.

## Usage

``` r
plot_equivalence(
  x = NULL,
  estimate = NULL,
  lower = NULL,
  upper = NULL,
  names = NULL,
  effect_lower = NULL,
  effect_upper = NULL,
  xlab = "Contrast",
  title = NULL,
  palette = "okabe_ito"
)
```

## Arguments

- x:

  Either a single result from
  [`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
  or a list of them (a named list supplies the row labels).
  Alternatively, supply `estimate`, `lower`, and `upper` directly.

- estimate, lower, upper:

  Numeric vectors of contrast estimates and their confidence limits,
  used when `x` is not supplied.

- names:

  Optional character vector of row labels.

- effect_lower, effect_upper:

  Equivalence bounds on the scale of the estimates, as positive
  magnitudes (the region drawn is \\(-\delta_L, +\delta_U)\\, with
  `effect_lower` as \\\delta_L\\ and `effect_upper` as \\\delta_U\\).
  Taken from `x` when it carries `equivalence_c` results; required
  otherwise. If only `effect_upper` is supplied, the bounds are
  symmetric.

- xlab:

  The horizontal axis label. Default `"Contrast"`.

- title:

  Optional plot title.

- palette:

  Character string naming the color palette for the verdict colors.
  Defaults to `"okabe_ito"`, base R's colorblind-safe Okabe-Ito palette;
  `"tableau"` is also available.

## Value

A `ggplot` object. Requires ggplot2 to be installed.

## Details

The shaded band is the equivalence region and the dashed vertical lines
are its bounds; the solid line at zero marks exact equality, which is
the null value of ordinary significance testing and is deliberately
*not* a decision boundary here. Verdicts are recomputed from the
supplied limits and bounds by the rule
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
applies, so the plot cannot disagree with it. The legend sets each
one-sided verdict beside the stronger one that implies it, with
inconclusive between them (equivalent; superior, noninferior only;
inconclusive; not relevantly superior only, inferior).

## References

Chattopadhyay, B., Bandyopadhyay, T., Kelley, K., & Padalunkal, J. J.
(2025). A sequential approach for noninferiority or equivalence of a
linear contrast under cost constraints. *Psychological Methods, 30*(2),
425–439. [doi:10.1037/met0000570](https://doi.org/10.1037/met0000570)

Schuirmann, D. J. (1987). A comparison of the two one-sided tests
procedure and the power approach for assessing the equivalence of
average bioavailability. *Journal of Pharmacokinetics and
Biopharmaceutics, 15*(6), 657–680.

## See also

[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`plot_ci`](https://yelleknek.github.io/DMAR/reference/plot_ci.md)

Other equivalence testing:
[`equivalence_c()`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md),
[`equivalence_proportion_paired()`](https://yelleknek.github.io/DMAR/reference/equivalence_proportion_paired.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`power_density_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_density_equivalence_md.md),
[`power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_c.md),
[`power_equivalence_md()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md)

Other plotting:
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/plot_fa_confirmatory.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md),
[`plot_irt_information()`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md),
[`plot_mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md),
[`plot_randomization_test()`](https://yelleknek.github.io/DMAR/reference/plot_randomization_test.md),
[`plot_regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
[`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Six constructed intervals, one per verdict, against bounds of 5
# (A equivalent, B superior, C noninferior only, D inconclusive,
#  E not relevantly superior only, F inferior):
plot_equivalence(estimate = c(-1.0, 7.0, 3.5, 0.5, -1.5, -8.0),
                 lower    = c(-3.2, 5.5, -1.4, -6.0, -6.6, -10.5),
                 upper    = c( 1.2, 8.5,  8.4,  7.0,  3.6,  -5.5),
                 names    = c("A", "B", "C", "D", "E", "F"),
                 effect_upper = 5)


# From equivalence_c() results; a named list supplies the labels.
res <- list(
  "Focal vs. reference" = equivalence_c(contrast = -5.28, se = 2.49,
                                 df_error = 399, contrast_upper = 5),
  "Within pipeline"     = equivalence_c(contrast = -0.53, se = 2.66,
                                 df_error = 399, contrast_upper = 5)
)
plot_equivalence(res)

```
