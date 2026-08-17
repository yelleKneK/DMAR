# Visualize a Standardized Mean Difference With Overlapping Distributions

Creates a publication-quality plot showing two normal distributions
separated by the standardized mean difference (*d*). The plot includes a
confidence interval for the population effect size and sample size
annotations, both shown by default.

## Usage

``` r
plot_smd(
  smd = NULL,
  n_1 = NULL,
  n_2 = NULL,
  group_1 = NULL,
  group_2 = NULL,
  conf_level = 0.95,
  show_ci = TRUE,
  show_n = TRUE,
  title = NULL,
  group_labels = c("Group 1", "Group 2"),
  palette = "okabe_ito",
  colors = NULL
)
```

## Arguments

- smd:

  The standardized mean difference (Cohen's *d*).

- n_1:

  Sample size for Group 1.

- n_2:

  Sample size for Group 2.

- group_1:

  Raw data for Group 1. When provided, `smd`, `n_1`, and `n_2` are
  computed from the data.

- group_2:

  Raw data for Group 2.

- conf_level:

  Confidence level for the confidence interval (default `0.95`).

- show_ci:

  Logical. If `TRUE` (the default), a confidence interval for the
  population standardized mean difference is displayed beneath the
  distributions. Requires both `n_1` and `n_2`.

- show_n:

  Logical. If `TRUE` (the default), per-group sample sizes are annotated
  on the plot.

- title:

  Optional character string for the plot title. Defaults to
  `"Standardized Mean Difference"`.

- group_labels:

  Character vector of length 2 giving labels for the two groups.
  Defaults to `c("Group 1", "Group 2")`.

- palette:

  Character string naming the color palette used when `colors` is
  `NULL`. Defaults to `"okabe_ito"`, base R's colorblind-safe Okabe-Ito
  palette; `"tableau"` is also available.

- colors:

  Optional character vector of length 2 giving fill colors for the two
  groups. When `NULL` (the default), the first two colors of `palette`
  are used.

## Value

A `ggplot2` object that can be further customized with standard ggplot2
layers, scales, and themes.

## Details

Two unit-variance normal distributions are drawn, centered at 0 (Group 2
/ reference) and *d* (Group 1 / focal). The semi-transparent fills make
the overlap visible, giving a direct visual impression of how much the
distributions differ.

When `show_ci = TRUE` and both `n_1` and `n_2` are available, the
function calls
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) to
compute the noncentral *t* based confidence interval and displays it as
a horizontal bar beneath the curves. A filled dot marks the point
estimate and vertical caps mark the confidence bounds.

## Note

Requires ggplot2 (listed in `Suggests`). Install it with
`install.packages("ggplot2")` if needed.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

## See also

[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`plot_ci`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_R2`](https://yelleknek.github.io/DMAR/reference/plot_R2.md)

Other plotting:
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
[`plot_cfa_k()`](https://yelleknek.github.io/DMAR/reference/plot_cfa_k.md),
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md),
[`plot_irt_information()`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md),
[`plot_mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md),
[`plot_randomization_test()`](https://yelleknek.github.io/DMAR/reference/plot_randomization_test.md),
[`plot_regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md),
[`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# From a known effect size and sample sizes.
plot_smd(smd = 0.50, n_1 = 50, n_2 = 50)


# The variations below are not run, since the call above already shows
# the default display and each additional figure has to be drawn. From
# raw data, where the standardized mean difference and both sample
# sizes are taken from the data:
# set.seed(113)
# g1 <- rnorm(40, mean = 0.6, sd = 1)
# g2 <- rnorm(40, mean = 0.0, sd = 1)
# plot_smd(group_1 = g1, group_2 = g2)

# Without the confidence interval or the sample size annotations:
# plot_smd(smd = 0.80, show_ci = FALSE, show_n = FALSE)

# Custom group labels and title:
# plot_smd(smd = 0.45, n_1 = 75, n_2 = 75,
#          group_labels = c("Treatment", "Control"),
#          title = "Treatment Effect on Reading Scores")
```
