# Visualize the Proportion of Variance Explained (\\R^2\\)

Creates a horizontal bar chart showing the observed \\R^2\\ as a
proportion of total variance, with an optional confidence interval
displayed beneath the bar and sample size / predictor-count annotations.

## Usage

``` r
plot_R2(
  R2,
  N = NULL,
  p = NULL,
  conf_level = 0.95,
  show_ci = TRUE,
  show_n = TRUE,
  random_predictors = TRUE,
  title = NULL,
  palette = "okabe_ito",
  colors = NULL
)
```

## Arguments

- R2:

  The observed squared multiple correlation coefficient (\\0 \le R^2 \le
  1\\).

- N:

  Total sample size.

- p:

  Number of predictors.

- conf_level:

  Confidence level for the confidence interval (default `0.95`).

- show_ci:

  Logical. If `TRUE` (the default), a confidence interval is shown
  beneath the proportion bar. Requires both `N` and `p`.

- show_n:

  Logical. If `TRUE` (the default), `N` and `p` are annotated on the
  plot.

- random_predictors:

  Logical. Whether the predictors are random (`TRUE`, the default) or
  fixed. Passed to
  [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md).

- title:

  Optional plot title.

- palette:

  Character string naming the color palette used for the “Explained”
  portion of the bar when `colors` is `NULL`. Defaults to `"okabe_ito"`,
  base R's colorblind-safe Okabe-Ito palette; `"tableau"` is also
  available.

- colors:

  Optional character vector of length 2: the first color fills the
  “Explained” portion of the bar, the second the “Unexplained” portion.
  When `NULL` (the default), the “Explained” portion uses the first
  color of `palette` and the “Unexplained” portion a neutral light gray.

## Value

A `ggplot2` object.

## Note

Requires ggplot2 (listed in `Suggests`).

## References

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

## See also

[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_R`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`plot_ci`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_smd`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

Other plotting:
[`plot_cfa_k()`](https://yelleknek.github.io/DMAR/reference/plot_cfa_k.md),
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
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
# Basic call.
plot_R2(R2 = 0.25, N = 100, p = 5)


# The variations below are not run, since the call above already shows
# the default display and each additional figure has to be drawn. With
# fixed predictors and a 90% confidence interval:
# plot_R2(R2 = 0.35, N = 200, p = 3, conf_level = 0.90,
#         random_predictors = FALSE)

# Without the annotations:
# plot_R2(R2 = 0.10, show_ci = FALSE, show_n = FALSE)
```
