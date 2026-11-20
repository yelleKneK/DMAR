# Forest-Plot-Style Confidence Interval Display

Plots one or more effect size estimates with their confidence intervals.

## Usage

``` r
plot_ci(
  ci = NULL,
  estimate = NULL,
  lower = NULL,
  upper = NULL,
  names = NULL,
  n = NULL,
  conf_level = 0.95,
  show_n = TRUE,
  reference_line = NULL,
  xlab = "Effect Size",
  title = NULL,
  palette = "okabe_ito"
)
```

## Arguments

- ci:

  A `data.frame` from a DMAR `ci_*` function. When supplied, the
  function auto-detects the format and extracts the point estimate(s),
  lower limit(s), and upper limit(s). Explicit `estimate`, `lower`, and
  `upper` arguments override values parsed from `ci`.

- estimate:

  Numeric vector of point estimates.

- lower:

  Numeric vector of lower confidence limits.

- upper:

  Numeric vector of upper confidence limits.

- names:

  Optional character vector of labels for each effect.

- n:

  Optional numeric vector (or scalar) of sample sizes. Recycled to match
  the number of effects.

- conf_level:

  Confidence level of the limits, used only for the title (default
  `0.95`). When `ci` records the level of its interval, as every DMAR
  table does, that level is used, and a `conf_level` that differs from
  it is refused, since the plot draws the limits it is given and cannot
  change their level.

- show_n:

  Logical. If `TRUE` (the default), the sample size is annotated above
  each estimate, with the estimate and its interval printed below.

- reference_line:

  Optional numeric value at which to draw a vertical reference line
  (e.g., `0` for mean differences, `1` for ratios).

- xlab:

  Label for the horizontal (effect size) axis. Defaults to
  `"Effect Size"`.

- title:

  Optional plot title.

- palette:

  Character string naming the color palette. The point estimates and
  interval bars are drawn in the palette's primary color. Defaults to
  `"okabe_ito"`, base R's colorblind-safe Okabe-Ito palette; `"tableau"`
  is also available.

## Value

A `ggplot2` object.

## Details

The function accepts either (a) a `data.frame` produced by an DMAR
`ci_*` function (e.g.,
[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)),
or (b) explicit numeric vectors for the estimate(s), lower bound(s), and
upper bound(s).

The function recognizes three DMAR output formats:

- **Long term/value with estimate row**:

  Output from
  [`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
  which includes a row for the point estimate (e.g., `term = "smd"`) in
  addition to `"ci_lower"` and `"ci_upper"`.

- **Long term/value without estimate**:

  Output from
  [`ci_R`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md)
  or [`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
  which contains only `"ci_lower"` and `"ci_upper"`. Supply the point
  estimate via the `estimate` argument.

- **Wide per-effect format**:

  Output from
  [`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md)
  or
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
  which has one row per effect with columns for each estimate and its
  limits (`omega_squared_partial`, `omega_squared_partial_ci_lower`,
  `omega_squared_partial_ci_upper`, and the same for the total) and `N`,
  or from
  [`ci_eta_squared_generalized`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
  with `ci_lower` and `ci_upper`. Each estimate that has limits is drawn
  (the total of an effect of a factorial design has none), and when a
  table has limits for two, each row is labeled with the effect and the
  estimate.

The default title names the level the limits were formed at. An interval
formed from unequal tail areas (`alpha_lower` and `alpha_upper` of a
`ci_*` function) is drawn with a subtitle that gives the two areas, or,
when one of them is 0, says which limit is the one-sided confidence
bound, as
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
does.

## Note

Requires ggplot2 (listed in `Suggests`).

## See also

[`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_R`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`plot_smd`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
[`plot_R2`](https://yelleknek.github.io/DMAR/reference/plot_R2.md)

Other plotting:
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
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
# From explicit values.
plot_ci(estimate = 0.45, lower = 0.15, upper = 0.75,
        names = "Cohen's d", n = 60, reference_line = 0)


# From ci_smd() output.
ci_result <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
plot_ci(ci_result, n = 100, reference_line = 0)


# Multiple effects from ci_omega_squared(): the expectancy treatment
# and the grade classification in the pygmalion data.
pyg <- pygmalion
pyg$grade <- factor(pyg$grade)
fit <- aov(iq_8 ~ treatment + grade, data = pyg)
omega_result <- ci_omega_squared(fit, estimand = "partial")
plot_ci(omega_result, reference_line = 0,
        xlab = expression(omega^2))

```
