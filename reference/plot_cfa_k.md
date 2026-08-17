# Plot the Estimates of a Multiple-Factor CFA

Displays the item-level estimates of a
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) fit, one
panel per factor, with each estimate's confidence interval. The display
is built to make the equality questions behind the classical measurement
structures visible: a dashed vertical line marks, per factor, either the
common (equated) estimate when the plotted parameter was constrained
equal, or the mean of the free estimates as an informal anchor for the
question "could these plausibly be one value?". Confidence intervals
that all cover the anchor are what equal loadings (or equal error
variances, or equal intercepts) would look like; an interval far from it
shows which item resists the constraint, and the likelihood ratio test
of the two nested
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) fits is
the formal companion (see the examples in
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md)).

## Usage

``` r
plot_cfa_k(
  x,
  what = c("loadings", "errors", "intercepts"),
  show_equal_reference = TRUE,
  xlab = NULL,
  title = NULL,
  palette = "okabe_ito"
)
```

## Arguments

- x:

  A `dmar_cfa_k` object from
  [`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) with
  the default `output = "verbose"`.

- what:

  Which parameter to display: `"loadings"` (default, the `lambda`
  terms), `"errors"` (the `psi` terms), or `"intercepts"` (the `nu`
  terms; requires a fit with the mean structure).

- show_equal_reference:

  Logical. If `TRUE` (default), draw the dashed per-factor reference
  line described above. When the parameter was constrained equal the
  line is the common estimate and is always drawn.

- xlab:

  Label for the horizontal axis. Defaults to a description of the
  plotted parameter.

- title:

  Optional plot title.

- palette:

  Character string naming the color palette. Defaults to `"okabe_ito"`,
  base R's colorblind-safe Okabe-Ito palette; `"tableau"` is also
  available.

## Value

A `ggplot2` object.

## Note

Requires ggplot2 (listed in `Suggests`).

## See also

[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) for the
fit; [`plot_ci`](https://yelleknek.github.io/DMAR/reference/plot_ci.md)
for the general forest-style confidence interval display.

Other plotting:
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
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
data(holzinger_swineford)
hs_factors <- list(
  verbal = c("t6_paragraph_comprehension",
             "t7_sentence", "t9_word_meaning"),
  deduction = c("t20_deduction", "t22_problem_reasoning",
             "t23_series_completion"))
res <- cfa_k(holzinger_swineford, hs_factors)

# Are equal loadings plausible? Compare each interval with the anchor.
plot_cfa_k(res)


# Two further displays are shown but not run here, since each draws
# another figure and the second refits the model as well. The same
# question for the error variances, the additional constraint that
# separates essentially parallel from essentially tau-equivalent:
# plot_cfa_k(res, what = "errors")
#
# After imposing the constraint, every item in a factor sits at the
# common estimate and the dashed line is that estimate rather than
# the mean of the free ones:
# res_equal <- cfa_k(holzinger_swineford, hs_factors,
#                    equal_loading = TRUE)
# plot_cfa_k(res_equal)
```
