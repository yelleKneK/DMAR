# Plot the Randomization Distribution Behind a Randomization Test

Displays the reference distribution that
[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md)
built by reassigning the observed scores to the two groups, with the
observed statistic marked and every reassignment at least as extreme as
the observed one shaded. The shaded proportion is the *p*-value, so the
figure shows where that number came from instead of only reporting it.

## Usage

``` r
plot_randomization_test(object, bins = 40L, palette = "okabe_ito", ...)
```

## Arguments

- object:

  A result of
  [`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md).

- bins:

  Number of histogram bins used to display the reference distribution.
  Defaults to `40`.

- palette:

  Character; the color palette. Defaults to `"okabe_ito"`, base R's
  colorblind-safe Okabe-Ito palette; `"tableau"` is also available.

- ...:

  Currently unused; present so the signature can grow without breaking
  existing calls.

## Value

A `ggplot` object, which can be printed or further modified with the
usual ggplot2 verbs.

## Details

Reading the figure is the point of it. The spread of the distribution is
what the reassignments alone can produce when the grouping is
irrelevant, which is the null hypothesis of the test. If the observed
statistic sits inside that spread, reassignment alone explains it. If it
sits out in a tail, few reassignments reproduce it, and that scarcity is
the evidence. No normal or *t* distribution appears anywhere in the
construction. This is the display Chapter 1 of Maxwell, Delaney, and
Kelley (2027) uses to introduce the logic of the randomization test.

## References

Fisher, R. A. (1935). *The design of experiments*. Oliver & Boyd.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 1 on the logic of the randomization test.)

## See also

[`randomization_test`](https://yelleknek.github.io/DMAR/reference/randomization_test.md)
for the test itself.

Other plotting:
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
[`plot_cfa_k()`](https://yelleknek.github.io/DMAR/reference/plot_cfa_k.md),
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md),
[`plot_irt_information()`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md),
[`plot_mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md),
[`plot_regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
[`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)

## Author

Ken Kelley

## Examples

``` r
treatment <- c(80, 84, 79, 88, 83)
control   <- c(72, 75, 68, 81, 74)
rt <- randomization_test(group_1 = treatment, group_2 = control)
plot_randomization_test(rt)

```
