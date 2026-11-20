# Visualize Observed Individual Trajectories in a Longitudinal Data Set

Plots one trajectory per subject from a long-format data frame,
optionally colored by a grouping variable, and optionally faceted into
one panel per subject. Returns a ggplot2 object that can be further
customized.

## Usage

``` r
plot_trajectories(
  data,
  id,
  time,
  outcome,
  group = NULL,
  ids = NULL,
  n_random = NULL,
  pct_random = NULL,
  facet = FALSE,
  nrow = NULL,
  ncol = NULL,
  show_points = TRUE,
  point_size = 1.5,
  linewidth = 0.5,
  alpha = 0.7,
  palette = "okabe_ito",
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  seed = NULL
)
```

## Arguments

- data:

  A long-format `data.frame` (one row per subject-occasion).

- id:

  Character. Column name in `data` identifying the subject.

- time:

  Character. Column name for the time / occasion variable (the *x*
  axis).

- outcome:

  Character. Column name for the outcome / score variable (the *y*
  axis).

- group:

  Optional character. Column name for a grouping variable used to color
  the trajectories (and panels, if faceted).

- ids:

  Optional vector of subject IDs to plot.

- n_random:

  Optional integer; randomly sample this many subjects.

- pct_random:

  Optional numeric; sample this percentage of subjects. Values \\\leq
  1\\ are interpreted as proportions; values \\\> 1\\ as percentages. At
  most one of `ids`, `n_random`, or `pct_random` may be supplied.

- facet:

  Logical. If `TRUE`, draw one panel per subject via
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
  If `FALSE` (default), overlay all trajectories in one panel.

- nrow, ncol:

  Optional integers passed to `facet_wrap()` when `facet = TRUE`.

- show_points:

  Logical. If `TRUE` (default), draw the observed points as well as the
  connecting lines.

- point_size:

  Size of the observed points (default `1.5`).

- linewidth:

  Line width for the connecting segments (default `0.5`).

- alpha:

  Transparency for points and lines (default `0.7`).

- palette:

  Character string naming the color palette used to color the
  trajectories when `group` is a discrete (factor, character, or
  logical) variable. Defaults to `"okabe_ito"`, base R's colorblind-safe
  Okabe-Ito palette; `"tableau"` is also available. Ignored when `group`
  is `NULL` or numeric.

- title, xlab, ylab:

  Optional plot labels. Sensible defaults are taken from `outcome` and
  `time` when these are `NULL`.

- seed:

  Optional integer random seed used when `n_random` or `pct_random` is
  supplied. Defaults to `NULL`, under which the subjects are sampled
  from, and advance, the current random number stream; supply an integer
  for reproducible subject sampling, which is set for the call only, the
  caller's random number state being restored when it returns.

## Value

A `ggplot` object.

## Details

The function reworks the original `vit()` (visualize individual
trajectories) function to return a single ggplot2 object instead of
producing graphical side effects. Saving is handled by the user via
[`ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html);
multi-page output via faceting and
[`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)'s
`nrow`/`ncol`.

## Note

Requires ggplot2 (a `Suggests` dependency).

## See also

[`plot_trajectories_fitted`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md)
for plotting observed trajectories together with a fitted multilevel
model's predictions.

Other plotting:
[`plot_R2()`](https://yelleknek.github.io/DMAR/reference/plot_R2.md),
[`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md),
[`plot_equivalence()`](https://yelleknek.github.io/DMAR/reference/plot_equivalence.md),
[`plot_fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/plot_fa_confirmatory.md),
[`plot_forest()`](https://yelleknek.github.io/DMAR/reference/plot_forest.md),
[`plot_irt_information()`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md),
[`plot_mediation_mbco()`](https://yelleknek.github.io/DMAR/reference/plot_mediation_mbco.md),
[`plot_randomization_test()`](https://yelleknek.github.io/DMAR/reference/plot_randomization_test.md),
[`plot_regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/plot_regions_of_significance.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md),
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md),
[`power_equivalence_md_plot()`](https://yelleknek.github.io/DMAR/reference/power_equivalence_md_plot.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The Orthodont data from nlme: 27 children, 4 measurements each.
d <- nlme::Orthodont

# Overlay all trajectories, colored by sex.
plot_trajectories(d, id = "Subject", time = "age",
                  outcome = "distance", group = "Sex")


# One panel per child, for six children drawn at random, few enough
# that the figure is quick to draw; drop n_random to get a panel for
# every child. The seed makes the draw reproducible, and the session's
# generator state is left as it was.
plot_trajectories(d, id = "Subject", time = "age",
                  outcome = "distance",
                  n_random = 6, facet = TRUE, ncol = 3,
                  seed = 113)

```
