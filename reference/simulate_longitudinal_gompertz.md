# Simulate Data From a Gompertz Change (Growth) Model

Generates longitudinal data from a random-coefficients Gompertz change
model: each unit (a person, an animal, a tree) follows an S-shaped curve
that, unlike the logistic, is not symmetric about its point of
inflection, the parameters vary randomly across units, and each
measurement adds level-one error. The deterministic part of the four
parameter Gompertz curve is \$\$\mu(t) = \alpha \exp\left(-\exp(-\gamma
(t - \beta))\right) + \zeta,\$\$ the parameterization of Kelley (2005,
2008), which generalizes the three parameter Gompertz of the literature
(Winsor, 1932; Ratkowsky, 1983) by adding \\\zeta\\ so the lower
asymptote is itself a modeled quantity rather than fixed at zero.

## Usage

``` r
simulate_longitudinal_gompertz(
  n,
  target_times = NULL,
  fixed_parameters,
  time_range = NULL,
  occasions = NULL,
  time_distribution = "uniform",
  random_variances = 0,
  random_correlation = NULL,
  error_variance = NULL,
  reliability = NULL,
  error_structure = c("independent", "ar1", "compound_symmetry", "toeplitz"),
  error_correlation = NULL,
  timing_sd = 0
)
```

## Arguments

- n:

  A single positive integer, the number of units (persons, animals,
  trees, classrooms) whose trajectories are drawn, or a vector giving
  the number of units for each population, one entry per parameter
  vector in `fixed_parameters`.

- target_times:

  Numeric vector of the nominal measurement times, one shared schedule
  for every unit. Give either this or `time_range`.

- fixed_parameters:

  The population parameters `c(alpha = , beta = , gamma = , zeta = )`
  (an unnamed length-4 vector is taken in that order), or a list of such
  vectors, one per population (distinct data generating parameter
  vectors):

  `alpha`

  :   The total change: the curve travels `alpha` units from the lower
      asymptote `zeta` to the upper asymptote \\\alpha + \zeta\\ (for
      \\\gamma \> 0\\).

  `beta`

  :   The point of inflection on the time axis. The Gompertz inflection
      is early and asymmetric: it occurs where the curve passes through
      \\\alpha / e + \zeta\\, about 36.8% of the total change, so growth
      accelerates briefly and decelerates over a long approach to the
      ceiling.

  `gamma`

  :   The curvature: how sharply the curve rises through its inflection.
      Negative values flip the curve to decreasing.

  `zeta`

  :   The lower asymptote (for \\\gamma \> 0\\). The intercept is \\\phi
      = \alpha \exp(-\exp(\gamma \beta)) + \zeta\\.

- time_range:

  Alternative to `target_times`: `c(lower, upper)` bounds from which
  each unit draws its own measurement times, so no two units share a
  schedule (e.g., age in weeks at testing rather than a fixed grade).
  Requires `occasions`; with `time_range`, the level-one error must be a
  single `error_variance` with the default independent structure, and
  `timing_sd` does not apply.

- occasions:

  With `time_range`: a single positive integer (every unit measured the
  same number of times) or `c(min, max)`, from which each unit's number
  of measurement times is drawn uniformly.

- time_distribution:

  Distribution of the unit-specific times over `time_range`; currently
  `"uniform"`.

- random_variances:

  Between-unit variances of `(alpha, beta, gamma, zeta)`, a single
  number recycled to all four or a length-4 vector. Default `0`. A named
  vector over any subset of the parameter names (e.g., `c(beta = 1.2)`)
  varies only those named and leaves the rest fixed.

- random_correlation:

  Optional 4-by-4 correlation matrix among the random parameters;
  default uncorrelated.

- error_variance, reliability, error_structure, error_correlation,
  timing_sd:

  The level-one error and assessment-time machinery, with the same
  meaning as in
  [`simulate_longitudinal_polynomial`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md):
  specify exactly one of `error_variance` or `reliability` (solved here
  through the first-order delta method true-score variance; because the
  curve is nonlinear in its parameters, that approximation, and the
  `reliability_by_occasion` attribute with it, can drift from the
  realized variance ratio when the random variances are large relative
  to the mean curve); `error_structure` and `error_correlation` set the
  across-occasion error correlation; `timing_sd` jitters the actual
  assessment times.

## Value

A long-format `data.frame` with columns `id`, `population`, `occasion`,
`target_time`, `time`, `true_score`, and `y`, directly usable with
[`plot_trajectories`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md)
and nonlinear mixed-model fitters such as
[`nlme::nlme()`](https://rdrr.io/pkg/nlme/man/nlme.html). Attributes
carry the `model`, `fixed_parameters`, `random_covariance`,
`error_variance`, `error_covariance`, `reliability_by_occasion`, and
`schedule` (`"shared"` or `"unit_specific"`).

## Details

The choice between the Gompertz and the logistic is substantive, not
cosmetic: both are S-shaped, but the logistic spends equal time
approaching floor and ceiling while the Gompertz commits to an early
inflection (36.8% of total change) followed by a long deceleration.
Processes with rapid early gains and slow consolidation, common in
learning and development, are natural Gompertz candidates. The Gompertz
is the \\\delta \rightarrow 0\\ limit of the Richards curve
([`simulate_longitudinal_richards`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md)),
which frees the inflection entirely (Kelley, 2005, 2008).

## References

Kelley, K. (2005). *Estimating nonlinear change models in heterogeneous
populations when class membership is unknown: Defining and developing
the latent classification differential change model* (Doctoral
dissertation). University of Notre Dame.

Kelley, K. (2008). Nonlinear change models in populations with
unobserved heterogeneity. *Methodology, 4*(3), 97–112.

Ratkowsky, D. A. (1983). *Nonlinear regression modeling: A unified
practical approach*. Marcel Dekker.

Winsor, C. P. (1932). The Gompertz curve as a growth curve. *Proceedings
of the National Academy of Sciences, 18*(1), 1–8.

## See also

[`simulate_longitudinal_logistic`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_logistic.md)
for the symmetric sibling;
[`simulate_longitudinal_richards`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md)
for the family that subsumes both;
[`simulate_longitudinal_negative_exponential`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_negative_exponential.md);
[`simulate_longitudinal_polynomial`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md);
[`plot_trajectories`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md).

Other data simulators:
[`simulate_ancova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
[`simulate_ancova_factorial_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_factorial_data.md),
[`simulate_anova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
[`simulate_longitudinal_logistic()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_logistic.md),
[`simulate_longitudinal_negative_exponential()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_negative_exponential.md),
[`simulate_longitudinal_polynomial()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_polynomial.md),
[`simulate_longitudinal_richards()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md),
[`simulate_regression_data()`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The six-curve illustration from Kelley (2005): alpha = 0.75 and
# zeta = 0.25 throughout, so every curve crosses its inflection at
# the same height, 0.75 / exp(1) + 0.25, about 0.526. The three
# rising curves share the inflection time beta = 2 and differ only
# in curvature; the three falling curves share beta = 3. Each curve
# is its own population of size one, so a single call draws the whole panel.
panel <- simulate_longitudinal_gompertz(
  n = 1, target_times = seq(0, 6, by = 0.1),
  fixed_parameters = list(
    c(alpha = 0.75, beta = 2, gamma =  1.75, zeta = 0.25),
    c(alpha = 0.75, beta = 2, gamma =  1.00, zeta = 0.25),
    c(alpha = 0.75, beta = 2, gamma =  0.45, zeta = 0.25),
    c(alpha = 0.75, beta = 3, gamma = -0.35, zeta = 0.25),
    c(alpha = 0.75, beta = 3, gamma = -0.60, zeta = 0.25),
    c(alpha = 0.75, beta = 3, gamma = -2.00, zeta = 0.25)),
  error_variance = 0
)
plot_trajectories(panel, id = "id", time = "time",
                  outcome = "true_score", group = "population")


# Individual differences in a single parameter: only the inflection
# time varies (a named entry leaves every other variance at zero),
# so every trajectory shares the floor and the ceiling but reaches
# its fastest growth at its own moment.
set.seed(113)
d_beta <- simulate_longitudinal_gompertz(
  n = 25, target_times = seq(0, 8, by = 0.5),
  fixed_parameters = c(alpha = 75, beta = 3, gamma = 0.55, zeta = 10),
  random_variances = c(beta = 0.8), error_variance = 0
)
plot_trajectories(d_beta, id = "id", time = "time",
                  outcome = "true_score")


# Individual differences in every parameter at once, plus level-one
# error: the realistic sampling model.
set.seed(113)
d <- simulate_longitudinal_gompertz(
  n = 30, target_times = 0:12,
  fixed_parameters = c(alpha = 75, beta = 3, gamma = 0.55, zeta = 10),
  random_variances = c(alpha = 36, beta = 0.8, gamma = 0.01, zeta = 9),
  error_variance = 16
)
plot_trajectories(d, id = "id", time = "time", outcome = "y")

```
