# Simulate Data From a Polynomial Change (Growth) Model

Generates longitudinal data from a random-coefficients polynomial change
model: each subject follows a degree-\\P\\ polynomial in time whose
coefficients vary randomly across subjects, and each measurement adds
independent level-one error. The polynomial order is general (order 0 is
a flat line, 1 a straight line, 2 a quadratic, and so on), one or
several populations may differ in their mean trajectories, the level-one
error can be set directly or pinned to a target measurement reliability,
and the actual time of each assessment may jitter away from its nominal
target. This is the Monte Carlo companion to
[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
which plans power for the same model under the closed-form assumptions
of Raudenbush and Liu (2001); the simulator can relax those assumptions
(notably the assumption of fixed, error-free, equally reliable
assessment times) and study what happens.

## Usage

``` r
simulate_longitudinal_polynomial(
  n,
  target_times = NULL,
  fixed_coefficients,
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

  A single positive integer (equal number of units in every population)
  or a numeric vector of length `G` giving the number of subjects in
  each of the `G` populations, where `G` is the number of mean
  trajectories supplied through `fixed_coefficients`.

- target_times:

  A numeric vector of the nominal (planned) measurement times, length
  \\M\\, one shared schedule for every unit. They need not be equally
  spaced. A degree-\\P\\ model requires \\M \ge P + 1\\ occasions. Give
  either this or `time_range`.

- fixed_coefficients:

  The population mean trajectory, as the coefficients of a polynomial in
  time, ordered from the intercept upward: `c(b0, b1, ..., bP)` encodes
  \\b_0 + b_1 t + b_2 t^2 + \dots + b_P t^P\\. The length sets the
  polynomial order \\P\\ (length 1 is order 0, a flat line at `b0`). For
  several populations, pass a *list* of equal-length coefficient
  vectors, one per population; the populations then differ in their mean
  trajectories but share the variance components below (the
  Raudenbush-Liu two-group setup).

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
  of measurement times is drawn uniformly. Every value must be at least
  \\P + 1\\ so each unit's trajectory identifies the polynomial.

- time_distribution:

  Distribution of the unit-specific times over `time_range`; currently
  `"uniform"`.

- random_variances:

  The between-subject variances of the polynomial coefficients (the
  diagonal of the level-two covariance matrix), as a single number
  recycled to all \\P + 1\\ coefficients or a vector of length \\P +
  1\\. Default `0`, a fixed common trajectory with no between-subject
  heterogeneity. Entries may be zero coefficient by coefficient (e.g., a
  random intercept with a fixed slope).

- random_correlation:

  Optional \\(P + 1) \times (P + 1)\\ correlation matrix among the
  random coefficients. Default `NULL` treats them as uncorrelated.
  Combined with `random_variances` it forms the level-two covariance \\T
  = D R D\\, with \\D\\ the diagonal matrix of coefficient standard
  deviations.

- error_variance:

  The level-one (within-subject) measurement error variance
  \\\sigma^2_e\\. Most simply a single non-negative number (the same
  error variance at every occasion). It may instead be a length-\\M\\
  vector of per-occasion (heteroscedastic) error variances, or a full
  \\M \times M\\ error covariance matrix when the errors are correlated
  across occasions in an arbitrary (unstructured) way. For the common
  structured cases, give a scalar or vector here and set
  `error_structure` and `error_correlation` instead of building the
  matrix by hand. Specify exactly one of `error_variance` or
  `reliability`.

- reliability:

  A target measurement reliability in \\(0, 1)\\ from which
  \\\sigma^2_e\\ is derived (see Details). The value is the *average*
  per-occasion reliability across `target_times`; the per-occasion
  reliabilities, which generally differ from one another, are returned
  in the `"reliability_by_occasion"` attribute. Requires at least one
  positive entry in `random_variances`. The solved error variance is
  homoscedastic and may still be given an across-occasion correlation
  through `error_structure`. Specify exactly one of `error_variance` or
  `reliability`.

- error_structure:

  The correlation pattern of the level-one errors across occasions:
  `"independent"` (the default, uncorrelated errors), `"ar1"` (a
  first-order autoregressive decay \\\rho^{\|j-k\|}\\, so errors at
  occasions closer in time are more alike), `"compound_symmetry"` (a
  constant correlation \\\rho\\ between every pair of occasions), or
  `"toeplitz"` (a banded structure set by the lag-1 through
  lag-\\(M-1)\\ correlations). Ignored when `error_variance` is a full
  covariance matrix, which already fixes the structure.

- error_correlation:

  The correlation parameter(s) for `error_structure`: a single number
  \\\rho\\ for `"ar1"` and `"compound_symmetry"`, or a vector of the
  lag-1 to lag-\\(M-1)\\ correlations for `"toeplitz"`. Left `NULL` for
  `"independent"`. The implied correlation matrix must be positive
  semidefinite.

- timing_sd:

  The standard deviation of the difference between a subject's actual
  and nominal assessment time, as a single number recycled to all
  occasions or a vector of length \\M\\. Default `0`, every subject
  measured exactly on schedule. A positive value draws each subject's
  actual time at occasion \\m\\ as \\\tau_m + N(0,
  \text{timing\\sd}\_m^2)\\ and evaluates that subject's true score at
  the actual time, while the nominal target is retained in a separate
  column (see Details).

## Value

A long-format `data.frame` with one row per subject-occasion and the
columns

- `id`:

  Factor uniquely identifying each subject.

- `population`:

  Factor with `G` levels (`"1"`, ...) giving each unit's population (its
  data generating parameter vector). With one parameter vector there is
  one level.

- `occasion`:

  Integer occasion index, `1` to \\M\\.

- `target_time`:

  The nominal (planned) measurement time.

- `time`:

  The actual measurement time (equal to `target_time` when
  `timing_sd = 0`, otherwise jittered).

- `true_score`:

  The subject's latent trajectory value at the actual time, before
  level-one error.

- `y`:

  The observed score, `true_score` plus level-one error.

The returned object carries attributes `"error_variance"` (the
\\\sigma^2_e\\ used, a scalar when the errors are homoscedastic and
independent, otherwise the vector of per-occasion error variances),
`"error_covariance"` (the full \\M \times M\\ level-one error covariance
actually used), `"reliability_by_occasion"` (the per-occasion
reliabilities at the nominal times), `"random_covariance"` (the
level-two covariance \\T\\), `"polynomial_order"` (\\P\\), and
`"schedule"` (`"shared"` or `"unit_specific"`). With `time_range` there
is no shared occasion grid, so `"error_covariance"` is `NA` and
`"reliability_by_occasion"` is `NA`. The format is directly usable with
[`plot_trajectories`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md)
and with mixed-model fitters such as
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) or
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).

## Details

**The model.** Subject \\i\\ in population \\g\\ has a random
coefficient vector \\\pi_i = (\pi\_{i0}, \dots, \pi\_{iP})\\ drawn from
a multivariate normal with mean the population's `fixed_coefficients`
\\\beta_g\\ and covariance \\T\\. The latent trajectory is the
polynomial \\\mu_i(t) = \sum\_{k=0}^{P} \pi\_{ik}\\ t^k\\, and the
observed score at a measurement time \\t\\ is \\y = \mu_i(t) + e\\, with
\\e \sim N(0, \sigma^2_e)\\ independent across occasions. Order 0
collapses to a flat line \\\mu_i(t) = \pi\_{i0}\\; order 1 is the
straight-line growth model underlying Raudenbush and Liu (2001) and
[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md).

**Coefficient metric.** The coefficients here are the ordinary (raw)
polynomial coefficients on \\t^k\\, which is the most transparent metric
for specifying a trajectory. The derivative-scaled change coefficient
used by
[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)
and the Raudenbush-Liu power formulas is \\P!\\ times the leading
(highest-order) coefficient supplied here, so a quadratic with
`fixed_coefficients = c(b0, b1, b2)` corresponds to a Raudenbush-Liu
quadratic change coefficient of \\2! \\ b_2 = 2 b_2\\.

**Reliability varies by occasion.** At a measurement time \\t\\ the
implied between-subject (true-score) variance is the quadratic form
\\c(t)^\top T\\ c(t)\\ with \\c(t) = (1, t, t^2, \dots, t^P)^\top\\, so
the classical reliability of the observed score, \$\$\rho\_{XX}(t) =
\frac{c(t)^\top T\\ c(t)}{c(t)^\top T\\ c(t) + \sigma^2_e},\$\$
generally *changes from occasion to occasion*: a growth measurement is
not equally reliable everywhere, because the spread of true scores
depends on where in time you measure relative to the centering of the
polynomial and the random-effect structure. When `reliability` is
supplied, \\\sigma^2_e\\ is solved (by
[`uniroot`](https://rdrr.io/r/stats/uniroot.html)) so that the average
of \\\rho\_{XX}(t)\\ over the nominal `target_times` equals the
requested value; the occasion-by-occasion reliabilities are returned in
the `"reliability_by_occasion"` attribute so the variation is visible
rather than hidden behind a single number. Reliability is only
meaningful when there is true-score variance to detect, so this route
requires `random_variances` to be positive for at least one coefficient.

**Measurement errors need not be independent or equal.** The simplest
model adds an independent, equal-variance error at every occasion, but
repeated measurements of the same person are often correlated (an
unmodeled state, a rater, or an instrument carries over from one wave to
the next) and may be more or less variable at different waves. The
level-one errors are drawn from \\N(0, \Sigma_e)\\, and \\\Sigma_e\\ can
be set three ways: a scalar or per-occasion `error_variance` combined
with an `error_structure` (`"ar1"` for autoregressive decay, the natural
choice when occasions are ordered in time; `"compound_symmetry"` for an
equicorrelated error; `"toeplitz"` for a general banded pattern), or a
full covariance matrix passed directly as `error_variance`. Because
classical reliability at an occasion is a marginal quantity, it depends
only on the diagonal of \\\Sigma_e\\; the across-occasion error
correlation leaves `"reliability_by_occasion"` unchanged but does affect
how a mixed model that assumes independent errors performs, which is
exactly the kind of misspecification this simulator is meant to let a
user study.

**Assessment timing is rarely exact.** Designs are written as if every
subject is measured at the same fixed times (“the 7-day follow-up”), but
in practice people arrive early or late, so the actual time differs from
the nominal target. Setting `timing_sd > 0` draws each subject's actual
time per occasion and evaluates the true score *at the time the
measurement really happened*, while `target_time` keeps the nominal
value an analyst would typically use. Analyzing on the nominal time when
the data were in fact collected on jittered times biases estimates of
the change coefficients, and the bias grows with the order of the trend
and with the size of the timing variability. The two time columns let a
user quantify that bias by fitting the same model on `time` versus
`target_time`.

**Why the closed-form Raudenbush-Liu planner does not cover all of
this.** The power formulas in Raudenbush and Liu (2001), carried by
[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
are exact under three assumptions this simulator can relax: every
subject is measured at the *same*, equally spaced, *error-free* occasion
times; the level-one error variance is a single constant (so the closed
form needs no notion of an occasion-varying reliability); and the
within-subject sampling variance of the change coefficient has the known
form \\V = \sigma^2_e f^{2p} (M - p - 1)! / \[K_p (M + p)!\]\\. Those
assumptions buy a clean formula, but real designs violate them:
assessments drift in time, and reliability is not the same at every
wave. This function is the Monte Carlo complement that lets a researcher
generate data under the messier reality and check how far the
closed-form power and the fitted estimates can be trusted.

## References

Kelley, K., & Rausch, J. R. (2011). Sample size planning for
longitudinal models: Accuracy in parameter estimation for polynomial
change parameters. *Psychological Methods, 16*(4), 391–405.
[doi:10.1037/a0023352](https://doi.org/10.1037/a0023352)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 15 on the analysis of repeated measures
and growth.)

Raudenbush, S. W., & Liu, X.-F. (2001). Effects of study duration,
frequency of observation, and sample size on power in studies of group
differences in polynomial change. *Psychological Methods, 6*(4),
387–401.
[doi:10.1037/1082-989X.6.4.387](https://doi.org/10.1037/1082-989X.6.4.387)

## See also

[`ss_power_pcm`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md)
for closed-form power planning on the same model,
[`plot_trajectories`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md)
to visualize the simulated curves, and
[`mvrnorm`](https://rdrr.io/pkg/MASS/man/mvrnorm.html) for the
random-coefficient draw.

Other data simulators:
[`simulate_ancova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_data.md),
[`simulate_ancova_factorial_data()`](https://yelleknek.github.io/DMAR/reference/simulate_ancova_factorial_data.md),
[`simulate_anova_data()`](https://yelleknek.github.io/DMAR/reference/simulate_anova_data.md),
[`simulate_longitudinal_gompertz()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_gompertz.md),
[`simulate_longitudinal_logistic()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_logistic.md),
[`simulate_longitudinal_negative_exponential()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_negative_exponential.md),
[`simulate_longitudinal_richards()`](https://yelleknek.github.io/DMAR/reference/simulate_longitudinal_richards.md),
[`simulate_regression_data()`](https://yelleknek.github.io/DMAR/reference/simulate_regression_data.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. One population of linear growers with a random intercept and a random slope,
#    measured yearly for four years (five occasions), with the level-one
#    error set directly.
set.seed(113)
d <- simulate_longitudinal_polynomial(
  n                  = 50,
  target_times       = 0:4,
  fixed_coefficients = c(10, 1.5),          # intercept 10, slope 1.5 per year
  random_variances   = c(4, 0.25),          # var(intercept) = 4, var(slope) = .25
  error_variance     = 1
)
head(d)
#>   id population occasion target_time time true_score         y
#> 1  1          1        1           0    0   9.733291  9.084685
#> 2  1          1        2           1    1  11.773957 11.949644
#> 3  1          1        3           2    2  13.814622 13.610660
#> 4  1          1        4           3    3  15.855288 15.918206
#> 5  1          1        5           4    4  17.895953 17.911464
#> 6  2          1        1           0    0   7.249557  8.428662

# 2. Two populations that differ only in their slope (a treatment that changes
#    the rate of growth). Pass a list of coefficient vectors, one per population.
set.seed(113)
two <- simulate_longitudinal_polynomial(
  n                  = c(40, 40),
  target_times       = 0:4,
  fixed_coefficients = list(control = c(10, 1.0), treatment = c(10, 1.8)),
  random_variances   = c(4, 0.25),
  error_variance     = 1
)
aggregate(y ~ population + occasion, data = two, FUN = mean)
#>    population occasion        y
#> 1           1        1  9.77354
#> 2           2        1 10.34420
#> 3           1        2 10.66252
#> 4           2        2 12.11865
#> 5           1        3 11.38262
#> 6           2        3 14.02571
#> 7           1        4 12.39497
#> 8           2        4 15.66385
#> 9           1        5 13.34494
#> 10          2        5 17.31383

# 3. Pin the level-one error to a target reliability instead of setting it
#    directly. The single number is the average reliability across occasions;
#    the per-occasion values differ and are returned as an attribute.
set.seed(113)
rel <- simulate_longitudinal_polynomial(
  n                  = 100,
  target_times       = 0:4,
  fixed_coefficients = c(10, 1.5),
  random_variances   = c(4, 0.25),
  reliability        = 0.80
)
attr(rel, "error_variance")
#> [1] 1.303557
round(attr(rel, "reliability_by_occasion"), 3)   # not constant across waves
#> occasion_1 occasion_2 occasion_3 occasion_4 occasion_5 
#>      0.754      0.765      0.793      0.827      0.860 

# 4. Assessment-time jitter: the nominal "yearly" schedule, but subjects
#    actually arrive a little early or late (SD of about six weeks on a
#    one-year scale). The nominal and actual times are kept in separate
#    columns so the consequences of analyzing on the nominal time can be
#    studied.
set.seed(113)
jit <- simulate_longitudinal_polynomial(
  n                  = 30,
  target_times       = 0:4,
  fixed_coefficients = c(10, 1.5),
  random_variances   = c(4, 0.25),
  error_variance     = 1,
  timing_sd          = 0.12
)
head(jit[, c("id", "occasion", "target_time", "time")])
#>   id occasion target_time        time
#> 1  1        1           0  0.16252055
#> 2  1        2           1  0.93632062
#> 3  1        3           2  2.08847793
#> 4  1        4           3  2.70577279
#> 5  1        5           4  4.00090495
#> 6  2        1           0 -0.02811366

# 4b. Unit-specific measurement times: each of 12 children is tested
#     between 40 and 90 weeks of age, five to nine times, no two on the
#     same schedule. The level-one error is a single variance; the
#     "schedule" attribute records the design.
set.seed(113)
ages <- simulate_longitudinal_polynomial(
  n                  = 12,
  time_range         = c(40, 90),
  occasions          = c(5, 9),
  fixed_coefficients = c(10, 0.5),
  random_variances   = c(4, 0.01),
  error_variance     = 2
)
attr(ages, "schedule")
#> [1] "unit_specific"
table(table(ages$id))   # units per occasion count
#> 
#> 5 6 7 8 9 
#> 4 1 2 4 1 

# 5. A flat line (order 0): no growth, only a random subject level and
#    measurement error. The coefficient vector has length one.
set.seed(113)
flat <- simulate_longitudinal_polynomial(
  n                  = 20,
  target_times       = 0:4,
  fixed_coefficients = 5,
  random_variances   = 2,
  error_variance     = 1
)
head(flat)
#>   id population occasion target_time time true_score        y
#> 1  1          1        1           0    0   5.188592 5.471883
#> 2  1          1        2           1    1   5.188592 6.936332
#> 3  1          1        3           2    2   5.188592 6.385949
#> 4  1          1        4           3    3   5.188592 4.852617
#> 5  1          1        5           4    4   5.188592 6.085924
#> 6  2          1        1           0    0   6.944857 8.731318

# 6. Autocorrelated measurement error: the same error variance at each wave,
#    but the level-one errors decay as an AR(1) process (errors at adjacent
#    occasions correlate 0.5), the kind of dependence a model assuming
#    independent errors would miss. The full error covariance is returned.
set.seed(113)
ar <- simulate_longitudinal_polynomial(
  n                  = 40,
  target_times       = 0:4,
  fixed_coefficients = c(10, 1.5),
  random_variances   = c(4, 0.25),
  error_variance     = 1,
  error_structure    = "ar1",
  error_correlation  = 0.5
)
round(attr(ar, "error_covariance"), 3)
#>       [,1]  [,2] [,3]  [,4]  [,5]
#> [1,] 1.000 0.500 0.25 0.125 0.062
#> [2,] 0.500 1.000 0.50 0.250 0.125
#> [3,] 0.250 0.500 1.00 0.500 0.250
#> [4,] 0.125 0.250 0.50 1.000 0.500
#> [5,] 0.062 0.125 0.25 0.500 1.000
```
