# Sensitivity Analysis for Sample Size Planning From the Accuracy in Parameter Estimation Perspective for the Coefficient of Variation

Simulates a sample size plan for the coefficient of variation to show
what it delivers when the population value differs from the planning
value, or at a sample size chosen in advance. It takes the population
coefficient of variation, either a planning value or a sample size, and
a target width, and returns summaries of the estimate, the interval
width, the proportion of intervals no wider than the target, and
coverage over `G` simulated normal samples (Kelley, 2007; Maxwell,
Delaney, & Kelley, 2027, Section 3.11 on power-based sample size
planning).

## Usage

``` r
ss_aipe_cv_sensitivity(
  true_cv = NULL,
  supposed_cv = NULL,
  width = NULL,
  assurance = NULL,
  true_mean = 100,
  specified_N = NULL,
  conf_level = 0.95,
  G = 1000,
  print_iter = FALSE,
  filename = NULL,
  seed = NULL
)
```

## Arguments

- true_cv:

  Population coefficient of variation (the data generating value); a
  single positive number

- supposed_cv:

  The supposed coefficient of variation, the planning value from which
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
  plans the study (a value the researcher supposes, not an estimate
  computed from data); must be positive. Supply this or `specified_N`
  but not both.

- width:

  Desired (full) width of the two-sided confidence interval for the
  population coefficient of variation; a single positive number. With
  `specified_N` it may be left `NULL`, and `pct_ci_less_width` is then
  `NaN`.

- assurance:

  Probability with which the realized interval should be no wider than
  `width`, passed to
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md):
  `NULL` (plan for the expected width) or a single number strictly
  between 0 and 1; see
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
  for how the procedure behaves near .50. With `specified_N` it is only
  echoed.

- true_mean:

  Population mean used by the simulator to generate data, a single
  positive number (the standard deviation is determined by `true_mean`
  and `true_cv`, since CV = sigma/mu). Default 100.

- specified_N:

  Pre-specified sample size to evaluate, a whole number of at least 4
  (use this when you want the sensitivity results at a fixed *N* rather
  than at the *N* that
  [`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
  would recommend); incompatible with `supposed_cv`.

- conf_level:

  Desired confidence level (i.e., 1 - Type I error rate); default 0.95

- G:

  Number of Monte Carlo replications; defaults to 1000. Increase (e.g.,
  5000 or 10000) for stable Type I error estimates.

- print_iter:

  Logical. If `TRUE` the simulation prints the iteration index after
  each replication (helpful for long runs); default `FALSE`.

- filename:

  An optional path for a comma separated file recording every
  replication (the two confidence limits, the realized coefficient of
  variation, a coverage indicator, and the interval width): nothing is
  written when `filename` is `NULL` (the default), a new file with a
  header row is created otherwise, an existing file at that path is
  appended to, and a throwaway run should point it at
  `tempfile(fileext = ".csv")`.

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the run reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the Monte
Carlo results across the `G` replications. The `term` entries are:
`"mean_cv"`, `"median_cv"`, `"sd_cv"` (mean / median / SD of the
observed sample coefficients of variation); `"mean_ci_width"`,
`"median_ci_width"`, `"sd_ci_width"` (corresponding summaries of the
realized interval widths); `"pct_ci_less_width"` (proportion of
intervals at or below the planning width `width`); `"pct_ci_miss_low"`
and `"pct_ci_miss_high"` (tail-specific non-coverage);
`"total_type_I_error"` (overall empirical non-coverage of `true_cv`);
plus the input echoes `"total_N"` (the sample size evaluated),
`"true_cv"`, `"supposed_cv"` (NA when `specified_N` was supplied
instead), `"width"`, `"conf_level"`, and `"assurance"` (present only
when an assurance was supplied). The proportion rows are on the 0 to 1
scale, not percentages, so `total_type_I_error` is the sum of
`pct_ci_miss_low` and `pct_ci_miss_high`.

**Which replications each row describes.** Three rules decide it.

1.  A replication whose sample mean is not positive has no confidence
    interval for the coefficient of variation
    ([`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)
    requires a positive mean), and its sample coefficient is negative or
    undefined; this can happen when `true_cv` is large and the sample is
    small. All ten statistic rows, `mean_cv`, `median_cv`, and `sd_cv`
    included, describe the replications with a positive mean, so that
    every row is about the same samples. Such a replication is recorded
    with its coefficient and missing limits in the file `filename`
    names. One warning of class `dmar_replications_without_interval`
    states how many there were, and the count is
    `attr(, "simulation")$failed`.

2.  An interval can have an infinite upper limit (see
    [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)).
    Such an interval has infinite width and counts as not meeting the
    width target in `pct_ci_less_width`. `mean_ci_width` and
    `sd_ci_width` are over the finite widths, since one infinite width
    would make both uninformative, and `median_ci_width` is over all
    widths, infinite ones included. The number of infinite widths is
    `attr(, "simulation")$n_infinite_width`.

3.  An interval with an infinite upper limit is a valid interval of
    [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
    which signals nothing for it. One warning of class
    `dmar_replications_with_infinite_limit` states how many replications
    had one and how they enter the rows; the condition carries the count
    as `n_replications` and the number of replications as `G`.

The attribute `simulation` (see the section on the record of the
simulation) has one element more here, `n_infinite_width`.

## Details

Sample size planning for the coefficient of variation under the Accuracy
in Parameter Estimation framework chooses *N* so that the expected (or,
with assurance, the high-probability) confidence interval width is no
larger than `width` (Kelley, 2007). Because the procedure assumes the
planning value `supposed_cv` matches the population value `true_cv`, in
practice the realized width will deviate from the planned width whenever
the planning value is wrong. This sensitivity analysis quantifies the
deviation by Monte Carlo simulation: the planned *N* is obtained from
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
with `supposed_cv`, then samples are drawn from the *true* population
(with coefficient of variation `true_cv`) and the realized confidence
interval widths are summarized.

For planning a sample size so that a confidence interval is sufficiently
narrow, see Maxwell, Delaney, & Kelley (2027, Section 3.9.1.5, and
Section 4.9 on assurance).

The samples are drawn from a normal population, which the interval and
the plan assume, so the function does not show what a plan delivers for
a population that is not normal. With a skewed population the coverage
of the interval departs from the stated level: in 40,000 samples at each
of the three plans for a coefficient of variation of .50, a width of
.20, and 95 percent confidence (79 without `assurance`, 97 with an
assurance of .80, and 142 with .99), coverage was .97 for a gamma
population with that coefficient of variation, and for a lognormal one
it was .915, .913, and .905, lowest at the largest plan. The proportion
of intervals no wider than `width` was .55, .88, and 1.00 for the gamma
population and .59, .84, and .99 for the lognormal one. For a normal
population that probability, computed from the noncentral *t*
distribution, is .53, .84, and 1.00 (see
[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md)
and [`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md)).

## The record of the simulation

Every sensitivity function of the package summarizes its replications
over one denominator: the replications that produced a confidence
interval (for a power sensitivity function, a test). The proportion of
intervals that met the width target (`pct_ci_less_width`), the
proportions of intervals that missed the true value, and the summaries
of the widths are all taken over those replications. A replication
without an interval is in none of these rows. It is counted, and the
returned table carries the counts in its attribute `simulation`, a list
read with `attr(x, "simulation")`:

- `attempted`:

  the number of replications asked for (`G`).

- `returned`:

  the number that produced an interval (a test).

- `failed`:

  the number that did not (`attempted - returned`).

- `denominator`:

  the number of replications the rows named above are taken over, which
  is `returned`.

- `width_target_probability`:

  the proportion of the replications with an interval whose interval met
  the width target, the value of the `pct_ci_less_width` row; `NA` for a
  function or a call without a width target, and when no replication
  produced an interval.

- `width_target_mc_ci`:

  the exact binomial (Clopper-Pearson) limits for that proportion, at
  the confidence level of the table. They describe the error of the
  simulation, which shrinks as `G` grows, and nothing else.

A function that records more about its run adds elements to the list,
and its help page names them. The line "Replications" beneath the
printed table gives `G`, and the three counts when a replication was
without an interval. With many replications without an interval the rows
describe only the studies in which the interval could be formed, so
`failed` is to be reported beside them.

## References

Chattopadhyay, B., & Kelley, K. (2016). Estimation of the coefficient of
variation with minimum risk: A sequential method for minimizing sampling
error and study cost. *Multivariate Behavioral Research, 51*(5),
627–648.
[doi:10.1080/00273171.2016.1203279](https://doi.org/10.1080/00273171.2016.1203279)

Kelley, K. (2007). Sample size planning for the coefficient of variation
from the accuracy in parameter estimation approach. *Behavior Research
Methods, 39*(4), 755–766.
[doi:10.3758/BF03192966](https://doi.org/10.3758/BF03192966)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on confidence intervals for effect
sizes.)

## See also

[`ss_aipe_cv`](https://yelleknek.github.io/DMAR/reference/ss_aipe_cv.md),
[`ci_cv`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`cv`](https://yelleknek.github.io/DMAR/reference/cv.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Textbook scenario (Kelley, 2007). A researcher plans a study to estimate the
# coefficient of variation of reaction times with a 95% confidence interval no
# wider than .10. They guess from prior work that the population coefficient of
# variation is around .25 and apply ss_aipe_cv() to obtain a planned N.

# Question 1: how does the realized interval width behave when the planning
# value is correct (well-specified case)? G = 100 keeps the example quick;
# a reported sensitivity analysis deserves the default G = 1000 or more.
set.seed(113)
ss_aipe_cv_sensitivity(
  true_cv          = .25,
  supposed_cv      = .25,
  width            = .10,
  assurance        = NULL,
  conf_level       = .95,
  G                = 100,
  print_iter       = FALSE
)
#>  term               value
#>  mean_cv            0.246
#>  median_cv          0.245
#>  sd_cv              0.021
#>  mean_ci_width      0.098
#>  median_ci_width    0.098
#>  sd_ci_width        0.009
#>  pct_ci_less_width  0.590
#>  pct_ci_miss_low    0.010
#>  pct_ci_miss_high   0.010
#>  total_type_I_error 0.020
#>  total_N            59   
#>  true_cv            0.250
#>  supposed_cv        0.250
#>  width              0.100
#>  conf_level         0.950
#> 
#> Replications: 100
#> 
#> Confidence level: 95%

# Question 2: what happens if the planning value is materially smaller than
# the true coefficient of variation (a common direction of misspecification,
# since planning values are often optimistic)? The intervals will be wider on
# average than the target and the pct_ci_less_width will fall.
set.seed(113)
ss_aipe_cv_sensitivity(
  true_cv          = .35,
  supposed_cv      = .25,
  width            = .10,
  assurance        = NULL,
  conf_level       = .95,
  G                = 100,
  print_iter       = FALSE
)
#>  term               value
#>  mean_cv            0.343
#>  median_cv          0.343
#>  sd_cv              0.030
#>  mean_ci_width      0.145
#>  median_ci_width    0.144
#>  sd_ci_width        0.015
#>  pct_ci_less_width  0.000
#>  pct_ci_miss_low    0.010
#>  pct_ci_miss_high   0.010
#>  total_type_I_error 0.020
#>  total_N            59   
#>  true_cv            0.350
#>  supposed_cv        0.250
#>  width              0.100
#>  conf_level         0.950
#> 
#> Replications: 100
#> 
#> Confidence level: 95%
```
