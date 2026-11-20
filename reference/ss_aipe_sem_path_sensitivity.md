# A Priori Monte Carlo Simulation for Sample Size Planning for SEM Targeted Effects

Conduct a priori Monte Carlo simulation to empirically study the effects
of (mis)specifications of input information on the calculated sample
size. Random data are generated from the true covariance matrix but fit
to the proposed model, whereas sample size is calculated based on the
input covariance matrix and proposed model.

## Usage

``` r
ss_aipe_sem_path_sensitivity(
  model,
  supposed_cov_matrix,
  true_cov_matrix = supposed_cov_matrix,
  which_path,
  width,
  specified_N = NULL,
  conf_level = 0.95,
  assurance = NULL,
  G = 100,
  filename = NULL,
  ...,
  seed = NULL
)
```

## Arguments

- model:

  A single character string giving the free analysis model in lavaan
  model syntax (see
  [`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)),
  the model that would be fit to the data. The target path must carry a
  parameter label so it can be referred to by name, for example
  `"f2 ~ b*f1"` labels the structural path `b`. The model may or may not
  be the true data generating model

- supposed_cov_matrix:

  the population covariance matrix supposed when the sample size is
  planned (passed to
  [`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)
  as its `true_cov_matrix`). It is a planning assumption, not a matrix
  estimated from data, and it may or may not be the true covariance
  matrix. The row names and column names of `supposed_cov_matrix` should
  be the same as the observed variables in `model`

- true_cov_matrix:

  the true population covariance matrix, which will be used to generate
  random data for the simulation study. The row names and column names
  of `true_cov_matrix` should be the same as the observed variables in
  `model`. Defaults to `supposed_cov_matrix`

- which_path:

  the parameter label of the targeted path, given as a character string,
  for example `"b"` for the path labeled `f2 ~ b*f1` in `model`

- width:

  desired confidence interval width for the model parameter of interest

- specified_N:

  the total sample size *N* of each set of random data. If it is `NULL`,
  it is the sample size
  [`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)
  plans from `supposed_cov_matrix`

- conf_level:

  confidence level (i.e., 1- Type I error rate)

- assurance:

  the assurance that the confidence interval obtained in a particular
  study will be no wider than desired (must be `NULL` or a value between
  0.50 and 1)

- G:

  number of replications in the Monte Carlo simulation

- filename:

  an optional path for a comma separated file recording every
  replication (the estimate of the targeted path, its standard error,
  the two confidence limits, and the interval width): nothing is written
  when `filename` is `NULL` (the default), a new file with a header row
  is created otherwise, an existing file at that path is appended to,
  and a throwaway run should point it at `tempfile(fileext = ".csv")`

- ...:

  further arguments passed to
  [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) for every fit: the
  planning fit to `supposed_cov_matrix` (through
  [`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)),
  the fit to `true_cov_matrix`, and the fit to each simulated sample
  (for example `std.lv = TRUE`). Every name given here must be an
  argument or an option of lavaan (see
  [`lavOptions`](https://rdrr.io/pkg/lavaan/man/lavOptions.html)); any
  other name, and any argument without a name, is refused before a model
  is fitted, with a condition of class `dmar_argument_invalid` that
  gives the name as it was typed

- seed:

  `NULL` (the default) to draw from the current random number stream, or
  a single number that makes the call reproducible; a supplied seed is
  set for this call only, and the caller's random number state is
  restored when the function returns.

## Value

A `data.frame` with columns `term` and `value` summarizing the a priori
Monte Carlo study. The `term` entries are: `"mean_path"`,
`"median_path"`, `"sd_path"` (summaries of the realized estimates of the
targeted path across the converged replications); `"mean_ci_width"`,
`"median_ci_width"`, `"sd_ci_width"` (summaries of the realized interval
widths); `"pct_ci_less_width"` (proportion of realized widths at or
below `width`); `"pct_ci_miss_low"` and `"pct_ci_miss_high"`
(tail-specific empirical non-coverage of the population path);
`"total_type_I_error"` (overall empirical non-coverage, the sum of the
two tails); and the echoes `"total_N"` (the *N* evaluated),
`"true_path"` (the population value of the targeted path under
`true_cov_matrix`), `"width"`, `"conf_level"`, and `"assurance"`
(present only when an assurance was supplied). The proportion rows are
on the 0 to 1 scale, not percentages. The number of replications that
gave an interval is `attr(x, "simulation")$returned`, and the number of
data sets drawn again is `attr(x, "simulation")$n_redrawn` (see the
sections below).

## Details

This function implements the sample size planning methods proposed in
Lai and Kelley (2011). It calls
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)
to plan the sample size and to identify the targeted path, then fits the
analysis model to data simulated from `true_cov_matrix` with
[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html). The analysis model is
written in lavaan model syntax with the targeted path given a parameter
label; see
[`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html) for
the syntax and [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) for the
fitting machinery. The population covariance matrices are most naturally
produced by
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md) from
a fully fixed population model. This function requires lavaan and MASS
to be installed.

## Note

Occasionally a replication fails to converge when the analysis model is
fit to a simulated data set. Such replications are not counted toward
the `G` converged replications; the simulation draws fresh data and
continues. A replication whose likelihood has no maximum, an error
variance running off below zero while the likelihood is flat or still
rising (see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)),
is treated the same way, and an analysis model with no maximum at
`supposed_cov_matrix` or `true_cov_matrix` stops the function, before
any replication is drawn, with an error of class `dmar_no_maximum` that
prints the lavaan code of that fit. A safety cap stops the loop after
`20 * G` attempts, and a single warning is issued if fewer than `G`
replications converged.

## Replications drawn again

A simulated data set whose model fit does not converge, or whose
estimate is not a maximum of the likelihood, is replaced by a new one,
so such a data set is not among the replications without an interval and
is in no summary row. At most `20 * G` data sets are drawn in all; if
that many do not give `G` converged fits, the function warns and
summarizes the fits it has. The number of data sets drawn and set aside
is the element `n_redrawn` of the `simulation` attribute, an integer
read with `attr(x, "simulation")$n_redrawn`: 0 when every fit converged,
and the total number of data sets drawn is `n_redrawn` plus the number
of converged fits. The other elements of the record keep their meaning:
`attempted` is `G`, `returned` the replications with an interval, and
`failed` those without one (a converged fit whose variance for the path
was negative or missing, or a replication never obtained). A large
`n_redrawn` beside `G` says the model is hard to fit at this sample
size, and the summaries then describe only the samples in which it could
be fitted.

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

Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
targeted effects in structural equation modeling: Sample size planning
for narrow confidence intervals. *Psychological Methods, 16*(2),
127–148. [doi:10.1037/a0021764](https://doi.org/10.1037/a0021764)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

## See also

[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md),
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The planning values a researcher would bring to ss_aipe_sem_path(): two
# factors, each measured by three indicators, joined by a structural path
# of 0.5, with every residual variance at 0.5.
planning_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f2 ~ 0.5*f1
  f1 ~~ 1*f1
  f2 ~~ 0.75*f2
  y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
  y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
"

# The population the study will actually sample from: the same structural
# path, but noisier indicators than the planning values assumed, with
# every residual variance at 0.8.
true_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f2 ~ 0.5*f1
  f1 ~~ 1*f1
  f2 ~~ 0.75*f2
  y1 ~~ 0.8*y1; y2 ~~ 0.8*y2; y3 ~~ 0.8*y3
  y4 ~~ 0.8*y4; y5 ~~ 0.8*y5; y6 ~~ 0.8*y6
"

# The free analysis model, with the targeted path labeled "b".
analysis_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f2 ~ b*f1
"

supposed_cov_matrix <- cov_sem(planning_model)$sigma_theta
true_cov_matrix <- cov_sem(true_model)$sigma_theta

# The sample size planned from the optimistic measurement quality is
# evaluated against the population that actually holds. Notice that the
# mean realized width (mean_ci_width) exceeds the desired 0.30, and that
# pct_ci_less_width, the proportion of intervals no wider than desired, falls
# well below the one half that the planned sample size would deliver if
# the planning values were right. G = 10 keeps the example quick; a
# reported sensitivity study deserves the default G = 100 or more.
set.seed(113)
ss_aipe_sem_path_sensitivity(model = analysis_model,
                             supposed_cov_matrix = supposed_cov_matrix,
                             true_cov_matrix = true_cov_matrix,
                             which_path = "b", width = 0.30, G = 10)
#>  term               value
#>  mean_path          0.497
#>  median_path        0.502
#>  sd_path            0.053
#>  mean_ci_width      0.334
#>  median_ci_width    0.338
#>  sd_ci_width        0.032
#>  pct_ci_less_width  0.200
#>  pct_ci_miss_low    0.000
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.000
#>  total_N            264  
#>  true_path          0.500
#>  width              0.300
#>  conf_level         0.950
#> 
#> Replications: 10
#> 
#> Confidence level: 95%
```
