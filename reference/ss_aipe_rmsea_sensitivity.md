# A Priori Monte Carlo Simulation for Sample Size Planning for RMSEA in SEM

Conduct a priori Monte Carlo simulation to empirically study the effects
of (mis)specifications of input information on the calculated sample
size. The sample size is planned so that the expected width of a
confidence interval for the population RMSEA is no larger than desired.
Random data are generated from the true covariance matrix but fit to the
proposed model, whereas the sample size is calculated based on the input
covariance matrix and proposed model.

## Usage

``` r
ss_aipe_rmsea_sensitivity(
  width,
  model,
  true_cov_matrix,
  specified_N = NULL,
  conf_level = 0.9,
  G = 200,
  filename = NULL,
  ...,
  seed = NULL
)
```

## Arguments

- width:

  desired confidence interval width for the population RMSEA.

- model:

  the model the researcher proposes, which may or may not be the true
  model, written in lavaan model syntax (see
  [`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)).
  The observed variable names in the model must match the row and column
  names of `true_cov_matrix`.

- true_cov_matrix:

  the true population covariance matrix, which is used to generate
  random data for the simulation study. The row and column names of
  `true_cov_matrix` must match the observed variables in `model`.

- specified_N:

  the total sample size *N*. If `specified_N` is given, random samples
  of that size are generated. Otherwise the sample size is calculated
  with the sample size planning method so that the expected width of a
  confidence interval for the population RMSEA is no larger than
  `width`, for the interval formed as lavaan forms it (see Details).

- conf_level:

  confidence level of the planned and the simulated intervals (i.e., 1 -
  the Type I error rate). The default, 0.90, is the conventional level
  of the RMSEA interval and the default of
  [`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md)
  and
  [`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md).

- G:

  number of replications in the Monte Carlo simulation.

- filename:

  an optional path for a comma separated file recording every converged
  replication (its index, the RMSEA estimate, the two confidence limits,
  and the interval width): nothing is written when `filename` is `NULL`
  (the default), a new file with a header row is created otherwise, an
  existing file at that path is appended to, and a throwaway run should
  point it at `tempfile(fileext = ".csv")`.

- ...:

  additional arguments passed to
  [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) when fitting the
  model (for example `estimator` or `missing`). Every name given here
  must be an argument or an option of lavaan (see
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
Monte Carlo study. The `term` entries are: `"mean_rmsea"`,
`"median_rmsea"`, `"sd_rmsea"` (summaries of the realized RMSEA
estimates across the converged replications); `"mean_ci_width"`,
`"median_ci_width"`, `"sd_ci_width"` (summaries of the realized interval
widths); `"pct_ci_less_width"` (proportion of intervals narrower than
the target width); `"pct_ci_miss_low"` and `"pct_ci_miss_high"`
(tail-specific empirical non-coverage of the population RMSEA);
`"total_type_I_error"` (overall empirical non-coverage, the sum of the
two tails); and the echoes `"total_N"` (the *N* evaluated), `"df"`
(model degrees of freedom), `"true_rmsea"` (the population RMSEA
recovered from fitting `model` to `true_cov_matrix`), `"width"`, and
`"conf_level"`. The proportion rows are on the 0 to 1 scale, not
percentages. The number of converged replications is
`attr(x, "simulation")$returned` (see the section on the record of the
simulation).

## Details

This function implements the sample size planning method proposed in
Kelley and Lai (2011). It uses
[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) to fit the proposed
model to the population covariance matrix, which recovers the population
RMSEA (the model misspecification) and the model degrees of freedom, and
to fit the model to each simulated sample, and it forms the confidence
interval for the population RMSEA in each replication from that sample's
chi square statistic, inverted with
[`ci_nc_chisq`](https://yelleknek.github.io/DMAR/reference/ci_nc_chisq.md)
as
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
inverts it. The model is specified in lavaan syntax, so lavaan must be
installed.

The population RMSEA is read from
[`lavaan::fitMeasures()`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html),
which is computed reliably for the large-sample population fit.

Each replication's interval is formed from the statistic \\T\\ of that
replication's fit and the sample size multiplier \\a\\ lavaan formed
\\T\\ with (`sample_size_convention` in
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)):
\\a = N\\ for maximum likelihood under lavaan's default normal
likelihood, and \\a = N - 1\\ under `likelihood = "wishart"` and for the
least squares estimators, either of which `...` can request. That is the
interval lavaan itself reports for the sample, whatever the size of
\\T\\. When \\T\\ is at or below its degrees of freedom, lavaan reports
an RMSEA of 0, from which
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)
could rebuild only \\T = df\\; the interval formed from \\T\\ itself is
narrower the smaller \\T\\ is, and its upper limit is 0 when \\T\\ lies
at or below the `(1 - conf_level) / 2` quantile of the central chi
square, the convention lavaan and
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
share. The replications of a model that fits well often have such
statistics, so the distinction matters most when the proposed model is
close to the true one. The RMSEA summaries are those of the RMSEA lavaan
reports.

The width of the interval depends on the sample size only through \\a\\,
and
[`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md)
plans with \\a = N - 1\\, the convention of Kelley and Lai (2011). The
simulated samples have the planner's `necessary_N` cases, the size a
study would collect, so under lavaan's default \\a = N\\ each realized
interval is the one lavaan reports for a sample of that size, whose
multiplier is one larger than the plan's, and the realized widths sit
below the planned width by an amount of order \\1/N\\.

The RMSEA is conventionally reported with a 90 percent interval, because
a \\100(1 - 2\alpha)\\ percent interval decides the one-sided tests of
close and not close fit at level \\\alpha\\ (MacCallum, Browne, &
Sugawara, 1996, p. 136), and that is the default here, as in
[`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md).

## Note

Replications in which lavaan fails to converge, or for which the RMSEA
is undefined, are skipped; the number of converged replications is
reported as `attr(x, "simulation")$returned`. A fit whose likelihood has
no maximum, an error variance running off below zero while the
likelihood is flat or still rising (see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)),
is used as it is, at `true_cov_matrix` and in a replication: its
estimates are wherever lavaan's optimizer stopped on that ridge, but the
RMSEA is formed from the chi square alone, which is not, since the
likelihood levels off along the ridge and the optimizer stops only where
its steps no longer gain. Increase `G` if many replications fail to
converge.

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

Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix
that yields a specified minimizer and a specified minimum discrepancy
function value. *Psychometrika, 57*, 357–369.
[doi:10.1007/BF02295424](https://doi.org/10.1007/BF02295424)

Kelley, K., & Lai, K. (2011). Accuracy in parameter estimation for the
root mean square error of approximation: Sample size planning for narrow
confidence intervals. *Multivariate Behavioral Research, 46*, 1–32.
[doi:10.1080/00273171.2011.543027](https://doi.org/10.1080/00273171.2011.543027)

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power
analysis and determination of sample size for covariance structure
modeling. *Psychological Methods, 1*(2), 130–149.
[doi:10.1037/1082-989X.1.2.130](https://doi.org/10.1037/1082-989X.1.2.130)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

## See also

[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`ss_aipe_rmsea`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea.md),
[`ci_rmsea`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# True data generating model: two correlated factors, each measured by
# three standardized indicators. The factor correlation is 0.5 and every
# loading is 0.7. The implied population covariance matrix is assembled
# from the loading matrix, the factor correlation matrix, and the
# residual variances.
Lambda <- matrix(0, 6, 2)
Lambda[1:3, 1] <- 0.7
Lambda[4:6, 2] <- 0.7
Phi   <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
true_cov_matrix <- Lambda %*% Phi %*% t(Lambda) + diag(1 - 0.7^2, 6)
dimnames(true_cov_matrix) <- list(paste0("x", 1:6), paste0("x", 1:6))

# Proposed (misspecified) model: a single common factor.
proposed <- "g =~ x1 + x2 + x3 + x4 + x5 + x6"

# The proposed model is fit once at a very large N to recover the
# population RMSEA, the sample size is planned so that the expected width
# of the 90 percent interval is 0.05, and a fresh sample of that size is
# drawn and fit on every replication. Notice that true_rmsea is about
# 0.20, since a single factor is a poor description of two-factor data,
# that the realized widths sit close to the target, and that
# pct_ci_less_width is well away from 0 and 1: planning for the expected
# width aims the average width at the target, not the share of intervals
# narrower than it, which is what assurance controls. G = 10 keeps
# the example quick and estimates that share only roughly; a reported
# sensitivity study deserves the default G = 200 or more.
set.seed(113)
ss_aipe_rmsea_sensitivity(width = 0.05, model = proposed,
                          true_cov_matrix = true_cov_matrix, G = 10)
#>  term               value
#>  mean_rmsea         0.205
#>  median_rmsea       0.208
#>  sd_rmsea           0.017
#>  mean_ci_width      0.050
#>  median_ci_width    0.050
#>  sd_ci_width        0.000
#>  pct_ci_less_width  0.700
#>  pct_ci_miss_low    0.000
#>  pct_ci_miss_high   0.000
#>  total_type_I_error 0.000
#>  total_N            493  
#>  df                 9    
#>  true_rmsea         0.204
#>  width              0.050
#>  conf_level         0.900
#> 
#> Replications: 10
#> 
#> Confidence level: 90%
```
