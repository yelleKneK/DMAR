# Bootstrap Test of Exact Fit for a Structural Equation Model With Complete or Incomplete Data

Gives the chi square test of exact fit of a structural equation model a
reference distribution that does not assume normal data: the model based
bootstrap of Bollen and Stine (1992), in which the data are first
transformed so that the model holds in them exactly, as estimated, and
the statistic is then recomputed on `B` samples drawn with replacement
from the transformed data. The share of bootstrap statistics at least as
large as the observed one is the bootstrap *p*-value. For incomplete
data fitted by full information maximum likelihood (FIML) the
transformation is extended so that the FIML estimates of the saturated
model on the transformed data are the model implied moments (Enders,
2002; Savalei & Yuan, 2009), and each bootstrap sample keeps the missing
values of the rows it draws, so that the model is fitted to it by FIML
as it was to the data. Requires lavaan (in Suggests).

## Usage

``` r
sem_exact_fit_bootstrap(
  fit,
  B = 1000,
  transformation = c("common", "pattern"),
  improper = c("keep", "drop"),
  conf_level = 0.95,
  seed = NULL,
  ncpus = 1L
)
```

## Arguments

- fit:

  A fitted lavaan model: a maximum likelihood fit (estimator ML, MLM,
  MLMV, or MLR) of raw continuous data, complete or incomplete under
  `missing = "ml"` (with auxiliary variables through lavaan's `aux` or
  without), in one or several groups, that converged to an admissible
  solution. See Details for what is refused.

- B:

  Number of bootstrap samples, a single whole number of at least 100.
  Defaults to `1000`: every sample refits the model. A Monte Carlo
  *p*-value is run only from a request of at least 100 samples, at which
  its standard error is at most 0.05; a smaller `B` is refused before
  anything is computed (class `dmar_count_invalid`). See
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- transformation:

  For incomplete data, how the transformation makes the saturated FIML
  estimates equal the model implied moments: `"common"` (the default)
  finds one matrix that transforms every row on its observed variables,
  symmetric on the scale of the variables divided by their saturated
  standard deviations, so that the FIML estimating equations hold summed
  over the missing data patterns; `"pattern"` transforms the rows of
  each missing data pattern on their own, with the pattern's own
  covariance matrix and means, so that each pattern's covariance matrix
  and means equal the model implied ones exactly, which needs every
  pattern to have more rows than observed variables. Both are solved on
  that standardized scale, so the result does not depend on the units of
  the variables. See Details. Ignored with complete data.

- improper:

  Whether a bootstrap sample whose fit converged to an inadmissible
  solution counts in the *p*-value: `"keep"` (the default) or `"drop"`.
  A solution is inadmissible when lavaan's post-estimation check finds a
  negative variance estimate or a covariance matrix that is not positive
  definite, of the latent variables, of the residuals of the observed
  variables, or of the residuals of the structural part of the model;
  with auxiliary variables the check is made on the model's own
  parameters (see What Is Counted). Either way the row `B_improper`
  reports how many there were.

- conf_level:

  Confidence level of the interval for the Monte Carlo error of the
  bootstrap *p*-value. Defaults to 0.95.

- seed:

  Random number seed for the bootstrap samples. Defaults to `NULL`,
  under which the samples are drawn from, and advance, the current
  random number stream; supply an integer for a reproducible result,
  which is set for the call only, the caller's random number state being
  restored when it returns.

- ncpus:

  The number of processes over which the bootstrap refits are spread.
  The default `1` runs them in this R session; a larger number uses
  forked workers or, on Windows, a socket cluster. Every bootstrap
  sample is drawn before any refit runs, so a given `seed` gives the
  same result whatever `ncpus` is. A session that forks workers and also
  starts background processes through callr or processx (a Shiny app, a
  job runner) needs the environment variable
  `PROCESSX_NOTIFY_OLD_SIGCHLD` set to a nonempty value before processx
  is loaded: processx otherwise replaces the signal handler through
  which parallel learns that a forked worker has finished, and parallel
  then reports at exit that it could not terminate its workers (processx
  documents the interaction and the switch). The results are not
  affected.

## Value

A `dmar_tbl` (class `dmar_sem_fit`) with columns `term` and `value` and,
always in this order, the rows `chi_square`, `df`, `p_chi_square` (the
asymptotic *p*-value of the standard statistic), `p_bootstrap`,
`p_bootstrap_ci_lower`, `p_bootstrap_ci_upper`, `bootstrap_mean`, `B`,
`B_used`, `B_improper`, `B_failed`, `N`, and `n_groups`. The attributes
`B`, `B_used`, `n_dropped`, `B_required` (100), and `bootstrap_status`
record the samples in the names every bootstrap result of the package
carries. The attribute `bootstrap` is a data frame with one row per
bootstrap sample, its `chi_square` (`NA` when it was set aside) and its
`status` (`proper`, `improper`, `not_converged`, or `error`); the
attribute `transformation` records, per group, the transformation
applied, the number of missing data patterns, the Newton steps and the
largest estimating equation at the solution (on the standardized scale),
and the largest difference between the transformed means and the model
implied ones; and the `sem_fit` attribute records the statistic and its
degrees of freedom, the transformation requested and the ones applied,
the counts, and the lines of the printed footer.

## Details

### The Model Based Bootstrap

The chi square statistic \\T\\ of a maximum likelihood fit follows its
chi square distribution under the model only for normal data in large
samples. Bollen and Stine (1992) proposed to find its distribution under
the model by bootstrap instead: resampling the data as they are would
reproduce the misfit of the model as well as the sampling variation, so
the data are first transformed to \$\$Z = (Y - \bar{Y}) S^{-1/2}
\hat{\Sigma}^{1/2}\$\$ (Enders, 2002, Eq. 1), where \\S\\ is the sample
covariance matrix, \\\hat{\Sigma}\\ the covariance matrix the fitted
model implies, and the square roots are symmetric. The transformed data
have covariance matrix exactly \\\hat{\Sigma}\\, so the model fits them
exactly, and their distribution otherwise keeps the shape of the data.
Each bootstrap sample draws rows of \\Z\\ with replacement, within each
group as many as the group has (each group is transformed with its own
moments), the model is fitted to it, and its statistic \\T^\*\\ is
recorded. The bootstrap *p*-value is the share of the \\T^\*\\ at least
as large as the observed \\T\\ (row `p_bootstrap`); a small value says
that the observed statistic is unusual for data in which the model holds
and have this distribution. Under a mean structure the model implied
means are added back, so the means hold too. The statistic is the
standard (uncorrected) likelihood ratio statistic in every case: the
bootstrap replaces the correction a scaled statistic applies, and a fit
by MLM or MLR contributes its estimates and its standard statistic. This
is the bootstrap `lavaan::bootstrapLavaan(type = "bollen.stine")` runs
for complete data, with the same transformation and resamples: the
`seed` given here, passed to lavaan as its `iseed` (with lavaan run
serially), gives the same \\T^\*\\. lavaan called without `iseed` draws
its own seed from the caller's random number stream, so a
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the two calls
does not make them agree.

### Incomplete Data

With missing values the transformation cannot be applied row by row,
since a row observes only some of the variables. Enders (2002, Eq. 4)
first extended it by transforming each row with the submatrices, for the
variables the row observes, of the saturated and the model implied
covariance matrices estimated by FIML from the whole sample, and adding
back the model implied means; the model then does not hold exactly in
the transformed data, since the rows of a pattern keep their own
covariance matrix, and that transformation is not offered here. Savalei
and Yuan (2009) gave two transformations under which the FIML estimates
of the saturated model on the transformed data are exactly \\\hat{\mu}\\
and \\\hat{\Sigma}\\, so that the model fits the transformed data
exactly under FIML, as the complete data transformation makes it fit
exactly under complete data. The `"pattern"` transformation (their
first) transforms the rows of each missing data pattern with the
pattern's own covariance matrix and means, so each pattern's rows take
exactly the model implied moments of the variables they observe; it
needs every pattern to have more rows than observed variables, since the
pattern's covariance matrix is inverted, and a pattern that falls short
is refused with its index and counts. The `"common"` transformation
(their second, and the default) centers each pattern at its own means,
transforms every row with one matrix, symmetric on the standardized
scale described below, chosen so that the covariance estimating
equations of the saturated FIML model, summed over the patterns, hold at
\\\hat{\Sigma}\\, and shifts each pattern to \\\hat{\mu}\\; it does not
need large patterns. DMAR finds the matrix by Newton's method on the
standardized variables, starting from the complete data transformation
applied to the EM moments, and records the number of steps and the
largest estimating equation at the solution in the `transformation`
attribute, where the pattern transformation records the same check
computed afresh. Both transformations are solved on the variables
divided by their saturated standard deviations, and the transformed
values are multiplied back: the matrix of the common transformation and
the square roots of the pattern transformation are symmetric on that
standardized scale. The bootstrap therefore does not depend on the units
of the variables (rescaling an item rescales its transformed values and
changes nothing else). Symmetry belongs to the scale on which it is
imposed, so the same transformations solved in raw units, as
`semTools::bsBootMiss()` solves the common one (its transformation 2),
give transformed data that differ element by element from these and
depend on the units; both versions make the saturated FIML estimates
exactly \\\hat{\mu}\\ and \\\hat{\Sigma}\\, so the model holds exactly
in either. Each bootstrap sample draws rows of the transformed data with
their missing values, so it has the missing data patterns of the data in
bootstrap proportions, and the model is fitted to it by FIML. In
Enders's (2002) simulation of his transformation, with 10 to 20 percent
of the values missing completely at random from nine nonnormal
indicators, the bootstrap held the Type I error rate near .05 at \\N\\ =
250 and 500, where the uncorrected FIML statistic rejected far too
often, and was conservative at \\N\\ = 100, the more so the more
nonnormal the data and the more missing values; he suggests a sample of
about 200 or more with incomplete data. The calibration of the two
transformations offered here on the current code is planned, not done.

### What Is Counted

A bootstrap sample whose fit did not converge, raised an error, or
returned no finite statistic is set aside, and one that converged to an
inadmissible solution is counted or set aside as `improper` says. That
includes a sample whose likelihood has no maximum, an error variance
running off below zero while the likelihood is flat or still rising (see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)):
its estimates are wherever lavaan's optimizer stopped on that ridge, but
its chi square is not, since the likelihood levels off along the ridge
and the optimizer stops only where its steps no longer gain, which puts
the statistic within the optimizer's tolerance of its limit. A sample's
solution is judged as the fit's is (see What Is Refused); with auxiliary
variables that is the check of the model's own parameters, since
lavaan's own check would also count as improper the samples in which
only the residual covariances of the auxiliary variables are indefinite,
which with an auxiliary variable related to a factor can be many samples
whose solution is proper. The rows `B`, `B_used`, `B_improper`, and
`B_failed` count them, and a single warning (class
`dmar_bootstrap_replications_dropped`) reports the samples set aside.
The *p*-value is formed from the `B_used` statistics; when fewer than
half of the attempted samples were used, it is `NA` with a warning of
class `dmar_bootstrap_insufficient`, since a *p*-value from the samples
that happened to converge would describe a selected distribution. The
100 required are held to the number requested, `B`, and not to `B_used`:
a request that met the count and then had samples set aside still
reports its *p*-value, the warning says how many were used, and the
Wilson interval shows the Monte Carlo error from that number. Where
refits can fail, request more than the minimum.
`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration; every such call gives a
warning of class `dmar_bootstrap_exploratory`, and the line under the
table says so (see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)).
The rows `p_bootstrap_ci_lower` and `p_bootstrap_ci_upper` are the
Wilson interval for the Monte Carlo error of the *p*-value at
`conf_level`, the uncertainty that comes from using `B_used` samples
rather than infinitely many
([`ci_proportion`](https://yelleknek.github.io/DMAR/reference/ci_proportion.md));
it is not an interval for anything about the model. The row
`bootstrap_mean` is the mean of the used statistics, the bootstrap
estimate of the expected statistic under the model, which the
Bollen-Stine adjustments of fit indices use in place of the degrees of
freedom.

### What Is Refused

A fit is refused when it did not converge (an error of class
`dmar_not_converged`) or has an inadmissible solution (a refusal of
class `dmar_inadmissible` that names the negative variance or the matrix
at fault), as for every fit measure of the SEM family (see
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
which explains the check and how it treats auxiliary variables). It is
refused with an error of class `dmar_exact_fit_bootstrap_scope`, whose
`reason` element names the case (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)),
when it was made from summary statistics rather than raw data (there are
no rows to transform); when it has ordered indicators or was fitted by
an estimator other than maximum likelihood (the transformation is
derived for the normal theory fit function); when it is a two-level
model, has a cluster variable (the resamples draw rows independently,
which the refit could not undo) or sampling weights, uses two-stage
estimation (whose statistic is not the likelihood ratio statistic the
bootstrap refits), or conditions on observed exogenous variables
(`fixed.x = TRUE` with such variables, or `conditional.x = TRUE`: the
transformation moves every variable, and the conditional likelihood
holds some fixed; refit with `fixed.x = FALSE`, and with
`conditional.x = FALSE` too, since lavaan does not take the first beside
`conditional.x = TRUE`); and when its analysis data keep missing values
without FIML (as under pairwise deletion) or it has no standard test
statistic. A missing data pattern with no more rows than observed
variables is refused for the pattern transformation with an error of
class `dmar_pattern_too_small`, which carries the pattern, its counts,
and its group, and a transformation that is undefined for the data (a
moment that is not positive definite, a Newton iteration that does not
converge) with one of class `dmar_transformation_undefined`. These three
classes, and `dmar_inadmissible`, inherit `dmar_refusal`;
`dmar_not_converged` is caught by its own name. The checks of `B`,
`conf_level`, and `ncpus` are argument checks and give plain errors. A
fit whose cases were removed by listwise deletion is bootstrapped on its
complete cases, as it was fitted, with a warning. A fit by MLM, MLMV, or
MLR is accepted: its estimates are the ML estimates, and the bootstrap
refits and resamples its standard statistic.

## References

Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit
measures in structural equation models. *Sociological Methods &
Research, 21*(2), 205–229.
[doi:10.1177/0049124192021002004](https://doi.org/10.1177/0049124192021002004)

Enders, C. K. (2002). Applying the Bollen-Stine bootstrap for
goodness-of-fit measures to structural equation models with missing
data. *Multivariate Behavioral Research, 37*(3), 359–377.
[doi:10.1207/S15327906MBR3703_3](https://doi.org/10.1207/S15327906MBR3703_3)

Savalei, V., & Yuan, K.-H. (2009). On the model-based bootstrap with
missing data: Obtaining a *p*-value for a test of exact fit.
*Multivariate Behavioral Research, 44*(6), 741–763.
[doi:10.1080/00273170903333590](https://doi.org/10.1080/00273170903333590)

## See also

[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
for the asymptotic test of exact fit beside the other fit measures, and
[`equivalence_rmsea`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md)
for the test of not close fit.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_adjusted_r()`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`fa_bifactor()`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
[`fa_second_order()`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
[`factor_congruence()`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
[`factor_scores()`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`sem_fit_measures()`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
[`sem_information_criteria()`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Two correlated factors for six of the holzinger_swineford tests, with
# the cube and sentence scores of about a third of the students missing
# at random given the general information score, which the FIML fit
# takes as an auxiliary variable. lavaan notes that it adds the
# auxiliary variable and warns that the residual covariance matrix,
# which includes it, is not positive definite; the model's own solution
# is proper (see What Is Refused). B = 10 keeps the example quick and
# is below the 100 a reported p-value needs; the option lets it run
# and the warning says so. Use the default of 1000 or more in an
# analysis.
d <- holzinger_swineford[c("t1_visual_perception", "t2_cubes", "t4_lozenges",
                           "t6_paragraph_comprehension", "t7_sentence",
                           "t9_word_meaning", "t5_general_information")]
set.seed(113)
drop <- runif(nrow(d)) < plogis((d$t5_general_information - 30) / 10) - 0.35
d[drop, c("t2_cubes", "t7_sentence")] <- NA
fit <- lavaan::cfa(
  "visual  =~ t1_visual_perception + t2_cubes + t4_lozenges
   textual =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning",
  data = d, missing = "ml", std.lv = TRUE, aux = "t5_general_information")
#> lavaan NOTE:  
#>    auxiliary (aux) variable(s) added to the model as saturated correlates 
#>    (estimated with FIML and hidden from the output): t5_general_information
#> Warning: lavaan->lav_object_post_check():  
#>    the covariance matrix of the residuals of the observed variables (theta) 
#>    is not positive definite; use lavInspect(fit, "theta") to investigate.
old <- options(dmar.bootstrap_exploratory = TRUE)
sem_exact_fit_bootstrap(fit, B = 10, seed = 113)
#> Warning: sem_exact_fit_bootstrap(): B = 10 is below the 100 replications a Monte Carlo p-value is reported from (100 for a Monte Carlo p-value). The p-value is shown because options(dmar.bootstrap_exploratory = TRUE) is set. With 10 replications it will move from seed to seed; use the result to explore or to demonstrate, and report a p-value only from B of at least 100. B, B_used, and B_required are recorded in the result.
#>  term                 value   
#>  chi_square           23.859  
#>  df                   8       
#>  p_chi_square         0.0024  
#>  p_bootstrap          < 0.0001
#>  p_bootstrap_ci_lower 0.0000  
#>  p_bootstrap_ci_upper 0.2775  
#>  bootstrap_mean       6.604   
#>  B                    10      
#>  B_used               10      
#>  B_improper           0       
#>  B_failed             0       
#>  N                    301     
#>  n_groups             1       
#> 
#> Confidence level: 95%
#> 
#> Transformation: incomplete data, one matrix for every pattern (Savalei & Yuan, 2009); the statistic is the standard likelihood ratio statistic under FIML
#> Bootstrap: B = 10 samples, 10 used; p_bootstrap is the share of used statistics at least the observed one, with the Wilson interval for its Monte Carlo error at 95 percent. Exploratory: below the 100 samples a reported p-value needs
#> Sample size: N = 301 in 1 group
options(old)
```
