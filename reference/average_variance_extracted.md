# Average Variance Extracted (AVE)

The average variance extracted (AVE) is the mean proportion of indicator
variance a factor accounts for, a standard convergent validity summary
for a reflective measurement block in confirmatory factor analysis and
structural equation modeling. With standardized loadings \\\ell_j\\,
\$\$\mathrm{AVE} = \frac{1}{J} \sum_j \ell_j^2.\$\$ The quantity itself
is elementary. In a model with cross-loadings, an item contributes its
loading to the AVE of every factor it loads on; the same shared variance
then counts toward each factor's summary, so compare AVE values across
factors of such a model with that overlap in mind. A loading the model
holds at zero, fixed (`0*item`) or constrained (`a*item` with `a == 0`),
is not a loading: the model written with it is the model written without
it, and the item is not one of the factor's indicators. Fornell and
Larcker (1981) are credited for establishing it as a validity criterion:
a construct shows convergent validity when its AVE reaches the
conventional 0.50 (the construct explains at least half its indicators'
variance), and the Fornell-Larcker discriminant criterion compares each
construct's AVE with its squared correlations with the other constructs.
AVE is closely related to composite reliability (omega); for
discriminant validity,
[`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) gives the
heterotrait-monotrait ratio of correlations, a criterion that does not
use the AVE.

## Usage

``` r
average_variance_extracted(
  fit = NULL,
  loadings = NULL,
  conf_level = 0.95,
  ci_method = c("none", "percentile"),
  B = 1000L,
  seed = NULL
)
```

## Arguments

- fit:

  Optional lavaan fit of one group and one level (for example from
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) or
  [`lavaan::cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html)); the
  standardized loadings of every factor are extracted and one AVE is
  reported per factor. A fit whose optimizer did not converge, and a
  multiple-group or multilevel fit, are refused (see *The fit* under
  Details).

- loadings:

  Optional numeric vector of standardized loadings for a single block,
  as an alternative to `fit`. Supply exactly one of `fit` and
  `loadings`. No interval can be constructed from loadings alone (their
  sampling variability is not carried by the numbers), so
  `ci_method = "percentile"` requires `fit`.

- conf_level:

  Confidence level for the bootstrap interval (default `0.95`); used
  when `ci_method = "percentile"`.

- ci_method:

  Interval method: `"none"` (the default) or `"percentile"`. No
  closed-form interval for the AVE is in common use; the percentile
  bootstrap is the standard route in the validity literature. The cases
  in the fitted data are resampled with replacement `B` times, the model
  is refit to each resample, and each factor's interval is the pair of
  empirical quantiles of its `B` AVE values (Efron & Tibshirani, 1993).
  Replications whose refit fails or does not converge are set aside, and
  the interval is computed from those that return a value; a single
  warning reports how many were set aside, and gives the caution
  described under *Details* when it holds. A replication whose solution
  is improper for a factor is kept and counted (see *Details*). An
  interval requires that no more than half of the attempted replications
  were set aside and that at least the number the level requires
  returned an estimate (400 at the 95 percent level; see `B`). Otherwise
  the point estimates and counts are retained, with `NA` limits and a
  warning of class `dmar_bootstrap_insufficient`. These floors limit the
  error of drawing a finite number of resamples; they are not coverage
  guarantees.

- B:

  Number of bootstrap replications when `ci_method = "percentile"`, a
  single whole number. Default `1000`, smaller than the package's usual
  `10000` because every replication refits the model; raise it for a
  reported analysis when time allows. The fewest accepted is 400 at the
  default `conf_level` (200 at 0.90, and 2000 at 0.99, which the default
  does not meet); a smaller or fractional `B` is refused before anything
  is computed (class `dmar_count_invalid`). See *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- seed:

  Optional integer seed for the bootstrap. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  internally and the prior state restored on exit.

## Value

A `data.frame` (classes `dmar_ave` and `dmar_tbl`) with one row per
factor: `factor` (label), `ave`, and `ci_lower` / `ci_upper` (the
percentile bootstrap limits; `NA` when `ci_method = "none"`). All three
are `NA` for a factor whose solution is improper (see Details). With
`ci_method = "percentile"` and at least one factor whose solution is
proper, so that an interval is formed, the table carries the attributes
`conf_level`, `B` (requested), `B_used` (the replications the intervals
rest on), `n_dropped`, `B_required`, and `bootstrap_status` (see
*Bootstrap replications* under Details), and `bootstrap`, the counts in
the form the reliability family records them (see
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)):
`attempted` (`B`), `returned`, `failed` (the replications whose refit
failed or did not converge, which lavaan's bootstrap does not tell
apart, counted in `reasons` as non-finite), `reasons`, `error_messages`,
`undetermined`, `improper`, and `no_solution` (all zero here),
`improper_kept` (the replications used whose solution was improper for
at least one factor, including when its interval is withheld),
`not_converged` (zero), `refit_failed` (the same count as `failed`,
under the reason that the line under the table and the warning give for
both causes, that the refit failed or did not converge), `improper_rule`
(`"keep"`), `caution`, and `improper_kept_by_factor`, the improper
replications kept for each factor, with `factors_with_intervals`
recording which intervals were eligible.

## Details

**The fit.** The AVE is read from a solution, so a fit whose optimizer
did not converge (`lavaan::lavInspect(fit, "converged")` is `FALSE`, as
when it reached its iteration limit) is refused with an error of class
`dmar_not_converged`: its estimates mark where the optimizer stopped. A
fit that converged to an improper solution is not refused; it is
described under *Improper solutions* below.

A factor's standardized loadings are those of one group at one level. A
factor of a multiple-group fit therefore has an AVE in each group, and a
factor of a two-level fit has one at each level, even when it has the
same name at both. Neither the unweighted mean of those values nor a
mean weighted by group size is the AVE of any group, so the function
describes a fit of one group and one level and refuses a multiple-group
or multilevel fit with an error of class `dmar_argument_invalid` that
names the groups (as `groups`) or the levels (as `levels`). The AVE of
one group or level of such a model is
`average_variance_extracted(loadings = )` given the factor's
standardized loadings in that group or at that level: the factor's
column of that group's or level's `lambda` matrix in
`lavaan::lavInspect(fit, "std")`, without the zeros of the items that do
not load on it. For a group, a fit of the model to that group's data
alone also gives it, without any equality constraints across groups.

**Improper solutions.** When the error variance of an indicator is
estimated below zero (a Heywood case; for an ordered item, the residual
variance of its latent response), its communality, one minus its
standardized error variance, exceeds one, and so does its squared
standardized loading in a model without cross-loadings. The AVE averages
those communalities, so for a factor with such an indicator it is not a
proportion of variance: its `ave` (and its interval) is `NA`, with a
warning naming the indicators and their communalities, as in the
measurement table of
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).
On `holzinger_swineford`, a factor measured by visual perception and
lozenges alone, beside one measured by general information, paragraph
comprehension, and sentence completion, gives visual perception an error
variance below zero and a communality of 1.30, so the first factor's AVE
is `NA` (the mean of its squared standardized loadings would be .726),
while the verbal factor's, .705, is reported. A larger sample, or
reconsidering the model for those indicators, are the usual remedies.

No interval is formed for such a factor, so when every factor of the fit
is improper the bootstrap is not run, and the table records no
confidence level.

The bootstrap meets the same problem one resample at a time. A
replication whose solution is improper for a factor that is proper in
the observed fit is kept, as the reliability family keeps one under its
default `improper = "keep"` (see
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)),
and counted. Its AVE for that factor averages a communality above one
and can itself exceed one, so such replications gather in the upper tail
of the factor's bootstrap distribution. The line printed under the table
says how many replications were attempted, how many returned a value,
and how many of those came from an improper solution, with, for two or
more factors, the count for each; when the replications set aside and
the improper ones kept for a factor reach the tail area beyond a limit,
\\(1 - C)/2\\ of `B` for confidence level \\C\\, it adds a caution,
since the limits then depend on how those replications were treated;
with two or more factors the caution gives the share of each factor it
holds for, counted for that factor alone. The warning the call gives
says the same, as the reliability family's does, naming each factor the
caution holds for with its count. On 80 pupils drawn from
`holzinger_swineford` (`set.seed(113)`), a factor of visual perception,
cubes, and lozenges beside a verbal factor is proper in the sample, with
an AVE of .449, but 101 of the 198 replications that returned a value
were improper for it (`B = 200`, `seed = 113`), and its upper limit,
1.19, is not a proportion of variance; the verbal factor had none. With
the 2 replications whose refit failed or did not converge, 103 of the
200 are in question for the spatial factor, against the 5 that lie
beyond each limit of a .95 interval, and the warning names it. The
counts are those the console and a knitted report print, and they are
kept in the `bootstrap` attribute.

**Bootstrap replications.** The percentile interval reads its two limits
from the quantiles of the `B` replicated AVE values of each factor, so
it is run only when `B` puts at least 10 replications beyond each limit
(400 at the 95 percent level). A replication in which an AVE is
undefined (the refit failed or did not converge) is dropped and counted,
never replaced, and the interval is computed from the rest. The result
records `B` (requested), `B_used` (the number the interval rests on),
`n_dropped`, and `B_required` (the fewest replications a request may ask
for; that count is held to `B`, not to `B_used`) as attributes, beside
the fuller `bootstrap` list, with `bootstrap_status` (`"complete"`,
`"dropped"`, `"exploratory"`, `"insufficient_count"`, or
`"insufficient_share"`). One warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the 2.5 percent of `B` that lies beyond one limit at the 95
percent level, because they could then decide by themselves where that
limit falls. When more than half were dropped, the limits of every
factor are `NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`.
`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration; every such call gives a
warning of class `dmar_bootstrap_exploratory` and the result is not one
to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

The percentile bootstrap interval resamples the cases behind `fit` and
refits the model once per replication, so its cost is `B` model fits.
That refitting is why the examples below stop at the point estimates:
the fewest replications a 95 percent interval is reported from,
`B = 400`, run well beyond the time an example is given on the
two-factor model there. To obtain the interval, pass
`ci_method = "percentile"` together with a `seed`, as in
`average_variance_extracted(fit, ci_method = "percentile", seed = 113)`,
and keep the default `B = 1000` or more for a reported analysis.

A supplied two-stage fit requires eligible unrestricted first-stage
moments. Bootstrap intervals currently refuse separately stored
two-stage auxiliary variables because the refits do not carry their
resampled values. Resample and refit items and auxiliaries together
explicitly when that analysis is needed.

## References

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation
models with unobservable variables and measurement error. *Journal of
Marketing Research, 18*(1), 39–50.

## See also

[`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) for
discriminant validity;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for the composite reliability of the same block (coefficient omega is
what the composite-reliability literature computes);
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) to obtain
the fit.

Other multivariate and latent variable methods:
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
[`sem_exact_fit_bootstrap()`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
[`sem_fit_measures()`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
[`sem_information_criteria()`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Directly from the standardized loadings a paper reports:
average_variance_extracted(loadings = c(.8, .7, .6))
#>  factor ave   ci_lower ci_upper
#>  f      0.497 <NA>     <NA>    

# From a fitted model, one AVE per factor (requires lavaan).
data(holzinger_swineford)
fit <- lavaan::cfa(
  "verbal    =~ t6_paragraph_comprehension + t7_sentence +
                t9_word_meaning
   deduction =~ t20_deduction + t22_problem_reasoning +
                t23_series_completion",
  data = holzinger_swineford)
ave_tbl <- average_variance_extracted(fit)
ave_tbl
#>  factor    ave   ci_lower ci_upper
#>  deduction 0.469 <NA>     <NA>    
#>  verbal    0.719 <NA>     <NA>    

# The Fornell and Larcker (1981) discriminant criterion compares each
# factor's AVE with the squared correlation between the factors: a
# factor should account for more of its own indicators' variance than
# it shares with the other factor. Here the comparison favors verbal
# and goes against deduction, whose AVE falls below the shared
# variance. Fitting with fa_confirmatory(..., output = "measurement") puts the
# AVE values and the latent correlations in one table.
lavaan::lavInspect(fit, "cor.lv")["verbal", "deduction"]^2
#> [1] 0.5327517

# The broom verbs: one row per factor.
generics::tidy(ave_tbl)
#>        term  estimate estimand ci_lower ci_upper conf_level
#> 1 deduction 0.4694312      ave       NA       NA         NA
#> 2    verbal 0.7193173      ave       NA       NA         NA
generics::glance(ave_tbl)
#>   n_terms conf_level B_used
#> 1       2         NA     NA
```
