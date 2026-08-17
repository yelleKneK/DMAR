# Average Variance Extracted (AVE)

The average variance extracted (AVE) is the mean proportion of indicator
variance a factor accounts for, a standard convergent validity summary
for a reflective measurement block in confirmatory factor analysis and
structural equation modeling. With standardized loadings \\\ell_j\\,
\$\$\mathrm{AVE} = \frac{1}{J} \sum_j \ell_j^2.\$\$ The quantity itself
is elementary. In a model with cross-loadings, an item contributes its
loading to the AVE of every factor it loads on; the same shared variance
then counts toward each factor's summary, so compare AVE values across
factors of such a model with that overlap in mind. Fornell and Larcker
(1981) are credited for establishing it as a validity criterion: a
construct shows convergent validity when its AVE reaches the
conventional 0.50 (the construct explains at least half its indicators'
variance), and the Fornell-Larcker discriminant criterion compares each
construct's AVE with its squared correlations with the other constructs.
AVE is closely related to composite reliability (omega); the modern
complement on the discriminant side is
[`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md).

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

  Optional lavaan fit (for example from
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) or
  [`lavaan::cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html)); the
  standardized loadings of every factor are extracted and one AVE is
  reported per factor.

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
  Replications whose refit fails or does not converge are dropped, and
  the interval is computed from those that return a value; a single
  warning reports how many were dropped.

- B:

  Number of bootstrap replications when `ci_method = "percentile"`
  (default `1000`). The default is smaller than the package's usual
  `10000` because every replication refits the model; raise it for a
  reported analysis when time allows.

- seed:

  Optional integer seed for the bootstrap. The default `NULL` uses the
  current state of the random number generator; a supplied seed is set
  internally and the prior state restored on exit.

## Value

A `data.frame` (class `dmar_tbl`) with one row per factor: `factor`
(label), `ave`, and `ci_lower` / `ci_upper` (the percentile bootstrap
limits; `NA` when `ci_method = "none"`).

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
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
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
# variance. Fitting with cfa_k(..., output = "measurement") puts the
# AVE values and the latent correlations in one table.
lavaan::lavInspect(fit, "cor.lv")["verbal", "deduction"]^2
#> [1] 0.5327517

# An interval comes from ci_method = "percentile", which resamples the
# cases and refits the model once per replication. That refitting is
# why it is not run here; the call is
#   average_variance_extracted(fit, ci_method = "percentile",
#                              B = 1000, seed = 113)
# and a reported interval deserves the default B = 1000 or more.

# The broom verbs: one row per factor.
generics::tidy(ave_tbl)
#>        term  estimate ci_lower ci_upper
#> 1 deduction 0.4694312       NA       NA
#> 2    verbal 0.7193173       NA       NA
generics::glance(ave_tbl)
#>   n_terms conf_level B_used
#> 1       2         NA     NA
```
