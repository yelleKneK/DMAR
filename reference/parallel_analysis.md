# Parallel Analysis of Common-Factor Eigenvalues

Compares factor-analytic eigenvalues with independently generated null
references. Each correlation matrix is reduced by replacing its diagonal
with the communalities of a fixed one-factor ML solution. The procedure
concerns common factors, not principal components. The reference
percentile is a null threshold; it is not a sampling confidence interval
for an observed eigenvalue or a confidence statement about a factor
count.

## Usage

``` r
parallel_analysis(
  data = NULL,
  S = NULL,
  N = NULL,
  items = NULL,
  missing = "fiml",
  ordered = NULL,
  reference = "normal",
  B = 1000L,
  percentile = 0.95,
  conf_level = 0.95,
  seed = NULL
)
```

## Arguments

- data:

  Numeric raw observations, with `items` selecting columns. Supply
  exactly one of `data` or `S`.

- S:

  Named positive definite covariance matrix, with `N`. Available only
  for a normal-reference analysis of continuous items.

- N:

  Reference sample size for covariance input. For raw incomplete data it
  may be supplied as an explicit benchmark. With complete data the
  default is the number of cases; with FIML it is the complete-case
  count, and the output explicitly identifies this benchmark.

- items:

  The items to analyze, by name: columns of `data`, or rows and columns
  of `S`, which is then reduced to them, in the order given; only the
  reduced matrix need be finite and positive definite. `NULL`, the
  default, analyzes every column of `data` or of `S`. Names that cannot
  be read as at least three distinct items are refused with an error of
  class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

- missing:

  `"fiml"` (default), `"ml"`, `"listwise"`, or `"pairwise"`. Continuous
  FIML obtains the unrestricted correlation matrix from
  [`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md).
  Ordered items follow the rule of
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
  and
  [`factor_count`](https://yelleknek.github.io/DMAR/reference/factor_count.md):
  with no value missing, leaving `missing` out analyzes the items as
  they stand (listwise deletion with nothing to delete, recorded as
  `"listwise"`); with a value missing, the choice of `"pairwise"` or
  `"listwise"` is the user's, and leaving `missing` out is refused;
  giving `"fiml"` or `"ml"` for ordered items is refused whether or not
  a value is missing. Each refusal has class `dmar_argument_invalid`, as
  in those functions, and comes before anything is computed.

- ordered:

  `NULL`, `TRUE`, or names of ordered items. The current ordered
  interface requires all analysis items to be ordered.

- reference:

  `"normal"` (default) generates independent normal responses. With
  ordered items these are categorized using the observed univariate
  category probabilities, then polychoric correlations are estimated.
  `"resample"` samples each item's observed responses independently with
  replacement, preserving its empirical marginal distribution in
  expectation; ordered correlations remain polychoric.

- B:

  Number of independent reference datasets, a single whole number,
  default 1000 (each one is generated and fit by a one-factor model).
  The fewest accepted depends on `percentile`: the reference quantile is
  reported only when at least 10 reference datasets lie beyond it and
  never from fewer than 100, which takes 200 at the default 95th
  percentile, 1000 at the 99th, and 100 at the 90th or at the median. A
  smaller `B` is refused before anything is computed (class
  `dmar_count_invalid`); see *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- percentile:

  Null reference percentile, default 0.95.

- conf_level:

  Confidence level for Monte Carlo uncertainty in the reference
  percentile, default 0.95. This is not observed-data sampling
  uncertainty.

- seed:

  Optional seed, default `NULL`. A supplied seed restores the caller's
  random stream on exit.

## Value

A `dmar_tbl` with ranks, observed `estimate`, `reference_quantile`,
`reference_mean`, `reference_quantile_ci_lower`,
`reference_quantile_ci_upper` (the Monte Carlo confidence limits for the
reference quantile), and `exceeds_reference`. Attributes retain
`suggested_factors`, `reference_draws`, `counts`, `observed_boundary`,
`method`, and an explicit `analysis_note`, with the counts every result
of the package that ran replications carries: `B`, `B_used`,
`n_dropped`, `B_required`, and `bootstrap_status` (see
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)).
Observed sampling confidence limits are unavailable and are not inferred
from the null distribution.

## Details

The fixed communality reduction is computed by
[`stats::factanal()`](https://rdrr.io/r/stats/factanal.html) with its
lower uniqueness bound of 0.005. Observed and null matrices use exactly
the same reduction. Defined boundary solutions are retained and counted;
failed fits are counted separately.

**Reference datasets.** The threshold of each rank is the `percentile`
quantile of the eigenvalues of the `B` reference datasets, so it is
located by the reference datasets that lie beyond it. The function runs
only when `B` puts at least 10 of them beyond the quantile, and never
with fewer than 100 (200 at the 95th percentile). A reference dataset
whose fit fails is dropped and counted, never replaced. The count is
held to the number requested. Thresholds and a factor count are reported
only when no more than half of the attempted reference datasets failed;
otherwise they are `NA`, the observed eigenvalues are returned, and the
warning has class `dmar_bootstrap_insufficient`. When some failed and
the thresholds were formed the warning has class
`dmar_bootstrap_replications_dropped`. The result records `B`, `B_used`,
`n_dropped`, `B_required`, and `bootstrap_status` as attributes beside
the fuller `counts` list. `options(dmar.bootstrap_exploratory = TRUE)`
lets a smaller `B` (2 or more) run for a first look or a demonstration;
every such call gives a warning of class `dmar_bootstrap_exploratory`
and the line under the table says so. These rules limit the error of
drawing a finite number of reference datasets; they are not a coverage
guarantee.

The suggested factor count is the number of consecutive ranks, starting
at one, whose observed eigenvalues exceed their reference thresholds. It
is evidence to combine with theory, model fit, and interpretability.
Selection itself has no confidence interval here. Monte Carlo limits for
the reference quantile use binomial order-statistic bounds; the quantile
estimate uses R's type 7 interpolation. When reference fits fail, these
quantities describe the null distribution conditional on returning a
defined fit. Discarding failures can select a different distribution;
the half-return rule does not repair that selection or establish the
intended unconditional null quantile.

An incomplete-data FIML correlation matrix uses all observed values,
whereas the complete-case reference size is only an explicitly labeled
benchmark. It is not an effective sample size and is not guaranteed to
be conservative. These references do not reproduce a study's missingness
mechanism. They do not establish general MAR calibration. If fewer than
three complete cases remain, choose a justified reference `N` explicitly
rather than silently substituting the total sample size.

## References

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. *Psychometrika, 30*, 179–185.
[doi:10.1007/BF02289447](https://doi.org/10.1007/BF02289447)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
           "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
# B = 200 is the fewest reference data sets the 95th percentile is
# reported from and keeps the example quick; the default is B = 1000.
parallel_analysis(holzinger_swineford, items = items, B = 200, seed = 113)
#>  term   rank estimate reference_quantile reference_mean
#>  rank_1 1    2.444    1.032              0.675         
#>  rank_2 2    0.568    0.170              0.122         
#>  rank_3 3    0.014    0.093              0.053         
#>  rank_4 4    -0.029   0.029              0.004         
#>  rank_5 5    -0.208   -0.019             -0.052        
#>  rank_6 6    -0.377   -0.074             -0.127        
#>  reference_quantile_ci_lower reference_quantile_ci_upper exceeds_reference
#>  1.031                       1.037                        TRUE            
#>  0.164                       0.208                        TRUE            
#>  0.083                       0.107                       FALSE            
#>  0.025                       0.039                       FALSE            
#>  -0.023                      -0.011                      FALSE            
#>  -0.085                      -0.063                      FALSE            
#> 
#> Common-factor parallel analysis: fixed one-factor ML communality reduction; normal reference; N = 301; 200 of 200 reference datasets usable. Reference percentile: 95%. Its Monte Carlo limits have nominal 95% coverage conditional on defined reference fits; discarded failures can select a different null distribution. They are not sampling CIs for the observed eigenvalues or the suggested factor count. Defined boundary reference fits: 92.
```
