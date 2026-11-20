# Factor Count Comparison: Model Fit, Parallel Analysis, and MAP

Fits the exploratory common-factor model with one factor, two, and so
on, and sets beside each count what bears on it: the fit of that model
and the test of the previous count against it, its information criteria,
the reduced eigenvalue at that rank against its parallel analysis
reference, and Velicer's average squared partial correlation after that
many components. The table applies no cutoff and chooses nothing: the
lines beneath it state what each criterion shows, and the count is a
judgment that combines them with theory and the interpretability of the
rotated solutions, which
[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
gives at a chosen count.

## Usage

``` r
factor_count(
  data = NULL,
  S = NULL,
  N = NULL,
  items = NULL,
  max_factors = NULL,
  estimator = "ML",
  missing = "fiml",
  ordered = NULL,
  criterion = "BIC",
  conf_level = 0.9,
  B = 1000L,
  percentile = 0.95,
  seed = NULL
)

# S3 method for class 'dmar_factor_count'
print(x, digits = getOption("dmar.digits", 3L), ...)
```

## Arguments

- data:

  Numeric raw observations. Supply exactly one of `data` or `S`.

- S:

  Named positive definite covariance matrix, with `N`.

- N:

  Sample size for covariance input, or a reference sample size for the
  parallel analysis of incomplete raw data (see
  [`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)).

- items:

  The items to analyze, by name: columns of `data`, or rows and columns
  of `S`, which is then reduced to them, in the order given; only the
  reduced matrix need be finite and positive definite. `NULL`, the
  default, analyzes every column of `data` or of `S`. Names that cannot
  be read as at least three distinct items are refused with an error of
  class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)).

- max_factors:

  The largest count fitted. The default is the largest count whose model
  keeps at least one degree of freedom, the largest \\k\\ with \\(p -
  k)^2 - (p + k) \ge 2\\ for \\p\\ items; a larger count is accepted up
  to the one with zero degrees of freedom, whose model reproduces the
  covariances exactly.

- estimator, missing, ordered:

  As for
  [`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md),
  which fits each count with `rotation = "none"`, since the fit of a
  count does not depend on the rotation. The missing data method is
  resolved once, as
  [`fa_exploratory()`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
  resolves it, and the same method serves every fit, the correlations of
  the eigenvalues and the MAP, and the parallel analysis. With ordered
  items and no item value missing, leaving `missing` out fits the items
  as they stand (listwise deletion with nothing to delete), exactly as
  `missing = "listwise"` does, and the `method` attribute records
  `"listwise"`; with a value missing, ordered items require an explicit
  `"pairwise"` or `"listwise"`, and leaving `missing` out is refused
  before any correlation, fit, or reference data set is computed.

- criterion:

  The information criterion by which the proper counts are compared, as
  [`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
  compares a set of models: `"BIC"` (the default), `"SABIC"`, `"HBIC"`,
  or `"AIC"`, in either case. IBIC and SPBIC, the information matrix
  based criteria of Bollen, Ray, Zavisca, and Harden (2012), are not
  offered: lavaan's exploratory fits carry singular information matrices
  (their parameterization counts the rotational constraints), which
  [`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md)
  refuses for them.

- conf_level:

  Confidence level of the RMSEA intervals, default 0.90, the
  conventional level, whose limits decide the tests of close and not
  close fit at .05 (see
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md));
  `conf_level = 0.95` gives 95 percent intervals. It sets nothing else:
  the Monte Carlo limits of the parallel analysis reference in the
  `parallel_analysis` attribute keep the 95 percent default of
  [`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md).

- B:

  Number of reference datasets for the parallel analysis, a single whole
  number, default 1000; `0` omits the reference, and the table then
  carries the reduced eigenvalues alone. A positive `B` must be at least
  the number
  [`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
  requires at `percentile`: 200 at the default 95th percentile, 1000 at
  the 99th, 100 at the 90th or at the median (at least 10 reference
  datasets beyond the quantile and never fewer than 100). A smaller
  positive `B` is refused before any model is fit (class
  `dmar_count_invalid`); see
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- percentile:

  Reference percentile of the parallel analysis, default 0.95.

- seed:

  Optional seed for the parallel analysis; default `NULL` uses the
  current random stream. A supplied seed restores the caller's stream.

- x:

  For [`print()`](https://rdrr.io/r/base/print.html), a table returned
  by `factor_count()`.

- digits:

  For [`print()`](https://rdrr.io/r/base/print.html), the decimal places
  of the display (default `getOption("dmar.digits", 3)`; see
  [`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)).

- ...:

  For [`print()`](https://rdrr.io/r/base/print.html), passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).

## Value

A `dmar_tbl` data frame (class `dmar_factor_count`) with one row per
count: `term`, `factors`, `solution` (`"proper"`, `"improper"`,
`"no maximum"`, or `"not converged"`), `chi_square`, `df`, `p_value`,
`rmsea`, `rmsea_ci_lower`, `rmsea_ci_upper`, `cfi`, `tli`, `srmr`,
`aic`, `bic`, `sabic`, `hbic`, `criterion_delta`, `criterion_weight`,
`criterion_log_evidence` (the comparison by `criterion`),
`delta_chi_square`, `delta_df`, `p_delta`, `eigenvalue`,
`reference_quantile`, `exceeds_reference`, and `map`. Its `tidy()` keeps
the count as `estimate` and adds a `conf_level` column, the level of the
RMSEA intervals (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)),
which `glance()` also reports. Attributes retain `fits` (the lavaan
object of each proper count, `NULL` otherwise), `fit_measures` (each
count's table, as
[`fa_exploratory()`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
returns it), `parallel_analysis` (the
[`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
table, when `B > 0`), `comparison` (the
[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
table of the proper counts, when two or more were compared), `map` (the
average squared partial correlation for \\m = 0, 1, \dots, p - 2\\),
`reasons` (why a count has no proper solution), `method`, and
`analysis_note`.

## Details

**The fit of each count.** Each count is fitted by
[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
with the same estimator and missing data method, and its `fit_measures`
supply the chi square statistic with its degrees of freedom and p-value,
the RMSEA with its interval, the CFI, TLI, and SRMR, under the same
version policy (see
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md);
after FIML the indices are the corrected ones of Zhang and Savalei,
2023). The model with one factor fewer is nested in the model with one
more, and `delta_chi_square`, `delta_df`, and `p_delta` give
[`lavaan::lavTestLRT()`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html)'s
test of the previous count against this one (the Satorra-Bentler scaled
difference under `estimator = "MLR"`). The AIC, BIC, and sample size
adjusted BIC are those of
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md);
smaller is preferred, and only the differences between rows carry
meaning. Under `estimator = "ULS"` the statistic has no reference
distribution, so the p-values, the RMSEA, CFI, TLI, difference tests,
and information criteria are `NA` and the SRMR remains (see
[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)).

**Differences, weights, and evidence.** The values of a criterion mean
nothing on their own; what carries the evidence is the difference
between two counts. The proper counts are therefore compared by
`criterion` as
[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
compares a set of models fitted to the same data: `criterion_delta` is
each count's gap from the lowest value, `criterion_weight` is
\\\exp(-\Delta/2)\\ normalized over the proper counts, and, for BIC and
HBIC, which approximate \\-2\\ times the log marginal likelihood (Bollen
et al., 2012), `criterion_log_evidence` is \\\Delta/2\\, the approximate
log Bayes factor of the best count over this one, and the weights
approximate posterior probabilities under equal prior odds (Burnham &
Anderson, 2004); for AIC the weights are Akaike weights and for SABIC
they carry no further interpretation, so the log evidence column is `NA`
for both. The weights are conditional on the set of proper counts
compared, and the line under the table gives the best count, its weight,
and how far the nearest other count stands from it. Every count is
fitted to the same cases under the same missing data method and mean
structure, which
[`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
checks. The comparison needs two proper likelihood fits, so under ULS,
or with one proper count, these columns are `NA`.

A count whose model does not converge, or converges to an improper
solution (lavaan's post-estimation check fails: a residual variance
below zero, or a covariance matrix that is not positive definite), or
whose likelihood has no maximum (a residual variance runs off below zero
while the likelihood is flat or still rising, and the optimizer stopped
somewhere along the way; see
[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)),
has no fit to report: its `solution` column says so, its fit columns are
`NA`, and a line under the table names it. Such a count is informative
in its own right, since the unrestricted model with that many factors
drives an item's residual variance to or below zero on these data, but
it is not a stopping rule: the counts beyond it are fitted and reported
in the same way.
[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
refuses an improper solution for the same reason.

**Eigenvalues and their reference.** The `eigenvalue` column holds the
eigenvalues of the reduced correlation matrix, the item correlations
with the communalities of a fixed one-factor maximum likelihood solution
on the diagonal, as
[`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
forms them; the `reference_quantile` is that function's `percentile` of
the same eigenvalue in `B` independent normal reference datasets, and
`exceeds_reference` compares the two (Horn, 1965). The correlation
matrix is the one
[`parallel_analysis()`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
uses: the full information maximum likelihood estimate of
[`cov_fiml`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md) for
incomplete continuous items, the polychoric correlations for ordered
items, the sample correlations otherwise. The reference percentile is a
null threshold, not a confidence interval for an eigenvalue or a
confidence statement about the count. The parallel analysis is run when
every item is ordered or none is, which is the interface
[`parallel_analysis()`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
supports.

**The minimum average partial.** Velicer's (1976) rule extracts
principal components of the correlation matrix one at a time, partials
the first \\m\\ out of it, and averages the squared partial correlations
that remain; no component is retained after the average reaches its
minimum. The `map` column holds the average after as many components as
the row has factors; the `map` attribute holds it for every \\m\\ from 0
(the average squared correlation itself) to \\p - 2\\, and the line
under the table names the \\m\\ at which it is smallest. The rule
concerns components used as a first-stage description of the common
variance; it is reported beside the fit of the factor models, not as a
substitute for it.

**Display.** The printed table, and `knit_print()` and
[`as_kable`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
of it, show one row per count with its solution status, the chi square
statistic, degrees of freedom, and p-value, the RMSEA with its interval
(whose level the line beneath the table names), the CFI, TLI, and SRMR,
the chosen criterion with its `delta` and `weight`, the p-value of the
difference test against the previous count, the eigenvalue and its
reference, and the MAP; columns with nothing to show (the reference when
`B = 0`, the likelihood quantities under ULS, the comparison with one
proper count) are left out. The table itself keeps every column. Beneath
the table: the counts with the lowest AIC, BIC, and SABIC among the
proper solutions, the criterion comparison, how far the eigenvalues
exceed their reference, where the MAP is smallest, and the counts
without a proper solution.

## References

Bollen, K. A., Ray, S., Zavisca, J., & Harden, J. J. (2012). A
comparison of Bayes factor approximation methods including two new
methods. *Sociological Methods & Research, 41*(2), 294–324.
[doi:10.1177/0049124112452393](https://doi.org/10.1177/0049124112452393)

Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
Understanding AIC and BIC in model selection. *Sociological Methods &
Research, 33*(2), 261–304.
[doi:10.1177/0049124104268644](https://doi.org/10.1177/0049124104268644)

Horn, J. L. (1965). A rationale and test for the number of factors in
factor analysis. *Psychometrika, 30*, 179–185.
[doi:10.1007/BF02289447](https://doi.org/10.1007/BF02289447)

Velicer, W. F. (1976). Determining the number of components from the
matrix of partial correlations. *Psychometrika, 41*(3), 321–327.
[doi:10.1007/BF02293557](https://doi.org/10.1007/BF02293557)

Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI
following FIML and TS estimation with missing data. *Psychological
Methods, 28*(2), 263–283.
[doi:10.1037/met0000445](https://doi.org/10.1037/met0000445)

## See also

[`fa_exploratory`](https://yelleknek.github.io/DMAR/reference/fa_exploratory.md)
for the solution at a chosen count,
[`parallel_analysis`](https://yelleknek.github.io/DMAR/reference/parallel_analysis.md)
for the eigenvalue reference with its Monte Carlo limits,
[`compare_sem`](https://yelleknek.github.io/DMAR/reference/compare_sem.md)
for the comparison the criterion columns reproduce, and
[`sem_information_criteria`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md).

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
items <- c("t1_visual_perception", "t2_cubes", "t4_lozenges",
           "t6_paragraph_comprehension", "t7_sentence", "t9_word_meaning")
# B = 200 is the fewest reference data sets the 95th percentile of the
# parallel analysis is reported from; the default is B = 1000.
factor_count(holzinger_swineford, items = items, B = 200, seed = 113)
#> Factor counts 1 to 2 for 6 items, N = 301
#> Estimator: ML; missing data: FIML; unrotated solutions (the fit of a count does not depend on the rotation)
#> 
#>  factors solution chi_square df  p_value                rmsea   cfi   tli  srmr
#>        1   proper    103.230  9 < 0.0001 0.187 [0.155, 0.220] 0.856 0.760 0.114
#>        2   proper      2.109  4   0.7157 0.000 [0.000, 0.064] 1.000 1.011 0.008
#>      bic delta weight  p_delta eigenvalue reference   map
#>  11054.5  72.6  0.000       NA      2.444     1.032 0.110
#>  10981.9   0.0  1.000 < 0.0001      0.568     0.170 0.102
#> (each count's model: chi square test, RMSEA with 90% interval, CFI, TLI, SRMR, BIC with delta (its gap from the lowest) and weight (its relative weight among the proper solutions); p_delta: the test of the previous count against this one; eigenvalue: of the reduced correlation matrix, with its parallel analysis 95th percentile; map: Velicer's average squared partial correlation after this many components)
#> 
#> Information criteria over the proper solutions: lowest AIC at 2 factors, lowest BIC at 2 factors, lowest SABIC at 2 factors.
#> BIC comparison over the proper solutions (compare_sem()): lowest at 2 factors, weight 1.000; the nearest other count, 1 factor, is 72.6 higher (weight 0.000, log Bayes factor 36.3 for 2 factors over it).
#> Parallel analysis (normal reference, B = 200, 95th percentile, reference N = 301): the reduced eigenvalues exceed their reference through rank 2.
#> Velicer's MAP: the average squared partial correlation is smallest after 2 components (0.102; before any component is removed it is 0.164).
#> 
#> Confidence level: 90%
#> 
#> No cutoff is applied and no count is chosen: the count is a judgment that combines the fit of each model, the difference tests, the information criteria, the eigenvalues against their reference, and the MAP with theory and the interpretability of the rotated solutions, which fa_exploratory() gives at a chosen count. Each model's interval and p-values condition on its own count, and the criterion weights are conditional on the set of proper counts compared.
```
