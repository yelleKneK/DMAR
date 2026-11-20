# Publication-Ready Display of DMAR Result Tables

The [`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
print layer formats a result table for the console. These helpers carry
the same formatting into the two places a researcher writes up an
analysis: a knitted report and a results sentence.

## Usage

``` r
# S3 method for class 'dmar_tbl'
knit_print(x, ...)

as_kable(x, ...)

# S3 method for class 'dmar_tbl'
as_kable(x, format = NULL, ...)

results_sentence(x, label = NULL, digits = 2, term = NULL, conf_level = NULL)
```

## Arguments

- x:

  A `dmar_tbl` object (a `data.frame` returned by a DMAR function). For
  `results_sentence`, a `dmar_tbl` that carries a confidence interval,
  either as `ci_lower` / `ci_upper` rows of a long table or as
  `ci_lower` / `ci_upper` columns of a one-row wide table, such as one
  row of
  [`common_method_adjusted_r`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md)
  or of a `tidy()` result, or as `<estimate>_ci_lower` /
  `<estimate>_ci_upper` rows beside the estimate's row, as in a
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
  table, or as `<estimate>_ci_lower` / `<estimate>_ci_upper` columns
  beside the estimate's column of a one-row wide table, as in one row of
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md)
  or
  [`ci_omega_squared`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md).

- ...:

  Additional arguments. For `knit_print.dmar_tbl` and
  `as_kable.dmar_tbl`, passed to
  [`kable`](https://rdrr.io/pkg/knitr/man/kable.html) (other than
  `caption`, which `as_kable` writes itself).

- format:

  Passed to [`kable`](https://rdrr.io/pkg/knitr/man/kable.html) as its
  `format` argument (for example `"html"`, `"latex"`, `"pipe"`). The
  default `NULL` lets knitr choose based on the output context.

- label:

  For `results_sentence`, the label that leads the sentence (for example
  `"Cohen's d"`). The default `NULL` uses the name of the point-estimate
  term: the estimate's row of a long table, the quantity a `tidy()` row
  names in its `estimand` column (see
  [`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)),
  the row's `term` beside a wide table's `estimate` or `value` column,
  and otherwise the estimate column's name.

- digits:

  For `results_sentence`, the number of decimal places for the estimate
  and the interval limits. Default 2. Like the `digits` of the table
  display (see
  [`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)),
  it counts decimal places, but it keeps a default of its own and does
  not follow `options(dmar.digits = )`, so a sentence written for a
  manuscript does not change with the console display. As in the table,
  a value that rounds to zero is written without a minus sign (0.00, not
  -0.00).

- term:

  For `results_sentence`, the estimate whose interval a table holds in
  `<term>_ci_lower` and `<term>_ci_upper` rows, for example
  `"gamma_hat"` in a
  [`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
  table or `"rmsea"` in a
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
  table, or in columns of those names, for example `"eta_squared_total"`
  in a row of
  [`ci_eta_squared`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md).
  The default `NULL` reads the table's interval as described above,
  taking the first estimate with such rows when the table holds its
  interval in no other form, and, of such columns, the one estimate
  whose limits are not `NA` (a row with several is refused unless `term`
  names one). Either way the interval for the standardized mean
  difference that
  [`cles`](https://yelleknek.github.io/DMAR/reference/cles.md) and
  [`nnt_from_smd`](https://yelleknek.github.io/DMAR/reference/nnt_from_smd.md)
  echo beside their own, which theirs are computed from, is taken only
  when `term` is `"smd"`.

- conf_level:

  For `results_sentence`, the confidence level of the table's interval,
  for a table that records none (one assembled by hand, say). The
  default `NULL` reads the level the table records and refuses a table
  that records none, rather than assuming one; a value that differs from
  the level the table records is refused.

## Value

`knit_print.dmar_tbl` returns a `knit_asis` object (the rendered table)
for knitr to place in the document. `as_kable` returns a `knitr_kable`
object. `results_sentence` returns a length-one character string.

## Details

`knit_print.dmar_tbl` is the knitr print method, so a `dmar_tbl` dropped
into an R Markdown chunk renders as a formatted
[`kable`](https://rdrr.io/pkg/knitr/man/kable.html) (values to a fixed
number of decimal places, whole-number sample sizes without decimals,
*p*-values to fixed decimals) rather than as unrounded doubles.
`as_kable` is the explicit form of the same rendering: it returns the
`knitr_kable` object so the caller can pipe it into further styling or
embed it in a larger document. Both show the table the console shows,
formatted by
[`format.dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
(for most results the whole table; for a few, such as
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
which leaves out the columns that carry nothing for the call, and
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
which shows the first rows of each factor, the same view the console
prints), together with every line the console writes with it: the
sum-of-squares type, the confidence level, how an omega's denominator
was estimated and its largest residual correlation, the bootstrap
replication counts, the missing data a sample size plan assumed, and the
lines a result of its own adds (which version of a fit measure was
reported and why, what a model comparison checked, how a fraction of
missing information was estimated, which rule set a table of degrees of
freedom). What a reader sees in a report therefore matches what they saw
at the console, and neither rounds the stored numbers.
`knit_print.dmar_tbl` sets each of those lines as its own paragraph
beneath the table (a line that introduces the table, as the summary of
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
the measurement structure of
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
and the model summary of
[`mediation_mbco`](https://yelleknek.github.io/DMAR/reference/mediation_mbco.md)
do, comes before it, and a small table among the lines is set as a code
block so its columns stay aligned), and sets the whole off with a blank
line before and after it, so that two tables printed in one chunk stay
apart; `as_kable` joins them into the table's caption. Both write those
lines in the markup of the format knitr renders the table in, escaped
for LaTeX (a document knitted from `.Rnw`, or `format = "latex"`) and
for HTML, so that the percent sign of a confidence level or an
underscore in an item name prints as itself; knitr escapes a table's
cells and column names but writes a caption as it is given.

`results_sentence` turns a table that carries a confidence interval into
the one sentence an author puts in a results section, for example “smd =
0.50, 95% CI \[0.10, 0.90\]”. It reads the estimate and its limits from
the numeric columns at full precision and formats them for the sentence,
so the reported numbers are exact to the requested decimals rather than
transcribed from the rounded console display. In a wide table the
estimate is the column the row names in an `estimand` column, as the
`glance()` row of a long table names the column of the estimate its
interval belongs to (see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)),
else the `estimate` (or `value`) column when there is one, and otherwise
the first numeric column that is not a limit, a standard error, a test
statistic with its degrees of freedom and *p*-value, or a count; a
*p*-value the table carries is appended, to the four decimals of
[`format_p`](https://yelleknek.github.io/DMAR/reference/format_p.md),
and one below that floor as the bound it is (“p \< 0.0001”). In a long
table, whose interval is held in `ci_lower` and `ci_upper` rows, the
estimate is the row the table names as the one the interval belongs to,
when it names one (the corrected correlation of
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md),
whose table lists the observed correlation first, and the proportion of
superiority of
[`proportion_of_superiority`](https://yelleknek.github.io/DMAR/reference/proportion_of_superiority.md),
listed after the standardized mean difference it is computed from), and
otherwise its first row that is not a limit. Its *p*-value is likewise
the row the table names as the *p*-value of its test, when it names one,
and otherwise its first *p*-value row: the equivalence tables of
[`equivalence_r`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_smd`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
and
[`equivalence_c`](https://yelleknek.github.io/DMAR/reference/equivalence_c.md)
list the *p*-values of the two one-sided tests first and name the TOST
*p*-value, the test their \\1 - 2\alpha\\ interval decides. A long table
may instead hold an estimate's interval in rows of its own, named for
the estimate:
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
reports `rmsea` beside `rmsea_ci_lower` and `rmsea_ci_upper`, and
`gamma_hat` and `adjusted_gamma_hat` the same way. The sentence then
reports the estimate named by `term`, or the first such estimate in the
table (the RMSEA there), with no *p*-value, since the *p*-values of such
a table belong to other tests (of exact, close, and not close fit).
`term` also reaches the RMSEA rows of a
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
table. The confidence level is the interval's own level when the table
records it
([`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
records the level of its RMSEA interval, its `rmsea_conf_level`, in the
`rmsea_ci_level` row and in an attribute of the same name, which is read
when a selection of rows leaves that row out; a table stacked with
[`rbind()`](https://rdrr.io/r/base/cbind.html) keeps the attributes of
its first table only, so keep the `rmsea_ci_level` row in each selection
when stacking tables fit at different levels), and otherwise the table's
`conf_level` attribute, or its `conf_level` row or column (as the
reliability coefficients and the `tidy()` and `glance()` verbs carry
it). The `tidy()` of a long table holds the level in a `conf_level`
column on the rows of its interval, the estimate and the two limits, and
`NA` on the others; the sentence reads the level and the estimate from
those rows, so that rows of several `tidy()` results stacked with
[`rbind()`](https://rdrr.io/r/base/cbind.html) give the sentence of
their own table (a frame's attributes are those of the first frame of a
stack). Every DMAR table that forms an interval records its level, and a
sentence states no level the table does not record: a table that records
none, such as one assembled by hand, is refused unless `conf_level`
gives the level its interval was formed at, and a `conf_level` that
differs from the level the table records is refused.

A table whose two limits are not a confidence interval gets a sentence
that says so.
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md),
[`cliff_delta`](https://yelleknek.github.io/DMAR/reference/cliff_delta.md),
and
[`probability_of_superiority_paired`](https://yelleknek.github.io/DMAR/reference/probability_of_superiority_paired.md)
report both limits equal to the estimate when the estimated variance is
0 (completely separated groups, or no variation) and record the case in
the table's `interval_degenerate` attribute. The sentence then gives the
estimate and the reason and no interval, “A = 1.00 (no confidence
interval: the groups are completely separated, so the estimated variance
is 0)” or “A = 0.50 (no confidence interval: every score is equal, so
the estimated variance is 0)”, with “every pair favors the same member”
and “every pair is tied” for the paired function. A *p*-value is
appended only when the table has one (the permutation test of
[`vargha_delaney_A`](https://yelleknek.github.io/DMAR/reference/vargha_delaney_A.md)).
The `tidy()` and `glance()` rows of these tables carry the mark in a
logical column `interval_degenerate`, which says that the variance is 0
and not which case gave it, so their sentence reads “A = 1.00 (no
confidence interval: the estimated variance is 0)”.

The level a sentence states is the interval's coverage, which for an
interval formed from tail areas `alpha_lower` and `alpha_upper` is \\1 -
\alpha\_{lower} - \alpha\_{upper}\\ whatever the split. A sentence in
the form “95% CI \[L, U\]” is read as the equal-tailed interval, so an
interval whose tail areas differ says so, “smd = 0.50, 95% CI \[0.03,
0.85\] (unequal tail areas: 0.01 below, 0.04 above)”, and an interval
with a tail area of 0 is written as the one-sided confidence bound it
is, “smd = 0.50, one-sided 95% upper confidence bound 0.83”, since its
other limit is infinite or the bound of the parameter (0 for an RMSEA or
a squared multiple correlation) and carries no information. The split is
read from the table's record of it, never from a limit: an equal-tailed
interval whose lower limit is 0 because the statistic fell below its
critical value (Steiger, 2004) is written as the two-sided interval it
is. The line beneath a printed table states the same split, for example
“Confidence level: 95% (unequal tail areas: 0.01 below the lower limit,
0.04 above the upper limit)”, and the `tidy()` and `glance()` rows of
such an interval carry it in `alpha_lower` and `alpha_upper` columns
(see
[`dmar_tidiers`](https://yelleknek.github.io/DMAR/reference/dmar_tidiers.md)).

A confidence set that is not an interval is written as what it is. The
Fieller set of
[`ci_ratio_of_means`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means.md)
and
[`ci_ratio_of_means_paired`](https://yelleknek.github.io/DMAR/reference/ci_ratio_of_means_paired.md)
can be the complement of an interval or the whole line, and with unequal
tail areas an interval unbounded on one side or the union of two
intervals. The sentence reads the form from the table's `confidence_set`
attribute (or the `confidence_set` column of its `tidy()` row) and the
ends of the excluded gap from its `excluded_lower` and `excluded_upper`
rows (or columns), and writes, for example, “ratio = 0.32, 99.5%
confidence set: every value except those between 0.90 and 147.05”,
“ratio = 0.32, 99.9% confidence set: the whole line, which excludes no
value”, or “ratio = 0.32, 97.8% confidence set: the union of \[-0.61,
0.93\] and \[9.55, Inf) (unequal tail areas: 0.02 below, 0.002 above)”,
since “99.5% CI \[-Inf, Inf\]” would read as the whole line and hide the
gap. An infinite end is written open and a finite one closed. A bounded
Fieller interval is written as any other interval.

A table whose limits are `NA`, because no interval was formed (as from
[`krippendorff_alpha`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md)
without a bootstrap), has no sentence to write and gives an error.

## References

Steiger, J. H. (2004). Beyond the *F* test: Effect size confidence
intervals and tests of close fit in the analysis of variance and
contrast analysis. *Psychological Methods, 9*(2), 164–182.
[doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)

## See also

[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md) for
the console print layer and
[`format_p`](https://yelleknek.github.io/DMAR/reference/format_p.md) for
the *p*-value convention these helpers reuse.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)

# A knitr_kable that keeps every column, ready for a report.
as_kable(x)
#> 
#> 
#> Table: Confidence level: 95%
#> 
#> |term     |value |
#> |:--------|:-----|
#> |ci_lower |0.101 |
#> |smd      |0.500 |
#> |ci_upper |0.897 |

# The sentence an author writes in a results section.
results_sentence(x, label = "Cohen's d")
#> [1] "Cohen's d = 0.50, 95% CI [0.10, 0.90]"

# Wide tables (one interval per row) work the same way.
results_sentence(ci_R2(R2 = 0.25, N = 100, p = 5))
#> [1] "R2 = 0.25, 95% CI [0.08, 0.37]"

# So does one row of a tidy() result, here the same interval.
results_sentence(generics::tidy(x))
#> [1] "smd = 0.50, 95% CI [0.10, 0.90]"

# A table that holds an interval in rows named for its estimate: the
# RMSEA of a model with chi square 85.3 on 24 df and N = 301, with its
# 90 percent interval (the default level of every RMSEA interval), and
# then gamma hat.
fm <- sem_fit_measures(chisq = 85.3, df = 24, N = 301, p = 9)
results_sentence(fm, label = "RMSEA", digits = 3)
#> [1] "RMSEA = 0.092, 90% CI [0.072, 0.114]"
results_sentence(fm, term = "gamma_hat", digits = 3)
#> [1] "gamma_hat = 0.957, 90% CI [0.935, 0.973]"
```
