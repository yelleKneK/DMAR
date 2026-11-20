# Printing for DMAR Result Tables

Most DMAR estimation and testing functions return a tidy `data.frame`
with a `term` column and one or more numeric columns. Because a single
numeric column often holds quantities on very different scales (for
example, whole-number degrees of freedom alongside an *F* statistic, an
effect size, and a small *p*-value), the base
[`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html) method
formats the whole column with one common format and is easily pushed
into scientific notation with many trailing digits. The `dmar_tbl` class
supplies `print` and `format` methods that format each value on its own
terms: whole numbers (such as degrees of freedom and sample sizes) print
without a decimal part, other values print to a fixed number of decimal
places, and no value is shown in scientific notation.

## Usage

``` r
# S3 method for class 'dmar_tbl'
format(
  x,
  digits = getOption("dmar.digits", 3L),
  digits_p = 4L,
  digits_fixed = 3L,
  ...
)

# S3 method for class 'dmar_tbl'
print(
  x,
  digits = getOption("dmar.digits", 3L),
  digits_p = 4L,
  digits_fixed = 3L,
  ...
)

# S3 method for class 'dmar_tbl'
rbind(..., deparse.level = 1)
```

## Arguments

- x:

  A `dmar_tbl` object (a `data.frame` returned by a DMAR function).

- digits:

  Number of decimal places for values that are not whole numbers, and
  for measured values that are (see Details): a single whole number of 0
  or more. Defaults to `getOption("dmar.digits", 3L)`. A class whose
  display lays out cells of its own, or writes numbers in the lines
  beneath the table, uses the same `digits` there.

- digits_p:

  Number of decimal places for *p*-values. Defaults to 4. A *p*-value
  below `10^(-digits_p)` prints as “\< 0.0001” (with the threshold
  tracking `digits_p`).

- digits_fixed:

  Number of decimal places for `fixed_terms` rows (information criteria
  such as AIC and BIC, and log-likelihoods). Defaults to 3.

- ...:

  For [`print()`](https://rdrr.io/r/base/print.html), additional
  arguments passed to
  [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html). For
  [`rbind()`](https://rdrr.io/r/base/cbind.html), the tables to stack,
  and any of the named arguments of
  [`rbind.data.frame`](https://rdrr.io/r/base/cbind.html)
  (`make.row.names`, `stringsAsFactors`, `factor.exclude`).

- deparse.level:

  For [`rbind()`](https://rdrr.io/r/base/cbind.html): passed to
  [`rbind.data.frame`](https://rdrr.io/r/base/cbind.html), which does
  not use it for data frames. It keeps its base R spelling because
  [`rbind()`](https://rdrr.io/r/base/cbind.html) supplies it by that
  name.

## Value

`print.dmar_tbl` returns `x` invisibly. `format.dmar_tbl` returns a
`data.frame` whose numeric columns have been formatted to character for
display. [`rbind()`](https://rdrr.io/r/base/cbind.html) of DMAR tables
returns a `dmar_tbl` when every table stacked has the same classes and
attributes, and an ordinary `data.frame` otherwise (see the section on
stacking).

## Details

The stored numeric values are never rounded; only their display changes,
so downstream arithmetic on the returned object (confidence interval
widths, further calculations) uses full precision.

The same formatting applies to every `dmar_tbl`, whether the table is
long (a `term` column beside a single `value` or `estimate` column) or
wide (a leading label column such as `term`, `effect`, or `sample_type`
beside several typed numeric columns). Each numeric column is formatted
on its own terms, so the shape of the table does not matter.

Display precision is a number of decimal places, three by default, set
by the `digits` argument or, globally, by `options(dmar.digits = )`.
(The option is dot-named because that is the R convention for package
options, for example `dplyr.width` and `knitr.table.format`; it is not a
function or argument name and so is outside the package's snake_case
rule.) A value that is not a whole number prints to that many decimal
places with its trailing zeros, so every such value in a column shows
the same number of decimal places: 13.018 prints as 13.018, 0.56 as
0.560, 233.4107 as 233.411, and 1.5 as 1.500. A value too small for
those decimals prints as 0.000 rather than in scientific notation
(0.000123 and 3.35e-52 both print as 0.000), and a value that rounds to
zero prints as 0.000 whatever its sign, never as -0.000.

A whole number prints without a decimal part, as a count (degrees of
freedom, a sample size, a number of replications) does. A measured
quantity can land on a whole number as well, at the bound of its range
or where it is fixed, and it then prints to the decimal places of the
values of its kind beside it. The display recognizes a measured value in
four ways.

- *A bound and the estimate it belongs to.* A row or column named
  `ci_lower` or `ci_upper`, or with a name ending in `_lower` or
  `_upper`, is a bound: a confidence limit (`tau2_ci_lower`,
  `rmsea_ci_upper`, `cles_ci_lower`, `smd_ci_lower`), and in the same
  way an equivalence bound (`smd_lower`, `contrast_upper`) or a
  one-sided test statistic (`t_upper`). The estimate an interval belongs
  to is the row or column its limits are named for (`tau2`, `rmsea`),
  the row of a long table that its `ci_lower` and `ci_upper` rows belong
  to (the estimate
  [`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
  reports), the column
  [`tidy`](https://generics.r-lib.org/reference/tidy.html) takes as the
  estimate of a wide table with limit columns, and a row whose own
  limits sit beside it in limit columns (the `ci_lower` and `ci_upper`
  of a coefficient table). So a lower limit of exactly 0 reads 0.000
  beside an upper limit of 0.083, and a loading fixed at 1 reads 1.000
  beside its interval.

- *A measured kind, wherever it appears.* A fit index (CFI, TLI, RMSEA,
  SRMR, and the others, with their differences between models); a power;
  a proportion, probability, or fraction (the `pct_` rows of a
  sensitivity analysis, a coverage, an assurance, a fraction of missing
  information); a mean, median, standard deviation, variance, or higher
  moment (the `mean_`, `median_`, and `sd_` rows of a sensitivity
  analysis); a test statistic, a standard error, or a noncentrality
  parameter; the delta and the weights of a model comparison; and the
  heterogeneity indices of a meta-analysis. So a CFI of 1 reads 1.000
  beside a TLI of 1.012, a power of 1 reads 1.000 beside a power of
  0.840, and a standard deviation of the interval widths of 0 reads
  0.000 beside their mean of 1.176. A name that also names a count, such
  as `delta_df` or `median_N`, is a count.

- *A column of measured values.* A column that holds values that are not
  whole is measured throughout, so an *F* of 4 reads 4.000 beside 4.270
  and a standard error of 0 reads 0.000 beside 0.050. The exception is a
  column named for a count, with `n`, `N`, `df`, or `B` as its name or
  as a part of its name (`df_2`, `N_effective`), which counts its whole
  values even when others are not whole (degrees of freedom adjusted for
  sphericity read 2 beside 1.830), and a row that records the table's
  sum-of-squares type. The `value` or `estimate` column beside a `term`
  column is not read this way, since it can hold a different quantity on
  each row.

- *A row or column the producing function names.* A function names the
  rows or columns that hold measured values of its own kind through a
  `measured_terms` attribute (the congruence coefficients of
  [`factor_congruence`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
  the loadings and intercepts of
  [`dmacs`](https://yelleknek.github.io/DMAR/reference/dmacs.md)).

Elsewhere a whole number prints without a decimal part. That covers the
counts, the codes and indicators (the sum-of-squares type, whether an
equivalence test concluded equivalence), and a value the call was given
as a whole number, since nothing in a table says which of its other
whole numbers are counts: in a long table the `value` column holds a
sample size of 301 on one row and, on another, a planning value given as
15, which prints as 15.

*p*-values are shown to a number of decimal places of their own (four by
default, set by `digits_p`), which is the conventional way to report
them. A *p*-value smaller than the smallest magnitude those decimals can
represent prints as “\< 0.0001” instead of rounding to `0.0000`. A
column is treated as holding *p*-values when it is named `p_value`,
`p.value`, `p_adjusted` (the multiplicity-adjusted case, as in
[`ci_dunnett`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md)),
or `p_chi_square` (the exact-fit test of a fitted model, as in
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md));
in a long-format table whose quantities share a single `value` column,
the rows to format this way are named by the producing function through
a `p_terms` attribute. The limits of an interval for a *p*-value, rows
named for it with `_ci_lower` and `_ci_upper` (the Monte Carlo interval
of the bootstrap *p*-value of
[`sem_exact_fit_bootstrap`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md)),
print to the same decimal places, without the floor, since a limit of 0
is a limit and not a *p*-value below the smallest that can be shown.

Information criteria such as AIC and BIC, and log-likelihoods, keep a
number of decimal places of their own, `digits_fixed` (three by
default), which `digits` does not change, so a model comparison
difference of a few points is not rounded away when the rest of a table
is shown to fewer decimals; such a value prints its decimals even when
it is a whole number. The producing function names these rows through a
`fixed_terms` attribute.

A table that carries a `bootstrap` attribute (the reliability family's
bootstrap intervals; see
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md))
prints one more line stating how many replications were attempted, how
many returned a value, and why the rest were set aside, so the
replications the interval conditions on are always on the page.

A table that carries a `denominator` attribute (the omega coefficients
of
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
and
[`bifactor_indices`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md))
prints a line naming how the total variance in the denominator was
estimated, and a table that carries a `residual_cor` attribute prints
the largest residual correlation of the single-factor fit with the item
pair it belongs to (or, for a three-item fit, says that the saturated
model leaves no residual correlation to read), so the coefficient on the
page is identified without reference to a default and the local
dependence beside it is in view. The knitr methods
([`knit_print.dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
and
[`as_kable`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md))
write the same lines with a knitted table, as they do the lines a result
adds of its own (for example which version of a fit measure
[`sem_fit_measures`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md)
reported and why), so a report shows what the console shows.

A sample size plan made for incomplete data
([`ss_aipe_reliability`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability.md)
with `fmi` above 0) prints a line naming the fraction of missing
information it assumed and the sample size complete data would need, and
a sensitivity analysis that planned for or simulated missing data
([`ss_aipe_reliability_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reliability_sensitivity.md)
with `fmi` or `missing_prop` above 0) prints a line naming both, so the
assumptions behind the number stay beside it. In the same way, an
assurance plan of
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)
made with a nonzero `true_kurtosis` prints a line naming the excess
kurtosis of the errors its assurance allows for.

To see more precision than the display shows, raise `digits`
(`print(x, digits = 8)` shows eight decimal places) or read the columns
directly, since the stored values are never rounded: `x$value` or
`x[["p_value"]]` returns the numbers at full precision.

## Using the result in your own code

You do not need to know anything about S3 classes to use a `dmar_tbl`.
It is an ordinary `data.frame` with a print method, so everything you
already do with a data frame works: `x$value` pulls the numeric column,
`x[x$term == "smd", ]` selects a row, and the full-precision numbers are
right there for any further calculation. Three common needs:

- *Read one number.* Index it like any data frame, for example
  `x$value[x$term == "ci_upper"]`. The display rounds; the stored value
  does not, so this returns the number at full precision.

- *See more (or fewer) decimal places.* Use `print(x, digits = 6)` for a
  single table, or `options(dmar.digits = 6)` for the rest of the
  session.

- *Hand the result to other tools.* `tidy(x)` returns a one-row-per-term
  table in DMAR's own column vocabulary (`term`, `estimate`, `se`,
  `statistic`, `p_value`, `ci_lower`, `ci_upper`, and so on), which is
  the convenient “wide” view for plotting or joining; `glance(x)`
  returns a one-row model-level summary. Both come from the generics
  package and need no extra setup. For a single-estimand result such as
  [`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) the
  table is already one row, so `glance()` coincides with `tidy()` (there
  are no extra model-level statistics to report); for a multi-row result
  such as [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)
  they differ.

## Stacking tables with rbind()

A `dmar_tbl` carries facts about its numbers as attributes: the
confidence level, which rows are *p*-values, and whatever else its
function recorded
([`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)
lists them). The lines printed beneath a table are written from those
attributes, so a stacked table can keep them only when they are true of
every row. [`rbind()`](https://rdrr.io/r/base/cbind.html) of DMAR tables
therefore returns

- a `dmar_tbl` of the same class and with the same attributes, when
  every table stacked is a `dmar_tbl` and all of them have identical
  classes and attributes (two calls of one function at one confidence
  level, for example);

- an ordinary `data.frame` with the same rows and columns and none of
  the DMAR attributes, in every other case: tables at different
  confidence levels, tables from different functions, or a DMAR table
  stacked with a plain data frame. No warning is given. The numbers are
  all there; what is dropped is the class and the attributes, so that no
  printed line states a confidence level (or anything else) that does
  not hold for every row.

To compare results across confidence levels, stack the `tidy()` rows
instead, which carry the level in a `conf_level` column:
`rbind(tidy(x_90), tidy(x_95))`. `NULL` arguments are ignored, as
[`rbind()`](https://rdrr.io/r/base/cbind.html) ignores them.
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
does not call this method, so what is said here does not apply to it;
use [`rbind()`](https://rdrr.io/r/base/cbind.html), or stack `tidy()`
rows.

## See also

[`tidy`](https://generics.r-lib.org/reference/tidy.html) and
[`glance`](https://generics.r-lib.org/reference/glance.html) for the
wide one-row-per-term and the one-row summary views. For a gentle,
non-technical tour of how to read and use DMAR result tables, see the
“Reading DMAR result tables” vignette:
[`vignette("dmar_output", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/dmar_output.md).

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Every DMAR estimation function returns a table that prints this way.
x <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50)
x                       # rounded for reading, to three decimal places
#>  term     value
#>  ci_lower 0.101
#>  smd      0.500
#>  ci_upper 0.897
#> 
#> Confidence level: 95%

# The stored numbers keep full precision; only the display rounds.
x$value[x$term == "smd"]
#> [1] 0.5
print(x, digits = 8)    # ask the display for more decimal places
#>  term     value     
#>  ci_lower 0.10058571
#>  smd      0.50000000
#>  ci_upper 0.89694143
#> 
#> Confidence level: 95%

# Pull a single number out, exactly as you would from a data frame.
x$value[x$term == "ci_upper"]
#> [1] 0.8969414

# The broom verbs give the programmer-friendly wide and summary views.
generics::tidy(x)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95
generics::glance(x)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95

# The same display rules apply to a wide table, one with several typed
# columns, for example an effect size with its confidence interval per
# effect.
ci_eta_squared(aov(iq_8 ~ treatment, data = pygmalion))
#>  effect    eta_squared_total eta_squared_total_ci_lower
#>  treatment 0.020             0.001                     
#>  eta_squared_total_ci_upper eta_squared_partial eta_squared_partial_ci_lower
#>  0.061                      0.020               0.001                       
#>  eta_squared_partial_ci_upper F_value df_effect df_error N  
#>  0.061                        6.335   1         308      310
#> 
#> Confidence level: 95%

# Stacking two tables computed at different confidence levels gives an
# ordinary data frame, because no one confidence level describes every
# row. The tidy() rows keep each level in a column.
x_90 <- ci_smd(smd = 0.5, n_1 = 50, n_2 = 50, conf_level = 0.90)
class(rbind(x, x_90))
#> [1] "data.frame"
rbind(generics::tidy(x), generics::tidy(x_90))
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95
#> 2  smd      0.5 0.1645872 0.8329100       0.90
```
