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
without a decimal part, other values print to a small number of
significant figures, and scientific notation is reserved for magnitudes
where it is the clearer choice (for example, a very small *p*-value).

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
```

## Arguments

- x:

  A `dmar_tbl` object (a `data.frame` returned by a DMAR function).

- digits:

  Number of significant figures for non-integer values. Defaults to
  `getOption("dmar.digits", 3L)`.

- digits_p:

  Number of decimal places for *p*-values. Defaults to 4. A *p*-value
  below `10^(-digits_p)` prints as “\< 0.0001” (with the threshold
  tracking `digits_p`).

- digits_fixed:

  Number of decimal places for `fixed_terms` rows (information criteria
  such as AIC and BIC, and log-likelihoods). Defaults to 3.

- ...:

  Additional arguments passed to
  [`print.data.frame`](https://rdrr.io/r/base/print.dataframe.html).

## Value

`print.dmar_tbl` returns `x` invisibly. `format.dmar_tbl` returns a
`data.frame` whose numeric columns have been formatted to character for
display.

## Details

The stored numeric values are never rounded; only their display changes,
so downstream arithmetic on the returned object (confidence interval
widths, further calculations) uses full precision.

The same formatting applies to every `dmar_tbl`, whether the table is
long (a `term` column beside a single `value` or `estimate` column) or
wide (a leading label column such as `term`, `effect`, or `sample_type`
beside several typed numeric columns). Each numeric column is formatted
on its own terms, so the shape of the table does not matter.

Display precision is controlled by the `digits` argument or, globally,
by `options(dmar.digits = )`. The default is 3 significant figures. (The
option is dot-named because that is the R convention for package
options, for example `dplyr.width` and `knitr.table.format`; it is not a
function or argument name and so is outside the package's snake_case
rule.)

*p*-values are shown to a fixed number of decimal places (four by
default, set by `digits_p`) rather than to significant figures, which is
the conventional way to report them. A *p*-value smaller than the
smallest magnitude those decimals can represent prints as “\< 0.0001”
instead of rounding to `0.0000`. A column is treated as holding
*p*-values when it is named `p_value`, `p.value`, `p_adjusted` (the
multiplicity-adjusted case, as in
[`dunnett_ci`](https://yelleknek.github.io/DMAR/reference/dunnett_ci.md)),
or `p_chi_square` (the exact-fit test of a fitted model, as in
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md));
in a long-format table whose quantities share a single `value` column,
the rows to format this way are named by the producing function through
a `p_terms` attribute.

A few quantities read better at a fixed number of decimal places than at
significant figures even though they are not *p*-values: information
criteria such as AIC and BIC, and log-likelihoods, where a model
comparison difference of a few points would be rounded away by three
significant figures (an AIC of 2284.830 would otherwise print as 2280).
The producing function names these rows through a `fixed_terms`
attribute, and they print to `digits_fixed` decimal places (three by
default).

To see more precision than the display shows, raise `digits` (for
example `print(x, digits = 8)`) or read the columns directly, since the
stored values are never rounded: `x$value` or `x[["p_value"]]` returns
the numbers at full precision.

## Using the result in your own code

You do not need to know anything about S3 classes to use a `dmar_tbl`.
It is an ordinary `data.frame` with a print method, so everything you
already do with a data frame works: `x$value` pulls the numeric column,
`x[x$term == "smd", ]` selects a row, and the full-precision numbers are
right there for any further calculation. Three common needs:

- *Read one number.* Index it like any data frame, for example
  `x$value[x$term == "upper_limit"]`. The display rounds; the stored
  value does not, so this returns the number at full precision.

- *See more (or fewer) digits.* Use `print(x, digits = 6)` for a single
  table, or `options(dmar.digits = 6)` for the rest of the session.

- *Hand the result to other tools.* `tidy(x)` returns a one-row-per-term
  table with broom-style column names (`term`, `estimate`, `conf.low`,
  `conf.high`, and so on), which is the convenient “wide” view for
  plotting or joining; `glance(x)` returns a one-row model-level
  summary. Both come from the generics package and need no extra setup.
  For a single-estimand result such as
  [`ci_smd`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) the
  table is already one row, so `glance()` coincides with `tidy()` (there
  are no extra model-level statistics to report); for a multi-row result
  such as [`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)
  they differ.

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
x                       # rounded for reading; sample sizes have no decimals
#>  term        value
#>  lower_limit 0.101
#>  smd         0.5  
#>  upper_limit 0.897
#> 
#> Confidence level: 95%

# The stored numbers keep full precision; only the display rounds.
x$value[x$term == "smd"]
#> [1] 0.5
print(x, digits = 8)    # ask the display for more digits
#>  term        value     
#>  lower_limit 0.10058571
#>  smd         0.5       
#>  upper_limit 0.89694143
#> 
#> Confidence level: 95%

# Pull a single number out, exactly as you would from a data frame.
x$value[x$term == "upper_limit"]
#> [1] 0.8969414

# The broom verbs give the programmer-friendly wide and summary views.
generics::tidy(x)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95
generics::glance(x)
#>   term estimate  ci_lower  ci_upper conf_level
#> 1  smd      0.5 0.1005857 0.8969414       0.95

# The same display rules apply to wide tables (several typed columns),
# for example an effect size with its confidence interval per effect.
ci_eta_squared(aov(len ~ supp * factor(dose), data = ToothGrowth))
#>  effect            eta_squared lower_limit upper_limit F_value df_effect
#>  supp              0.224       0.0531      0.377       15.6    1        
#>  factor(dose)      0.773       0.638       0.823       92      2        
#>  supp:factor(dose) 0.132       0.00132     0.274       4.11    2        
#>  df_error N 
#>  54       60
#>  54       60
#>  54       60
```
