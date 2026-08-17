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

results_sentence(x, label = NULL, digits = 2)
```

## Arguments

- x:

  A `dmar_tbl` object (a `data.frame` returned by a DMAR function). For
  `results_sentence`, a `dmar_tbl` that carries a confidence interval,
  either as `lower_limit` / `upper_limit` rows of a long table or as
  `lower_limit` / `upper_limit` columns of a wide table.

- ...:

  Additional arguments. For `knit_print.dmar_tbl` and
  `as_kable.dmar_tbl`, passed to
  [`kable`](https://rdrr.io/pkg/knitr/man/kable.html).

- format:

  Passed to [`kable`](https://rdrr.io/pkg/knitr/man/kable.html) as its
  `format` argument (for example `"html"`, `"latex"`, `"pipe"`). The
  default `NULL` lets knitr choose based on the output context.

- label:

  For `results_sentence`, the label that leads the sentence (for example
  `"Cohen's d"`). The default `NULL` uses the name of the point-estimate
  term.

- digits:

  For `results_sentence`, the number of decimal places for the estimate
  and the interval limits. Default 2.

## Value

`knit_print.dmar_tbl` returns a `knit_asis` object (the rendered table)
for knitr to place in the document. `as_kable` returns a `knitr_kable`
object. `results_sentence` returns a length-one character string.

## Details

`knit_print.dmar_tbl` is the knitr print method, so a `dmar_tbl` dropped
into an R Markdown chunk renders as a formatted
[`kable`](https://rdrr.io/pkg/knitr/man/kable.html) (sensible rounding,
whole-number sample sizes, *p*-values to fixed decimals) rather than as
a raw dump of doubles. `as_kable` is the explicit form of the same
rendering: it returns the `knitr_kable` object so the caller can pipe it
into further styling or embed it in a larger document. Both reuse
[`format.dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md),
so what a reader sees in a report matches what they saw at the console,
and neither rounds the stored numbers.

`results_sentence` turns a table that carries a confidence interval into
the one sentence an author puts in a results section, for example “smd =
0.50, 95% CI \[0.10, 0.90\]”. It reads the estimate and its limits from
the numeric columns at full precision and formats them for the sentence,
so the reported numbers are exact to the requested decimals rather than
transcribed from the rounded console display.

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
#> |term        |value |
#> |:-----------|:-----|
#> |lower_limit |0.101 |
#> |smd         |0.5   |
#> |upper_limit |0.897 |

# The sentence an author writes in a results section.
results_sentence(x, label = "Cohen's d")
#> [1] "Cohen's d = 0.50, 95% CI [0.10, 0.90]"

# Wide tables (one interval per row) work the same way.
results_sentence(ci_R2(R2 = 0.25, N = 100, p = 5))
#> [1] "R2 = 0.25, 95% CI [0.08, 0.37]"
```
