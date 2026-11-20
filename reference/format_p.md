# Format *p*-Values for Display at a Fixed Number of Decimal Places

Render a vector of *p*-values as character strings at a fixed number of
decimal places (default 4) with a floor label “\< 10^(-digits_p)” for
values too small to express. Never uses scientific notation. This is the
package-wide convention for displaying *p*-values used by the
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md)
print layer, by
[`correlations_test`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
and by the helper display functions
[`print_anova`](https://yelleknek.github.io/DMAR/reference/print_anova.md)
and
[`print_summary`](https://yelleknek.github.io/DMAR/reference/print_summary.md).
Exposing it as a user-facing function lets analysts apply the same
convention to ad hoc *p*-values that they want to report in prose or in
a manually constructed table.

## Usage

``` r
format_p(p, digits_p = 4L)
```

## Arguments

- p:

  Numeric vector of *p*-values. `NA` values pass through as
  `NA_character_`. The function only formats: it shows each value as
  given and does not check that it is a probability. A value below 0
  (and `-Inf`) is below the floor and shows as the floor label, a value
  above 1 is shown at `digits_p` decimal places (1.1 as “1.1000”), `Inf`
  is shown as Inf (with the leading spaces
  [`formatC`](https://rdrr.io/r/base/formatc.html) gives it), and `NaN`
  is returned as `NA_character_`.

- digits_p:

  Number of decimal places, a single whole number from 1 to 324, the
  most decimal places [`formatC`](https://rdrr.io/r/base/formatc.html)
  writes (`3` and `3L` are accepted; `3.9` is refused with a condition
  of class `dmar_argument_invalid`, rather than cut to 3). Default `4L`.
  Values strictly below \\10^{-\mathrm{digits\\p}}\\ render as the floor
  label “\< 10^(-digits_p)” (for example “\< 0.0001” when
  `digits_p = 4`).

## Value

A character vector of the same length as `p`.

## Details

The function does not modify the input value; the underlying numeric
*p*-value retains full precision and can still be indexed out of
whatever object holds it. Only the returned display string is rounded.

## See also

[`print_anova`](https://yelleknek.github.io/DMAR/reference/print_anova.md),
[`print_summary`](https://yelleknek.github.io/DMAR/reference/print_summary.md),
[`dmar_tbl`](https://yelleknek.github.io/DMAR/reference/dmar_tbl.md).

## Author

Ken Kelley

## Examples

``` r
# Round to four decimals with a "< 0.0001" floor:
format_p(c(0.5, 0.0234, 0.0001234, 1e-10, NA))
#> [1] "0.5000"   "0.0234"   "0.0001"   "< 0.0001" NA        

# Six decimals when more precision is wanted:
format_p(0.0001234, digits_p = 6)
#> [1] "0.000123"

# Use inline in prose for a publication-style summary:
fit <- lm(weight ~ Time + Diet, data = ChickWeight)
p_time <- summary(fit)$coefficients["Time", "Pr(>|t|)"]
paste0("The Time coefficient was significant (p = ", format_p(p_time), ").")
#> [1] "The Time coefficient was significant (p = < 0.0001)."
```
