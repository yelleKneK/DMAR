# Print a Model Comparison or ANOVA Table With DMAR *p*-value Formatting

Pretty-print an ANOVA-like object (the output of
[`stats::anova`](https://rdrr.io/r/stats/anova.html),
[`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html),
`lmerTest::anova`, etc.) with *p*-values formatted at a fixed number of
decimal places (default 4) and with a “\< 10^(-digits_p)” floor for
values too small to express. The default behavior of `print.anova`
routes *p*-values through `stats::format.pval`, which applies its own
digit rule (`max(1L, getOption("digits") - 2L)`) and switches to
scientific notation for tiny values. `print_anova()` sidesteps that by
converting the *p*-value columns to character strings up front and
printing as a data frame.

## Usage

``` r
print_anova(x, digits_p = 4L)
```

## Arguments

- x:

  An ANOVA-like data frame with one or more `Pr(...)` columns. Accepts
  `anova` objects from
  [`stats::anova`](https://rdrr.io/r/stats/anova.html),
  [`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html),
  [`car::Manova`](https://rdrr.io/pkg/car/man/Anova.html), and
  `lmerTest::anova`.

- digits_p:

  Integer number of decimal places for the *p*-value column(s). Default
  `4L`.

## Value

The input `x`, invisibly and unchanged.

## Details

The returned object is the input `x` invisibly, unchanged: the
underlying numeric *p*-values retain full precision and can still be
indexed (for example as `x[["Pr(>F)"]]`).

Any column whose name starts with `Pr(` is formatted as a *p*-value
column. Other columns print at whatever `getOption("digits")` dictates
(so set `options(digits = 4)` for a uniformly compact display).

## See also

[`format_p`](https://yelleknek.github.io/DMAR/reference/format_p.md),
[`print_summary`](https://yelleknek.github.io/DMAR/reference/print_summary.md).

## Author

Ken Kelley

## Examples

``` r
fit <- lm(weight ~ Time + Diet, data = ChickWeight)
print_anova(anova(fit))
#> Analysis of Variance Table
#> 
#> Response: weight
#> 
#>            Df    Sum Sq     Mean Sq    F value   Pr(>F)
#> Time        1 2042343.7 2042343.749 1576.45969 < 0.0001
#> Diet        3  129876.1   43292.019   33.41657 < 0.0001
#> Residuals 573  742336.1    1295.526         NA     <NA>

print_anova(car::Anova(fit, type = "III"))
#> Anova Table (Type III tests)
#> 
#> Response: weight
#> 
#>                 Sum Sq  Df    F value   Pr(>F)
#> (Intercept)   13689.64   1   10.56687   0.0012
#> Time        2016357.15   1 1556.40096 < 0.0001
#> Diet         129876.06   3   33.41657 < 0.0001
#> Residuals    742336.12 573         NA     <NA>

# Underlying numeric p-values are untouched:
a <- anova(fit)
print_anova(a)
#> Analysis of Variance Table
#> 
#> Response: weight
#> 
#>            Df    Sum Sq     Mean Sq    F value   Pr(>F)
#> Time        1 2042343.7 2042343.749 1576.45969 < 0.0001
#> Diet        3  129876.1   43292.019   33.41657 < 0.0001
#> Residuals 573  742336.1    1295.526         NA     <NA>
a[["Pr(>F)"]]   # full-precision doubles
#> [1] 1.226523e-166  6.473189e-20            NA
```
