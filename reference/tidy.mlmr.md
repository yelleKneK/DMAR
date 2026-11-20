# An Mlmr Fit

Returns a one-row-per-coefficient `data.frame` in the column convention
used by the broom ecosystem (`term`, `estimate`, `se`, `statistic`,
`p_value`, `wald_df`, and optionally `ci_lower`, `ci_upper`). `wald_df`
is the degrees of freedom of the *t* reference of the row's statistic
and *p*-value (the `wald_df` argument of
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md)), `Inf` for
the normal reference, so that rows on the two references can be told
apart after tables are stacked. Use
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on the
fit for the DMAR-style table (snake_case columns) stored at
`fit$coef_table`.

## Usage

``` r
# S3 method for class 'mlmr'
tidy(x, conf.int = FALSE, conf_level = NULL, standardized = FALSE, ...)
```

## Arguments

- x:

  An object of class `"mlmr"`.

- conf.int:

  Logical; if `TRUE`, append `ci_lower` and `ci_upper` columns using the
  confidence intervals already computed at fit time, and a `conf_level`
  column holding the level they were computed at. Defaults to `FALSE`.

- conf_level:

  Ignored; the confidence interval comes from the fit object at
  `x$conf_level`. Present for compatibility with the broom generic.

- standardized:

  Logical; if `TRUE` and effect sizes were computed at fit time, append
  a `std_estimate` column. Defaults to `FALSE`.

- ...:

  Unused.

## Value

A `data.frame`.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
fit <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
              t9_word_meaning,
            data = holzinger_swineford, ci_method = "wald")
generics::tidy(fit)
#>                     term   estimate         se statistic      p_value wald_df
#> 1            (Intercept) 2.38038350 0.47731165  4.987064 6.130392e-07     Inf
#> 2 t5_general_information 0.08482737 0.01642376  5.164917 2.405460e-07     Inf
#> 3        t9_word_meaning 0.21956217 0.02651353  8.281136 1.220057e-16     Inf
generics::tidy(fit, conf.int = TRUE)
#>                     term   estimate         se statistic      p_value wald_df
#> 1            (Intercept) 2.38038350 0.47731165  4.987064 6.130392e-07     Inf
#> 2 t5_general_information 0.08482737 0.01642376  5.164917 2.405460e-07     Inf
#> 3        t9_word_meaning 0.21956217 0.02651353  8.281136 1.220057e-16     Inf
#>     ci_lower  ci_upper conf_level
#> 1 1.44486986 3.3158971       0.95
#> 2 0.05263738 0.1170174       0.95
#> 3 0.16759660 0.2715277       0.95
generics::tidy(fit, conf.int = TRUE, standardized = TRUE)
#>                     term   estimate         se statistic      p_value wald_df
#> 1            (Intercept) 2.38038350 0.47731165  4.987064 6.130392e-07     Inf
#> 2 t5_general_information 0.08482737 0.01642376  5.164917 2.405460e-07     Inf
#> 3        t9_word_meaning 0.21956217 0.02651353  8.281136 1.220057e-16     Inf
#>     ci_lower  ci_upper conf_level std_estimate
#> 1 1.44486986 3.3158971       0.95           NA
#> 2 0.05263738 0.1170174       0.95    0.3007216
#> 3 0.16759660 0.2715277       0.95    0.4821600
```
