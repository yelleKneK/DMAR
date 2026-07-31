# An Mlmr Fit

Returns a one-row-per-coefficient `data.frame` in the column convention
used by the broom ecosystem (`term`, `estimate`, `se`, `statistic`,
`p_value`, and optionally `ci_lower`, `ci_upper`). Use
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
  confidence intervals already computed at fit time. Defaults to
  `FALSE`.

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
fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
generics::tidy(fit)
#>          term    estimate          se statistic       p_value
#> 1 (Intercept) 37.22727012 1.522000401 24.459435 3.993931e-132
#> 2          wt -3.87783074 0.602344344 -6.437897  1.211403e-10
#> 3          hp -0.03177295 0.008596028 -3.696236  2.188195e-04
generics::tidy(fit, conf.int = TRUE)
#>          term    estimate          se statistic       p_value    ci_lower
#> 1 (Intercept) 37.22727012 1.522000401 24.459435 3.993931e-132 34.24420415
#> 2          wt -3.87783074 0.602344344 -6.437897  1.211403e-10 -5.05840396
#> 3          hp -0.03177295 0.008596028 -3.696236  2.188195e-04 -0.04862085
#>      ci_upper
#> 1 40.21033609
#> 2 -2.69725752
#> 3 -0.01492504
generics::tidy(fit, conf.int = TRUE, standardized = TRUE)
#>          term    estimate          se statistic       p_value    ci_lower
#> 1 (Intercept) 37.22727012 1.522000401 24.459435 3.993931e-132 34.24420415
#> 2          wt -3.87783074 0.602344344 -6.437897  1.211403e-10 -5.05840396
#> 3          hp -0.03177295 0.008596028 -3.696236  2.188195e-04 -0.04862085
#>      ci_upper std_estimate
#> 1 40.21033609           NA
#> 2 -2.69725752   -0.6295545
#> 3 -0.01492504   -0.3614507
```
