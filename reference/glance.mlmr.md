# Glance at an Mlmr Fit

Returns a one-row `data.frame` of model-level summaries in the column
convention used by the broom ecosystem (`R2`, `adj_R2`, `sigma`,
`statistic`, `p_value`, `df`, `logLik`, `AIC`, `BIC`, `deviance`,
`df_residual`, `nobs`). The `statistic` and `p_value` columns report the
omnibus likelihood ratio test of all slopes equal to zero (the FIML
analog of the `lm` omnibus *F*-test).

## Usage

``` r
# S3 method for class 'mlmr'
glance(x, ...)
```

## Arguments

- x:

  An object of class `"mlmr"`.

- ...:

  Unused.

## Value

A one-row `data.frame`.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
fit <- mlmr(mpg ~ wt + hp, data = mtcars, ci_method = "wald")
generics::glance(fit)
#>          R2    adj_R2    sigma       f2 statistic     p_value df    logLik
#> 1 0.8267855 0.8148396 2.468854 4.773187  56.10318 6.56674e-13  9 -289.6083
#>        AIC      BIC deviance df_residual nobs
#> 1 597.2166 610.4082 579.2166          29   32
```
