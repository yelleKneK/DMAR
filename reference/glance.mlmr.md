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
fit <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
              t9_word_meaning,
            data = holzinger_swineford, ci_method = "wald")
generics::glance(fit)
#>          R2    adj_R2    sigma       f2 statistic      p_value df   logLik
#> 1 0.5372996 0.5341942 2.371619 1.161226  231.9733 4.242586e-51  9 -2791.77
#>       AIC      BIC deviance df_residual nobs
#> 1 5601.54 5634.904  5583.54         298  301
```
