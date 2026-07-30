# Likelihood Ratio Test for Nested Mlmr Fits

Compares two or more nested
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) fits with
the likelihood ratio test, delegating the chi square computation to
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html). The
models must be nested (every parameter in the more restricted model is
also in the more general one) and must be fit to the same data with the
same missing data handling.

## Usage

``` r
# S3 method for class 'mlmr'
anova(object, ...)
```

## Arguments

- object:

  An `mlmr` fit.

- ...:

  Additional `mlmr` fits, nested with respect to `object`.

## Value

A `data.frame` of class `anova` reporting the degrees of freedom, AIC,
BIC, log-likelihood, chi square test statistic, and *p*-value for each
consecutive pairwise comparison.

## Details

"The same data" means the same observations, with the variables the fits
share holding the same values. It does not mean the same number of
complete cases. Under `missing = "fiml"` a legitimate nested comparison
routinely has different complete-case counts: adding a predictor that is
itself incompletely observed lowers the number of rows complete on every
modeled variable, yet the observations, and the information each of them
contributes to the likelihood, are unchanged. Comparing `y ~ x1` with
`y ~ x1 + x2` when `x2` has missing values is exactly the comparison
full information maximum likelihood exists to support, and it is
accepted. Fits run on genuinely different data are refused.

Each smaller model is refit as a constrained version of the largest one
on the largest fit's data, with the slopes of its absent predictors
fixed at zero. That keeps the comparison on a single joint observed
variable set, which is what makes the chi square difference
interpretable when the predictors are modeled as random
(`fixed_x = FALSE`, the
[`mlmr`](https://yelleknek.github.io/DMAR/reference/mlmr.md) default).

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
fit1 <- mlmr(mpg ~ wt,       data = mtcars, ci_method = "wald",
             effect_sizes = FALSE)
fit2 <- mlmr(mpg ~ wt + hp,  data = mtcars, ci_method = "wald",
             effect_sizes = FALSE)
anova(fit1, fit2)
#> Likelihood ratio test for nested mlmr fits
#> Model 2: mpg ~ wt + hp
#> Model 1: mpg ~ wt
#>         Df    AIC    BIC  Chisq Chisq diff Df diff Pr(>Chisq)    
#> Model 2  0 597.22 610.41  0.000                                  
#> Model 1  1 606.59 618.32 11.377     11.377       1  0.0007436 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
