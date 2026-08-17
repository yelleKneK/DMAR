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
# Both fits ask for the Wald interval and skip the effect size block,
# since the likelihood ratio test needs neither and each costs refits.
fit1 <- mlmr(t6_paragraph_comprehension ~ t5_general_information,
             data = holzinger_swineford, ci_method = "wald",
             effect_sizes = FALSE)
fit2 <- mlmr(t6_paragraph_comprehension ~ t5_general_information +
               t9_word_meaning,
             data = holzinger_swineford, ci_method = "wald",
             effect_sizes = FALSE)
anova(fit1, fit2)
#> Likelihood ratio test for nested mlmr fits
#> Model 2: t6_paragraph_comprehension ~ t5_general_information + t9_word_meaning
#> Model 1: t6_paragraph_comprehension ~ t5_general_information
#>         Df    AIC    BIC Chisq Chisq diff Df diff Pr(>Chisq)    
#> Model 2  0 5601.5 5634.9  0.00                                  
#> Model 1  1 5661.3 5691.0 61.78      61.78       1   3.84e-15 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
