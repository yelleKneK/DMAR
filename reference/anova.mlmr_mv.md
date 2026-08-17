# Compare Nested Multivariate FIML Regression Fits

Compares two or more nested
[`mlmr_mv`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md) fits
with the likelihood ratio test, delegating the chi square computation to
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html). The
models must be nested (every predictor in the more restricted model is
also in the more general one) and must be fit to the same outcomes and
the same data with the same missing data handling. This is the
multivariate counterpart of
[`anova.mlmr`](https://yelleknek.github.io/DMAR/reference/anova.mlmr.md):
because every outcome is regressed on the shared predictor set, dropping
a predictor drops its slope on every outcome, so the nesting constraint
sets that slope to zero across all outcomes at once.

## Usage

``` r
# S3 method for class 'mlmr_mv'
anova(object, ...)
```

## Arguments

- object:

  An `mlmr_mv` fit.

- ...:

  Additional `mlmr_mv` fits, nested with respect to `object`.

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
contributes to the likelihood, are unchanged. That comparison is
accepted. Fits run on genuinely different data are refused.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Both fits ask for the Wald interval and skip the effect size block,
# since the likelihood ratio test needs neither and each costs refits.
fit1 <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
                  t5_general_information,
                data = holzinger_swineford,
                ci_method = "wald", effect_sizes = FALSE)
fit2 <- mlmr_mv(cbind(t6_paragraph_comprehension, t9_word_meaning) ~
                  t5_general_information + t7_sentence,
                data = holzinger_swineford,
                ci_method = "wald", effect_sizes = FALSE)
anova(fit1, fit2)
#> Likelihood ratio test for nested mlmr_mv fits
#> Model 2: cbind(t6_paragraph_comprehension, t9_word_meaning) ~ t5_general_information +      t7_sentence
#> Model 1: cbind(t6_paragraph_comprehension, t9_word_meaning) ~ t5_general_information
#>         Df    AIC    BIC  Chisq Chisq diff Df diff Pr(>Chisq)    
#> Model 2  0 7131.8 7183.7   0.00                                  
#> Model 1  2 7232.4 7276.9 104.65     104.65       2  < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
