# Variance of the Squared Multiple Correlation Coefficient

Computes the sampling variance of the squared multiple correlation
coefficient from the population value, the sample size, and the number
of predictors, the quantity that governs how precisely \\R^2\\ is
estimated at a given design size.

## Usage

``` r
var_R2(population_R2, N, p)
```

## Arguments

- population_R2:

  Population squared multiple correlation coefficient

- N:

  Sample size

- p:

  The number of predictor variables

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` value
is `"var_R2"` and `value` is the asymptotic variance of \\R^2\\.

## Details

Uses the hypergeometric function as discussed in and section 28 of
Stuart, Ord, and Arnold (1999) in order to obtain the *correct* value
for the variance of the squared multiple correlation coefficient.

## Note

The Gauss hypergeometric function \\{}\_2F_1\\ is computed in base R
(see the internal `.hyperg_2F1`); no GSL system library is required.

## References

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Stuart, A., Ord, J. K., & Arnold, S. (1999). *Kendall's advanced theory
of statistics, volume 2A: Classical inference and the linear model* (6th
ed.). Arnold.

## See also

[`expected_R2`](https://yelleknek.github.io/DMAR/reference/expected_R2.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
var_R2(.5, 10, 5)
#>  term   value 
#>  var_R2 0.0268
var_R2(.5, 25, 5)
#>  term   value 
#>  var_R2 0.0169
var_R2(.5, 50, 5)
#>  term   value  
#>  var_R2 0.00926
var_R2(.5, 100, 5)
#>  term   value  
#>  var_R2 0.00482
```
