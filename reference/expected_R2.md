# Expected Value of the Squared Multiple Correlation Coefficient

Computes the expected value of the observed squared multiple correlation
coefficient given the population squared multiple correlation
coefficient, the sample size, and the number of predictors. The sample
\\R^2\\ is a positively biased estimator of its population value, and
the expected value quantifies how large that bias is for a particular
design.

## Usage

``` r
expected_R2(population_R2, N, p)
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
is `"expected_value_population_R2"` and `value` is the expected value of
\\R^2\\ under random sampling.

## Details

Uses the hypergeometric function as discussed in section 28 of Stuart,
Ord, and Arnold (1999) in order to obtain the *correct* value for the
squared multiple correlation coefficient. Many times an exact value is
given that ignores the hypergeometric function. This function yields the
correct value.

## References

Kelley, K. (2008). Sample size planning for the squared multiple
correlation coefficient: Accuracy in parameter estimation via narrow
confidence intervals. *Multivariate Behavioral Research, 43*, 524–555.
[doi:10.1080/00273170802490632](https://doi.org/10.1080/00273170802490632)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 3 on \\R^2\\ as a model comparison effect
size.)

Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain
correlation coefficients. *The Annals of Mathematical Statistics,
29*(1), 201–211.

Stuart, A., Ord, J. K., & Arnold, S. (1999). *Kendall's advanced theory
of statistics, volume 2A: Classical inference and the linear model* (6th
ed.). Arnold.

## See also

[`ss_aipe_R2`](https://yelleknek.github.io/DMAR/reference/ss_aipe_R2.md),
[`ci_R2`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`var_R2`](https://yelleknek.github.io/DMAR/reference/var_R2.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
expected_R2(.5, 10, 5)
#>  term                         value
#>  expected_value_population_R2 0.754
expected_R2(.5, 25, 5)
#>  term                         value
#>  expected_value_population_R2 0.588
expected_R2(.5, 50, 5)
#>  term                         value
#>  expected_value_population_R2 0.542
expected_R2(.5, 100, 5)
#>  term                         value
#>  expected_value_population_R2 0.52 
expected_R2(.5, 1000, 5)
#>  term                         value
#>  expected_value_population_R2 0.502
expected_R2(.5, 10000, 5)
#>  term                         value
#>  expected_value_population_R2 0.5  
```
