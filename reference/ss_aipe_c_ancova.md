# Sample Size Planning for a Contrast in Randomized ANCOVA From the Accuracy in Parameter Estimation (AIPE) Perspective

Plans the sample size per group so that the confidence interval for an
unstandardized contrast in a one-covariate randomized ANCOVA is
sufficiently narrow, following the accuracy in parameter estimation
(AIPE) approach. To the extent the covariate correlates with the
response, the covariate adjustment shrinks the error variance, so the
desired precision is reached with a smaller sample size than the
corresponding ANOVA design requires.

## Usage

``` r
ss_aipe_c_ancova(
  error_var_ancova = NULL,
  error_var_anova = NULL,
  rho = NULL,
  c_weights,
  width,
  conf_level = 0.95,
  assurance = NULL
)
```

## Arguments

- error_var_ancova:

  The population error variance of the ANCOVA model (i.e., the mean
  square within of the ANCOVA model)

- error_var_anova:

  The population error variance of the ANOVA model (i.e., the mean
  square within of the ANOVA model)

- rho:

  The population correlation coefficient of the response and the
  covariate

- c_weights:

  The contrast weights

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage, (i.e., 1 - Type I error
  rate)

- assurance:

  Parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be NULL or between zero and unity)

## Value

- sample_size:

  The necessary sample size *per group*

## Details

Either the error variance of the ANCOVA model or of the ANOVA model can
be used to plan the appropriate sample size per group. When using the
error variance of the ANOVA model to plan sample size, the correlation
coefficient of the response and the covariate is also needed.

## References

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or
obtaining precision: Delineating methods of sample size planning.
*Evaluation and the Health Professions, 26*(3), 258–287.
[doi:10.1177/0163278703255242](https://doi.org/10.1177/0163278703255242)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ss_aipe_c`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose the population error variance of some three-group ANOVA model
# is believed to be 40, and the population correlation coefficient
# of the response and the covariate is 0.22. The researcher is
# interested in the difference between the mean of group 1 and
# the average of means of group 2 and 3. To plan the sample size so
# that, with 90 percent certainty, the obtained 95 percent full
# confidence interval width is no wider than 3:

ss_aipe_c_ancova(error_var_anova = 40, rho = .22, c_weights = c(1, -0.5, -0.5),
                 width = 3, assurance = .90)
#>  term                  value
#>  necessary_n_per_group 109  
#> 
#> Confidence level: 95%
```
