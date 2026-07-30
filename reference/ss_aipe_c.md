# Sample Size Planning for an ANOVA Contrast From the Accuracy in Parameter Estimation (AIPE) Perspective

Plans the sample size *per group* so that the confidence interval for an
unstandardized contrast of means in a fixed effects analysis of variance
is sufficiently narrow, following the accuracy in parameter estimation
(AIPE) approach: the design goal is a contrast estimated with the
precision the research question requires, not merely one detected as
nonzero. AIPE sample size planning for ANOVA and ANCOVA contrasts is
developed in Lai and Kelley (2012).

## Usage

``` r
ss_aipe_c(
  error_variance = NULL,
  c_weights,
  width,
  conf_level = 0.95,
  assurance = NULL,
  MSwithin = NULL,
  SD = NULL,
  ...
)
```

## Arguments

- error_variance:

  The common error variance; i.e., the mean square error

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

- MSwithin:

  An alias for `error_variance`

- SD:

  The standard deviation of the common error in ANOVA model

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

- sample_size:

  the necessary sample size *per group*

## Note

Be sure to use the error variance and not its square root (i.e., the
standard deviation of the errors).

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
ed.). Routledge.

## See also

[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
[`ss_aipe_c_ancova`](https://yelleknek.github.io/DMAR/reference/ss_aipe_c_ancova.md),
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Suppose the population error variance of some three-group ANOVA model
# is believed to be 40. The researcher is interested in the difference
# between the mean of group 1 and the average of means of group 2 and 3.
# To plan the sample size so that, with 90 percent certainty, the
# obtained 95 percent full confidence interval width is no wider than 3:

ss_aipe_c(error_variance = 40, c_weights = c(1, -0.5, -0.5),
          width = 3, assurance = .90)
#>  term                  value
#>  necessary_n_per_group 114  
#> 
#> Confidence level: 95%
```
