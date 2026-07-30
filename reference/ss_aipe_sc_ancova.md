# Sample Size Planning From the AIPE Perspective for Standardized ANCOVA Contrasts

Sample size planning from the accuracy in parameter estimation (AIPE)
perspective for standardized ANCOVA contrasts.

## Usage

``` r
ss_aipe_sc_ancova(
  psi = NULL,
  sigma_anova = NULL,
  sigma_ancova = NULL,
  psi_standardized = NULL,
  ratio = NULL,
  rho = NULL,
  divisor = "s_ancova",
  c_weights,
  width,
  conf_level = 0.95,
  alpha_lower = NULL,
  alpha_upper = NULL,
  assurance = NULL,
  ...
)
```

## Arguments

- psi:

  The population unstandardized ANCOVA (adjusted) contrast

- sigma_anova:

  The population error standard deviation of the ANOVA model

- sigma_ancova:

  The population error standard deviation of the ANCOVA model

- psi_standardized:

  The population standardized ANCOVA (adjusted) contrast

- ratio:

  The ratio of `sigma_ancova` over `sigma_anova`

- rho:

  The population correlation coefficient between the response and the
  covariate

- divisor:

  Which error standard deviation to be used in standardizing the
  contrast; the value can be either `"s_ancova"` or `"s_anova"`

- c_weights:

  Contrast weights

- width:

  The desired full width of the obtained confidence interval

- conf_level:

  The desired confidence interval coverage (i.e., 1 - Type I error
  rate). Default is `.95`, which gives a symmetric two-sided interval.
  Specify either `conf_level` or both of `alpha_lower` and
  `alpha_upper`, not both.

- alpha_lower:

  Lower-tail Type I error rate, used to plan an asymmetric confidence
  interval. When supplied together with `alpha_upper`, the planned
  interval has lower-tail probability `alpha_lower` and upper-tail
  probability `alpha_upper`. Set `conf_level = NULL` when supplying
  these.

- alpha_upper:

  Upper-tail Type I error rate, used together with `alpha_lower` to plan
  an asymmetric confidence interval.

- assurance:

  Parameter to ensure that the obtained confidence interval width is
  narrower than the desired width with a specified degree of certainty
  (must be `NULL` or between zero and unity)

- ...:

  Allows one to potentially include parameter values for inner functions

## Value

A 1-row `data.frame` with columns `term` and `value`. The `term` is
`"necessary_n_per_group"` and `value` is the per-group sample size
needed for the planned ANCOVA contrast.

## Details

The sample size planning method this function is based on is developed
in the context of simple (i.e., one-response-one-covariate) ANCOVA model
and randomized design (i.e., same population covariate mean across
groups).

An ANCOVA contrast can be standardized in at least two ways: (a) divided
by the error standard deviation of the ANOVA model, (b) divided by the
error standard deviation of the ANCOVA model. This function can be used
to analyze both types of standardized ANCOVA contrasts.

Not all of the effect size arguments need to be specified. When
`divisor="s_ancova"` the input is either (a) `psi_standardized`, or (b)
`psi` (the unstandardized ANCOVA contrast) and `sigma_ancova`. When
`divisor="s_anova"`, the valid input combinations are (a)
`psi_standardized` and `ratio`; (b) `psi_standardized` and `rho`; or (c)
`psi`, `sigma_anova`, and `sigma_ancova`.

## Note

When `divisor="s_anova"` and the argument `assurance` is specified, the
necessary sample size *per group* returned by the function with
`assurance` specified is slightly underestimated. The method to obtain
exact sample size in the above situation has not been developed yet. A
practical solution is to use the sample size returned as the starting
value to conduct a priori Monte Carlo simulations with function
[`ss_aipe_sc_ancova_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md),
as discussed in Lai & Kelley (2012).

## References

Kelley, K. (2007). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24.
[doi:10.18637/jss.v020.i08](https://doi.org/10.18637/jss.v020.i08)

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.
[doi:10.1037/1082-989X.11.4.363](https://doi.org/10.1037/1082-989X.11.4.363)

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for
ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence
intervals. *British Journal of Mathematical and Statistical Psychology,
65*, 350–370.
[doi:10.1111/j.2044-8317.2011.02029.x](https://doi.org/10.1111/j.2044-8317.2011.02029.x)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval
estimation and the evaluation of statistical methods. In L. L. Harlow,
S. A. Mulaik, & J. H. Steiger (Eds.), *What if there were no
significance tests?* (pp. 221–257). Mahwah, NJ: Lawrence Erlbaum.

## See also

[`ss_aipe_sc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc.md),
[`ss_aipe_sc_ancova_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sc_ancova_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
ss_aipe_sc_ancova(psi_standardized = .8, width = .5, c_weights = c(.5, .5, 0, -1))
#>  term                  value
#>  necessary_n_per_group 98   
#> 
#> Confidence level: 95%

ss_aipe_sc_ancova(psi_standardized = .8, ratio = .6, width = .5,
                  c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#>  term                  value
#>  necessary_n_per_group 39   
#> 
#> Confidence level: 95%

ss_aipe_sc_ancova(psi_standardized = .5, rho = .4, width = .3,
               c_weights = c(.5, .5, 0, -1), divisor = "s_anova")
#>  term                  value
#>  necessary_n_per_group 221  
#> 
#> Confidence level: 95%
```
