# Kuder-Richardson Formula 20 (KR-20) With a Confidence Interval

Estimates Kuder-Richardson formula 20 (Kuder & Richardson, 1937) for a
homogeneous composite scored on dichotomous (0/1) items, and returns a
confidence interval for the population coefficient.

## Usage

``` r
reliability_kr20(
  data,
  ci_method = c("feldt", "bonett", "fisher", "hakstian_whalen", "ml", "ml_logistic",
    "adf", "adf_logistic", "bootstrap_se", "bootstrap_se_logistic", "percentile", "bca",
    "none"),
  conf_level = 0.95,
  B = 10000,
  seed = NULL
)
```

## Arguments

- data:

  A numeric matrix or data frame of 0/1 item scores (rows are
  respondents, columns are items). Rows with any missing values are
  listwise-deleted. Non-binary values trigger an error.

- ci_method:

  Method for constructing the confidence interval; see
  [`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md)
  for the full list. Defaults to `"feldt"`.

- conf_level:

  Confidence level for the interval (1 - Type I error rate). Defaults to
  `0.95`.

- B:

  Number of bootstrap replications when a bootstrap method is selected.
  Defaults to `10000`.

- seed:

  Random number seed used for bootstrap reproducibility. Defaults to
  `NULL`, which leaves the user's current RNG state intact; supply an
  integer for reproducibility.

## Value

A `data.frame` with columns `term` and `value` and rows `"estimate"`
(sample KR-20), `"se"` (standard error, `NA` for methods that do not
produce one), `"lower_limit"` and `"upper_limit"` (clamped to \[0, 1\]),
`"conf_level"`, `"N"` (effective sample size after listwise deletion),
`"N_complete"` (the complete cases; equal to `"N"` here, and carried so
the whole reliability family returns one shape), and `"J"` (number of
items). Attributes `coefficient` (`"kr20"`) and `ci_method` record the
computation; bootstrap calls also record `B`.

## Details

Kuder and Richardson's (1937) formula 20 for a \\J\\-item composite of
binary items is \$\$KR\_{20} = \frac{J}{J-1}\left(1 - \frac{\sum\_{j}
p\_{j} q\_{j}}{s\_{Y}^{2}}\right),\$\$ where \\p\_{j}\\ is the
proportion of respondents endorsing item *j* (i.e., scoring 1), \\q\_{j}
= 1 - p\_{j}\\, and \\s\_{Y}^{2}\\ is the variance of the composite
score. For dichotomous items \\p\_{j} q\_{j}\\ is the item variance, so
KR-20 is algebraically identical to coefficient \\\alpha\\ (Guttman,
1945; Cronbach, 1951) computed from the item covariance matrix. KR-20
predates coefficient \\\alpha\\ by 14 years and Feldt's (1965)
*F*-distribution interval was derived specifically for the sampling
distribution of KR-20.

Because KR-20 is a special case of coefficient \\\alpha\\ (Guttman,
1945; Cronbach, 1951), the same considerations apply: KR-20 equals the
population reliability of the composite under essential
\\\tau\\-equivalence (i.e., equal factor loadings) and serves as a lower
bound otherwise. For dichotomous items see also
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
(categorical omega), which treats the items via a probit-link
single-factor model and does not assume equal loadings.

Available confidence interval methods (set via `ci_method`) are the same
as for
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md);
see that function's *Details* for full descriptions. The default for
KR-20 is `"feldt"`, the *F*-distribution interval originally developed
for KR-20.

**Comparison with other packages.** The psych package computes the same
quantity via [`alpha`](https://rdrr.io/pkg/psych/man/alpha.html) (since
\\\alpha\\ on 0/1 data *is* KR-20). `reliability_kr20` restricts input
to raw 0/1 data so the dichotomous-items assumption cannot be quietly
violated, presents the historical formula in the documentation, and
accompanies the point estimate with a confidence interval drawn from the
methods compared in Kelley and Pornprasertmanit (2016).

## References

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. *Psychometrika, 16*(3), 297–334.

Feldt, L. S. (1965). The approximate sampling distribution of
Kuder-Richardson reliability coefficient twenty. *Psychometrika, 30*,
357–370.

Guttman, L. (1945). A basis for analyzing test-retest reliability.
*Psychometrika, 10*(4), 255–282.

Kelley, K., & Cheng, Y. (2012). Estimation of and confidence interval
formation for reliability coefficients of homogeneous measurement
instruments. *Methodology, 8*, 39–50.
[doi:10.1027/1614-2241/a000036](https://doi.org/10.1027/1614-2241/a000036)

Kelley, K., & Pornprasertmanit, S. (2016). Confidence intervals for
population reliability coefficients: Evaluation of methods,
recommendations, and software for composite measures. *Psychological
Methods, 21*, 69–92.
[doi:10.1037/a0040086](https://doi.org/10.1037/a0040086)

Kuder, G. F., & Richardson, M. W. (1937). The theory of the estimation
of test reliability. *Psychometrika, 2*, 151–160.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Terry, L. J., & Kelley, K. (2012). Sample size planning for composite
reliability coefficients: Accuracy in parameter estimation via narrow
confidence intervals. *British Journal of Mathematical and Statistical
Psychology, 65*, 371–401.
[doi:10.1111/j.2044-8317.2011.02030.x](https://doi.org/10.1111/j.2044-8317.2011.02030.x)

## See also

[`reliability`](https://yelleknek.github.io/DMAR/reference/reliability.md)
(general wrapper),
[`reliability_alpha`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md),
[`alpha`](https://rdrr.io/pkg/psych/man/alpha.html).

Other reliability:
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`diagnosis_agreement`](https://yelleknek.github.io/DMAR/reference/diagnosis_agreement.md),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`reliability()`](https://yelleknek.github.io/DMAR/reference/reliability.md),
[`reliability_H()`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
[`reliability_alpha()`](https://yelleknek.github.io/DMAR/reference/reliability_alpha.md),
[`reliability_omega()`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_omega_categorical()`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
set.seed(113)
# Ten dichotomous items with a single underlying ability.
N <- 300
J <- 10
ability <- rnorm(N)
loadings <- rep(0.6, J)
latent <- outer(ability, loadings) +
          matrix(rnorm(N * J, sd = sqrt(1 - 0.6^2)), N, J)
items <- (latent > 0) * 1
colnames(items) <- paste0("y", seq_len(J))

reliability_kr20(data = items)
#>  term        value
#>  estimate    0.774
#>  se          <NA> 
#>  lower_limit 0.734
#>  upper_limit 0.81 
#>  conf_level  0.95 
#>  N           300  
#>  N_complete  300  
#>  J           10   
reliability_kr20(data = items, ci_method = "bonett")
#>  term        value 
#>  estimate    0.774 
#>  se          0.0864
#>  lower_limit 0.732 
#>  upper_limit 0.809 
#>  conf_level  0.95  
#>  N           300   
#>  N_complete  300   
#>  J           10    
reliability_kr20(data = items, ci_method = "percentile", B = 200)
#>  term        value 
#>  estimate    0.774 
#>  se          0.0193
#>  lower_limit 0.728 
#>  upper_limit 0.805 
#>  conf_level  0.95  
#>  N           300   
#>  N_complete  300   
#>  J           10    
```
