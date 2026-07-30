# Tidy a Multiple-Factor CFA Fit

Returns the parameter rows of a
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) table
(loadings, error variances, intercepts, latent correlations, and the
defined measurement properties) in the column convention used by the
broom ecosystem. Fit-information rows belong in
[`glance.dmar_cfa_k`](https://yelleknek.github.io/DMAR/reference/glance.dmar_cfa_k.md).

## Usage

``` r
# S3 method for class 'dmar_cfa_k'
tidy(x, ...)
```

## Arguments

- x:

  A `dmar_cfa_k` object returned by
  [`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md).

- ...:

  Unused.

## Value

A `data.frame` with columns `term`, `estimate`, `se`, `statistic`,
`p_value`, `ci_lower`, `ci_upper`.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
data(holzinger_swineford)
res <- cfa_k(holzinger_swineford,
             list(verbal = c("t6_paragraph_comprehension",
                             "t7_sentence", "t9_word_meaning"),
                  deduction = c("t20_deduction",
                                "t22_problem_reasoning",
                                "t23_series_completion")))
generics::tidy(res)
#>                     term    estimate          se statistic      p_value
#> 1        lambda_verbal_1   2.9466143  0.16934224 17.400350 0.000000e+00
#> 2        lambda_verbal_2   4.3888203  0.24939123 17.598134 0.000000e+00
#> 3        lambda_verbal_3   6.4894296  0.37137393 17.474112 0.000000e+00
#> 4             phi_verbal   1.0000000  0.00000000        NA           NA
#> 5           psi_verbal_1   3.4734500  0.41765103  8.316632 0.000000e+00
#> 6           psi_verbal_2   7.2948357  0.90173426  8.089784 6.661338e-16
#> 7           psi_verbal_3  16.5088761  2.00512280  8.233349 2.220446e-16
#> 8     lambda_deduction_1  11.1491106  1.14813971  9.710587 0.000000e+00
#> 9     lambda_deduction_2   6.7383698  0.52382197 12.863855 0.000000e+00
#> 10    lambda_deduction_3   6.6808810  0.52073318 12.829759 0.000000e+00
#> 11         phi_deduction   1.0000000  0.00000000        NA           NA
#> 12       psi_deduction_1 248.2536995 23.49872202 10.564562 0.000000e+00
#> 13       psi_deduction_2  38.9026007  4.76641082  8.161823 2.220446e-16
#> 14       psi_deduction_3  38.6263555  4.70933491  8.202083 2.220446e-16
#> 15  phi_verbal_deduction   0.7298985  0.04282311 17.044501 0.000000e+00
#> 16    loading_sum_verbal  13.8248642  0.64703777 21.366394 0.000000e+00
#> 17      error_sum_verbal  27.2771619  2.00065797 13.634096 0.000000e+00
#> 18          omega_verbal   0.8751069  0.01361438 64.278149 0.000000e+00
#> 19            ave_verbal   0.7193173  0.02277669 31.581292 0.000000e+00
#> 20              H_verbal   0.8849393  0.01149083 77.012649 0.000000e+00
#> 21 loading_sum_deduction  24.5683615  1.54417976 15.910299 0.000000e+00
#> 22   error_sum_deduction 325.7826557 23.76883069 13.706297 0.000000e+00
#> 23       omega_deduction   0.6494650  0.03543796 18.326816 0.000000e+00
#> 24         ave_deduction   0.4694311  0.03323371 14.125148 0.000000e+00
#> 25           H_deduction   0.7384531  0.02725357 27.095650 0.000000e+00
#>       ci_lower    ci_upper
#> 1    2.6147096   3.2785190
#> 2    3.9000225   4.8776181
#> 3    5.7615500   7.2173091
#> 4    1.0000000   1.0000000
#> 5    2.6548691   4.2920310
#> 6    5.5274690   9.0622024
#> 7   12.5789076  20.4388446
#> 8    8.8987982  13.3994231
#> 9    5.7116976   7.7650420
#> 10   5.6602628   7.7014993
#> 11   1.0000000   1.0000000
#> 12 202.1970507 294.3103484
#> 13  29.5606072  48.2445942
#> 14  29.3962287  47.8564823
#> 15   0.6459668   0.8138303
#> 16  12.5566934  15.0930349
#> 17  23.3559443  31.1983794
#> 18   0.8484232   0.9017906
#> 19   0.6746758   0.7639588
#> 20   0.8624176   0.9074609
#> 21  21.5418248  27.5948982
#> 22 279.1966036 372.3687078
#> 23   0.5800079   0.7189222
#> 24   0.4042942   0.5345680
#> 25   0.6850371   0.7918692
generics::glance(res)
#>   chi_square df    p_value       cfi       tli      rmsea rmsea_low rmsea_high
#> 1   13.76907  8 0.08798665 0.9928478 0.9865897 0.04894684         0 0.09151619
#>         srmr     AIC      BIC    logLik
#> 1 0.02341405 11755.5 11803.69 -5864.748
```
