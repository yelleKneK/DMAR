# Tidy a Multiple-Factor CFA Fit

Returns the parameter rows of a
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
table (loadings, error variances, intercepts, latent correlations, and
the defined measurement properties) in the column convention used by the
broom ecosystem. Fit-information rows belong in
[`glance.dmar_fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/glance.dmar_fa_confirmatory.md).

## Usage

``` r
# S3 method for class 'dmar_fa_confirmatory'
tidy(x, ...)
```

## Arguments

- x:

  A `dmar_fa_confirmatory` object returned by
  [`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md).

- ...:

  Unused.

## Value

A `data.frame` with columns `term`, `estimate`, `se`, `statistic`,
`p_value`, `ci_lower`, `ci_upper`, `conf_level`, the level of the
parameter intervals, and `wald_df`, the degrees of freedom of the *t*
reference of the row's interval and *p*-value (`wald_df` in
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)),
`Inf` for the normal reference. Under `wald_df = "missing_information"`
only the loadings and the omega rows take a *t* reference, so the column
tells them from the other rows, as it tells the rows of one fit from
another's after tables are stacked.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
data(holzinger_swineford)
res <- fa_confirmatory(holzinger_swineford,
             list(verbal = c("t6_paragraph_comprehension",
                             "t7_sentence", "t9_word_meaning"),
                  deduction = c("t20_deduction",
                                "t22_problem_reasoning",
                                "t23_series_completion")))
generics::tidy(res)
#>                     term    estimate          se statistic      p_value
#> 1        lambda_verbal_1   2.9466143  0.16947344 17.386880 0.000000e+00
#> 2        lambda_verbal_2   4.3888203  0.24924860 17.608205 0.000000e+00
#> 3        lambda_verbal_3   6.4894296  0.37258272 17.417420 0.000000e+00
#> 4             phi_verbal   1.0000000  0.00000000        NA           NA
#> 5           psi_verbal_1   3.4734500  0.41949480  8.280079 2.220446e-16
#> 6           psi_verbal_2   7.2948357  0.89868934  8.117194 4.440892e-16
#> 7           psi_verbal_3  16.5088761  2.04254621  8.082498 6.661338e-16
#> 8            nu_verbal_1   9.1827243  0.20096103 45.694055 0.000000e+00
#> 9            nu_verbal_2  17.3621262  0.29703171 58.452096 0.000000e+00
#> 10           nu_verbal_3  15.2990033  0.44131174 34.667111 0.000000e+00
#> 11    lambda_deduction_1  11.1491106  1.15197228  9.678280 0.000000e+00
#> 12    lambda_deduction_2   6.7383698  0.52462060 12.844272 0.000000e+00
#> 13    lambda_deduction_3   6.6808810  0.52396112 12.750719 0.000000e+00
#> 14         phi_deduction   1.0000000  0.00000000        NA           NA
#> 15       psi_deduction_1 248.2536995 23.59181333 10.522875 0.000000e+00
#> 16       psi_deduction_2  38.9026007  4.78233925  8.134638 4.440892e-16
#> 17       psi_deduction_3  38.6263555  4.77282996  8.092967 6.661338e-16
#> 18        nu_deduction_1  26.8903654  1.11253260 24.170407 0.000000e+00
#> 19        nu_deduction_2  26.2392027  0.52923887 49.579130 0.000000e+00
#> 20        nu_deduction_3  18.1362126  0.52594015 34.483415 0.000000e+00
#> 21  phi_verbal_deduction   0.7298985  0.04322880 16.884543 0.000000e+00
#> 22    loading_sum_verbal  13.8248642  0.64674568 21.376044 0.000000e+00
#> 23      error_sum_verbal  27.2771619  1.99834709 13.649862 0.000000e+00
#> 24          omega_verbal   0.8751069  0.01358534 64.415537 0.000000e+00
#> 25            ave_verbal   0.7193173  0.02279408 31.557206 0.000000e+00
#> 26              H_verbal   0.8849393  0.01150594 76.911537 0.000000e+00
#> 27 loading_sum_deduction  24.5683615  1.54842271 15.866702 0.000000e+00
#> 28   error_sum_deduction 325.7826557 23.89015376 13.636691 0.000000e+00
#> 29       omega_deduction   0.6494650  0.03564292 18.221433 0.000000e+00
#> 30         ave_deduction   0.4694311  0.03330481 14.094994 0.000000e+00
#> 31           H_deduction   0.7384531  0.02722139 27.127684 0.000000e+00
#>       ci_lower    ci_upper conf_level wald_df h1_diagnostic
#> 1    2.6144525   3.2787761       0.95     Inf          <NA>
#> 2    3.9003020   4.8773386       0.95     Inf          <NA>
#> 3    5.7591809   7.2196783       0.95     Inf          <NA>
#> 4    1.0000000   1.0000000       0.95     Inf          <NA>
#> 5    2.6512553   4.2956447       0.95     Inf          <NA>
#> 6    5.5334370   9.0562344       0.95     Inf          <NA>
#> 7   12.5055591  20.5121931       0.95     Inf          <NA>
#> 8    8.7888479   9.5766006       0.95     Inf          <NA>
#> 9   16.7799548  17.9442977       0.95     Inf          <NA>
#> 10  14.4340482  16.1639584       0.95     Inf          <NA>
#> 11   8.8912865  13.4069348       0.95     Inf          <NA>
#> 12   5.7101323   7.7666073       0.95     Inf          <NA>
#> 13   5.6539361   7.7078259       0.95     Inf          <NA>
#> 14   1.0000000   1.0000000       0.95     Inf          <NA>
#> 15 202.0145951 294.4928040       0.95     Inf          <NA>
#> 16  29.5293880  48.2758134       0.95     Inf          <NA>
#> 17  29.2717807  47.9809303       0.95     Inf          <NA>
#> 18  24.7098416  29.0708893       0.95     Inf          <NA>
#> 19  25.2019135  27.2764918       0.95     Inf          <NA>
#> 20  17.1053889  19.1670364       0.95     Inf          <NA>
#> 21   0.6451716   0.8146254       0.95     Inf          <NA>
#> 22  12.5572659  15.0924624       0.95     Inf          <NA>
#> 23  23.3604735  31.1938502       0.95     Inf          <NA>
#> 24   0.8459606   0.8993942       0.95     Inf          <NA>
#> 25   0.6725634   0.7617601       0.95     Inf          <NA>
#> 26   0.8603958   0.9056410       0.95     Inf          <NA>
#> 27  21.5335087  27.6032142       0.95     Inf          <NA>
#> 28 278.9588148 372.6064967       0.95     Inf          <NA>
#> 29   0.5768445   0.7157637       0.95     Inf          <NA>
#> 30   0.4050389   0.5348574       0.95     Inf          <NA>
#> 31   0.6817237   0.7882124       0.95     Inf          <NA>
#>    wald_calibration_note
#> 1                   <NA>
#> 2                   <NA>
#> 3                   <NA>
#> 4                   <NA>
#> 5                   <NA>
#> 6                   <NA>
#> 7                   <NA>
#> 8                   <NA>
#> 9                   <NA>
#> 10                  <NA>
#> 11                  <NA>
#> 12                  <NA>
#> 13                  <NA>
#> 14                  <NA>
#> 15                  <NA>
#> 16                  <NA>
#> 17                  <NA>
#> 18                  <NA>
#> 19                  <NA>
#> 20                  <NA>
#> 21                  <NA>
#> 22                  <NA>
#> 23                  <NA>
#> 24                  <NA>
#> 25                  <NA>
#> 26                  <NA>
#> 27                  <NA>
#> 28                  <NA>
#> 29                  <NA>
#> 30                  <NA>
#> 31                  <NA>
generics::glance(res)
#>   chi_square df    p_value       cfi       tli      rmsea rmsea_low rmsea_high
#> 1   13.76907  8 0.08798665 0.9928478 0.9865897 0.04894684         0 0.09151619
#>         srmr     AIC      BIC    logLik estimator missing     test fit_indices
#> 1 0.02341405 11767.5 11837.93 -5864.748        ML    fiml standard      robust
#>   rmsea_ci_level h1_diagnostic wald_calibration_note
#> 1            0.9          <NA>                  <NA>
```
