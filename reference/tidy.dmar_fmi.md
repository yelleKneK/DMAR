# Tidy and Glance Methods for the Fraction of Missing Information

`tidy()` returns one row per parameter of a
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
table, with the FMI as the estimate, its percentile limits, and the
effective sample size; `glance()` returns a one-row summary of the
table.

## Usage

``` r
# S3 method for class 'dmar_fmi'
tidy(x, ...)

# S3 method for class 'dmar_fmi'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_fmi` object returned by
  [`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md).

- ...:

  Unused.

## Value

`tidy()` returns a `data.frame` with columns `term` (the parameter),
`group`, `estimate` (its FMI), `estimand` (`"fmi"`, which names the
estimate, so that
[`results_sentence`](https://yelleknek.github.io/DMAR/reference/dmar_output_helpers.md)
on a row reports `"fmi = 0.24, 90% CI [0.18, 0.32]"` and not the FMI as
the parameter's value), `ci_lower`, `ci_upper`, `conf_level` (all three
`NA` without the percentile interval), and `N_effective`. `glance()`
returns a one-row `data.frame` with `n_parameters`, `fmi_max`,
`fmi_median`, `N` (the cases the fit used, summed over groups), and
`information`.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
fit <- fa_confirmatory(data = holzinger_swineford,
             factors = list(spatial = c("t1_visual_perception",
               "t2_cubes", "t3_paper_form_board", "t4_lozenges",
               "t25_paper_form_board_r", "t26_flags")),
             missing = "fiml", output = "fit")
fmi <- fraction_missing_information(fit)
generics::tidy(fmi)
#>                                                term group      estimate
#> 1                   spatial =~ t1_visual_perception  <NA>  1.125059e-01
#> 2                               spatial =~ t2_cubes  <NA>  6.309454e-02
#> 3                    spatial =~ t3_paper_form_board  <NA>  7.458424e-02
#> 4                            spatial =~ t4_lozenges  <NA>  1.112695e-01
#> 5                 spatial =~ t25_paper_form_board_r  <NA>  5.116258e-01
#> 6                              spatial =~ t26_flags  <NA>  5.061975e-01
#> 7      t1_visual_perception ~~ t1_visual_perception  <NA>  2.140064e-01
#> 8                              t2_cubes ~~ t2_cubes  <NA>  5.769930e-02
#> 9        t3_paper_form_board ~~ t3_paper_form_board  <NA>  1.197454e-01
#> 10                       t4_lozenges ~~ t4_lozenges  <NA>  2.029862e-01
#> 11 t25_paper_form_board_r ~~ t25_paper_form_board_r  <NA>  5.165528e-01
#> 12                           t26_flags ~~ t26_flags  <NA>  5.109202e-01
#> 13                          t1_visual_perception ~1  <NA>  0.000000e+00
#> 14                                      t2_cubes ~1  <NA> -4.440892e-16
#> 15                           t3_paper_form_board ~1  <NA>  4.440892e-16
#> 16                                   t4_lozenges ~1  <NA> -8.881784e-16
#> 17                        t25_paper_form_board_r ~1  <NA>  4.436849e-01
#> 18                                     t26_flags ~1  <NA>  4.304065e-01
#> 19                              loading_sum_spatial  <NA>  1.837593e-01
#> 20                                error_sum_spatial  <NA>  3.479011e-01
#> 21                                    omega_spatial  <NA>  2.442136e-01
#> 22                                      ave_spatial  <NA>  2.213080e-01
#> 23                                        H_spatial  <NA>  2.163195e-01
#>    estimand ci_lower ci_upper conf_level N_effective
#> 1       fmi       NA       NA         NA    267.1357
#> 2       fmi       NA       NA         NA    282.0085
#> 3       fmi       NA       NA         NA    278.5501
#> 4       fmi       NA       NA         NA    267.5079
#> 5       fmi       NA       NA         NA    147.0006
#> 6       fmi       NA       NA         NA    148.6345
#> 7       fmi       NA       NA         NA    236.5841
#> 8       fmi       NA       NA         NA    283.6325
#> 9       fmi       NA       NA         NA    264.9566
#> 10      fmi       NA       NA         NA    239.9012
#> 11      fmi       NA       NA         NA    145.5176
#> 12      fmi       NA       NA         NA    147.2130
#> 13      fmi       NA       NA         NA    301.0000
#> 14      fmi       NA       NA         NA    301.0000
#> 15      fmi       NA       NA         NA    301.0000
#> 16      fmi       NA       NA         NA    301.0000
#> 17      fmi       NA       NA         NA    167.4509
#> 18      fmi       NA       NA         NA    171.4476
#> 19      fmi       NA       NA         NA    245.6884
#> 20      fmi       NA       NA         NA    196.2818
#> 21      fmi       NA       NA         NA    227.4917
#> 22      fmi       NA       NA         NA    234.3863
#> 23      fmi       NA       NA         NA    235.8878
generics::glance(fmi)
#>   n_parameters   fmi_max fmi_median   N  information
#> 1           23 0.5165528  0.2029862 301 unstructured
```
