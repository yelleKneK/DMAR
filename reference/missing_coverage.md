# Covariance Coverage of Incomplete Data

Reports, for every pair of items, the number and proportion of cases on
which both were observed (the covariance coverage), with each item's own
response rate on the diagonal. A covariance or a regression slope is
informed most directly by the cases that observed both of its variables,
so the coverage of a pair is a first, model-free look at how much data
stand behind the parameters that involve it. Full information maximum
likelihood (FIML) uses every case, and the other variables can carry
information about a pair with low coverage, but a pair with no jointly
observed case leaves its covariance with no direct information at all:
its fraction of missing information is then one (Savalei & Rhemtulla,
2012, p. 488).
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
gives that fraction for each parameter of a fitted model.

## Usage

``` r
missing_coverage(data, items = NULL)
```

## Arguments

- data:

  A `data.frame`.

- items:

  Optional character vector naming the columns of `data` to describe;
  `NULL` (the default) uses every column. Any column type is allowed,
  since only whether a value is `NA` matters.

## Value

A `data.frame` (class `dmar_tbl`) with one row per pair of items, the
item itself included, in the order of `items` (the pairs \\(i, j)\\ with
\\i \le j\\), and the columns `var_1`, `var_2`, `n_both` (the cases that
observed both; for `var_1 == var_2`, the cases that observed the item),
and `coverage` (`n_both` as a proportion of the cases with at least one
observed item). The attributes are `matrix` (the symmetric coverage
matrix), `N_used` (the denominator), `N_all_missing`, and
`min_coverage`.

## Details

The denominator of every proportion is the number of cases with at least
one observed item, the cases a FIML fit uses. That is the denominator of
`lavaan::lavInspect(fit, "coverage")`, which the `matrix` attribute
reproduces for a single group fit to the same items; cases with every
item missing are counted in the `N_all_missing` attribute instead.

## References

Savalei, V., & Rhemtulla, M. (2012). On obtaining estimates of the
fraction of missing information from full information maximum
likelihood. *Structural Equation Modeling, 19*(3), 477–494.
[doi:10.1080/10705511.2012.687669](https://doi.org/10.1080/10705511.2012.687669)

## See also

[`missing_patterns`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)
for the patterns behind the coverage;
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
for what the missing values cost each parameter of a FIML fit.

Other missing data:
[`cov_fiml()`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md),
[`df_barnard_rubin()`](https://yelleknek.github.io/DMAR/reference/df_barnard_rubin.md),
[`df_missing()`](https://yelleknek.github.io/DMAR/reference/df_missing.md),
[`fraction_missing_information()`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
[`missing_patterns()`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The revised paper form board (t25) and the flags test (t26) were
# given only at Grant-White, so every pair involving either is covered
# by the 145 Grant-White pupils, 48 percent of the 301.
spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
cov_tbl <- missing_coverage(holzinger_swineford, items = spatial)
cov_tbl
#>  var_1                  var_2                  n_both coverage
#>  t1_visual_perception   t1_visual_perception   301    1.000   
#>  t1_visual_perception   t2_cubes               301    1.000   
#>  t1_visual_perception   t3_paper_form_board    301    1.000   
#>  t1_visual_perception   t4_lozenges            301    1.000   
#>  t1_visual_perception   t25_paper_form_board_r 145    0.482   
#>  t1_visual_perception   t26_flags              145    0.482   
#>  t2_cubes               t2_cubes               301    1.000   
#>  t2_cubes               t3_paper_form_board    301    1.000   
#>  t2_cubes               t4_lozenges            301    1.000   
#>  t2_cubes               t25_paper_form_board_r 145    0.482   
#>  t2_cubes               t26_flags              145    0.482   
#>  t3_paper_form_board    t3_paper_form_board    301    1.000   
#>  t3_paper_form_board    t4_lozenges            301    1.000   
#>  t3_paper_form_board    t25_paper_form_board_r 145    0.482   
#>  t3_paper_form_board    t26_flags              145    0.482   
#>  t4_lozenges            t4_lozenges            301    1.000   
#>  t4_lozenges            t25_paper_form_board_r 145    0.482   
#>  t4_lozenges            t26_flags              145    0.482   
#>  t25_paper_form_board_r t25_paper_form_board_r 145    0.482   
#>  t25_paper_form_board_r t26_flags              145    0.482   
#>  t26_flags              t26_flags              145    0.482   
round(attr(cov_tbl, "matrix"), 2)
#>                        t1_visual_perception t2_cubes t3_paper_form_board
#> t1_visual_perception                   1.00     1.00                1.00
#> t2_cubes                               1.00     1.00                1.00
#> t3_paper_form_board                    1.00     1.00                1.00
#> t4_lozenges                            1.00     1.00                1.00
#> t25_paper_form_board_r                 0.48     0.48                0.48
#> t26_flags                              0.48     0.48                0.48
#>                        t4_lozenges t25_paper_form_board_r t26_flags
#> t1_visual_perception          1.00                   0.48      0.48
#> t2_cubes                      1.00                   0.48      0.48
#> t3_paper_form_board           1.00                   0.48      0.48
#> t4_lozenges                   1.00                   0.48      0.48
#> t25_paper_form_board_r        0.48                   0.48      0.48
#> t26_flags                     0.48                   0.48      0.48
```
