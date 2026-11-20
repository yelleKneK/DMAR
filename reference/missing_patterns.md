# Missing Data Patterns

Tabulates the distinct patterns of observed and missing values across a
set of items: which items each pattern leaves missing, and how many
cases show it. The table is the first description of incomplete data,
the one a reader needs before any model is fit. It shows whether the
missing values are scattered across many patterns or concentrated in a
few, whether they arrive together (two tests given to only one school
share one pattern), and how many cases are complete, which is all
listwise deletion would keep. Full information maximum likelihood (FIML)
uses every pattern, and the fraction of missing information
([`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md))
then says what the patterns cost each parameter.

## Usage

``` r
missing_patterns(data, items = NULL)
```

## Arguments

- data:

  A `data.frame`.

- items:

  Optional character vector naming the columns of `data` to describe;
  `NULL` (the default) uses every column. Any column type is allowed,
  since only whether a value is `NA` matters.

## Value

A `data.frame` (class `dmar_tbl`) with one row per distinct pattern and
the same columns whatever the number of items: `pattern` (the string of
1s and 0s), `missing_items` (the names of the items the pattern leaves
missing, separated by commas; `""` for the complete pattern), `n` (the
number of cases), `proportion` (of all rows of `data`), `n_observed`,
and `n_missing` (the numbers of items observed and missing in the
pattern). The attributes are `matrix` (the logical pattern matrix,
`TRUE` for observed, one row per pattern in the order of the table, the
layout of `lavaan::lavInspect(fit, "patterns")`), `N` (the rows of
`data`), `N_complete`, `N_all_missing`, `n_patterns`, and `items`.

## Details

A pattern is written as a string of 1s and 0s, one character per item in
the order of `items`, with 1 for an observed value and 0 for a missing
one (the convention of
[`mice::md.pattern()`](https://amices.org/mice/reference/md.pattern.html)).
The rows are sorted by the number of cases, largest first, then by the
number of items missing, fewest first, so the complete cases usually
lead and the rare patterns trail.

The pattern in which every item is missing is kept as a row
(`n_observed = 0`) when it occurs. lavaan sets such cases aside before
fitting, without reporting them among its own patterns, so the number of
cases a FIML fit uses is `N` minus `N_all_missing`.

## See also

[`missing_coverage`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md)
for the proportion of cases observed on each pair of items;
[`fraction_missing_information`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md)
for what the missing values cost each parameter of a FIML fit.

Other missing data:
[`cov_fiml()`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md),
[`df_barnard_rubin()`](https://yelleknek.github.io/DMAR/reference/df_barnard_rubin.md),
[`df_missing()`](https://yelleknek.github.io/DMAR/reference/df_missing.md),
[`fraction_missing_information()`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
[`missing_coverage()`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The six spatial tests of the Holzinger and Swineford data: two
# patterns, the 145 Grant-White pupils with every test and the 156
# Pasteur pupils, who were not given the revised paper form board
# (t25) or the flags test (t26).
spatial <- c("t1_visual_perception", "t2_cubes", "t3_paper_form_board",
             "t4_lozenges", "t25_paper_form_board_r", "t26_flags")
missing_patterns(holzinger_swineford, items = spatial)
#>  pattern missing_items                     n   proportion n_observed n_missing
#>  111100  t25_paper_form_board_r, t26_flags 156 0.518      4          2        
#>  111111                                    145 0.482      6          0        

# Values missing completely at random on two more tests add patterns.
set.seed(113)
d <- holzinger_swineford[spatial]
d$t1_visual_perception[runif(nrow(d)) < .10] <- NA
d$t2_cubes[runif(nrow(d)) < .10] <- NA
pat <- missing_patterns(d)
pat
#>  pattern missing_items                                           n   proportion
#>  111100  t25_paper_form_board_r, t26_flags                       128 0.425     
#>  111111                                                          115 0.382     
#>  011111  t1_visual_perception                                    16  0.053     
#>  101100  t2_cubes, t25_paper_form_board_r, t26_flags             14  0.047     
#>  011100  t1_visual_perception, t25_paper_form_board_r, t26_flags 14  0.047     
#>  101111  t2_cubes                                                13  0.043     
#>  001111  t1_visual_perception, t2_cubes                          1   0.003     
#>  n_observed n_missing
#>  4          2        
#>  6          0        
#>  5          1        
#>  3          3        
#>  3          3        
#>  5          1        
#>  4          2        
attr(pat, "N_complete")
#> [1] 115
```
