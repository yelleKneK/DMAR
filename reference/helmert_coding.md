# Helmert-Coding Contrast Matrix for a Factor

Builds the Helmert-coding contrast matrix for a factor with \\a\\
levels. The \\k\\-th column contrasts the \\(k + 1)\\-th level against
the average of all preceding levels, giving a fully orthogonal set under
equal sample sizes. The returned matrix has columns named after the
contrasted level rather than the numeric column names produced by
[`stats::contr.helmert()`](https://rdrr.io/r/stats/contrast.html).

## Usage

``` r
helmert_coding(levels)
```

## Arguments

- levels:

  Either an integer giving the number of levels or a character / factor
  vector giving the level labels. If integer, the labels default to
  `"L1"`, `"L2"`, ...

## Value

A numeric \\a \times (a - 1)\\ matrix with row names = the factor levels
and column names of the form `"L2_vs_prior"`, `"L3_vs_prior"`, ...

## Details

**Why Helmert.** Helmert contrasts are the canonical "sequential"
orthogonal contrast set: under equal-\\n\\, every column is orthogonal
to every other column and to the intercept. They are useful when the
factor has a natural ordering and the research questions are "does the
\\k\\-th level differ from the average of the preceding levels?"

**Equivalent to.**
[`stats::contr.helmert()`](https://rdrr.io/r/stats/contrast.html) but
with interpretable column names.

## References

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum.

## See also

[`contr.helmert`](https://rdrr.io/r/stats/contrast.html),
[`effects_coding`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`is_orthogonal_set`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md)

Other design utilities:
[`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md),
[`design_effect()`](https://yelleknek.github.io/DMAR/reference/design_effect.md),
[`effects_coding()`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`is_orthogonal_set()`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md),
[`orthogonal_polynomial()`](https://yelleknek.github.io/DMAR/reference/orthogonal_polynomial.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Helmert coding for a 4-level factor:
helmert_coding(c("baseline", "week1", "week2", "week3"))
#>          week1_vs_prior week2_vs_prior week3_vs_prior
#> baseline             -1             -1             -1
#> week1                 1             -1             -1
#> week2                 0              2             -1
#> week3                 0              0              3

# 2. Confirm orthogonality:
M <- helmert_coding(4)
is_orthogonal_set(M)
#>  term                           value
#>  all_orthogonal                 1    
#>  all_contrasts_sum_to_zero      1    
#>  n_contrasts                    3    
#>  dot[L2_vs_prior . L3_vs_prior] 0    
#>  dot[L2_vs_prior . L4_vs_prior] 0    
#>  dot[L3_vs_prior . L4_vs_prior] 0    
```
