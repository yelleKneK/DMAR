# Effects-Coding Contrast Matrix for a Factor

Builds the effects-coding (also called *deviation coding* or
*sum-to-zero*) contrast matrix for a factor with \\a\\ levels. Each
non-reference level contrasts with the grand mean (rather than with a
reference category as in dummy coding). The returned matrix has rows =
levels and columns named after the levels, replacing the numeric column
names produced by
[`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html).

## Usage

``` r
effects_coding(levels, reference = NULL)
```

## Arguments

- levels:

  Either an integer giving the number of levels or a character / factor
  vector giving the level labels. If integer, the labels default to
  `"L1"`, `"L2"`, ...

- reference:

  Optional character name of the reference level (whose coefficients are
  all \\-1\\). Default is the last level.

## Value

A numeric \\a \times (a - 1)\\ matrix with row names = the factor levels
and column names = the non-reference levels. Suitable for assignment to
`contrasts(factor)` or use in manual contrast construction.

## Details

**Why effects coding.** Effects coding gives the regression intercept
the interpretation of the grand mean (rather than the reference-category
mean), and each slope coefficient becomes the deviation of that level's
mean from the grand mean (rather than the difference vs the reference
category). For balanced designs the effect coefficients are orthogonal
to the intercept.

**Equivalent to.**
[`stats::contr.sum()`](https://rdrr.io/r/stats/contrast.html) but with
meaningful column names (the level labels), which is what is lost in the
base-R implementation.

## References

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). *Applied
multiple regression/correlation analysis for the behavioral sciences*
(3rd ed.). Lawrence Erlbaum. (See Chapter 8.)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 4, 7.)

## See also

[`contr.sum`](https://rdrr.io/r/stats/contrast.html),
[`helmert_coding`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`is_orthogonal_set`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md)

Other design utilities:
[`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md),
[`design_effect()`](https://yelleknek.github.io/DMAR/reference/design_effect.md),
[`helmert_coding()`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`is_orthogonal_set()`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md),
[`orthogonal_polynomial()`](https://yelleknek.github.io/DMAR/reference/orthogonal_polynomial.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Effects coding for a 4-level factor:
effects_coding(c("low", "med", "high", "very_high"))
#>           low med high
#> low         1   0    0
#> med         0   1    0
#> high        0   0    1
#> very_high  -1  -1   -1

# 2. With "low" as the reference category:
effects_coding(c("low", "med", "high", "very_high"),
               reference = "low")
#>           med high very_high
#> low        -1   -1        -1
#> med         1    0         0
#> high        0    1         0
#> very_high   0    0         1

# 3. Assign to a factor for modeling:
f <- factor(c("a", "b", "c", "a", "b", "c"))
contrasts(f) <- effects_coding(levels(f))
model.matrix(~ f)
#>   (Intercept) fa fb
#> 1           1  1  0
#> 2           1  0  1
#> 3           1 -1 -1
#> 4           1  1  0
#> 5           1  0  1
#> 6           1 -1 -1
#> attr(,"assign")
#> [1] 0 1 1
#> attr(,"contrasts")
#> attr(,"contrasts")$f
#>    a  b
#> a  1  0
#> b  0  1
#> c -1 -1
#> 
```
