# Check Whether a Set of Contrasts Is Mutually Orthogonal

Tests whether every pair of columns in a contrast-coefficient matrix is
orthogonal under either the equal-\\n\\ convention \\\sum_i c\_{ik}
c\_{ij} = 0\\ or the unequal-\\n\\ convention \\\sum_i c\_{ik} c\_{ij} /
n_i = 0\\ (Maxwell, Delaney, & Kelley, 2027, Sec. 4.10; Kirk, 2013).
Also checks that each column sums to zero (the contrast property).

## Usage

``` r
is_orthogonal_set(contrasts, n = NULL, tol = 1e-08)
```

## Arguments

- contrasts:

  A numeric \\a \times m\\ matrix or `data.frame`, where \\a\\ is the
  number of groups and \\m\\ is the number of contrasts.

- n:

  Optional integer vector of length \\a\\ giving the per- group sample
  sizes. If supplied, the unequal-\\n\\ convention is used; otherwise
  the equal-\\n\\ convention is assumed.

- tol:

  Numerical tolerance for declaring orthogonality. Default `1e-8`.

## Value

A `data.frame` with rows for the overall orthogonality flag (`1` = all
pairs orthogonal, `0` = not), the contrast-sum-to-zero flag, the number
of contrasts tested, and one row per pairwise dot-product, named by
contrast pair.

## Details

**Equal-\\n\\.** Two contrasts \\\mathbf c, \mathbf d\\ on \\a\\ groups
of equal size are orthogonal iff \\\sum\_{i=1}^{a} c_i d_i = 0\\.

**Unequal-\\n\\.** With sample sizes \\n_1, \ldots, n_a\\, the
orthogonality condition that yields uncorrelated sample contrasts is
\\\sum\_{i=1}^{a} c_i d_i / n_i = 0\\.

**Useful for design checks.** Before performing planned comparisons or
partitioning the omnibus sums of squares, the user typically wants
confirmation that the chosen contrast set is orthogonal so that its
component SS sum to the omnibus SS.

## References

Kirk, R. E. (2013). *Experimental design: Procedures for the behavioral
sciences* (4th ed.). Sage.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Sec. 4.10.)

## See also

[`effects_coding`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`helmert_coding`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`scheffe_ci`](https://yelleknek.github.io/DMAR/reference/scheffe_ci.md)

Other design utilities:
[`deft()`](https://yelleknek.github.io/DMAR/reference/deft.md),
[`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md),
[`effects_coding()`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`helmert_coding()`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`orthogonal_polynomial()`](https://yelleknek.github.io/DMAR/reference/orthogonal_polynomial.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Two orthogonal contrasts on a 4-group design (equal n):
cmat <- cbind(
  c_linear  = c(-3, -1,  1,  3),
  c_quad    = c( 1, -1, -1,  1)
)
is_orthogonal_set(cmat)
#>  term                      value
#>  all_orthogonal            1    
#>  all_contrasts_sum_to_zero 1    
#>  n_contrasts               2    
#>  dot[c_linear . c_quad]    0    

# 2. Same contrasts under unequal sample sizes:
is_orthogonal_set(cmat, n = c(10, 8, 12, 9))
#>  term                      value
#>  all_orthogonal            0    
#>  all_contrasts_sum_to_zero 1    
#>  n_contrasts               2    
#>  dot[c_linear . c_quad]    0.075

# 3. Non-orthogonal pair:
cmat_bad <- cbind(
  c_diff_1  = c( 1, -1,  0,  0),
  c_diff_2  = c( 1,  0, -1,  0)
)
is_orthogonal_set(cmat_bad)
#>  term                      value
#>  all_orthogonal            0    
#>  all_contrasts_sum_to_zero 1    
#>  n_contrasts               2    
#>  dot[c_diff_1 . c_diff_2]  1    
```
