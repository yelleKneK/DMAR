# Orthogonal-Polynomial (Trend) Contrast Coefficients

Builds the set of orthogonal-polynomial contrasts (linear, quadratic,
cubic, ...) for a factor whose \\a\\ levels are quantitative, so that
between-group variation can be decomposed into independent trend
components. The columns are mutually orthogonal and each sums to zero
(so each is a contrast on the group means). By default the coefficients
are returned in *orthonormal* form, meaning each column also has unit
length, \\\sum_i c_i^2 = 1\\; the alternative `type = "integer"`
rescales every column to the small whole numbers used in the published
orthogonal-polynomial tables, which are easier to read by eye and match
hand computation. The per-trend sum of squared coefficients \\\sum_i
c_i^2\\ is carried on the returned object and shown when it is printed.

## Usage

``` r
orthogonal_polynomial(
  levels,
  scores = NULL,
  type = c("orthonormal", "integer"),
  degree = NULL
)
```

## Arguments

- levels:

  One of: a single integer giving the number of levels \\a\\ (labels
  default to `"L1"`, `"L2"`, ..., and the levels are treated as equally
  spaced); a character or factor vector of level labels (equally spaced
  unless `scores` is given); or a numeric vector of the quantitative
  level values themselves, which are used both as the labels and as the
  spacing.

- scores:

  Optional numeric vector of length \\a\\ giving the quantitative
  position of each level. Supply this when the labels are non-numeric
  but the spacing is unequal (e.g., doses 0, 1, 2, and 4). Defaults to
  equally spaced positions `1:a`.

- type:

  Either `"orthonormal"` (default) for unit-length columns (\\\sum_i
  c_i^2 = 1\\, identical to
  [`contr.poly`](https://rdrr.io/r/stats/contrast.html)) or `"integer"`
  for the minimal whole-number coefficients of the classic
  trend-coefficient table. `"integer"` requires equally spaced `scores`.

- degree:

  Highest-order trend to return, an integer between `1` and \\a - 1\\.
  Defaults to \\a - 1\\ (the full set the design can support).

## Value

A numeric \\a \times \text{degree}\\ matrix with row names = the level
labels and column names = the trend names (`"linear"`, `"quadratic"`,
`"cubic"`, `"quartic"`, ...). The matrix can be assigned directly to
`contrasts(factor)` or passed to
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`is_orthogonal_set`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md),
or [`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md). The
per-column sum of squared coefficients is stored as `attr(*, "sum_sq")`
(a named numeric vector), and the level spacing and `type` are stored as
`attr(*, "scores")` and `attr(*, "type")`. The object carries class
`"orthogonal_polynomial"` so that printing shows the coefficients
alongside \\\sum_i c_i^2\\. When printed, the object is shown in the
textbook Table A.10 orientation (trends in rows, levels in columns) with
the \\\sum_i c_i^2\\ values as a final column; the stored matrix is the
transpose of that display (levels in rows) so it can be assigned
directly to [`contrasts()`](https://rdrr.io/r/stats/contrasts.html).

## Details

**What a trend contrast is.** When the levels of a factor are
quantitative and ordered (minutes of study, dose, day), the omnibus
between-group variation can be partitioned into a linear trend (does the
mean rise or fall steadily?), a quadratic trend (is there curvature?), a
cubic trend (an S-shape?), and so on, up to order \\a - 1\\. Each trend
is a single 1-*df* contrast on the group means, \\\hat\psi = \sum_i c_i
\bar Y_i\\, and because the contrasts are mutually orthogonal their sums
of squares add up to the omnibus between-group sum of squares exactly.

**Orthonormal versus integer scaling.** The two `type` values return the
*same* trends (same directions, same hypotheses, identical *F*, *t*, and
*p* for a given trend); they differ only in how each column is scaled.

- `"orthonormal"` normalizes each column to unit length, \\\sum_i c_i^2
  = 1\\. The trend portion of the design is then an orthonormal basis,
  the sum of squares for a trend is simply \\\hat\psi^2\\, and in a
  fitted [`lm`](https://rdrr.io/r/stats/lm.html) the `.linear` /
  `.quadratic` coefficients are the trend estimates on a common scale.
  This is the form [`lm`](https://rdrr.io/r/stats/lm.html) and
  [`aov`](https://rdrr.io/r/stats/aov.html) use internally and is
  defined for any spacing, equal or unequal.

- `"integer"` rescales each column to the smallest whole numbers with
  the same ratios, reproducing the published orthogonal-polynomial
  coefficient table (e.g., for \\a = 4\\ the linear contrast is \\-3,
  -1, 1, 3\\). These are easy to read and to compute with by hand, but
  the columns no longer share a common length: \\\sum_i c_i^2\\ varies
  by trend (for \\a = 4\\: 20, 4, and 20), so the sum of squares for a
  trend is \\n\\ \hat\psi^2 / \sum_i c_i^2\\. Integer coefficients exist
  only for equally spaced levels.

The reported \\\sum_i c_i^2\\ (shown on printing and stored in
`attr(*, "sum_sq")`) is exactly the divisor in that sum-of-squares
formula; for the orthonormal form it is 1 for every trend.

**Unequal spacing.** With unequally spaced `scores` the orthonormal
polynomials are still uniquely defined and are returned by the
`"orthonormal"` type. Whole-number coefficients generally do not exist
in that case, so `type = "integer"` is an error.

**Relation to base R.** For equally spaced levels the orthonormal output
is identical to
[`contr.poly`](https://rdrr.io/r/stats/contrast.html)`(a)`;
`orthogonal_polynomial` adds the integer table form, an explicit
`scores` argument for unequal spacing, meaningful trend names, and the
\\\sum_i c_i^2\\ report.

## References

Fisher, R. A., & Yates, F. (1953). *Statistical tables for biological,
agricultural and medical research* (4th ed.). Oliver and Boyd. (Origin
of the tabulated integer coefficients.)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 6 on trend analysis; Appendix Table A.10
reports these coefficients.)

## See also

[`contr.poly`](https://rdrr.io/r/stats/contrast.html),
[`effects_coding`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`helmert_coding`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`is_orthogonal_set`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md),
[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`ci_c`](https://yelleknek.github.io/DMAR/reference/ci_c.md)

Other design utilities:
[`design_consequences()`](https://yelleknek.github.io/DMAR/reference/design_consequences.md),
[`design_effect()`](https://yelleknek.github.io/DMAR/reference/design_effect.md),
[`effects_coding()`](https://yelleknek.github.io/DMAR/reference/effects_coding.md),
[`helmert_coding()`](https://yelleknek.github.io/DMAR/reference/helmert_coding.md),
[`is_orthogonal_set()`](https://yelleknek.github.io/DMAR/reference/is_orthogonal_set.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Trend (orthogonal-polynomial) analysis decomposes the omnibus effect of a
# quantitative factor (dose, time, trial block, stimulus intensity) into
# independent linear, quadratic, cubic, ... components, as in the trend
# analysis of Maxwell, Delaney, and Kelley (2027, Chapter 6).

# 1. The a = 4 coefficient table in the integer form of the textbook
#    appendix (Table A.10). Printing puts the trends in rows and appends a
#    final sum-of-squared-coefficients (sum c^2) column: linear -3 -1 1 3,
#    and per-trend sum c^2 of 20, 4, 20.
orthogonal_polynomial(4, type = "integer")
#> Orthogonal-polynomial (trend) contrasts (integer), 4 levels
#> Layout: trends in rows, levels in columns (MDK Table A.10); final column is sum c^2.
#> Stored object is the transpose (levels x trends) for contrasts().
#> 
#>           L1 L2 L3 L4 sum c^2
#> linear    -3 -1  1  3      20
#> quadratic  1 -1 -1  1       4
#> cubic     -1  3 -3  1      20

# 2. The default orthonormal form (identical to stats::contr.poly(4)); every
#    trend has sum c^2 = 1, the scale lm() and aov() use internally.
orthogonal_polynomial(4)
#> Orthogonal-polynomial (trend) contrasts (orthonormal), 4 levels
#> Layout: trends in rows, levels in columns (MDK Table A.10); final column is sum c^2.
#> Stored object is the transpose (levels x trends) for contrasts().
#> 
#>               L1     L2     L3    L4 sum c^2
#> linear    -0.671 -0.224  0.224 0.671       1
#> quadratic  0.500 -0.500 -0.500 0.500       1
#> cubic     -0.224  0.671 -0.671 0.224       1

# 3. A worked trend analysis in the style of MDK Chapter 6. An outcome is
#    measured at a = 5 equally spaced levels of a quantitative factor (say
#    stimulus intensity 1..5), n = 8 per level, from a population with a
#    strong linear and a mild quadratic trend. Assigning the trend contrasts
#    to the factor makes lm() report each trend test directly.
set.seed(113)
intensity <- factor(rep(1:5, each = 8))
contrasts(intensity) <- orthogonal_polynomial(levels(intensity))
y <- c(10, 9, 7, 6, 6)[as.integer(intensity)] +
       rnorm(length(intensity), sd = 1.2)
round(coef(summary(lm(y ~ intensity))), 3)
#>                    Estimate Std. Error t value Pr(>|t|)
#> (Intercept)           7.784      0.197  39.532    0.000
#> intensitylinear      -3.044      0.440  -6.913    0.000
#> intensityquadratic    0.273      0.440   0.620    0.539
#> intensitycubic        0.802      0.440   1.821    0.077
#> intensityquartic      0.205      0.440   0.465    0.645

# 4. The same trends through DMAR's contrast_test(), which reports each
#    trend's estimate, standard error, t, p, and confidence interval from a
#    fitted one-way model. The integer columns are the contrast weights.
op <- orthogonal_polynomial(5, type = "integer")
contrast_test(aov(y ~ intensity),
              contrasts = list(linear    = op[, "linear"],
                               quadratic = op[, "quadratic"],
                               cubic     = op[, "cubic"]))
#>  contrast  estimate se   t     df p_value  p_adjusted ci_lower ci_upper
#>  linear    -9.63    1.39 -6.91 35 < 0.0001 < 0.0001   -12.5    -6.8    
#>  quadratic 1.02     1.65 0.62  35 0.5390   0.5390     -2.32    4.37    
#>  cubic     2.54     1.39 1.82  35 0.0771   0.0771     -0.291   5.36    
#> 
#> Confidence level: 95%

# 5. Unequally spaced doses (0, 1, 2, 4 mg): integer coefficients no longer
#    exist, but the orthonormal trends remain uniquely defined.
orthogonal_polynomial(c("0 mg", "1 mg", "2 mg", "4 mg"),
                      scores = c(0, 1, 2, 4))
#> Orthogonal-polynomial (trend) contrasts (orthonormal), 4 levels
#> Layout: trends in rows, levels in columns (MDK Table A.10); final column is sum c^2.
#> Stored object is the transpose (levels x trends) for contrasts().
#> 
#>             0 mg   1 mg   2 mg  4 mg sum c^2
#> linear    -0.592 -0.254  0.085 0.761       1
#> quadratic  0.564 -0.322 -0.645 0.403       1
#> cubic     -0.286  0.763 -0.572 0.095       1

# 6. Confirm a returned set is mutually orthogonal.
is_orthogonal_set(orthogonal_polynomial(5, type = "integer"))
#>  term                      value
#>  all_orthogonal            1    
#>  all_contrasts_sum_to_zero 1    
#>  n_contrasts               4    
#>  dot[linear . quadratic]   0    
#>  dot[linear . cubic]       0    
#>  dot[linear . quartic]     0    
#>  dot[quadratic . cubic]    0    
#>  dot[quadratic . quartic]  0    
#>  dot[cubic . quartic]      0    
```
