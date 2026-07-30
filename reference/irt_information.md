# Item and Test Information for the Graded Response Model

Evaluates the item information functions and the test information
function of a graded response model on a grid of latent trait values,
together with the standard error of the latent trait estimate,
\\SE(\theta) = 1 / \sqrt{I(\theta)}\\. Reliability is a single number
that describes a scale at one place on the latent continuum; the
information function is the same idea expressed as a function of where
the respondent sits, so an item pool can be judged on where it measures
precisely rather than on one global summary. Because information is
additive across items, the curve also shows which items carry the
precision, and over what range, which is what makes it useful for
building and trimming a scale.

## Usage

``` r
irt_information(
  a,
  b = NULL,
  item = NULL,
  theta = seq(-4, 4, length.out = 81),
  grm = NULL
)
```

## Arguments

- a:

  Discriminations, a numeric vector of positive values. Supply either
  one value per item (named with the item names, or in the order the
  items first appear in `item`) or one value per boundary row (constant
  within an item). Not used when `grm` is supplied.

- b:

  Boundary locations (category thresholds), a numeric vector with one
  element per category boundary. An item with \\m\\ categories has \\m -
  1\\ boundaries, which the model requires to be in ascending order
  within the item. Not used when `grm` is supplied.

- item:

  Item labels, a character or factor vector the same length as `b`
  naming the item each boundary belongs to. When `NULL` (default), each
  element of `b` is treated as its own dichotomous item, named `item_1`,
  `item_2`, and so on, and `a` must then have the same length as `b`.
  Not used when `grm` is supplied.

- theta:

  Latent trait values at which to evaluate the information functions.
  Any finite numeric vector; the default, `seq(-4, 4, length.out = 81)`,
  covers the range in which almost all of a standard normal trait
  distribution falls, in steps of 0.1.

- grm:

  Optionally, the result of
  [`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md):
  a `data.frame` with one row per item and category boundary and columns
  `item`, `a`, and `b` (a `category` column, when present, orders the
  boundaries within an item). Supply this or the parameters, not both.

## Value

A `data.frame` (class `dmar_tbl`) with one row per value of `theta` and
columns:

- `theta`:

  The latent trait value, as supplied.

- `test_information`:

  Test information at that value, the sum of the item information
  functions.

- `se`:

  The standard error of the latent trait estimate, \\1 /
  \sqrt{I(\theta)}\\. It is `Inf` where test information is zero, which
  is the correct statement that the items carry no information there.

The result carries these attributes:

- `"item_information"`:

  A numeric matrix of item information with `theta` in the rows (row
  names are the `theta` values) and items in the columns (column names
  are the item names). Its row sums are `test_information`.

- `"item"`:

  The item names, in the order they appear in the columns of
  `"item_information"`.

- `"a"`:

  The discrimination used for each item, a numeric vector named by item.

- `"b"`:

  The boundary locations used, a numeric vector in item order and,
  within an item, in ascending order, named by the item each boundary
  belongs to.

- `"theta_max_information"`:

  The value of `theta` at which test information peaks on the supplied
  grid (the first such value if there are ties). It is a grid value, not
  the result of an optimization, so a finer `theta` locates the peak
  more sharply.

## Details

For item *i* with discrimination \\a_i\\ and ordered boundary locations
\\b\_{i1} \< b\_{i2} \< \cdots \< b\_{i,m-1}\\ for *m* categories, the
normal ogive graded response model of Samejima (1969) defines the
boundary response function \$\$P^\*\_{ik}(\theta) = \Phi\[a_i (\theta -
b\_{ik})\],\$\$ the probability of responding *above* boundary *k*, that
is, in any category higher than the *k*th, with the conventions
\\P^\*\_{i0} = 1\\ and \\P^\*\_{im} = 0\\. The category response
function is the difference of adjacent boundary functions,
\$\$P\_{ik}(\theta) = P^\*\_{i,k-1}(\theta) - P^\*\_{ik}(\theta),\$\$
and differentiating with respect to \\\theta\\ gives
\$\$P'\_{ik}(\theta) = a_i \\\phi\[a_i (\theta - b\_{i,k-1})\] -
\phi\[a_i (\theta - b\_{ik})\]\\,\$\$ where \\\phi\\ is the standard
normal density and the density terms vanish at the two extreme
categories (there is no \\b\_{i0}\\ and no \\b\_{im}\\). Item
information is \$\$I_i(\theta) = \sum\_{k=1}^{m}
\frac{\[P'\_{ik}(\theta)\]^2}{P\_{ik}(\theta)},\$\$ test information is
\\I(\theta) = \sum_i I_i(\theta)\\, and the standard error of the
maximum likelihood estimate of \\\theta\\ is \\SE(\theta) = 1 /
\sqrt{I(\theta)}\\.

Two properties make the curve worth reading. Information is additive
across items, so an item's contribution can be read off directly and a
pool can be assembled to cover a targeted range. And the reciprocal
relation to the squared standard error means the peak of the curve
locates where the scale estimates the trait most precisely, reported
here as the `"theta_max_information"` attribute.

For a dichotomous item the model reduces to the two parameter normal
ogive, whose information has the closed form \$\$I_i(\theta) =
\frac{a_i^2 \phi\[a_i(\theta - b_i)\]^2}{ \Phi\[a_i(\theta - b_i)\]
\\1 - \Phi\[a_i(\theta - b_i)\]\\},\$\$ which the general expression
above reproduces; that identity is one of the tests of this function.

The category probabilities underflow to zero for \\\theta\\ far from
every boundary, where the ratio \\(P')^2 / P\\ would be \\0/0\\. A
category whose probability is not strictly positive contributes zero to
the sum, which is the limit the ratio approaches, so the returned
information is finite and nonnegative on any grid, however extreme, and
is never `NaN`. In the regime where \\(P')^2\\ underflows but \\P\\ does
not, the ratio is formed as \\\exp\[2 \log \|P'\| - \log P\]\\ so the
contribution is kept rather than flushed to zero. Where two boundaries
of an item coincide, the category between them has probability zero
everywhere and, by the same guard, contributes nothing.

The parameters are in the normal ogive metric, which is what
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md)
returns by default. The logistic metric used by much of the item
response theory software scales the discrimination by approximately
1.702 (Camilli, 1994); a logistic \\a\\ is put on the normal ogive scale
by dividing by that constant. The two metrics give information functions
that are proportional in shape but not equal in value, so a
cross-software comparison is a comparison of curves, not of numbers.

## References

Baker, F. B., & Kim, S.-H. (2004). *Item response theory: Parameter
estimation techniques* (2nd ed.). Marcel Dekker.

Camilli, G. (1994). Teacher's corner: Origin of the scaling constant *d*
= 1.7 in item response theory. *Journal of Educational and Behavioral
Statistics, 19*(3), 293–295.
[doi:10.3102/10769986019003293](https://doi.org/10.3102/10769986019003293)

Embretson, S. E., & Reise, S. P. (2000). *Item response theory for
psychologists*. Lawrence Erlbaum.

Lord, F. M. (1980). *Applications of item response theory to practical
testing problems*. Lawrence Erlbaum.

Samejima, F. (1969). Estimation of latent ability using a response
pattern of graded scores. *Psychometrika Monograph Supplement, 34*(4,
Pt. 2), 1–97.

## See also

[`plot_irt_information`](https://yelleknek.github.io/DMAR/reference/plot_irt_information.md)
for the curve,
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for the single-number companion.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Three items: a five-category rating item and two dichotomous items.
# The discriminations are named, so they are matched to the item labels.
info <- irt_information(
  a = c(mood_1 = 1.4, mood_2 = 0.9, mood_3 = 1.1),
  b = c(-1.5, -0.5, 0.5, 1.5, 0.0, 0.8),
  item = c(rep("mood_1", 4), "mood_2", "mood_3")
)
head(info)
#>  theta test_information se  
#>  -4    0.00833          11  
#>  -3.9  0.0126           8.91
#>  -3.8  0.0187           7.31
#>  -3.7  0.0275           6.04
#>  -3.6  0.0396           5.03
#>  -3.5  0.056            4.23

# Where does this three-item set measure most precisely?
attr(info, "theta_max_information")
#> [1] 0.5

# Each item's contribution; the rows sum to the test information.
head(attr(info, "item_information"))
#>           mood_1      mood_2       mood_3
#> -4   0.006419777 0.001906444 2.328048e-06
#> -3.9 0.010014688 0.002567436 4.055379e-06
#> -3.8 0.015303102 0.003428216 6.976757e-06
#> -3.7 0.022905443 0.004538642 1.185366e-05
#> -3.6 0.033582642 0.005957571 1.988939e-05
#> -3.5 0.048230063 0.007753454 3.295750e-05

# A dichotomous item matches the two parameter normal ogive closed form.
one <- irt_information(a = 1.5, b = 0.25, theta = c(-1, 0, 1))
z <- 1.5 * (c(-1, 0, 1) - 0.25)
1.5^2 * dnorm(z)^2 / (pnorm(z) * (1 - pnorm(z)))
#> [1] 0.3612187 1.3607816 0.8913543
one$test_information
#> [1] 0.3612187 1.3607816 0.8913543
```
