# Graded Response Model for Ordered Categorical Items

Estimates Samejima's (1969) graded response model for a set of ordered
categorical items (Likert items, symptom severity ratings, rubric scored
performance items) and reports each item's discrimination and its
category boundary locations. The model is fitted as the single-factor
categorical factor analysis model it provably is (Takane and de Leeuw,
1987), using lavaan's categorical estimator on the polychoric
correlations, and the solution is then converted to the normal ogive (or
logistic) item response theory parameterization. A researcher who
already fits confirmatory factor analysis models therefore gets item
response theory item parameters without adopting a second estimation
engine, and the two analyses of the same items stay in one modeling
tradition.

## Usage

``` r
irt_grm(
  data,
  items = NULL,
  estimator = "WLSMV",
  metric = c("normal_ogive", "logistic")
)
```

## Arguments

- data:

  A `data.frame` or matrix of ordered item responses, one column per
  item, coded with integer category values (for example 1, 2, 3, 4, 5).
  Rows with any missing value on the analyzed items are
  listwise-deleted. Each item must have at least 2 and at most 20
  distinct observed categories; a column with more than 20 distinct
  values, or with non-integer values, is treated as continuous and
  rejected.

- items:

  Optional character vector naming the columns of `data` to analyze.
  Defaults to `NULL`, which uses every column. At least 3 items are
  required.

- estimator:

  Character; the lavaan estimator for categorical data. The choices
  differ in the weight matrix applied to the polychoric correlations and
  in whether the test statistic is corrected. `"WLSMV"` (the default) is
  diagonally weighted least squares with a mean- and variance-adjusted
  test statistic (Muthén, 1984; Muthén, du Toit, & Spisic, 1997), the
  standard estimator for ordered categorical items; `"WLSM"` applies the
  mean adjustment only. `"DWLS"` is the same diagonal-weight estimator
  with no correction. `"WLS"` uses the full weight matrix (the
  asymptotic distribution free approach; Browne, 1984), which is
  unstable unless the sample is large relative to the number of
  thresholds. `"ULS"`, unweighted least squares, uses an identity weight
  matrix, an option worth considering in small samples where even the
  diagonal weights are noisy; `"ULSMV"` and `"ULSM"` add the corrected
  test statistics. When in doubt keep the default.

- metric:

  Which discrimination metric to report in the `a` column:
  `"normal_ogive"` (default) or `"logistic"`. Both are always computed
  and the one not reported in `a` is attached as an attribute, so the
  returned table has the same columns and the same number of rows either
  way. The boundary locations `b` are identical under the two metrics.

## Value

A `data.frame` (class `dmar_tbl`) with one row per item and category
boundary and the columns

- `item`:

  Item name, taken from the column name.

- `factor`:

  Name of the latent variable, the same for every row in this
  unidimensional model.

- `category`:

  Boundary index \\k\\, running from 1 to one fewer than the item's
  number of categories.

- `lambda`:

  Standardized factor loading \\\lambda_i\\ of the item's latent
  response variate on the factor, repeated across the item's boundaries.

- `tau`:

  Standardized threshold \\\tau\_{ik}\\.

- `a`:

  Discrimination in the metric named by `metric`, repeated across the
  item's boundaries.

- `b`:

  Boundary location \\b\_{ik}\\ on the \\\theta\\ scale.

The attributes are `"fit"` (the fitted lavaan object), `"fit_measures"`
(the full named numeric vector from
[`lavaan::fitMeasures`](https://rdrr.io/pkg/lavaan/man/fitMeasures.html),
unrounded), `"metric"` (the reported discrimination metric),
`"estimator"`, `"n_categories"` (named integer vector of the number of
observed categories per item), `"N"` (the analyzed sample size),
`"factor_sign_flipped"` (a single logical recording whether the
direction of the latent variable was reversed to satisfy the sign
convention described in Details), and whichever of `"a_logistic"` or
`"a_normal_ogive"` was not reported in the `a` column (a named numeric
vector, one element per item).

## Details

**One model, two parameterizations.** Samejima's (1969) graded response
model and the single-factor categorical factor analysis model of Muthén
(1984) are the same model written in different parameterizations; Takane
and de Leeuw (1987) proved the equivalence, and Kamata and Bauer (2008)
give the algebra item by item. Each observed response \\X_i\\ is a
categorization of a latent continuous response variate \\X_i^{\*}\\ at
thresholds \\\tau\_{ik}\\, and \\X_i^{\*} = \lambda_i \theta +
\varepsilon_i\\ with \\\theta\\ standard normal and \\X_i^{\*}\\
standardized. Fitting that model on the polychoric correlations and
converting the solution gives the normal ogive graded response model
directly. For item \\i\\ with standardized loading \\\lambda_i\\ and
standardized thresholds \\\tau\_{ik}\\, \$\$a_i =
\frac{\lambda_i}{\sqrt{1 - \lambda_i^2}}, \qquad b\_{ik} =
\frac{\tau\_{ik}}{\lambda_i}.\$\$ The boundary response function is
\$\$P^{\*}\_{ik}(\theta) = \Phi\\\left\[a_i (\theta -
b\_{ik})\right\],\$\$ the probability of responding above boundary
\\k\\, with \\P^{\*}\_{i0}(\theta) \equiv 1\\ and \\P^{\*}\_{iK}(\theta)
\equiv 0\\; the probability of the individual category is the difference
of adjacent boundary functions, \\P\_{ik}(\theta) =
P^{\*}\_{i,k-1}(\theta) - P^{\*}\_{ik}(\theta)\\. When \\a_i \> 0\\,
that is, for an item keyed in the same direction as the rest of the
scale, \\P^{\*}\_{ik}\\ is monotone increasing in \\\theta\\, the
boundary locations of the item are ordered, \\b\_{i1} \< b\_{i2} \<
\cdots\\, and \\b\_{ik}\\ is the value of \\\theta\\ at which the
probability of responding above boundary \\k\\ reaches 0.50. An item
keyed in the opposite direction has \\\lambda_i \< 0\\, hence \\a_i \<
0\\ and boundary locations that run from high to low; see the two
paragraphs on direction below.

**The direction of the latent variable.** A single-factor model fixes
\\\theta\\ only up to its direction. Relabeling \\\theta\\ as
\\-\theta\\ changes the sign of every loading and leaves the fitted
model, the thresholds, and every fit measure exactly as they were, so it
is a renaming of the latent direction rather than a different model. The
thresholds are untouched because \\\tau\_{ik}\\ cuts the item's own
latent response variate \\X_i^{\*}\\, which the relabeling does not
move; the sign change therefore passes straight through to \\a_i =
\lambda_i / \sqrt{1 - \lambda_i^2}\\ and to \\b\_{ik} = \tau\_{ik} /
\lambda_i\\, both of which change sign. lavaan returns whichever
direction its starting values point toward, and for a scale that
contains a reverse-keyed item that direction can turn on something as
incidental as the order of the columns. The solution is therefore put in
a fixed direction before it is converted: if the standardized loadings
sum to a negative number the whole factor is flipped, so that \\\theta\\
runs in the direction the scale as a whole measures. The result is the
same table no matter how the columns are ordered. Whether the flip was
applied is recorded on the `"factor_sign_flipped"` attribute. The lavaan
object on the `"fit"` attribute is the fit as lavaan produced it, so
when a flip was applied its loadings carry the opposite sign to the
`lambda` column.

**Reverse-keyed items.** An item whose loading is still negative after
the direction is fixed is keyed opposite to the rest of the scale, which
is a property of the item rather than an artifact of the sign
indeterminacy. Its discrimination is negative and its boundary locations
run from high to low, so it does not satisfy the graded response model
as written above and its parameters do not belong on the same scale as
the others. Such items are named in a warning. Reverse score them (for
example `x <- (min(x) + max(x)) - x`) and refit; that puts the item in
the direction the rest of the scale measures and restores \\a_i \> 0\\
and the ordering \\b\_{i1} \< b\_{i2} \< \cdots\\.

**The two discrimination metrics and the constant 1.702.** The
conversion above puts \\a_i\\ in the normal ogive metric, where the
boundary function is a normal cumulative distribution function. The item
response theory literature more often writes the graded response model
with a logistic boundary function, and the two agree closely once the
logistic argument is stretched by a scaling constant: \\\|\Phi(x) -
\Psi(1.702 x)\| \< 0.01\\ for every \\x\\, where \\\Psi\\ is the
standard logistic cumulative distribution function. The value 1.702 is
the constant that minimizes that maximum discrepancy (Haley, 1952; see
Camilli, 1994, for the history), so \\a_i(\mathrm{logistic}) = 1.702 \\
a_i(\mathrm{normal\\ ogive})\\ and software that reports logistic slopes
(for example mirt and the classical two parameter logistic tradition)
gives values about 1.7 times larger for the same items. The scaling
multiplies the slope and leaves the location alone, so \\b\_{ik}\\ does
not depend on the metric.

**Estimation and what to expect.** lavaan estimates the thresholds and
the polychoric correlations, then fits the single-factor model to those
correlations by (diagonally) weighted least squares. This is limited
information estimation: it uses the univariate and bivariate margins of
the response table, whereas marginal maximum likelihood (the usual item
response theory approach, as in mirt) uses the full response pattern
likelihood. The two are consistent for the same population parameters
and agree closely in practice, but they are different estimators and
will not return identical numbers on a finite sample. Limited
information estimation scales well to many items and brings the whole
apparatus of factor analysis fit assessment (CFI, TLI, RMSEA) along with
it; the fit measures are returned on the `"fit_measures"` attribute and
the lavaan object itself on `"fit"`, so any lavaan accessor can be
applied to the result.

The model is unidimensional by construction. A standardized loading at
or beyond one is an improper (Heywood) solution: the implied
discrimination is infinite and the conversion is not interpretable. That
case is flagged with a warning rather than silently returned as a
number.

This function requires lavaan to be installed.

## References

Camilli, G. (1994). Teacher's corner: Origin of the scaling constant *d*
= 1.7 in item response theory. *Journal of Educational and Behavioral
Statistics, 19*(3), 293–295.
[doi:10.3102/10769986019003293](https://doi.org/10.3102/10769986019003293)

Haley, D. C. (1952). *Estimation of the dosage mortality relationship
when the dose is subject to error* (Technical Report No. 15). Applied
Mathematics and Statistics Laboratory, Stanford University.

Kamata, A., & Bauer, D. J. (2008). A note on the relation between factor
analytic and item response theory models. *Structural Equation Modeling,
15*(1), 136–153.
[doi:10.1080/10705510701758406](https://doi.org/10.1080/10705510701758406)

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Samejima, F. (1969). Estimation of latent ability using a response
pattern of graded scores. *Psychometrika Monograph Supplement, 34*(4,
Pt. 2), 1–97.

Takane, Y., & de Leeuw, J. (1987). On the relationship between item
response theory and factor analysis of discretized variables.
*Psychometrika, 52*(3), 393–408.

Wirth, R. J., & Edwards, M. C. (2007). Item factor analysis: Current
approaches and future directions. *Psychological Methods, 12*(1), 58–79.
[doi:10.1037/1082-989X.12.1.58](https://doi.org/10.1037/1082-989X.12.1.58)

## See also

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) (the same
single-factor model reported in the factor analysis parameterization),
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
(reliability for the same class of items),
[`cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html).

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
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Six five-category items generated from a known graded response model.
set.seed(113)
n <- 800
a_pop <- c(1.2, 0.9, 1.5, 1.0, 1.3, 1.1)
b_pop <- rbind(c(-1.6, -0.6, 0.3, 1.2), c(-1.4, -0.4, 0.5, 1.5),
               c(-1.8, -0.7, 0.2, 1.1), c(-1.2, -0.2, 0.7, 1.6),
               c(-1.5, -0.5, 0.4, 1.3), c(-1.3, -0.3, 0.6, 1.4))
theta <- rnorm(n)
responses <- vapply(seq_along(a_pop), function(i) {
  p_star <- outer(theta, b_pop[i, ], function(z, b) pnorm(a_pop[i] * (z - b)))
  as.integer(1 + rowSums(runif(n) < p_star))
}, integer(n))
colnames(responses) <- paste0("item", seq_along(a_pop))
responses <- as.data.frame(responses)

# Item parameters in the normal ogive metric.
grm <- irt_grm(responses)
grm
#>  item  factor category lambda tau    a     b     
#>  item1 theta  1        0.779  -1.3   1.24  -1.66 
#>  item1 theta  2        0.779  -0.496 1.24  -0.636
#>  item1 theta  3        0.779  0.215  1.24  0.275 
#>  item1 theta  4        0.779  0.869  1.24  1.11  
#>  item2 theta  1        0.64   -0.979 0.832 -1.53 
#>  item2 theta  2        0.64   -0.283 0.832 -0.442
#>  item2 theta  3        0.64   0.279  0.832 0.437 
#>  item2 theta  4        0.64   0.954  0.832 1.49  
#>  item3 theta  1        0.856  -1.49  1.65  -1.74 
#>  item3 theta  2        0.856  -0.617 1.65  -0.721
#>  item3 theta  3        0.856  0.186  1.65  0.217 
#>  item3 theta  4        0.856  0.873  1.65  1.02  
#>  item4 theta  1        0.694  -0.954 0.964 -1.38 
#>  item4 theta  2        0.694  -0.208 0.964 -0.3  
#>  item4 theta  3        0.694  0.423  0.964 0.609 
#>  item4 theta  4        0.694  1.2    0.964 1.73  
#>  item5 theta  1        0.795  -1.21  1.31  -1.52 
#>  item5 theta  2        0.795  -0.392 1.31  -0.493
#>  item5 theta  3        0.795  0.322  1.31  0.405 
#>  item5 theta  4        0.795  1.04   1.31  1.31  
#>  item6 theta  1        0.742  -0.989 1.11  -1.33 
#>  item6 theta  2        0.742  -0.237 1.11  -0.32 
#>  item6 theta  3        0.742  0.419  1.11  0.565 
#>  item6 theta  4        0.742  0.989  1.11  1.33  

# The generating discriminations, for comparison.
a_pop
#> [1] 1.2 0.9 1.5 1.0 1.3 1.1

# Model fit travels with the item parameters.
attr(grm, "fit_measures")[c("cfi", "tli", "rmsea", "srmr")]
#>        cfi        tli      rmsea       srmr 
#> 1.00000000 1.00131440 0.00000000 0.01127189 

# Two further calls, each of which refits the model and so is not run
# here. The first reports the same fit with logistic slopes (about
# 1.702 times larger than the normal ogive slopes above), the second
# fits a subset of the items selected by name:
# irt_grm(responses, metric = "logistic")
# irt_grm(responses, items = c("item1", "item3", "item5"))

# The boundary response function of the first item at theta = 0.
first <- grm[grm$item == "item1", ]
pnorm(first$a * (0 - first$b))
#> [1] 0.98069860 0.78564848 0.36591990 0.08279093
```
