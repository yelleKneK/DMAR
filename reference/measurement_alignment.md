# Approximate Measurement Invariance by Factor Alignment

Estimates the group factor means and factor variances that make the
measurement parameters as nearly invariant as possible across groups,
following the alignment method of Asparouhov and Muthén (2014). The
invariance ladder in
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
asks whether loadings and intercepts are exactly equal across groups.
With many groups that hypothesis is essentially never true, the ladder
stalls at configural invariance, and the comparison of factor means the
researcher wanted never happens. Alignment takes the configural solution
as given and searches for the group factor means and variances that
concentrate the noninvariance in a few parameters instead of spreading
it thinly over many, so that factor means stay comparable without
claiming that exact invariance holds. Requires lavaan.

## Usage

``` r
measurement_alignment(
  data,
  items,
  group,
  model = NULL,
  alignment = c("fixed", "free"),
  estimator = "ML",
  n_starts = 10,
  seed = NULL,
  epsilon = 0.01
)
```

## Arguments

- data:

  A `data.frame` holding the items and the grouping variable.

- items:

  Character vector naming the indicator columns of `data`. Three or more
  items are required, and each must be numeric.

- group:

  Single character string naming the grouping column of `data`. Two or
  more groups are required.

- model:

  Optional lavaan model syntax for the measurement model, for the case
  where the default single factor over `items` is not what is wanted (a
  residual covariance, for example). The syntax must define exactly one
  latent variable and its indicators must be exactly `items`. The factor
  is standardized (mean 0, variance 1) in every group by this function,
  so do not set its scale in the syntax. When `NULL` (the default) the
  model `f =~ item1 + item2 + ...` is used.

- alignment:

  Identification rule for the metric of the factor, either `"fixed"`
  (the default) or `"free"`. Under `"fixed"` the first group's factor
  mean is held at 0 and its factor variance at 1. Under `"free"` the
  group factor means are constrained to average 0 and the group factor
  standard deviations to have a geometric mean of 1, so no group serves
  as the reference. The choice sets the metric the factor means and
  variances are reported on; see Details.

- estimator:

  Estimator passed to lavaan for the configural model. Defaults to
  `"ML"`; `"MLR"` supplies the robust corrections.

- n_starts:

  Number of starting values for the optimizer, a positive integer
  (default 10). The first start sets every factor mean to 0 and every
  factor variance to 1, the second uses a median based heuristic, and
  any remaining starts are random. The simplicity function has local
  minima, so several starts is the default rather than a refinement.

- seed:

  Optional integer seed for the random starts. Defaults to `NULL`, which
  uses the caller's random number state and leaves it alone. When an
  integer is supplied the seed is set locally and the caller's state is
  restored on exit.

- epsilon:

  Smoothing constant \\\epsilon \> 0\\ inside the component loss
  function, default 0.01. Smaller values push the loss closer to
  \\\sqrt{\|x\|}\\ and make the surface harder to optimize; larger
  values smooth it at the cost of blurring the distinction between a few
  large differences and many small ones.

## Value

A `data.frame` (class `dmar_tbl`) with one row per group and columns
`group` (the group label), `n` (the number of cases the configural model
used in that group), `factor_mean` (the estimated \\\alpha_g\\), and
`factor_variance` (the estimated \\\psi_g\\). The rows follow the group
ordering lavaan uses.

Attributes carry the rest of the solution:

- `aligned_loadings`, `aligned_intercepts`:

  Matrices of the aligned \\\lambda\_{gi}\\ and \\\nu\_{gi}\\, groups in
  rows and items in columns.

- `configural_loadings`, `configural_intercepts`:

  Matrices of the configural \\\lambda^{0}\_{gi}\\ and
  \\\nu^{0}\_{gi}\\, same layout.

- `simplicity_function`:

  The achieved minimum of \\F\\.

- `simplicity_starts`:

  The value of \\F\\ achieved from each starting value, in start order.
  A start the optimizer could not complete, or one that ended in the
  degenerate branch described in Details, is recorded as `Inf` and
  discarded.

- `converged`:

  `TRUE` when the optimizer reported convergence at the retained
  solution.

- `n_starts`:

  The number of starting values used.

- `n_optima`:

  The number of distinct local minima the converged starts reached.

- `alignment`:

  The identification rule used, `"fixed"` or `"free"`.

- `epsilon`:

  The smoothing constant used.

- `R2_loadings`, `R2_intercepts`:

  Named numeric vectors, one entry per item, holding the per item
  \\R^2\\ invariance measures defined in Details.

- `R2_total`:

  The same two measures pooled over items, as a named numeric vector
  with elements `loadings` and `intercepts`. This is the overall effect
  size of approximate invariance reported by Asparouhov and Muthén
  (2014).

- `item_loss`:

  Matrix with one row per item and columns `loadings`, `intercepts`, and
  `total`, splitting the achieved simplicity function into per item
  contributions.

- `fit`:

  The configural lavaan fit object.

## Details

Alignment starts from the *configural* solution: the multiple group
single factor model with every loading, intercept, and residual variance
free across groups and the factor standardized (mean 0 and variance 1)
in each group. Write \\\lambda^{0}\_{gi}\\ and \\\nu^{0}\_{gi}\\ for the
resulting loading and intercept of item \\i\\ in group \\g\\. That
solution is not unique: for any group factor means \\\alpha_g\\ and
factor variances \\\psi_g\\ the reparameterized measurement parameters
\$\$\lambda\_{gi} = \lambda^{0}\_{gi} / \sqrt{\psi_g}\$\$ \$\$\nu\_{gi}
= \nu^{0}\_{gi} - \alpha_g \lambda^{0}\_{gi} / \sqrt{\psi_g}\$\$
reproduce the observed means and covariances exactly and so fit the data
identically. Alignment picks the member of that family whose measurement
parameters are closest to invariant, by minimizing the total simplicity
function \$\$F = \sum_i \sum\_{g_1 \< g_2} w\_{g_1 g_2} f(\lambda\_{g_1
i} - \lambda\_{g_2 i}) + \sum_i \sum\_{g_1 \< g_2} w\_{g_1 g_2}
f(\nu\_{g_1 i} - \nu\_{g_2 i})\$\$ over \\\alpha_g\\ and \\\psi_g\\,
with weights \\w\_{g_1 g_2} = \sqrt{N\_{g_1} N\_{g_2}}\\ and the
component loss function \$\$f(x) = \sqrt{\sqrt{x^2 + \epsilon}} = (x^2 +
\epsilon)^{1/4}.\$\$

The fourth root is the substance of the method, not a technical detail.
A squared loss would spread a fixed amount of noninvariance evenly over
all the parameters, because halving two differences beats zeroing one.
The fourth root is concave in \\\|x\|\\, so the marginal penalty falls
as a difference grows: the criterion prefers a solution in which most
parameters agree closely and a few disagree substantially, which is
exactly the pattern of approximate invariance a researcher wants to find
and report. The constant \\\epsilon\\ only rounds off the kink of
\\\sqrt{\|x\|}\\ at zero so the criterion is differentiable. One
consequence is worth knowing: each of the \\2I\\ terms in a group pair
contributes at least \\\epsilon^{1/4}\\, so the smallest attainable
value of \\F\\ is \\2 I \epsilon^{1/4} \sum\_{g_1 \< g_2} w\_{g_1
g_2}\\, reached only when every aligned parameter is exactly invariant.
The achieved value is comparable across runs on the same items, groups,
and \\\epsilon\\, not across data sets.

Two constraints identify the solution. The scale of the factor is
genuinely undetermined by \\F\\: multiplying every \\\sqrt{\psi_g}\\ and
every \\\alpha_g\\ by the same constant leaves the aligned intercepts
alone and shrinks the aligned loadings toward each other, so without a
constraint the criterion is minimized by letting the factor variances
run away. The location constraint plays the same role in the limiting
case of exact invariance, where shifting all the factor means by a
constant leaves \\F\\ unchanged. Under `alignment = "fixed"` the
constraints are \\\alpha_1 = 0\\ and \\\psi_1 = 1\\. Under
`alignment = "free"` they are \\\sum_g \alpha_g = 0\\ and \\\prod_g
\sqrt{\psi_g} = 1\\, which treats the groups symmetrically and is the
more natural choice when no group is a meaningful reference. The two
rules give genuinely different solutions, not a relabeling of one
another, because \\F\\ is not invariant to shifting the factor means.

The rule also sets the metric the answer is reported on, which matters
when the estimates are compared with anything else. The constraint
\\\psi_1 = 1\\ makes the first group's factor the common metric, so a
reported factor mean of 0.4 says that group's factor mean is 0.4 of the
*first group's* factor standard deviations above the first group's, and
a reported factor variance of 1.3 says its factor variance is 1.3 times
the first group's. Under `"free"` the common metric is the one in which
the group factor standard deviations have a geometric mean of 1, and the
factor means are deviations from their own average on that metric.
Simulating data with known group factor means and then comparing them
with the estimates requires dividing the generating means by the
generating factor standard deviation of the reference group first.

The minimization runs [`optim`](https://rdrr.io/r/stats/optim.html) with
the BFGS method and the analytic gradient, over the factor means
\\\alpha_g\\ and the log factor standard deviations \\\log
\sqrt{\psi_g}\\, reduced to the coordinates the identification rule
leaves free. The optimizer is run from `n_starts` starting values and
the best solution is kept, because a single start is not safe. Two
things can go wrong, and they are different problems.

The first is ordinary multimodality. When several loadings and
intercepts are noninvariant the criterion has more than one local
minimum, typically within a percent or so of each other in value but at
different factor means. The number of distinct minima the starts reached
is returned in the `"n_optima"` attribute and the value achieved from
each start in `"simplicity_starts"`. More than one is a signal to raise
`n_starts` and to check whether the reported solution is stable.

The second is a degenerate branch, and it is the one that bites. Sending
the factor variances of every group but the reference off to infinity
drives those groups' aligned loadings to zero, which flattens the
loading half of the criterion; the optimizer stops there and reports
convergence with factor variances of \\10^{36}\\ or `Inf`. That is not
an estimate, and on real data a meaningful share of random starts finds
it. A start whose factor standard deviations or factor means leave a
very wide sanity range (a factor of \\10^4\\ either way) is therefore
discarded and its entry in `"simplicity_starts"` is `Inf`. If every
start is discarded the function stops rather than return the degenerate
solution.

Item level invariance is summarized with the \\R^2\\ measure of
Asparouhov and Muthén (2014), which asks how much of the group to group
variation in an item's configural parameter is accounted for by the
estimated factor means and variances alone. Let \\\bar{\lambda}\_i\\ and
\\\bar{\nu}\_i\\ be the aligned parameters of item \\i\\ averaged over
groups. The factor means and variances by themselves imply the
configural values \\\sqrt{\psi_g} \bar{\lambda}\_i\\ and
\\\bar{\nu}\_i + \alpha_g \bar{\lambda}\_i\\. Write \\T_i\\ for the sum
over groups of the squared configural parameter of item \\i\\ and
\\E_i\\ for the sum over groups of its squared departure from that
implied value. Then \$\$R^2_i = 1 - E_i / T_i,\$\$ computed separately
for the loadings and for the intercepts. A value of 1 for every item on
the loadings is metric invariance, and on both loadings and intercepts
is scalar invariance. The complementary view is the `"item_loss"`
attribute, which splits the achieved simplicity function into the
contribution of each item, so the items carrying the noninvariance can
be named.

The criterion adds up raw parameter differences, so items on very
different measurement scales do not contribute equally: an item whose
raw variance is a hundred times another's dominates the sum. When the
items are not already on a common metric, put them on one (standardizing
them is the simple choice) before aligning.

The method is defined for two groups and is computed here for two, but
it has little to offer there. With two groups the standard invariance
ladder is tractable and interpretable, and the alignment criterion has
only one group pair to work with. Alignment earns its keep when the
number of groups makes exact invariance implausible and the ladder
uninformative.

## References

Asparouhov, T., & Muthén, B. (2014). Multiple-group factor analysis
alignment. *Structural Equation Modeling, 21*(4), 495–508.
[doi:10.1080/10705511.2014.919210](https://doi.org/10.1080/10705511.2014.919210)

Marsh, H. W., Guo, J., Parker, P. D., Nagengast, B., Asparouhov, T.,
Muthén, B., & Dicke, T. (2018). What to do when scalar invariance fails:
The extended alignment method for multi-group factor analysis comparison
of latent means across many groups. *Psychological Methods, 23*(3),
524–545. [doi:10.1037/met0000113](https://doi.org/10.1037/met0000113)

Muthén, B., & Asparouhov, T. (2014). IRT studies of many groups: The
alignment method. *Frontiers in Psychology, 5*, Article 978.
[doi:10.3389/fpsyg.2014.00978](https://doi.org/10.3389/fpsyg.2014.00978)

Muthén, B., & Asparouhov, T. (2018). Recent methods for the study of
measurement invariance with many groups: Alignment and random effects.
*Sociological Methods & Research, 47*(4), 637–664.
[doi:10.1177/0049124117701488](https://doi.org/10.1177/0049124117701488)

Robitzsch, A. (2025). *sirt: Supplementary item response theory models*.
R package version 4.2-133. <https://CRAN.R-project.org/package=sirt>

## See also

[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
for the exact invariance ladder alignment is meant to rescue;
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) for the
single group measurement model;
[`htmt`](https://yelleknek.github.io/DMAR/reference/htmt.md) for
discriminant validity of the same items.

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
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Five groups that differ in factor mean and factor variance, with two
# deliberately noninvariant measurement parameters (the loading of item
# 3 in group 2 and the intercept of item 5 in group 4). The first group's
# factor standard deviation is 1, which is the metric the default "fixed"
# rule reports on, so the estimates compare directly with the generating
# values below.
set.seed(113)
n_g <- c(400, 450, 380, 500, 420)
factor_mean <- c(0, 0.30, -0.50, 0.80, 0.20)
factor_sd   <- c(1, 1.20,  0.80, 1.10, 0.90)
Lambda <- matrix(0.8, nrow = 5, ncol = 6)
Nu     <- matrix(1.0, nrow = 5, ncol = 6)
Lambda[2, 3] <- 0.3
Nu[4, 5]     <- 1.8
d <- do.call(rbind, lapply(1:5, function(g) {
  eta <- rnorm(n_g[g], factor_mean[g], factor_sd[g])
  x <- sapply(1:6, function(i)
    Nu[g, i] + Lambda[g, i] * eta + rnorm(n_g[g], 0, 0.6))
  data.frame(x, cohort = paste0("cohort_", g))
}))
names(d)[1:6] <- paste0("x", 1:6)

out <- measurement_alignment(d, items = paste0("x", 1:6),
                             group = "cohort", seed = 113)
out                                   # recovered means and variances
#>  group    n   factor_mean factor_variance
#>  cohort_1 400 0           1              
#>  cohort_2 450 0.328       1.53           
#>  cohort_3 380 -0.54       0.581          
#>  cohort_4 500 0.861       1.28           
#>  cohort_5 420 0.213       0.783          
attr(out, "R2_loadings")              # item 3 stands out
#>        x1        x2        x3        x4        x5        x6 
#> 0.9999162 0.9976032 0.9079726 0.9994612 0.9996394 0.9983087 
attr(out, "R2_intercepts")            # item 5 stands out
#>        x1        x2        x3        x4        x5        x6 
#> 0.9999097 0.9989362 0.9968054 0.9995481 0.9525477 0.9992106 
attr(out, "item_loss")
#>    loadings intercepts    total
#> x1 1359.250   1371.846 2731.096
#> x2 1482.405   1470.454 2952.859
#> x3 2028.494   1381.886 3410.381
#> x4 1389.880   1418.496 2808.376
#> x5 1374.485   2399.982 3774.468
#> x6 1420.217   1427.319 2847.536

# Holzinger and Swineford's verbal tests across four groups formed by
# crossing school with sex. The raw tests are on very different scales,
# so they are standardized first (see Details).
data(holzinger_swineford)
hs <- holzinger_swineford
hs$school_sex <- interaction(hs$school, hs$sex, sep = ", ")
verbal <- c("t5_general_information", "t6_paragraph_comprehension",
            "t7_sentence", "t8_word_classification", "t9_word_meaning")
hs[verbal] <- scale(hs[verbal])
ma <- measurement_alignment(hs, items = verbal, group = "school_sex",
                            seed = 113)
ma
#>  group               n  factor_mean factor_variance
#>  Pasteur, Male       74 0           1              
#>  Pasteur, Female     82 0.0429      0.904          
#>  Grant-White, Male   72 0.608       0.802          
#>  Grant-White, Female 73 0.803       1.09           

# The broom verbs: one row per group, and the alignment summary.
generics::tidy(ma)
#>                  term  n factor_mean factor_variance
#> 1       Pasteur, Male 74  0.00000000       1.0000000
#> 2     Pasteur, Female 82  0.04288633       0.9043518
#> 3   Grant-White, Male 72  0.60761608       0.8023580
#> 4 Grant-White, Female 73  0.80291885       1.0856178
generics::glance(ma)
#>   n_groups alignment epsilon simplicity_function R2_loadings_mean
#> 1        4     fixed    0.01            1757.193        0.9942458
#>   R2_intercepts_mean  R2_total converged n_starts n_optima
#> 1          0.8839346 0.9345257      TRUE       10        1
```
