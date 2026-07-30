# The dMACS Effect Size of Measurement Noninvariance

Quantifies how much a violation of measurement invariance actually
matters for an item, rather than only whether it is statistically
detectable. A likelihood ratio or score test can flag a loading or
intercept difference that is too small to change anyone's score in a
meaningful way, and with a large sample size it usually will. The dMACS
index of Nye and Drasgow (2011) answers the size question directly: it
is the expected difference between the reference group's and the focal
group's measurement equations for that item, averaged over the focal
group's latent distribution and standardized by the pooled item standard
deviation, so it reads on the familiar standardized mean difference
scale. Input is either a fitted multiple group lavaan model (requires
lavaan) or the loadings, intercepts, and pooled standard deviations a
paper reports.

## Usage

``` r
dmacs(
  fit = NULL,
  reference = NULL,
  focal = NULL,
  lambda_reference = NULL,
  lambda_focal = NULL,
  nu_reference = NULL,
  nu_focal = NULL,
  mean_focal = 0,
  sd_focal = 1,
  sd_pooled = NULL,
  item_names = NULL
)
```

## Arguments

- fit:

  A fitted multiple group lavaan object carrying a mean structure, with
  the two groups' loadings and intercepts on a common metric (see
  Details). Supply either `fit` or the parameter vectors, never both.

- reference, focal:

  Which groups play the reference and focal roles, each given as a
  single group label or a single group index in
  `lavInspect(fit, "group.label")`. When both are `NULL` (default) and
  the fit has exactly two groups, the first is the reference and the
  second the focal; with more than two groups they must be named.

- lambda_reference, lambda_focal:

  Numeric vectors of unstandardized loadings, one per item, in the
  reference and focal groups. Any finite values are admissible.

- nu_reference, nu_focal:

  Numeric vectors of unstandardized intercepts, one per item, in the
  reference and focal groups, in the same item order as the loadings.
  Any finite values are admissible.

- mean_focal:

  Focal group latent mean \\\mu_F\\ on the common metric, a single
  finite number (default `0`, the usual identification in which the
  reference group's latent mean is fixed at zero). Used only on the
  parameter vector path; with a `fit` it is read from the fit.

- sd_focal:

  Focal group latent standard deviation \\\sigma_F\\ on the common
  metric, a single positive number (default `1`). Used only on the
  parameter vector path; with a `fit` it is read from the fit.

- sd_pooled:

  Pooled observed standard deviation of each item across the two groups,
  either one positive number applied to every item or one per item.
  Required on the parameter vector path; with a `fit` it is computed
  from the group sample sizes and observed item variances.

- item_names:

  Optional character vector of item labels, one per item. Defaults to
  the names carried by the parameter vectors, to the indicator names in
  the fit, or to `item_1`, `item_2`, and so on.

## Value

A wide `data.frame` (class `dmar_tbl`) with one row per item and columns

- `item`:

  Item label.

- `lambda_reference`:

  Reference group unstandardized loading.

- `lambda_focal`:

  Focal group unstandardized loading.

- `nu_reference`:

  Reference group unstandardized intercept.

- `nu_focal`:

  Focal group unstandardized intercept.

- `sd_pooled`:

  Pooled observed standard deviation of the item.

- `dmacs`:

  The dMACS effect size, nonnegative.

The returned object carries four attributes: `"reference"` and
`"focal"`, the two group labels, and `"mean_focal"` and `"sd_focal"`,
the focal group's latent mean and standard deviation used in the
integral (named by latent variable on the fit path).

## Details

For item \\i\\, let \\\nu_R\\ and \\\lambda_R\\ be the reference group's
intercept and loading, let \\\nu_F\\ and \\\lambda_F\\ be the focal
group's, and let the focal group's latent variable be \\\eta \sim
N(\mu_F, \sigma_F^2)\\ with density \\f_F\\. Each group's measurement
equation gives an expected item score at every value of \\\eta\\, and
dMACS is the root mean squared vertical distance between those two lines
over the focal group's latent distribution, divided by the pooled item
standard deviation: \$\$d\_{MACS, i} = \frac{1}{SD_i} \sqrt{\int \left\[
(\nu_R + \lambda_R \eta) - (\nu_F + \lambda_F \eta) \right\]^2 f_F(\eta)
\\ d\eta}.\$\$

Writing \\a = \nu_R - \nu_F\\ for the intercept difference and \\b =
\lambda_R - \lambda_F\\ for the loading difference, the integrand is
\\(a + b\eta)^2\\ and the integral is the second moment of a linear
function of a normal variate, so it has the closed form \\a^2 +
2ab\mu_F + b^2(\sigma_F^2 + \mu_F^2)\\. No numerical integration is
performed. The standardizer is the pooled observed standard deviation of
the item, \$\$SD_i = \sqrt{\frac{(n_R - 1)s_R^2 + (n_F - 1)s_F^2}{n_R +
n_F - 2}},\$\$ the same pooling used by the standardized mean
difference, which puts dMACS on a scale a reader of
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) already
understands.

dMACS is a root mean square and so is nonnegative by construction; it
carries no sign and is not given one here. The direction of the
violation is read from the returned components: \\\nu_R - \nu_F\\ says
which group is scored higher at the mean of the latent variable, and
\\\lambda_R - \lambda_F\\ says in which group the item discriminates
more sharply. When only the intercepts differ, the integral collapses to
\\a^2\\ and dMACS reduces to \\\|a\| / SD_i\\, a plain standardized
intercept difference.

The index is interpretable only when the two groups' loadings and
intercepts are expressed on a common metric. In practice that means a
partial invariance model in which a set of anchor items is constrained
equal across groups while the suspect items are freed, and the focal
group's latent mean and variance are freely estimated. A configural
model, which sets each group's latent scale separately, does not put the
groups on a common metric, and dMACS computed from one is not
meaningful. The usual workflow is therefore
[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
to locate where the ladder breaks, a partial invariance refit that frees
the offending parameters, and then `dmacs()` on that refit to judge
whether the violation is large enough to matter.

On the fit path each group's item variance is computed from the data in
the fit with the usual \\n - 1\\ divisor, using the number of nonmissing
observations for that item in that group. When the model was fitted from
sample moments rather than raw data, the sample covariances in the fit
are used instead, rescaled to the \\n - 1\\ divisor when the fit's
likelihood option calls for it.

## References

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525–543.
[doi:10.1007/BF02294825](https://doi.org/10.1007/BF02294825)

Millsap, R. E. (2011). *Statistical approaches to measurement
invariance*. Routledge.

Nye, C. D., Bradburn, J., Olenick, J., Bialko, C., & Drasgow, F. (2019).
How big are my effects? Examining the magnitude of effect sizes in
studies of measurement equivalence. *Organizational Research Methods,
22*(3), 678–709.
[doi:10.1177/1094428118761122](https://doi.org/10.1177/1094428118761122)

Nye, C. D., & Drasgow, F. (2011). Effect size indices for analyses of
measurement equivalence: Understanding the practical importance of
differences between groups. *Journal of Applied Psychology, 96*(5),
966–980. [doi:10.1037/a0022955](https://doi.org/10.1037/a0022955)

## See also

[`measurement_invariance`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md)
for the invariance ladder that locates a violation;
[`smd`](https://yelleknek.github.io/DMAR/reference/smd.md) for the
standardized mean difference whose pooling and scale dMACS borrows;
[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) for the
single group measurement model.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Reported measurement equations: the two items share loadings, and the
# second item's intercept is 0.30 higher in the reference group. With no
# loading difference, dMACS is just 0.30 divided by the pooled SD.
dmacs(lambda_reference = c(0.80, 0.75), lambda_focal = c(0.80, 0.75),
      nu_reference = c(2.00, 2.30), nu_focal = c(2.00, 2.00),
      sd_pooled = c(1.20, 1.10), item_names = c("optimism", "worry"))
#>  item     lambda_reference lambda_focal nu_reference nu_focal sd_pooled dmacs
#>  optimism 0.8              0.8          2            2        1.2       0    
#>  worry    0.75             0.75         2.3          2        1.1       0.273

# A partial invariance model for the four spatial tests at the two
# Holzinger and Swineford schools. The anchors are constrained equal; the
# cubes and lozenges tests are freed, so only those two can move.
data(holzinger_swineford)
items <- c("t1_visual_perception", "t2_cubes",
           "t3_paper_form_board", "t4_lozenges")
model <- paste("spatial =~", paste(items, collapse = " + "))
fit <- lavaan::cfa(model, data = holzinger_swineford, group = "school",
                   group.equal = c("loadings", "intercepts"),
                   group.partial = c("spatial =~ t2_cubes", "t2_cubes ~ 1",
                                     "spatial =~ t4_lozenges",
                                     "t4_lozenges ~ 1"))
dmacs(fit)
#>  item                 lambda_reference lambda_focal nu_reference nu_focal
#>  t1_visual_perception 1                1            29.6         29.6    
#>  t2_cubes             0.495            0.52         23.9         24.8    
#>  t3_paper_form_board  0.308            0.308        14.2         14.2    
#>  t4_lozenges          1.28             1.37         19.9         15.8    
#>  sd_pooled dmacs   
#>  7.02      0       
#>  4.7       0.175   
#>  2.83      4.35e-16
#>  8.85      0.461   

# The broom verbs: one row per item, and the group metadata.
generics::tidy(dmacs(fit))
#>                   term     estimate lambda_reference lambda_focal nu_reference
#> 1 t1_visual_perception 0.000000e+00        1.0000000    1.0000000     29.56940
#> 2             t2_cubes 1.754320e-01        0.4953945    0.5199470     23.93590
#> 3  t3_paper_form_board 4.353960e-16        0.3084177    0.3084177     14.21526
#> 4          t4_lozenges 4.606217e-01        1.2796675    1.3720483     19.89744
#>   nu_focal sd_pooled
#> 1 29.56940  7.016213
#> 2 24.75043  4.697740
#> 3 14.21526  2.834122
#> 4 15.83472  8.845986
generics::glance(dmacs(fit))
#>         n_items reference       focal mean_focal sd_focal
#> spatial       4   Pasteur Grant-White 0.09533115 4.444806
```
