# Bifactor Model Dimensionality and Reliability Indices

Computes the indices that describe how much of the common variance of a
scale's items, and of the variance of its total and subscale scores, a
general factor accounts for beside the group factors of its subscales
(Rodriguez, Reise, and Haviland, 2016): the explained common variance
(ECV), coefficient omega and omega hierarchical (omega_H) for the
general factor, omega hierarchical subscale (omega_HS) for each group
factor, the percentage of uncontaminated correlations (PUC), and the
index those authors call coefficient H, the closed form of Hancock and
Mueller (2001) applied to one factor's loadings, which counts the other
factors' common variance as uncorrelated error and so is not the maximal
reliability of the factor under the fitted model; that maximal
reliability is returned beside the table. The input is a fitted bifactor
model in which one general factor loads on every item and each item
loads on at most one of two or more orthogonal group factors.
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
fits such a model and reports these indices with delta method intervals,
with omega, omega_H, and omega_HS for the sum of the items as scored
rather than of the standardized items (see *Details*).

## Usage

``` r
bifactor_indices(
  fit,
  general = NULL,
  denominator = c("model_implied", "observed")
)
```

## Arguments

- fit:

  A fitted bifactor lavaan model: one general factor on all items plus
  two or more orthogonal group factors, each item on one group factor at
  most. An item on no group factor is on the general factor alone: it
  counts in the general factor's indices and in PUC, and in no subscale.
  An item listed under two or more group factors is refused unless all
  but one of those loadings are held at zero, fixed (`0*item`) or
  constrained (`a*item` with `a == 0`); such a loading describes the
  same model as leaving the item off that factor. The factors must be
  orthogonal (the bifactor specification), and a fit that correlates any
  two of them is refused. Residual covariances among the items are
  allowed and enter the model implied denominators. A fit with fewer
  than two group factors, a multiple-group fit, a multilevel fit, and a
  fit whose optimizer did not converge are refused (see *Details*).

- general:

  Optional name of the general factor. When `NULL` (default) the general
  factor is detected as the one that loads on every item.

- denominator:

  How the total-score variance in the denominators of omega, omega_H,
  and omega_HS is estimated: `"model_implied"` (the default; the
  variance the fitted bifactor model reproduces, as in Rodriguez, Reise,
  and Haviland, 2016, with any residual covariances the model specifies)
  or `"observed"` (the variance of the composite estimated from the
  data, the denominator of robust omega in
  [`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)).
  For a fit to ordered items both describe the sum of the latent
  responses, not of the summed categories. See *Details*.

## Value

A `data.frame` (class `dmar_tbl`) with one row per factor (the general
factor first, then each group factor) and columns `factor`, `ECV`,
`omega`, `omega_H`, `omega_HS`, `PUC`, and `H`. Quantities that do not
apply to a row are `NA` (for example omega_H and PUC on a group factor).
On a group factor's row, ECV is the general factor's share of the common
variance of that factor's items and omega is the reliability of their
total score (see *Details*); H is the index of Rodriguez, Reise, and
Haviland (2016) on every row. The `"improper"` attribute flags a Heywood
solution, the `"maximal_reliability"` attribute holds the maximal
reliability the fitted model implies for each factor (a vector named by
factor, in the order of the rows; `NA` for an improper solution), the
`"denominator"` attribute records the `denominator` setting, and the
`"metric"` attribute records the metric of the items the indices
describe: `"observed"` for continuous items, `"latent_response"` for
ordered items, whose coefficients describe the sum of the latent
responses.

## Details

Let \\\lambda^g_i\\ be item \\i\\'s standardized loading on the general
factor, \\\lambda^s_i\\ its loading on its group factor (0 for an item
on the general factor alone), and \\\theta_i\\ its standardized residual
variance. With orthogonal factors the overall ECV is \\\sum_i
(\lambda^g_i)^2 / \sum_i \[(\lambda^g_i)^2 + (\lambda^s_i)^2\]\\; omega
and omega_H share the total-score variance \\(\sum_i \lambda^g_i)^2 +
\sum_g (\sum\_{i \in g} \lambda^s_i)^2 + \sum_i \theta_i\\ as
denominator, with the general part \\(\sum_i \lambda^g_i)^2\\ in the
numerator of omega_H. Each group factor's omega_HS uses the analogous
numerator and that subscale's own total variance. PUC is one minus the
share of item pairs that fall within the same group factor. Omega_H
against omega says how much of the total score's reliable variance is
due to the general factor; omega_HS against a subscale's omega says how
much of the subscale score's reliable variance is its own, beyond the
total; ECV and PUC describe how close the items come to measuring one
thing. The table reports the numbers; whether a total or subscale score
is worth reporting is the researcher's judgment, and the package draws
no line, as on the page of
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md).

**The metric of omega.** The loadings and residual variances above are
standardized, so omega, omega_H, and omega_HS describe the sum of the
standardized items, as Rodriguez, Reise, and Haviland (2016) write them.
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
reports them for the sum of the items as scored, the coefficients
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
and `semTools::compRelSEM()` compute. The two agree when the items have
equal variances and differ otherwise: in the first example below, whose
items have nearly equal variances, omega_H is 0.67348 here and 0.67305
for the sum of the items as scored, while on the 15 tests of the example
of
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
whose variances range from 8.0 to 153.3, it is 0.577 here and 0.552
there. ECV, *H*, the maximal reliability, and PUC do not differ.

**A group factor's row.** The row of group factor \\s\\ describes its
items, the subscale \\s\\, and its ECV and omega are not the general
row's quantities restricted to a group factor. Its ECV is the general
factor's share of the subscale's common variance, \$\$\frac{\sum\_{i \in
s} (\lambda^g_i)^2} {\sum\_{i \in s} \[(\lambda^g_i)^2 +
(\lambda^s_i)^2\]},\$\$ so one minus it is the group factor's own share
of that variance, and a high value says the subscale's common variance
is mostly general (it goes with a low omega_HS). It is not the group
factor's share of the common variance of all the items, which is
\\\sum\_{i \in s} (\lambda^s_i)^2\\ over the denominator of the overall
ECV. Its omega is the omega of the subscale's total score, \\(\sum\_{i
\in s} \lambda^g_i)^2 + (\sum\_{i \in s} \lambda^s_i)^2\\ over that
score's variance, the reliability of the subscale score with both
factors counted as true score; omega_HS keeps only the group factor's
part of the numerator. Its H is the index described next, computed from
the group factor's loadings. The printed table states what ECV and omega
mean on these rows.

**Coefficient H.** The column `H` is the index Rodriguez, Reise, and
Haviland (2016) report: the closed form of Hancock and Mueller (2001),
\$\$H = \frac{\sum_i \lambda_i^2 / (1 - \lambda_i^2)} {1 + \sum_i
\lambda_i^2 / (1 - \lambda_i^2)},\$\$ over the general loadings for the
general factor and over a group factor's loadings for that factor. It
counts \\1 - \lambda_i^2\\ as item \\i\\'s error variance, and so treats
the variance an item shares with the other factors as error that is
uncorrelated across items. In a bifactor model that variance is
correlated (the items of a subscale share their group factor, and every
item shares the general factor), so the index is not the maximal
reliability of the factor under the fitted model, which is what
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md)
computes for a block of items whose errors are correlated (Gabler &
Raykov, 2017). That maximal reliability is \\\lambda_f' R^{-1}
\lambda_f\\, where \\\lambda_f\\ holds the factor's standardized
loadings (zero for an item off the factor) and \\R\\ is the item
correlation matrix the model reproduces: the squared multiple
correlation of the factor with all of the items, the largest reliability
any weighted sum of them reaches as a measure of that factor, and the
square of the determinacy of its regression factor score in
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md).
It is returned as the `"maximal_reliability"` attribute and printed
under the table, and
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
gives it a confidence interval, as it does every index of the table but
PUC, which depends on the model's structure alone.

The two differ in either direction, and what sets the direction is not
the same on the two kinds of row. On the general factor's row the index
counts the group factors' variance as uncorrelated error, but the items
of a subscale share it. When the loadings on each group factor are
positive and similar, no weighting of the items removes that shared
variance without removing general variance with it, and the index is the
larger: in the first example below it is 0.769 for the general factor,
whose maximal reliability is 0.678, barely above its omega_H of 0.673,
so no weighting of the items reaches 0.769. When an item is nearly a
pure indicator of its group factor, or a group factor's loadings differ
in sign, a weighting can cancel group variance that the index counts as
error, and the maximal reliability is the larger. On a group factor's
row the index counts the general factor's variance, with the unique
variance, as uncorrelated error. When the subscale's own items have
proportional loadings on the two factors (equal loadings, for instance),
no weighting of them alone removes the general factor without removing
the group factor with it, but the optimal weights use every item, and
the items of the other subscales measure the general factor and nothing
of this group factor, so the weights use them to remove it. The better
those items measure the general factor, the closer the maximal
reliability comes to the closed form above with each item's unique
variance \\\theta_i\\ in place of \\1 - (\lambda^s_i)^2\\, which is
above the index. Which of the two prevails depends on the general
loadings, the group loadings, and the number of items outside the
subscale together, not on the sign and similarity of the group loadings
alone. In a population of nine items in three subscales of three, with
every general loading .7 and every group loading .3 (an ECV of .84), the
index is 0.229 for each group factor and its maximal reliability 0.279,
while for the general factor they are 0.896 and 0.865; with general and
group loadings of .5, as in the first example below, the index is the
larger on every row (0.500 against 0.471 for a group factor in the
population); and with two subscales of six items, general loadings .8,
and group loadings .4, the index is the larger for a group factor (0.533
against 0.468), although the general factor is stronger than in the
first population.

The composite that reaches a group factor's maximal reliability is the
factor's regression score from all of the items, and it weights the
items of the other subscales negatively: in the first population the
optimal weights are proportional to 1 for each of the subscale's three
items and -0.405 for each of the other six. No score of the subscale's
own items reaches that reliability. The best weighting of them alone
reaches \\\lambda_s' R\_{ss}^{-1} \lambda_s\\, with \\\lambda_s\\ their
group loadings and \\R\_{ss}\\ their correlations, which is omega_HS
when the items have equal loadings on both factors, and then below the
index as well: 0.125 in the first population, against an index of 0.229
and a maximal reliability of 0.279. So a maximal reliability above the
index on a group factor's row, which is common when the items outside
the subscale measure a strong general factor well (a profile of high
ECV, in which the subscale rows are commonly consulted), describes a
factor score computed from all of the items, not the subscale score; how
much of a subscale score's variance is due to its group factor is
omega_HS. Read the column as the published index, comparable with values
reported beside it, and the attribute as the reliability that optimal
weights of all of the items could reach in the population. The sample
values of both overstate their population values on average, the
attribute most of all for a weakly determined group factor, and a
composite weighted with sample estimates reaches less, as
Aguirre-Urreta, Rönkkö, and McIntosh (2019) found for H of a single
factor;
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
gives the bootstrap estimate of that bias with the intervals.

The H of a factor any of whose standardized loadings is at or beyond one
in absolute value is `NA`, since the index then gives that item an error
variance that is not positive.

Each denominator is the sum of the corresponding block of an item
correlation matrix on the standardized metric of the loadings: all items
for omega and omega_H, a group factor's items for its omega_HS. The
default (`denominator = "model_implied"`) takes the correlation matrix
the fitted model reproduces. For an orthogonal bifactor model without
residual covariances that block sum is the total-score variance written
above, the form given by Rodriguez, Reise, and Haviland (2016); when the
model also specifies residual covariances among the items, the block sum
carries each of them as well, since they are part of the variance of the
composite the model reproduces, although not of its common variance, so
the numerators leave them out. With `denominator = "observed"` the
matrix is instead the item correlation matrix estimated from the data
(lavaan's sample statistics, which under `missing = "ml"` are the full
information estimates). This is the denominator of robust omega
([`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md))
carried to the bifactor coefficients: it estimates the variance of the
composite consistently whether or not the bifactor model reproduces it,
while the numerators remain model based. When the model reproduces the
item correlations the two settings agree, and a gap between them shows
that it does not reproduce the composite variance. The printed table
names the denominator in a line under the table.

The indices assume uncorrelated factors: a correlation between the
general factor and a group factor, or between two group factors, is
common variance that belongs to neither, and the partition behind ECV,
omega_H, and omega_HS has no place for it. A fit that estimates such a
correlation is refused with an error naming the pair.

**A fit with one group factor.** A general factor beside a single group
factor is not a bifactor model with subscales: its indices would divide
the common variance between those two factors alone, every item would
belong to the one subscale, and PUC would be 0. Such a fit is refused,
as
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)
and
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
refuse it, with a message that names the group factor. A general factor
with one further factor can be fitted with
[`lavaan::cfa()`](https://rdrr.io/pkg/lavaan/man/cfa.html) and read from
its own parameter estimates.

**Ordered-categorical items.** For a fit with ordered indicators (WLSMV
or DWLS on polychoric correlations), lavaan's sample statistics are the
polychoric correlations, and the standardized loadings are on the metric
of the continuous latent responses behind the categories. Both settings
then give the reliability of the sum of those latent responses, not of
the summed category scores a user computes, and categorization lowers
reliability, often substantially (Green & Yang, 2009, Table 1). The
printed line says which metric the denominator is on, and the `"metric"`
attribute records it. For a single factor,
[`reliability_omega_categorical`](https://yelleknek.github.io/DMAR/reference/reliability_omega_categorical.md)
gives the reliability of the summed categories; DMAR does not compute
the bifactor indices on that metric.

A bifactor model that is over-parameterized for the data frequently
yields an improper (Heywood) solution, a standardized loading outside
\\\[-1, 1\]\\ or a negative residual variance, whose indices are not
trustworthy. When that happens the indices are still returned but a
warning is issued and the `"improper"` attribute is set to `TRUE`; H is
`NA` for a factor with a loading at or beyond one, as described above,
and the maximal reliability is `NA` for every factor, since it is not
defined for an improper solution.

A fit whose optimizer did not converge
(`lavaan::lavInspect(fit, "converged")` is `FALSE`, as when it reached
its iteration limit) is not a solution, proper or improper: its
estimates mark where the optimizer stopped. It is refused with an error
of class `dmar_not_converged`, as
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
refuses it.

**One group and one level.** The standardized loadings and residual
variances of a fit are those of one group at one level: each group of a
multiple-group fit has its own, standardized by that group's variances,
so the factors have indices in each group, and neither an unweighted
mean nor a mean weighted by group size is the index of any group. A
multiple-group fit is refused, as
[`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
refuses it, with an error of class `dmar_argument_invalid` that names
the groups (as `groups`); the indices of one group are those of the
model fitted to that group's data alone (which drops any equality
constraints across groups). A multilevel fit, whose loadings and
residual variances differ by level, is refused with an error of the same
class; the standardized solution of each level is in
`lavaan::lavInspect(fit, "std")`.

## References

Aguirre-Urreta, M. I., Rönkkö, M., & McIntosh, C. N. (2019). A
cautionary note on the finite sample behavior of maximal reliability.
*Psychological Methods, 24*(2), 236–252.
[doi:10.1037/met0000176](https://doi.org/10.1037/met0000176)

Gabler, S., & Raykov, T. (2017). Evaluation of maximal reliability for
unidimensional measuring instruments with correlated errors. *Structural
Equation Modeling, 24*(1), 104–111.
[doi:10.1080/10705511.2016.1159916](https://doi.org/10.1080/10705511.2016.1159916)

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
*Psychometrika, 74*, 155–167.
[doi:10.1007/s11336-008-9099-3](https://doi.org/10.1007/s11336-008-9099-3)

Hancock, G. R., & Mueller, R. O. (2001). Rethinking construct
reliability within latent variable systems. In R. Cudeck, S. du Toit, &
D. Sörbom (Eds.), *Structural equation modeling: Present and future*
(pp. 195–216). Scientific Software International.

Reise, S. P. (2012). The rediscovery of bifactor measurement models.
*Multivariate Behavioral Research, 47*(5), 667–696.
[doi:10.1080/00273171.2012.715555](https://doi.org/10.1080/00273171.2012.715555)

Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating
bifactor models: Calculating and interpreting statistical indices.
*Psychological Methods, 21*(2), 137–150.
[doi:10.1037/met0000045](https://doi.org/10.1037/met0000045)

## See also

[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
which fits the bifactor model and gives these indices delta method
intervals;
[`ci_bifactor_indices`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md)
for bootstrap confidence intervals for these indices and for the maximal
reliability of each factor;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md),
[`reliability_H`](https://yelleknek.github.io/DMAR/reference/reliability_H.md),
and
[`factor_scores`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
whose regression score determinacy is the square root of the maximal
reliability.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`cfa_1()`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
[`cfa_2()`](https://yelleknek.github.io/DMAR/reference/cfa_2.md),
[`ci_bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/ci_bifactor_indices.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_adjusted_r()`](https://yelleknek.github.io/DMAR/reference/common_method_adjusted_r.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`compare_sem()`](https://yelleknek.github.io/DMAR/reference/compare_sem.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
[`equivalence_rmsea()`](https://yelleknek.github.io/DMAR/reference/equivalence_rmsea.md),
[`fa_bifactor()`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md),
[`fa_confirmatory()`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md),
[`fa_second_order()`](https://yelleknek.github.io/DMAR/reference/fa_second_order.md),
[`factor_congruence()`](https://yelleknek.github.io/DMAR/reference/factor_congruence.md),
[`factor_scores()`](https://yelleknek.github.io/DMAR/reference/factor_scores.md),
[`htmt()`](https://yelleknek.github.io/DMAR/reference/htmt.md),
[`irt_grm()`](https://yelleknek.github.io/DMAR/reference/irt_grm.md),
[`irt_information()`](https://yelleknek.github.io/DMAR/reference/irt_information.md),
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`measurement_invariance()`](https://yelleknek.github.io/DMAR/reference/measurement_invariance.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`sem_exact_fit_bootstrap()`](https://yelleknek.github.io/DMAR/reference/sem_exact_fit_bootstrap.md),
[`sem_fit_measures()`](https://yelleknek.github.io/DMAR/reference/sem_fit_measures.md),
[`sem_information_criteria()`](https://yelleknek.github.io/DMAR/reference/sem_information_criteria.md),
[`sem_residual_fit()`](https://yelleknek.github.io/DMAR/reference/sem_residual_fit.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Nine items: one general factor and three orthogonal group factors.
set.seed(113)
n <- 600
g <- rnorm(n); grp <- list(rnorm(n), rnorm(n), rnorm(n))
X <- vapply(1:9, function(i)
  0.5 * g + 0.5 * grp[[ceiling(i / 3)]] + sqrt(0.5) * rnorm(n), numeric(n))
colnames(X) <- paste0("x", 1:9)
model <- "g  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
          f1 =~ x1 + x2 + x3
          f2 =~ x4 + x5 + x6
          f3 =~ x7 + x8 + x9"
fit <- lavaan::cfa(model, data = as.data.frame(X),
                   orthogonal = TRUE, std.lv = TRUE)
bifactor_indices(fit)
#>  factor ECV   omega omega_H omega_HS PUC   H    
#>  g      0.547 0.858 0.673   <NA>     0.750 0.769
#>  f1     0.609 0.727 <NA>    0.280    <NA>  0.412
#>  f2     0.543 0.756 <NA>    0.345    <NA>  0.477
#>  f3     0.490 0.743 <NA>    0.377    <NA>  0.506
#> 
#> Denominator: composite variance implied by the fitted model
#> 
#> For a group factor, ECV is the general factor's share of the common variance of that factor's items (one minus it is the group factor's share), and omega is the reliability of their total score
#> H is the index of Rodriguez, Reise, and Haviland (2016), which counts the other factors' common variance as uncorrelated error; the maximal reliability the fitted model implies is g 0.678, f1 0.396, f2 0.456, f3 0.494 (attribute "maximal_reliability"; ci_bifactor_indices() gives each an interval)

# The same indices with the observed composite variances in the
# denominators; when the model reproduces the item correlations the
# two settings agree.
bifactor_indices(fit, denominator = "observed")
#>  factor ECV   omega omega_H omega_HS PUC   H    
#>  g      0.547 0.858 0.673   <NA>     0.750 0.769
#>  f1     0.609 0.727 <NA>    0.280    <NA>  0.412
#>  f2     0.543 0.756 <NA>    0.345    <NA>  0.477
#>  f3     0.490 0.743 <NA>    0.377    <NA>  0.506
#> 
#> Denominator: observed composite variance
#> 
#> For a group factor, ECV is the general factor's share of the common variance of that factor's items (one minus it is the group factor's share), and omega is the reliability of their total score
#> H is the index of Rodriguez, Reise, and Haviland (2016), which counts the other factors' common variance as uncorrelated error; the maximal reliability the fitted model implies is g 0.678, f1 0.396, f2 0.456, f3 0.494 (attribute "maximal_reliability"; ci_bifactor_indices() gives each an interval)
```
