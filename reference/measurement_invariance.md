# Measurement Invariance Across Groups

Fits the standard ladder of multiple group invariance models for a
measurement model of any number of factors and reports the comparison
table researchers actually use: configural invariance (same pattern, all
parameters free by group), metric (equal loadings; required before
comparing relations involving the factor), scalar (equal loadings and
intercepts; required before comparing factor means), and strict (equal
residual variances as well; required before comparing observed-score
variances). With ordered indicators a thresholds rung comes first,
because thresholds rather than intercepts carry the location information
(Wu & Estabrook, 2016); with dichotomous indicators the two are not
separately identified, so that rung is folded into metric. Each rung is
tested against the previous with a likelihood ratio test, scaled when
the estimator in force carries a robust test, which is what declaring
ordered items arranges, and the practical-fit changes (delta CFI, delta
RMSEA) are reported alongside, since with large samples the chi square
will flag trivial differences (Cheung & Rensvold, 2002, suggest delta
CFI of about -.01 as a red flag). The fitted lavaan objects come back as
an attribute, so score tests and partial invariance refits do not
require refitting the ladder. Requires lavaan.

## Usage

``` r
measurement_invariance(
  data,
  model = NULL,
  group,
  items = NULL,
  levels = NULL,
  ordered = NULL,
  estimator = "ML",
  missing = "listwise",
  group_partial = NULL,
  parameterization = c("delta", "theta"),
  ...
)
```

## Arguments

- data:

  A `data.frame` with the items and the grouping variable.

- model:

  The measurement model, given either as lavaan model syntax (a single
  string, or a character vector of lines that is collapsed with
  newlines, such as
  `c("visual =~ x1 + x2 + x3", "verbal =~ x4 + x5 + x6")`) or as a named
  list mapping each factor name to a character vector of its item names,
  from which the syntax is built. Exactly one of `model` and `items`
  must be supplied.

- group:

  Single character string naming the grouping column (two or more
  groups).

- items:

  Character vector (three or more) naming the indicator columns of a
  single factor. A convenience for the one-factor (congeneric) case of
  [`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md),
  equivalent to `model = "f =~ item_1 + item_2 + ..."`. Exactly one of
  `model` and `items` must be supplied.

- levels:

  Which rungs of the ladder to fit, in order: any leading subset of the
  ladder in force. With continuous indicators the ladder is
  `c("configural", "metric", "scalar", "strict")`; with ordered
  indicators it is
  `c("configural", "thresholds", "metric", "scalar", "strict")`; when
  every indicator is ordered and dichotomous it is
  `c("configural", "metric", "scalar")`, because neither the thresholds
  nor the strict rung is identified for a two-category item (see
  `ordered` below). `NULL` (the default) fits the whole ladder in force.

- ordered:

  Ordered-categorical (including binary) items: `TRUE` for every
  indicator in the model, or a character vector naming the ordered ones.
  `NULL` (default) treats all indicators as continuous. Declaring
  ordered items switches the ladder and, unless the estimator already
  carries the mean and variance adjusted test, switches the estimator
  (see `estimator` below). The number of observed categories per
  declared item is counted before the ladder is built, because a
  dichotomous item has a single threshold that is not separately
  identified from the intercept of its underlying response (Millsap &
  Yun-Tein, 2004; Wu & Estabrook, 2016), and its residual variance is
  not a free parameter under either parameterization. When every
  declared item is dichotomous the thresholds rung is folded into
  metric, and the strict rung is dropped as well unless a continuous
  indicator is left to carry a residual variance; when only some
  declared items are dichotomous, both rungs are kept and a message
  names those items, whose thresholds and residual variances the rungs
  leave untested.

- estimator:

  Estimator passed to lavaan. Defaults to `"ML"`; use `"MLR"` for robust
  corrections, in which case the likelihood ratio tests use the scaled
  difference. When `ordered` is supplied, a maximum likelihood or
  generalized least squares (the `"GLS"` label; Browne, 1974) estimator
  is replaced by `"WLSMV"`, and `"DWLS"` and `"ULS"` are replaced by
  `"WLSMV"` and `"ULSMV"`. That second substitution changes only the
  test: in lavaan `"WLSMV"` is `"DWLS"` and `"ULSMV"` is `"ULS"` with
  the mean and variance adjusted statistic, which is what makes the
  rung-to-rung difference test asymptotically valid for ordered data.
  The one categorical estimator left alone is `"WLS"`, the full weight
  matrix estimator of Browne (1984) and Muthén (1984), whose statistic
  is asymptotically chi square already and whose weight matrix a
  substitution would silently change.

- missing:

  Missing data handling passed to lavaan: one of `"listwise"` (the
  default), `"ml"` or `"fiml"` (full information maximum likelihood, the
  reason to reach for this argument with continuous items and incomplete
  cases), or `"pairwise"`. Full information maximum likelihood is not
  available with the categorical estimator; asking for both is an error.

- group_partial:

  Character vector of parameters to leave free across groups at every
  rung, passed to lavaan's `group.partial`. Either user-supplied
  parameter labels or lavaan parameter specifications such as
  `"visual =~ x2"` or `"x2 ~1"`. This is the partial invariance case of
  Byrne, Shavelson, and Muthén (1989). `NULL` (default) constrains
  everything the rung calls for.

- parameterization:

  Identification of ordered indicators, passed to lavaan: `"delta"` (the
  default, lavaan's) or `"theta"`. Residual variances are free
  parameters only under `"theta"`, so a strict rung requested with
  ordered items switches to `"theta"` with a message. Ignored when no
  item is ordered.

- ...:

  Further arguments passed to
  [`cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html) at every rung (for
  example `std.lv`, `orthogonal`, `cluster`).

## Value

A tidy wide `data.frame` (class `dmar_tbl`) with one row per fitted
level and columns `level` (label), `chi_square`, `df`, `p_chi_square`
(exact-fit test), `cfi`, `rmsea`, and, from the second row on, the
step-comparison columns `delta_chi_square`, `delta_df`, `p_value` (the
likelihood ratio test against the previous rung), `delta_cfi`, and
`delta_rmsea`. `delta_chi_square` and `p_value` are `NA` wherever
`delta_df` is zero, since a difference test on zero degrees of freedom
tests nothing. Attributes carry the non-numeric information: `"fits"`,
the named list of fitted lavaan objects, one per level; `"estimator"`,
the estimator actually used; `"ordered"`, `TRUE` when any indicator was
declared ordered; `"test"`, naming the chi square difference test used
(`NA` when a single rung was fit and nothing was compared);
`"fit_indices"`, `"standard"` or `"robust"` according to which version
of the fit indices is tabled; and `"model"`, the lavaan syntax that was
fitted.

## Details

The models are nested by construction, each adding equality constraints
across groups to the previous. With continuous indicators the constraint
sets are none (beyond the configuration), then
`group.equal = "loadings"`, then `c("loadings", "intercepts")`, then
`c("loadings", "intercepts", "residuals")`.

The ordered ladder differs, and the difference is substantive rather
than cosmetic. For an ordered indicator the observed response is a
coarsening of an underlying continuous response at a set of thresholds,
and it is the thresholds, not an intercept, that locate the item on the
latent scale. Constraining loadings while leaving thresholds free across
groups therefore does not deliver what metric invariance is supposed to
deliver, and the accepted sequence (Millsap & Yun-Tein, 2004; Wu &
Estabrook, 2016) constrains thresholds first: `"thresholds"`, then
`c("thresholds", "loadings")` for metric, then
`c("thresholds", "loadings", "intercepts")` for scalar, then adding
`"residuals"` for strict. The intercept rung is not vacuous even when
every indicator is ordered: once thresholds and loadings are
constrained, lavaan frees the underlying-response intercepts in the
non-reference groups, and the scalar rung is what returns them to zero
and lets the latent means be estimated instead. Residual variances,
however, are free parameters only under the theta parameterization, so
`parameterization` switches to `"theta"` when a strict rung is requested
with ordered items.

Dichotomous items are the exception the ordered ladder has to make room
for. A two-category item contributes one threshold, and that threshold,
the intercept of the underlying response, and its residual variance are
not separately identified: the data give one proportion per group per
item, which pins down a single standardized location and nothing else
(Millsap & Yun-Tein, 2004; Wu & Estabrook, 2016). Constraining
thresholds alone across groups is then a reparameterization rather than
a restriction, since lavaan frees the underlying-response intercepts by
exactly as many parameters as the constraint removes, and the residual
variances stay fixed at one in every group whichever parameterization is
in force. So the function counts the observed categories of each
declared ordered indicator before building the ladder, and a message
explains what the count implies. When all of them are dichotomous the
thresholds rung is folded into metric, which constrains thresholds and
loadings together; the strict rung goes too when every indicator in the
model was declared ordered, leaving configural, metric, scalar. A
continuous indicator alongside dichotomous ones keeps the strict rung,
since its residual variance is a free parameter. When only some declared
items are dichotomous the full ordered ladder is kept, since the
polytomous items still carry testable thresholds and residual variances,
and the message names the dichotomous items so their contribution to
those two rungs is not overread. Every rung that is dropped is one that
would have cost zero degrees of freedom.

When the estimator carries a robust test, the difference between two chi
square statistics is not itself chi square distributed.
[`lavTestLRT`](https://rdrr.io/pkg/lavaan/man/lavTestLRT.html) then
returns the scaled difference test (the Satorra-Bentler or Satorra
correction, chosen by lavaan to match the test in force), and that is
what `delta_chi_square` and `p_value` report; the `"test"` attribute
names the test used. In that case the model-level `chi_square`,
`p_chi_square`, `cfi`, and `rmsea` columns are lavaan's scaled and
robust versions, so `delta_chi_square` will not equal the difference of
consecutive `chi_square` entries. That is a property of scaled
difference testing, not an inconsistency.

Which estimators carry such a test is worth being precise about, because
a naive chi square difference on ordered data is not asymptotically
valid. `"MLM"`, `"MLMV"`, `"MLR"`, `"WLSM"`, `"WLSMV"`, `"ULSM"`, and
`"ULSMV"` carry one. `"DWLS"` and `"ULS"` do not, which is why declaring
ordered items promotes them to `"WLSMV"` and `"ULSMV"`: the discrepancy
function and hence the parameter estimates are untouched, and only the
statistic changes. The remaining case is `"WLS"`, the full weight matrix
estimator, whose statistic is asymptotically chi square under the theory
of Browne (1984) and Muthén (1984) and so needs no correction; the
`"test"` attribute reports the standard difference test there. The
sample size that theory asks for is large, which is why full weighted
least squares is rarely the right choice in practice, but that is a
separate matter from whether the difference test is the right one.

A rung whose constraints cost no degrees of freedom tests nothing, so
`delta_chi_square` and `p_value` are `NA` when `delta_df` is zero rather
than reporting a difference in chi square that is numerical noise and
can come out negative. The situation is reported in a message. It arises
with dichotomous items (the ladder above avoids the two rungs where it
is structural) and, for instance, with `group_partial` specifications
that free every parameter a rung would constrain.

Failure at a rung does not end the conversation: partial invariance
(freeing the offending parameter) is the usual next step, for which
`group_partial` refits the whole ladder with named parameters free, and
the `"fits"` attribute gives the fitted objects to
[`lavTestScore`](https://rdrr.io/pkg/lavaan/man/lavTestScore.html) for a
score test of which constraint is doing the damage. This function
deliberately reports the standard ladder rather than automating
modification searches.

## References

Browne, M. W. (1974). Generalized least squares estimators in the
analysis of covariance structures. *South African Statistical Journal,
8*, 1–24.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. *British Journal of Mathematical and
Statistical Psychology, 37*(1), 62–83.

Byrne, B. M., Shavelson, R. J., & Muthén, B. (1989). Testing for the
equivalence of factor covariance and mean structures: The issue of
partial measurement invariance. *Psychological Bulletin, 105*(3),
456–466.

Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit
indexes for testing measurement invariance. *Structural Equation
Modeling, 9*(2), 233–255.
[doi:10.1207/S15328007SEM0902_5](https://doi.org/10.1207/S15328007SEM0902_5)

Meredith, W. (1993). Measurement invariance, factor analysis and
factorial invariance. *Psychometrika, 58*(4), 525–543.
[doi:10.1007/BF02294825](https://doi.org/10.1007/BF02294825)

Millsap, R. E. (2011). *Statistical approaches to measurement
invariance*. Routledge.

Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance in
ordered-categorical measures. *Multivariate Behavioral Research, 39*(3),
479–515.
[doi:10.1207/s15327906mbr3903_4](https://doi.org/10.1207/s15327906mbr3903_4)

Muthén, B. (1984). A general structural equation model with dichotomous,
ordered categorical, and continuous latent variable indicators.
*Psychometrika, 49*(1), 115–132.

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square
test statistic for moment structure analysis. *Psychometrika, 66*(4),
507–514. [doi:10.1007/BF02296192](https://doi.org/10.1007/BF02296192)

Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
analysis models of different levels of invariance for ordered
categorical outcomes. *Psychometrika, 81*(4), 1014–1045.
[doi:10.1007/s11336-016-9506-0](https://doi.org/10.1007/s11336-016-9506-0)

## See also

[`cfa_1`](https://yelleknek.github.io/DMAR/reference/cfa_1.md) and
[`cfa_k`](https://yelleknek.github.io/DMAR/reference/cfa_k.md) for the
single-group measurement models;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for the reliability of the composite the model justifies;
[`compare_cov_structures`](https://yelleknek.github.io/DMAR/reference/compare_cov_structures.md)
for covariance-structure comparisons outside the factor model;
[`lavTestScore`](https://rdrr.io/pkg/lavaan/man/lavTestScore.html) for
the score test that localizes a failed rung.

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
[`measurement_alignment()`](https://yelleknek.github.io/DMAR/reference/measurement_alignment.md),
[`procrustes_phi()`](https://yelleknek.github.io/DMAR/reference/procrustes_phi.md),
[`simple_structure()`](https://yelleknek.github.io/DMAR/reference/simple_structure.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Do the two schools measure the verbal and the reasoning construct
# in the same way? (Holzinger & Swineford, bundled.) The measurement
# model is a named list of factors; for a single factor, name its
# indicators with items = instead.
data(holzinger_swineford)
hs_factors <- list(
  verbal    = c("t6_paragraph_comprehension", "t7_sentence",
                "t9_word_meaning"),
  deduction = c("t20_deduction", "t22_problem_reasoning",
                "t23_series_completion"))

# The bottom two rungs. Configural invariance asks whether the same
# pattern of loadings holds at both schools; metric adds the
# constraint that the loadings are equal across schools, and it is
# the rung that has to hold before a relation involving one of these
# factors is compared across them. Naming those two in levels = stops
# the ladder there, which is what keeps this example quick.
mi <- measurement_invariance(holzinger_swineford, hs_factors,
                             group = "school",
                             levels = c("configural", "metric"))
mi
#>  level      chi_square df p_chi_square cfi   rmsea  delta_chi_square delta_df
#>  configural 22.5       16 0.1276       0.992 0.052  <NA>             <NA>    
#>  metric     36.3       20 0.0140       0.979 0.0737 13.8             4       
#>  p_value delta_cfi delta_rmsea
#>  <NA>    <NA>      <NA>       
#>  0.0078  -0.0127   0.0217     

# The fitted models travel with the table, so localizing a failed rung
# costs no refitting.
names(attr(mi, "fits"))
#> [1] "configural" "metric"    

# The broom verbs: one row per rung of the ladder, and the model-level
# summary (estimator, test flavor, fit index flavor).
generics::tidy(mi)
#>         term chi_square df p_chi_square       cfi      rmsea delta_chi_square
#> 1 configural   22.50471 16   0.12762979 0.9915883 0.05197395               NA
#> 2     metric   36.34593 20   0.01400184 0.9788619 0.07369221         13.84122
#>   delta_df    p_value  delta_cfi delta_rmsea
#> 1       NA         NA         NA          NA
#> 2        4 0.00781944 -0.0127264  0.02171826
generics::glance(mi)
#>   n_levels estimator ordered                                test fit_indices
#> 1        2        ML   FALSE standard chi square difference test    standard

# Each of the calls below refits the ladder from the bottom, so they
# are shown here rather than run. Leaving levels = at its default
# fits the whole ladder, configural through metric, scalar, and
# strict:
# measurement_invariance(holzinger_swineford, hs_factors,
#                        group = "school")
#
# Partial invariance frees one loading across the schools at every
# rung (Byrne, Shavelson, & Muthén, 1989), so each constrained rung
# loses one degree of freedom relative to full invariance:
# measurement_invariance(holzinger_swineford, hs_factors,
#                        group = "school",
#                        group_partial = "verbal =~ t7_sentence")
#
# With the indicators coded as ordered categories, ordered = TRUE puts
# the thresholds rung first, because thresholds rather than intercepts
# carry the location information there, and makes the rung-to-rung
# tests the scaled difference tests:
# hs_ordered <- holzinger_swineford
# for (item in unlist(hs_factors, use.names = FALSE)) {
#   hs_ordered[[item]] <- as.integer(cut(
#     holzinger_swineford[[item]],
#     breaks = quantile(holzinger_swineford[[item]],
#                       c(0, .25, .5, .75, 1)),
#     include.lowest = TRUE))
# }
# mi_ordered <- measurement_invariance(hs_ordered, hs_factors,
#                                      group = "school",
#                                      ordered = TRUE)
# mi_ordered
# attr(mi_ordered, "test")
#
# The measurement model can also be given as lavaan syntax, and
# missing = "fiml" fits the ladder by full information maximum
# likelihood when continuous items have incomplete cases:
# hs_missing <- holzinger_swineford
# set.seed(113)
# hs_missing$t7_sentence[sample(nrow(hs_missing), 20)] <- NA
# measurement_invariance(
#   hs_missing,
#   model = "verbal =~ t6_paragraph_comprehension + t7_sentence +
#            t9_word_meaning",
#   group = "school", missing = "fiml")
```
