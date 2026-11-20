# Heterotrait-Monotrait Ratio of Correlations (HTMT)

Computes the HTMT discriminant-validity index of Henseler, Ringle, and
Sarstedt (2015) for every pair of constructs: the average correlation
between items of *different* constructs, divided by the geometric mean
of the average correlations among items *within* each construct. Two
constructs whose HTMT approaches 1 are empirically indistinguishable
however cleanly the model draws them; the customary red flags are 0.85
(strict) or 0.90 (liberal). An optional bootstrap gives a one-sided
upper confidence bound, the quantity actually compared against the
cutoff in the validity literature.

## Usage

``` r
htmt(data, blocks, B = 0, conf_level = 0.95, seed = NULL)
```

## Arguments

- data:

  A `data.frame` of item responses.

- blocks:

  Named list of character vectors: each element names a construct and
  gives its item columns (two or more per construct, two or more
  constructs). The blocks are disjoint: each item is listed once, in one
  block, and each construct has a name of its own (see *Disjoint blocks*
  under *Details*).

- B:

  Number of bootstrap resamples for the upper confidence bound, a single
  whole number; `0` (default) skips the bootstrap and reports the point
  estimates only. When positive it must be at least 200 at the default
  `conf_level` (100 at 0.90, 1000 at 0.99): the bound is one quantile of
  the resamples, reported only when at least 10 of them lie beyond it
  and never from fewer than 100. A smaller positive `B` is refused
  before anything is computed (class `dmar_count_invalid`). See
  *Bootstrap replications* under *Details* and
  [`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md).

- conf_level:

  Confidence level for the one-sided upper bound. Defaults to 0.95.

- seed:

  Optional integer seed for the bootstrap, used locally (the caller's
  random number generator state is restored on exit).

## Value

A `data.frame` (class `dmar_tbl`) with one row per construct pair and
the same four columns whatever `B`: `construct_1`, `construct_2`,
`htmt`, and `ci_upper`, the one-sided `conf_level` bootstrap percentile
bound, which is `NA` when `B = 0` (no bootstrap was run). Both are `NA`
for a pair whose ratio is undefined on the data (see Details). With
`B > 0` the table also carries the attributes `conf_level`, `B`
(requested), `B_used` (the number of resamples the bound was computed
from), `n_dropped`, `B_required`, and `bootstrap_status` (see *Bootstrap
replications* under Details).

## Details

All correlations are Pearson, computed on pairwise-complete
observations, and every one of them enters its average as an absolute
value, the monotrait correlations among the items of a block as well as
the heterotrait correlations between blocks. This is the ratio of
arithmetic means that semTools' `htmt()` computes with `htmt2 = FALSE`
and its default `absolute = TRUE`; on incomplete data the two agree only
when semTools is given the same pairwise-complete correlations (through
its `sample.cov` argument), since by default it deletes incomplete rows
listwise, and its `missing = "pairwise"` option estimates the
correlations through lavaan, whose values differ from
`stats::cor(use = "pairwise.complete.obs")`. Absolute values make the
ratio indifferent to the direction in which an item is keyed, since
reversing an item changes the signs of its correlations and not their
sizes. With signed correlations, a reverse keyed item that was not
recoded lowers its block's average within-construct correlation, which
raises the ratio, while its negative correlations with the items of the
other block partly cancel the positive ones in the numerator, so the
ratio would depend on how the items happen to be scored. Values near or
above 1 indicate that the two item sets correlate across constructs
about as strongly as within them. HTMT is a correlation-based screen,
deliberately model-free; the confirmatory companion is the latent
correlation between the two factors (see
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md)
and its factor-model discussion, and the `phi` rows of
[`fa_confirmatory`](https://yelleknek.github.io/DMAR/reference/fa_confirmatory.md)
with `output = "measurement"`).

**Disjoint blocks.** The ratio compares two sets of items, so the blocks
are disjoint sets of distinct items under distinct names. An item listed
twice in one block would enter that block's average with its correlation
of 1 with itself, and its other correlations would be counted twice; an
item placed in two blocks would enter the average correlation between
them with its correlation of 1 with itself; and two blocks under one
name could not be told apart. Each of these, and an item whose name
belongs to two or more columns of `data`, is refused with an error of
class `dmar_htmt_blocks_invalid` (see
[`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md))
before any correlation is computed or any resample drawn, whatever `B`.
The condition's `reason` element names the check that failed
(`"construct_names"`, `"repeated"`, `"shared"`, or `"data_names"`), and
its `blocks` and `items` elements the names at issue. A ratio for item
sets that overlap would be a different quantity, which this function
does not compute.

A pair's ratio is undefined when a correlation it averages cannot be
computed: an item with no variance, or two items never observed
together, as in a planned missing data design whose forms give the items
of two blocks to different respondents. That pair's `htmt` (and
`ci_upper`) is `NA`, with a warning that names it, and the other pairs
are reported; the call stops only when no pair is defined.

The bootstrap, when requested (`B > 0`), resamples the rows of `data`
with replacement `B` times and recomputes every pairwise HTMT on each
resample; the reported `ci_upper` is the `conf_level` empirical quantile
of each pair's bootstrap distribution, a one-sided upper percentile
bound (Efron & Tibshirani, 1993). That bound is the only interval
offered, matching how the validity literature uses HTMT (the question is
whether the ratio credibly exceeds the cutoff); no two-sided or
bias-corrected and accelerated (BCa) variant is provided. The ratio of
each resample is the one reported, with every correlation an absolute
value, so a resample that draws a negative correlation within a block
enters with its size. Bootstrap results vary from run to run; supply
`seed` for reproducibility.

**Bootstrap replications.** The bound of each pair is the `conf_level`
quantile of its `B` resampled ratios, so it is run only when `B` puts at
least 10 resamples beyond it (200 at the 95 percent level). A
replication in which the HTMT of some pair defined on the data is
undefined (an item that is constant in the resample, or two items never
observed together) is dropped and counted, never replaced, and the bound
is computed from the rest. The result records `B` (requested), `B_used`
(the number the bound rests on), `n_dropped`, and `B_required` (the
fewest replications a request may ask for; that count is held to `B`,
not to `B_used`) as attributes, with `bootstrap_status` (`"complete"`,
`"dropped"`, `"exploratory"`, `"insufficient_count"`, or
`"insufficient_share"`). One warning of class
`dmar_bootstrap_replications_dropped` states the counts when any
replication was dropped, with a caution when the dropped replications
are at least the 5 percent of `B` that lies beyond one limit at the 95
percent level, because they could then decide by themselves where that
limit falls. When more than half were dropped, the bound of every pair
is `NA`, the estimate is returned, and the warning has class
`dmar_bootstrap_insufficient`.
`options(dmar.bootstrap_exploratory = TRUE)` lets a smaller `B` (2 or
more) run for a first look or a demonstration; every such call gives a
warning of class `dmar_bootstrap_exploratory` and the result is not one
to report. See
[`dmar_bootstrap`](https://yelleknek.github.io/DMAR/reference/dmar_bootstrap.md)
for the rule and its reasons.

## References

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the
bootstrap*. New York, NY: Chapman & Hall/CRC.

Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for
assessing discriminant validity in variance-based structural equation
modeling. *Journal of the Academy of Marketing Science, 43*(1), 115–135.
[doi:10.1007/s11747-014-0403-8](https://doi.org/10.1007/s11747-014-0403-8)

## See also

[`average_variance_extracted`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md)
for the convergent side of the validity ledger;
[`reliability_omega`](https://yelleknek.github.io/DMAR/reference/reliability_omega.md)
for the composite reliability of each block;
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md)
for the latent correlation route.

Other multivariate and latent variable methods:
[`average_variance_extracted()`](https://yelleknek.github.io/DMAR/reference/average_variance_extracted.md),
[`bifactor_indices()`](https://yelleknek.github.io/DMAR/reference/bifactor_indices.md),
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
# Two clean constructs and their six items.
set.seed(113)
n <- 300
f1 <- rnorm(n); f2 <- 0.3 * f1 + sqrt(1 - 0.09) * rnorm(n)
d <- data.frame(
  a1 = .8 * f1 + rnorm(n, 0, .6), a2 = .7 * f1 + rnorm(n, 0, .7),
  a3 = .6 * f1 + rnorm(n, 0, .8),
  b1 = .8 * f2 + rnorm(n, 0, .6), b2 = .7 * f2 + rnorm(n, 0, .7),
  b3 = .6 * f2 + rnorm(n, 0, .8))
h <- htmt(d, blocks = list(A = c("a1", "a2", "a3"),
                           B = c("b1", "b2", "b3")))
h
#>  construct_1 construct_2 htmt  ci_upper
#>  A           B           0.269 <NA>    

# The broom verbs: one row per construct pair.
generics::tidy(h)
#>   term  estimate ci_upper
#> 1  A:B 0.2691979       NA
generics::glance(h)
#>   n_terms conf_level B_used
#> 1       1         NA     NA

# The upper confidence bound, which is the quantity the validity
# literature compares against 0.85 or 0.90, comes from a bootstrap
# that recomputes every pairwise ratio on each of B resamples and
# fills the ci_upper column, NA above. B = 500 keeps the example
# quick; a claim about discriminant validity deserves the bound at
# B = 10000 rather than the point estimate alone.
htmt(d, blocks = list(A = c("a1", "a2", "a3"),
                      B = c("b1", "b2", "b3")),
     B = 500, seed = 113)
#>  construct_1 construct_2 htmt  ci_upper
#>  A           B           0.269 0.378   
#> 
#> Confidence level: 95%
#> 
#> Bootstrap: 500 replications, 500 used (0 dropped; at least 200 must be requested)
```
