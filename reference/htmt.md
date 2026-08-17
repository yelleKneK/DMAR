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
  constructs).

- B:

  Number of bootstrap resamples for the upper confidence bound; `0`
  (default) skips the bootstrap and reports the point estimates only.

- conf_level:

  Confidence level for the one-sided upper bound. Defaults to 0.95.

- seed:

  Optional integer seed for the bootstrap, used locally (the caller's
  random number generator state is restored on exit).

## Value

A `data.frame` (class `dmar_tbl`) with one row per construct pair:
`construct_1`, `construct_2`, `htmt`, and, when `B > 0`, `upper_limit`
(the one-sided `conf_level` bootstrap percentile bound).

## Details

All correlations are Pearson, computed on pairwise-complete
observations. The statistic uses absolute average heterotrait
correlations enter as absolute values, the convention of later
implementations (the 2015 proposal used the plain correlations, which
can cancel when signs mix; Roemer, Schuberth, & Henseler, 2021,
recommend the absolute form for exactly that reason), and values near or
above 1 indicating that the two item sets correlate across constructs
about as strongly as within them. HTMT is a correlation-based screen,
deliberately model-free; the confirmatory companion is the latent
correlation between the two factors (see
[`correction_for_attenuation`](https://yelleknek.github.io/DMAR/reference/correction_for_attenuation.md)
and its factor-model discussion).

The bootstrap, when requested (`B > 0`), resamples the rows of `data`
with replacement `B` times and recomputes every pairwise HTMT on each
resample; the reported `upper_limit` is the `conf_level` empirical
quantile of each pair's bootstrap distribution, a one-sided upper
percentile bound (Efron & Tibshirani, 1993). That bound is the only
interval offered, matching how the validity literature uses HTMT (the
question is whether the ratio credibly exceeds the cutoff); no two-sided
or bias-corrected and accelerated (BCa) variant is provided. A resample
in which some block's average within-construct correlation is not
positive leaves HTMT undefined there; such resamples are dropped, a
single warning reports how many, and the bound is computed from the
resamples that remained (the call stops only when fewer than 100
remain). Bootstrap results vary from run to run; supply `seed` for
reproducibility.

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
[`cfa_k()`](https://yelleknek.github.io/DMAR/reference/cfa_k.md),
[`ci_eigenvalue()`](https://yelleknek.github.io/DMAR/reference/ci_eigenvalue.md),
[`common_method_marker()`](https://yelleknek.github.io/DMAR/reference/common_method_marker.md),
[`common_method_single_factor()`](https://yelleknek.github.io/DMAR/reference/common_method_single_factor.md),
[`dmacs()`](https://yelleknek.github.io/DMAR/reference/dmacs.md),
[`ecvi()`](https://yelleknek.github.io/DMAR/reference/ecvi.md),
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
#>  construct_1 construct_2 htmt 
#>  A           B           0.269

# The broom verbs: one row per construct pair.
generics::tidy(h)
#>   term  estimate
#> 1  A:B 0.2691979
generics::glance(h)
#>   n_terms conf_level B_used
#> 1       1         NA     NA

# The upper confidence bound, which is the quantity the validity
# literature compares against 0.85 or 0.90, comes from a bootstrap
# that recomputes every pairwise ratio on each of B resamples. It is
# shown rather than run; the call is
#   htmt(d, blocks = list(A = c("a1", "a2", "a3"),
#                         B = c("b1", "b2", "b3")),
#        B = 10000, seed = 113)
# and a claim about discriminant validity deserves that bound rather
# than the point estimate alone.
```
