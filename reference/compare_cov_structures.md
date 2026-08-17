# Likelihood-Ratio Comparison of Covariance Structures

Fits a long-format within-subjects regression under a menu of
variance-covariance structures, from independence through the
unstructured form, and returns a comparison table of log-likelihood,
AIC, BIC, and pairwise likelihood-ratio tests against the most general
structure (UN). Wraps [`gls`](https://rdrr.io/pkg/nlme/man/gls.html).

## Usage

``` r
compare_cov_structures(
  data,
  outcome,
  subject,
  time,
  fixed_effects = NULL,
  structures = c("IND", "CS", "CSH", "AR1", "ARH1", "TOEP", "TOEPH", "UN")
)
```

## Arguments

- data:

  Long-format `data.frame` with one row per subject-by-condition
  observation.

- outcome:

  Character name of the response column.

- subject:

  Character name of the subject-id column.

- time:

  Character name of the time / within-subjects factor column.

- fixed_effects:

  Right-hand-side formula for the fixed effects (default: `~ time`).

- structures:

  Character vector of structures to fit. Any subset of
  `c("IND", "CS", "CSH", "AR1", "ARH1", "TOEP", "TOEPH", "UN")`
  (default: all eight). Matching is case insensitive, so lowercase
  aliases such as `"cs"`, `"ar1"`, `"csh"`, `"arh1"`, `"toep"`, and
  `"un"` are accepted and normalized to their canonical uppercase
  labels.

## Value

A `data.frame` with one row per structure. Columns: `structure`,
`log_lik`, `AIC`, `BIC`, `n_par`, `LRT_vs_UN_chisq`, `LRT_vs_UN_df`,
`LRT_vs_UN_p`.

## Details

**Structures.** Every structure below is nested in UN, so the
likelihood-ratio test against UN is well defined for each.

- `IND`: independent observations within subject (`correlation = NULL`
  in `gls`). Provided as a baseline.

- `CS`: compound symmetry, a constant correlation and a single variance
  across time points:
  [`nlme::corCompSymm()`](https://rdrr.io/pkg/nlme/man/corCompSymm.html).

- `CSH`: heterogeneous compound symmetry, a constant correlation with a
  separate variance at each time point:
  [`nlme::corCompSymm()`](https://rdrr.io/pkg/nlme/man/corCompSymm.html)
  with [`nlme::varIdent()`](https://rdrr.io/pkg/nlme/man/varIdent.html).

- `AR1`: first-order autoregressive correlation with a single variance:
  [`nlme::corAR1()`](https://rdrr.io/pkg/nlme/man/corAR1.html).

- `ARH1`: heterogeneous first-order autoregressive correlation with a
  separate variance at each time point:
  [`nlme::corAR1()`](https://rdrr.io/pkg/nlme/man/corAR1.html) with
  [`nlme::varIdent()`](https://rdrr.io/pkg/nlme/man/varIdent.html).

- `TOEP`: Toeplitz (banded), a separate correlation at each lag with a
  single variance:
  [`nlme::corARMA()`](https://rdrr.io/pkg/nlme/man/corARMA.html) with
  autoregressive order one less than the number of time points and no
  moving-average term.

- `TOEPH`: heterogeneous Toeplitz, the Toeplitz correlation with a
  separate variance at each time point:
  [`nlme::corARMA()`](https://rdrr.io/pkg/nlme/man/corARMA.html) with
  [`nlme::varIdent()`](https://rdrr.io/pkg/nlme/man/varIdent.html).

- `UN`: unstructured, every variance and covariance free:
  [`nlme::corSymm()`](https://rdrr.io/pkg/nlme/man/corSymm.html) with
  [`nlme::varIdent()`](https://rdrr.io/pkg/nlme/man/varIdent.html).

**LRT.** Each restricted structure is compared against UN by the
likelihood-ratio test. Both fits are re-estimated under ML (not REML)
for the LRT, following `nlme` convention. The chi square statistic is
\\-2 (\log L\_{\mathrm{restricted}} - \log L\_{\mathrm{UN}})\\ on
degrees of freedom equal to the difference in parameter count.

**Caveats.** The likelihood-ratio test against UN is valid because each
listed structure is a restriction of UN. Two structures that are not
nested in each other (for example CS and AR(1)) should be compared by
AIC or BIC rather than by an LRT.

## References

Littell, R. C., Milliken, G. A., Stroup, W. W., Wolfinger, R. D., &
Schabenberger, O. (2006). *SAS for mixed models* (2nd ed.). SAS
Institute.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 15.)

Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-effects models in S and
S-PLUS*. Springer.

## See also

[`gls`](https://rdrr.io/pkg/nlme/man/gls.html),
[`logLik`](https://rdrr.io/r/stats/logLik.html),
[`anova`](https://rdrr.io/r/stats/anova.html)

Other hypothesis tests:
[`adjusted_means()`](https://yelleknek.github.io/DMAR/reference/adjusted_means.md),
[`ancova()`](https://yelleknek.github.io/DMAR/reference/ancova.md),
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
[`ci_dunnett()`](https://yelleknek.github.io/DMAR/reference/ci_dunnett.md),
[`ci_scheffe()`](https://yelleknek.github.io/DMAR/reference/ci_scheffe.md),
[`ci_tukey_kramer()`](https://yelleknek.github.io/DMAR/reference/ci_tukey_kramer.md),
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md),
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md),
[`equivalence_r()`](https://yelleknek.github.io/DMAR/reference/equivalence_r.md),
[`equivalence_smd()`](https://yelleknek.github.io/DMAR/reference/equivalence_smd.md),
[`factorial_anova()`](https://yelleknek.github.io/DMAR/reference/factorial_anova.md),
[`manova_split_plot()`](https://yelleknek.github.io/DMAR/reference/manova_split_plot.md),
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md),
[`mixed_anova()`](https://yelleknek.github.io/DMAR/reference/mixed_anova.md),
[`obrien_test()`](https://yelleknek.github.io/DMAR/reference/obrien_test.md),
[`pairwise_within()`](https://yelleknek.github.io/DMAR/reference/pairwise_within.md),
[`randomization_test()`](https://yelleknek.github.io/DMAR/reference/randomization_test.md),
[`randomization_test_paired()`](https://yelleknek.github.io/DMAR/reference/randomization_test_paired.md),
[`regions_of_significance()`](https://yelleknek.github.io/DMAR/reference/regions_of_significance.md),
[`simple_effects_AB()`](https://yelleknek.github.io/DMAR/reference/simple_effects_AB.md),
[`summary_t_test()`](https://yelleknek.github.io/DMAR/reference/summary_t_test.md),
[`welch_t()`](https://yelleknek.github.io/DMAR/reference/welch_t.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Four repeated measures on each of 30 subjects.
set.seed(113)
n <- 30; k <- 4
subj <- factor(rep(1:n, each = k))
tm   <- factor(rep(1:k, times = n))
y    <- as.vector(t(matrix(rnorm(n * k), n, k) +
                     rep(rnorm(n, 0, 1), each = k)))
d <- data.frame(y, subj, tm)

# All eight structures at once. Read the table by comparing AIC and BIC
# across rows, and use the likelihood-ratio test only for the nested
# comparison it reports, each structure against UN.
compare_cov_structures(d, outcome = "y", subject = "subj",
                       time = "tm")
#>   structure   log_lik      AIC      BIC n_par LRT_vs_UN_chisq LRT_vs_UN_df
#> 1       IND -198.6529 407.3058 421.2432     5        9.753075            9
#> 2        CS -198.6389 409.2777 426.0027     6        9.725060            8
#> 3       CSH -194.9519 407.9038 432.9913     9        2.351142            5
#> 4       AR1 -198.5164 409.0328 425.7578     6        9.480155            8
#> 5      ARH1 -194.7986 407.5971 432.6846     9        2.044453            5
#> 6      TOEP -198.4499 412.8998 435.1997     8        9.347092            6
#> 7     TOEPH -194.7293 411.4585 442.1210    11        1.905857            3
#> 8        UN -193.7763 415.5527 454.5776    14              NA           NA
#>   LRT_vs_UN_p
#> 1   0.3708427
#> 2   0.2848558
#> 3   0.7987270
#> 4   0.3034218
#> 5   0.8429605
#> 6   0.1549782
#> 7   0.5921744
#> 8          NA

# A subset, requested with lowercase aliases (matching is case
# insensitive).
compare_cov_structures(d, outcome = "y", subject = "subj",
                       time = "tm",
                       structures = c("cs", "csh", "ar1", "arh1"))
#>   structure   log_lik      AIC      BIC n_par LRT_vs_UN_chisq LRT_vs_UN_df
#> 1        CS -198.6389 409.2777 426.0027     6              NA           NA
#> 2       CSH -194.9519 407.9038 432.9913     9              NA           NA
#> 3       AR1 -198.5164 409.0328 425.7578     6              NA           NA
#> 4      ARH1 -194.7986 407.5971 432.6846     9              NA           NA
#>   LRT_vs_UN_p
#> 1          NA
#> 2          NA
#> 3          NA
#> 4          NA
```
