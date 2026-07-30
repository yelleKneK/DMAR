# Bryant–Paulson Simultaneous Confidence Intervals for Contrasts of Adjusted Means in ANCOVA

Constructs Tukey–Kramer-type *simultaneous* confidence intervals on one
or more contrasts of covariate-adjusted means in the analysis of
covariance (ANCOVA) when the covariate(s) are *random*, using the
Bryant–Paulson generalized studentized range (see
[`bryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)).
Unlike the per-comparison interval of
[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
these intervals control the *familywise* error rate over the whole
family of comparisons and, through the Bryant–Paulson critical value,
correctly account for the extra sampling uncertainty that comes from
estimating the covariate adjustment from random covariates. Naively
applying Tukey's method to adjusted means ignores that uncertainty and
produces intervals that are too narrow (below-nominal coverage); see the
package vignette “Bryant–Paulson simultaneous intervals: a simulation
study.”

## Usage

``` r
ci_c_ancova_bp(
  adj_means,
  s_ancova,
  c_weights = NULL,
  n,
  num_covariates = 1,
  df = NULL,
  conf_level = 0.95,
  contrast_type = c("pairwise", "allowance"),
  ...
)
```

## Arguments

- adj_means:

  A numeric vector of the covariate-adjusted group means (one per
  group). The number of groups \\k\\ is inferred from its length.

- s_ancova:

  The standard deviation of the errors from the ANCOVA model, i.e., the
  square root of the ANCOVA error mean square (use the standard
  deviation, *not* the variance from the source table).

- c_weights:

  Optional contrast weights. May be (i) a numeric vector of length \\k\\
  giving a single contrast, or (ii) a matrix/data.frame with \\k\\
  columns, each row a contrast. If `NULL` (the default), all
  \\k(k-1)/2\\ pairwise comparisons are returned. For each contrast the
  weights should sum to zero.

- n:

  Either a single number giving the common per-group sample size or a
  numeric vector of per-group sample sizes. The Bryant–Paulson
  distribution is exact for balanced designs; for unequal \\n\\ a
  Tukey–Kramer harmonic adjustment is used (see Details).

- num_covariates:

  The number of random covariates, \\p\\. Default `1`.

- df:

  Optional error degrees of freedom \\\nu\\. If `NULL` (default) it is
  computed for a one-way ANCOVA as \\\nu = \sum n - k - p\\. Supply it
  directly for other designs (e.g., a randomized-block ANCOVA, where
  \\\nu\\ differs).

- conf_level:

  The simultaneous (familywise) confidence level. Default `0.95`.

- contrast_type:

  One of `"pairwise"` (default) or `"allowance"`. `"pairwise"` uses the
  quadratic standard error \\\sqrt{\sum c_i^2 / n_i}\\ (Tukey–Kramer),
  exact for pairwise comparisons. `"allowance"` uses Tukey's allowance
  \\\tfrac12 \sum \|c_i\|\\, which yields intervals that hold
  simultaneously over *all* contrasts, including complex ones (this is
  the form in Eq. (2.4) of Bryant and Bruvold, 1980). The two coincide
  for pairwise comparisons.

- ...:

  Additional arguments (currently unused).

## Value

A `data.frame` (class `dmar_tbl`) with one row per contrast and columns
`contrast` (a label), `estimate` (the contrast of adjusted means
\\\hat\psi\\), `lower_limit`, and `upper_limit`. The Bryant–Paulson
critical value used is stored in the `"critical_value"` attribute and
the confidence level in the `"conf_level"` attribute (printed beneath
the table).

## Details

**The interval.** For a contrast \\\psi = \sum_i c_i \theta_i\\ of
adjusted means, the simultaneous interval is \$\$\hat\psi \\\pm\\
q\_{\alpha;\\p,k,\nu}\\ \hat\sigma\_{y\mid x}\\ w(c),\$\$ where
\\q\_{\alpha;\\p,k,\nu}\\ is the upper-\\\alpha\\ Bryant–Paulson
critical value
([`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md))
and the width factor is \\w(c) = \tfrac{1}{\sqrt2}\sqrt{\sum_i
c_i^2/n_i}\\ for `contrast_type = "pairwise"` or \\w(c) = \tfrac12
\sum_i \|c_i\| \sqrt{1/n}\\ for `contrast_type = "allowance"` (balanced
\\n\\). For a pairwise difference with common \\n\\ both reduce to
\\q\_{\alpha;\\p,k,\nu}\\\hat\sigma\_{y\mid x}\sqrt{1/n}\\, reproducing
the critical difference of Bryant and Bruvold (1980).

**No per-comparison covariate term.** By design the standard error here
does *not* include the \\(\bar X_i - \bar X_j)^2 /
SS\_{\mathrm{within}(x)}\\ term that appears in a single-comparison
ANCOVA interval
([`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md)).
In the Bryant–Paulson framework the random-covariate uncertainty is
carried by the (larger) critical value, which holds *on average* over
the covariate distribution; adding the per-pair term as well would
double-count it.

**Unequal sample sizes.** For unbalanced designs the `"pairwise"`
standard error uses \\\sqrt{c_i^2/n_i}\\ directly (the Tukey–Kramer
generalization); coverage is then approximate but typically very close
to nominal and slightly conservative.

## References

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*, 631–638.

Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures
in the analysis of covariance. *Journal of the American Statistical
Association, 75*(372), 874–880.
[doi:10.2307/2287175](https://doi.org/10.2307/2287175)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`bryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md)
for the underlying distribution;
[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md)
for the per-comparison interval;
[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md) for a
ANCOVA fit.

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_cc()`](https://yelleknek.github.io/DMAR/reference/ci_cc.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
[`ci_r()`](https://yelleknek.github.io/DMAR/reference/ci_r.md),
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md),
[`ci_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ci_reg_coef.md),
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md),
[`ci_sc()`](https://yelleknek.github.io/DMAR/reference/ci_sc.md),
[`ci_sc_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md),
[`ci_sm()`](https://yelleknek.github.io/DMAR/reference/ci_sm.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md),
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md),
[`ci_src()`](https://yelleknek.github.io/DMAR/reference/ci_src.md),
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md),
[`contrast_adjusted()`](https://yelleknek.github.io/DMAR/reference/contrast_adjusted.md),
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Bryant & Bruvold (1980) worked example: 6 panels, 1 covariate, nu = 14,
# ANCOVA error MS = 0.01326. Here the design is a randomized block with
# s = 4 blocks, so the per-group "n" for the adjusted-mean SE is 4 and the
# error df (14) must be supplied directly.
adj <- c(3.595, 3.619, 4.102, 4.515, 4.618, 4.876)
ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326),
               n = 4, num_covariates = 1, df = 14)
#>  contrast          estimate lower_limit upper_limit
#>  group_1 - group_2 -0.024   -0.302      0.254      
#>  group_1 - group_3 -0.507   -0.785      -0.229     
#>  group_1 - group_4 -0.92    -1.2        -0.642     
#>  group_1 - group_5 -1.02    -1.3        -0.745     
#>  group_1 - group_6 -1.28    -1.56       -1         
#>  group_2 - group_3 -0.483   -0.761      -0.205     
#>  group_2 - group_4 -0.896   -1.17       -0.618     
#>  group_2 - group_5 -0.999   -1.28       -0.721     
#>  group_2 - group_6 -1.26    -1.54       -0.979     
#>  group_3 - group_4 -0.413   -0.691      -0.135     
#>  group_3 - group_5 -0.516   -0.794      -0.238     
#>  group_3 - group_6 -0.774   -1.05       -0.496     
#>  group_4 - group_5 -0.103   -0.381      0.175      
#>  group_4 - group_6 -0.361   -0.639      -0.0829    
#>  group_5 - group_6 -0.258   -0.536      0.0201     
#> 
#> Confidence level: 95%
# Each pairwise critical difference is 0.278, matching the paper.

# A single complex contrast (panels 1-2 vs. 3-6), simultaneous over all
# contrasts via the allowance form:
ci_c_ancova_bp(adj_means = adj, s_ancova = sqrt(0.01326), n = 4, df = 14,
               c_weights = c(0.5, 0.5, -0.25, -0.25, -0.25, -0.25),
               contrast_type = "allowance")
#>  contrast   estimate lower_limit upper_limit
#>  contrast_1 -0.921   -1.2        -0.643     
#> 
#> Confidence level: 95%
```
