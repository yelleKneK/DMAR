# Provides the Critical Value for the Bryant–Paulson ANCOVA Multiple-Comparison Procedure

Computes the critical value of the Bryant–Paulson generalized
studentized range, the reference distribution for multiple comparisons
of adjusted means in an analysis of covariance with random covariates.
The single-step value is the multiplier for simultaneous confidence
intervals on pairwise differences of adjusted means; the `"duncan"`
option instead returns the significant range of the stepwise Duncan
multiple-range procedure.

## Usage

``` r
cv_bryant_paulson(
  alpha_level,
  df,
  groups,
  covariates = 1,
  procedure = c("tukey", "duncan"),
  verbose = TRUE
)
```

## Arguments

- alpha_level:

  Type I error rate (i.e., the false positive rate). As with
  [`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
  the full `alpha_level` applies to the upper tail of the (non-negative)
  generalized studentized range distribution; it is not split between
  two tails (see Details).

- df:

  The ANCOVA error degrees of freedom (a positive number; in a one-way
  ANCOVA, `df = N - groups - covariates`).

- groups:

  The number of groups whose adjusted means are being compared (an
  integer of at least 2).

- covariates:

  The number of *random* covariates in the ANCOVA, the parameter \\p\\
  of the Bryant–Paulson distribution (a non-negative integer). With
  `covariates = 0` the critical value reduces to the ordinary
  studentized range and the result equals \\\sqrt2\\ times
  [`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
  (see Details).

- procedure:

  One of `"tukey"` (default) or `"duncan"`. `"tukey"` returns the
  single-step *simultaneous* critical value \\q\_{\alpha;p,k,\nu}\\ used
  for familywise (Tukey–Kramer-type) confidence intervals; `"duncan"`
  returns the stepwise Duncan multiple-range significant range
  \\r\_{\alpha;p,k,\nu}\\ tabled by Bryant and Bruvold (1980, Table 2)
  (see Details).

- verbose:

  Provides extra information (the tail areas) about the critical value.

## Value

Returns the critical value in a output style (a `data.frame` with class
`dmar_tbl` and one row per critical value, following the format used by
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md)
and [`cv_t`](https://yelleknek.github.io/DMAR/reference/cv_t.md)). The
`value` is on the *studentized-range scale* (the scale on which Bryant
and Paulson tabulate their critical values and on which
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)
uses them). When `verbose = TRUE` and `procedure = "tukey"`, the upper-
and lower-tail areas of the Bryant–Paulson distribution at the critical
value are also returned; for `procedure = "duncan"` the tail areas are
`NA` because the significant range is a stepwise quantity rather than a
single quantile.

## Details

The Bryant–Paulson procedure is the analysis-of-covariance
generalization of Tukey's method
([`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md))
for comparing adjusted means when the covariate is *random*. Because the
covariate adjustment must be estimated, the studentized range of
adjusted means is stochastically larger than the ordinary studentized
range, so the Bryant–Paulson critical value exceeds Tukey's; using the
latter would give intervals that are too narrow and a familywise error
rate above `alpha_level`. The single-step (`procedure = "tukey"`) value
is the multiplier for a family of simultaneous confidence intervals on
the pairwise differences of adjusted means that jointly hold at level
\\1 - \alpha\\. Maxwell, Delaney, and Kelley (2027, Chapter 9) develop
multiple comparisons of adjusted means in the analysis of covariance,
the setting this critical value serves.

The reference distribution is the Bryant–Paulson generalized studentized
range, implemented in
[`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md).
Its quantiles are not a standard base-R distribution and are not the
multivariate *t* quantiles that
[`cv_dunnett`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md)
and [`cv_smm`](https://yelleknek.github.io/DMAR/reference/cv_smm.md)
obtain from mvtnorm; they are computed directly by
[`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md),
so this function depends on neither base-R nor mvtnorm multiple-mean
machinery.

**Scale.** The returned `value` is on the studentized-range scale,
\\q\_{\alpha;p,k,\nu}\\, the scale of Bryant and Paulson's (1976) and
Bryant and Bruvold's (1980) tables and of Eq. (2.4) of the latter. A
pair of adjusted means is declared different when \\\|\hat\theta_i -
\hat\theta_j\| \> q\_{\alpha;p,k,\nu}\\\hat\sigma\_{y\mid
x}\sqrt{1/n}\\. This differs from
[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
which divides its value by \\\sqrt2\\ to report on the pairwise
mean-difference scale; divide the value here by \\\sqrt2\\ to obtain
that scale. With `covariates = 0`, `cv_bryant_paulson` returns exactly
\\\sqrt2 \times\\ `cv_tukey_hsd`.

**Duncan multiple-range.** For `procedure = "duncan"` the value is the
“significant range” of Duncan's stepwise test as extended to ANCOVA by
Bryant and Bruvold (1980, Section 4). With variable protection levels
\\\alpha_k = 1 - (1-\alpha)^{k-1}\\, \$\$r\_{\alpha;p,2,\nu} =
q\_{\alpha;p,2,\nu}, \qquad r\_{\alpha;p,k,\nu} = \max\\\\
r\_{\alpha;p,k-1,\nu},\\ q\_{\alpha_k;p,k,\nu} \\\\, \quad k \> 2.\$\$
These are the values in Bryant and Bruvold's Table 2, reproduced by this
function to the tabled two-decimal precision (see the package tests).

## References

Bryant, J. L., & Paulson, A. S. (1976). An extension of Tukey's method
of multiple comparisons to experimental designs with random concomitant
variables. *Biometrika, 63*, 631–638.

Bryant, J. L., & Bruvold, N. T. (1980). Multiple comparison procedures
in the analysis of covariance. *Journal of the American Statistical
Association, 75*(372), 874–880.
[doi:10.2307/2287175](https://doi.org/10.2307/2287175)

Duncan, D. B. (1955). Multiple range and multiple F tests. *Biometrics,
11*, 1–42.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9.)

## See also

[`cv_tukey_hsd`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`qbryant_paulson`](https://yelleknek.github.io/DMAR/reference/bryant_paulson.md),
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)

Other critical values:
[`cv_bonferroni_f()`](https://yelleknek.github.io/DMAR/reference/cv_bonferroni_f.md),
[`cv_chisq()`](https://yelleknek.github.io/DMAR/reference/cv_chisq.md),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
[`cv_f()`](https://yelleknek.github.io/DMAR/reference/cv_f.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md),
[`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Multiple comparisons of adjusted means in ANCOVA (the setting of Maxwell,
# Delaney, and Kelley, 2027, Chapter 9), using the worked example of Bryant
# & Bruvold (1980): 6 panels, 1 random covariate, 14 error df, alpha_level = .05.
# The single-step value is the multiplier for simultaneous confidence
# intervals on the pairwise differences of adjusted means; here it is 4.83.
cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 1)
#>  term     value area_less area_greater
#>  upper_cv 4.83  0.95      0.05        

# With no covariates it is sqrt(2) times the Tukey HSD critical value.
cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 0)$value
#> [1] 4.638538
sqrt(2) * cv_tukey_hsd(alpha_level = .05, df = 14, groups = 6)$value
#> [1] 4.638538

# \donttest{
# Duncan multiple-range significant range (Table 2 of the paper).
cv_bryant_paulson(alpha_level = .05, df = 14, groups = 6, covariates = 1,
                  procedure = "duncan")
#>  term     value area_less area_greater
#>  upper_cv 3.5   <NA>      <NA>        
# }
```
