# Gwet's AC1 and AC2 Chance-Corrected Agreement Coefficients

Computes Gwet's (2008, 2014) AC1 (nominal data) and AC2 (ordinal data
with user-supplied weights) chance-corrected agreement coefficients for
two or more raters. AC1/AC2 are more robust than Cohen's \\\kappa\\ to
extreme marginal-prevalence imbalance and the trait-distribution
paradox.

## Usage

``` r
gwet_ac(
  ratings,
  weights = c("unweighted", "linear", "quadratic"),
  conf_level = 0.95
)
```

## Arguments

- ratings:

  A units \\\times\\ raters matrix or `data.frame`. Rows = units;
  columns = raters. `NA` entries are allowed.

- weights:

  One of `"unweighted"` (default; computes AC1) or `"linear"` /
  `"quadratic"` (compute AC2 with weight matrices used by
  weighted-\\\kappa\\ conventions).

- conf_level:

  Confidence level. Default `0.95`.

## Value

A `data.frame` with rows for the point estimate
\\\widehat{\mathrm{AC}}\\, the standard error, the CI lower and upper
limits, the percent agreement \\p_a\\, and the chance-agreement term
\\p_e\\.

## Details

**Coefficient.** \$\$\widehat{\mathrm{AC}} \\=\\ \frac{p_a - p_e}{1 -
p_e},\$\$ identical to Cohen's \\\kappa\\ in structure but with a
different chance-correction \\p_e\\: \$\$p_e \\=\\ \frac{T_w}{Q (Q - 1)}
\sum\_{k = 1}^{Q} \pi_k (1 - \pi_k),\$\$ where \\\pi_k\\ is the mean
within-unit proportion of category \\k\\, \\Q\\ is the number of
categories, and \\T_w\\ is the sum of all entries of the weight matrix.
For nominal data with unit weights (AC1), \\T_w = Q\\ and \\p_e\\
reduces to \\(1 / (Q - 1)) \sum_k \pi_k (1 - \pi_k)\\.

**Why AC over \\\kappa\\.** \\\kappa\\ can be near zero even when raters
agree on almost every unit if the trait is rare or very common (the
"kappa paradox"; Feinstein & Cicchetti, 1990). Gwet's AC keeps the same
chance-correction logic but uses a less extreme reference distribution.

**Variance.** The SE is Gwet's (2008) linearization variance,
\$\$\mathrm{Var}(\widehat{\mathrm{AC}}) \\=\\ \frac{1 - f}{n (n - 1)}
\sum_i (\widehat{\mathrm{AC}}\_i^{\*} - \widehat{\mathrm{AC}})^2,\$\$
where \\\widehat{\mathrm{AC}}\_i^{\*}\\ is the \\i\\th unit's influence
value, combining its agreement and chance-term contributions, \\n\\ is
the number of units, and \\f\\ is the sampling fraction (\\0\\ for an
infinite target population). The interval is \\\widehat{\mathrm{AC}} \pm
t\_{1 - \alpha / 2,\\ n - 1} \mathit{SE}\\, with the upper limit
truncated at \\1\\. These quantities match Gwet's (2014) reference
software.

## References

Feinstein, A. R., & Cicchetti, D. V. (1990). High agreement but low
kappa: I. The problems of two paradoxes. *Journal of Clinical
Epidemiology, 43*(6), 543–549.
[doi:10.1016/0895-4356(90)90158-L](https://doi.org/10.1016/0895-4356%2890%2990158-L)

Gwet, K. L. (2008). Computing inter-rater reliability and its variance
in the presence of high agreement. *British Journal of Mathematical and
Statistical Psychology, 61*(1), 29–48.
[doi:10.1348/000711006X126600](https://doi.org/10.1348/000711006X126600)

Gwet, K. L. (2014). *Handbook of inter-rater reliability* (4th ed.).
Advanced Analytics, LLC.

## See also

[`cohen_kappa`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md),
[`fleiss_kappa`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md),
[`krippendorff_alpha`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md)

Other agreement and measurement:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md),
[`loa()`](https://yelleknek.github.io/DMAR/reference/loa.md),
[`variance_components_mls()`](https://yelleknek.github.io/DMAR/reference/variance_components_mls.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Unweighted AC1, two raters, nominal:
set.seed(113)
r1 <- sample(c("A", "B", "C"), 50, replace = TRUE)
r2 <- ifelse(runif(50) < 0.8, r1, sample(c("A", "B", "C"), 50, TRUE))
gwet_ac(cbind(r1, r2))
#>  term              value 
#>  gwet_ac           0.881 
#>  se                0.0578
#>  lower_limit       0.764 
#>  upper_limit       0.997 
#>  percent_agreement 0.92  
#>  chance_agreement  0.33  
#>  n_units           50    
#> 
#> Confidence level: 95%

# 2. AC2 with linear weights, ordinal scale 1-5:
set.seed(113)
r1 <- sample(1:5, 60, replace = TRUE)
r2 <- pmin(5, pmax(1, r1 + sample(-1:1, 60, replace = TRUE)))
gwet_ac(cbind(r1, r2), weights = "linear")
#>  term              value 
#>  gwet_ac           0.67  
#>  se                0.0405
#>  lower_limit       0.589 
#>  upper_limit       0.751 
#>  percent_agreement 0.867 
#>  chance_agreement  0.596 
#>  n_units           60    
#> 
#> Confidence level: 95%
```
