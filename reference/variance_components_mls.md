# Modified-Large-Sample Confidence Intervals on Variance Components

Computes modified-large-sample (MLS) confidence intervals on the
between-group and within-group variance components of a balanced one-way
random-effects ANOVA, following Burdick & Graybill (1992). MLS intervals
have substantially better coverage than the Satterthwaite or simple-Wald
intervals when the components are far from zero, and they are the
standard interval method in generalizability theory (Brennan, 2001).

## Usage

``` r
variance_components_mls(
  ms_between,
  ms_within,
  df_between,
  df_within,
  n,
  conf_level = 0.95
)
```

## Arguments

- ms_between:

  Mean square between groups (numerator of the ANOVA *F*).

- ms_within:

  Mean square within groups (denominator of the ANOVA *F*).

- df_between:

  Degrees of freedom for the between-group MS (typically \\a - 1\\ for
  \\a\\ groups).

- df_within:

  Degrees of freedom for the within-group MS (typically \\a (n - 1)\\
  for \\a\\ groups of size \\n\\).

- n:

  Number of observations per group (assumed balanced).

- conf_level:

  Confidence level for the CIs. Default `0.95`.

## Value

A `data.frame` with rows for the point estimates and MLS lower / upper
CIs of the between-group variance component (\\\sigma^2_b\\), the
within-group component (\\\sigma^2_w\\), and the implied intraclass
correlation (\\\rho = \sigma^2_b / (\sigma^2_b + \sigma^2_w)\\).

## Details

**Point estimates.** For a one-way random-effects ANOVA on \\a\\ groups
of size \\n\\, the method-of-moments estimators are \$\$\hat\sigma^2_b
\\=\\ \max(0,\\ (\mathit{MS}\_b - \mathit{MS}\_w)/n), \qquad
\hat\sigma^2_w \\=\\ \mathit{MS}\_w.\$\$

**Modified-large-sample CIs (Burdick-Graybill 1992).** The MLS interval
for \\\sigma^2_b\\ is \$\$\left\[\frac{\mathit{MS}\_b - \mathit{MS}\_w -
\sqrt{V_L}}{n},\\\\ \frac{\mathit{MS}\_b - \mathit{MS}\_w +
\sqrt{V_U}}{n}\right\],\$\$ with \$\$V_L \\=\\ G_1^2 \mathit{MS}\_b^2 +
H_2^2 \mathit{MS}\_w^2 + G\_{12} \mathit{MS}\_b \mathit{MS}\_w, \qquad
V_U \\=\\ H_1^2 \mathit{MS}\_b^2 + G_2^2 \mathit{MS}\_w^2 + H\_{12}
\mathit{MS}\_b \mathit{MS}\_w,\$\$ where the constants \\G_1\\, \\G_2\\,
\\H_1\\, \\H_2\\ and the cross-term constants \\G\_{12}\\, \\H\_{12}\\
depend on the degrees of freedom and on \\\chi^2\\ and *F* quantiles at
the chosen confidence level (Burdick & Graybill, 1992, equations
2.4.1–2.4.5 give the explicit formulas). The lower limit is truncated at
zero. For the within-group component, the standard \\\chi^2\\-based CI
on \\\mathit{MS}\_w\\ (Searle, Casella, & McCulloch, 1992) is used.

**Caveats.** MLS intervals assume balanced data and homogeneous
variances within groups. For unbalanced data the appropriate analog is
the Burdick-Graybill MLS extension to unequal sample sizes (Burdick &
Graybill, 1992, Section 2.5), which is not implemented here.

## References

Brennan, R. L. (2001). *Generalizability theory*. Springer.

Burdick, R. K., & Graybill, F. A. (1992). *Confidence intervals on
variance components*. Marcel Dekker.

Searle, S. R., Casella, G., & McCulloch, C. E. (1992). *Variance
components*. Wiley.

## See also

[`icc`](https://yelleknek.github.io/DMAR/reference/icc.md),
[`var_icc`](https://yelleknek.github.io/DMAR/reference/var_icc.md),
[`ss_aipe_icc`](https://yelleknek.github.io/DMAR/reference/ss_aipe_icc.md)

Other agreement and measurement:
[`R2_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md),
[`content_validity_index()`](https://yelleknek.github.io/DMAR/reference/content_validity_index.md),
[`gwet_ac()`](https://yelleknek.github.io/DMAR/reference/gwet_ac.md),
[`icc_lmer()`](https://yelleknek.github.io/DMAR/reference/icc_lmer.md),
[`krippendorff_alpha()`](https://yelleknek.github.io/DMAR/reference/krippendorff_alpha.md),
[`limits_of_agreement()`](https://yelleknek.github.io/DMAR/reference/limits_of_agreement.md),
[`lin_ccc()`](https://yelleknek.github.io/DMAR/reference/lin_ccc.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# 1. Balanced one-way random-effects ANOVA: a = 10 groups, n = 5.
#        Hypothetical MS_b = 6.0, MS_w = 1.5.
variance_components_mls(ms_between = 6.0, ms_within = 1.5,
                        df_between = 9, df_within = 40, n = 5)
#>  term           value 
#>  sigma2_between 0.9   
#>  sigma2_b_lower 0.236 
#>  sigma2_b_upper 3.69  
#>  sigma2_within  1.5   
#>  sigma2_w_lower 1.01  
#>  sigma2_w_upper 2.46  
#>  icc            0.375 
#>  icc_lower      0.0876
#>  icc_upper      0.785 
#> 
#> Confidence level: 95%
```
