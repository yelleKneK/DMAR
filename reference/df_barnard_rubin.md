# Barnard and Rubin's Small-Sample Degrees of Freedom for Multiple Imputation

Computes the degrees of freedom of the *t* reference distribution for an
estimate pooled over *m* imputed data sets by Rubin's (1987) rules,
using Barnard and Rubin's (1999) small-sample formula, together with the
quantities it is built from: the total variance, the fraction of missing
information the imputations estimate, the relative increase in variance
due to the missing values, and the two components of the degrees of
freedom. The function is vectorized over parameters.

## Usage

``` r
df_barnard_rubin(between, within, m, df_complete = Inf, floor = 3)
```

## Arguments

- between:

  The between-imputation variance \\B\\ of each parameter: the variance,
  with divisor \\m - 1\\, of its *m* estimates. A numeric vector, one
  value per parameter, whose names label the rows. (It is not called `B`
  because DMAR keeps `B` for the number of bootstrap replications.)

- within:

  The within-imputation variance \\\bar{U}\\ of each parameter: the mean
  of its *m* squared standard errors. A numeric vector of the same
  length as `between`.

- m:

  The number of imputations, a single whole number of 2 or more.

- df_complete:

  The complete data degrees of freedom, of length 1 or the length of
  `between`. Defaults to `Inf`, which gives Rubin's (1987) large-sample
  degrees of freedom; a small-sample analysis should supply its complete
  data error degrees of freedom.

- floor:

  The lower bound on the degrees of freedom: a single number of 0 or
  more. Defaults to 3, von Hippel's (2016) recommendation; `0`
  reproduces Barnard and Rubin (1999).

## Value

A `data.frame` (class `dmar_tbl`) with one row per parameter and the
columns `term`, `df` (bounded below by `floor`), `df_unbounded`
(\\\nu\_{BR}\\), `gamma`, `relative_increase` (\\r\\), `total_variance`
(\\T\\), `nu_m`, `nu_obs`, and `at_floor`. The attributes `m`,
`df_complete`, and `floor` record the call, and the printed table, at
the console or knitted, is followed by a line naming the rule, *m*,
`df_complete`, and the floor. A missing input value gives missing values
in its row.

## Details

**Rubin's rules.** Each of the *m* completed data sets gives an estimate
\\\hat{\theta}\_d\\ and a squared standard error \\U_d\\. The pooled
estimate \\\bar{\theta}\\ is their mean, the within-imputation variance
is \\\bar{U} = m^{-1} \sum_d U_d\\, the between-imputation variance is
\\B = (m - 1)^{-1} \sum_d (\hat{\theta}\_d - \bar{\theta})^2\\, and the
total variance of the pooled estimate is \$\$T = \bar{U} + (1 + 1/m)
B.\$\$ The share of the total variance due to the missing values,
\$\$\gamma = (1 + 1/m) B / T,\$\$ is the multiple imputation estimate of
the fraction of missing information, which mice reports as `lambda`, and
\\r = (1 + 1/m) B / \bar{U}\\ is the relative increase in variance due
to nonresponse (mice's `riv`). mice's `fmi` column is a different
quantity, \\(r + 2/(\nu + 3))/(r + 1)\\, which adjusts \\\gamma\\ for
the degrees of freedom. Von Hippel's (2016) Eq. 25 prints the
within-imputation variance in place of \\B\\ in the numerator of
\\\gamma\\, a typographical error: that ratio would exceed 1 with no
missing values at all.

**The degrees of freedom.** Rubin's (1987) degrees of freedom, \$\$\nu_m
= (m - 1) / \gamma^2,\$\$ reflect only the finite number of imputations:
they treat each completed data analysis as if it had infinitely many
degrees of freedom, and when \\\gamma\\ is small they can exceed the
degrees of freedom the analysis would have with no missing values at
all, which cannot be right in a small sample. Barnard and Rubin (1999)
combine them with observed data degrees of freedom, \$\$\nu\_{obs} =
\nu\_{com} (1 - \gamma) \frac{\nu\_{com} + 1}{\nu\_{com} + 3},\$\$ where
\\\nu\_{com}\\, given as `df_complete`, is the complete data degrees of
freedom (for a regression coefficient, \\n - k\\, with \\k\\ the number
of regression parameters), as \$\$\nu\_{BR} = \left(\frac{1}{\nu_m} +
\frac{1}{\nu\_{obs}}\right)^{-1},\$\$ which never exceeds \\\nu\_{obs}\\
and hence never exceeds \\\nu\_{com}\\. The formulas are as von Hippel
(2016, Eqs. 23 and 24) writes them, and they agree with mice's
`barnard.rubin()`. With `df_complete = Inf`, the default, \\\nu\_{obs}\\
is infinite and \\\nu\_{BR}\\ is Rubin's \\\nu_m\\. The observed data
degrees of freedom are
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)'s
`"complete_df"` rule, and von Hippel's degrees of freedom for maximum
likelihood are the limit of \\\nu\_{BR}\\ as the number of imputations
grows, where \\\nu_m\\ becomes infinite.

**The floor.** An estimate of \\\gamma\\ near 1 drives the degrees of
freedom toward 0, where the *t* quantile grows without bound. In von
Hippel's (2016, p. 431) simulation of 25 cases with half the values of
one variable missing at random, 24% of the Barnard and Rubin degrees of
freedom were below 3, 11% below 2, and 2% below 1, and some intervals
were orders of magnitude longer than coverage required. He bounds the
degrees of freedom below at 3 (his Eq. 27), which is the default
`floor = 3` here; `floor = 0` reproduces Barnard and Rubin's formula,
and mice, exactly. The unbounded value is kept in `df_unbounded`, and
`at_floor` marks the rows where the bound applied.

DMAR does not impute; this function serves analyses pooled elsewhere and
the comparison of maximum likelihood with multiple imputation. For a
FIML fit,
[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
gives the corresponding degrees of freedom from the fit's own fraction
of missing information, without imputations.

## References

Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika, 86*(4), 948–955.
[doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948)

Rubin, D. B. (1987). *Multiple imputation for nonresponse in surveys*.
Wiley.

van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
imputation by chained equations in R. *Journal of Statistical Software,
45*(3), 1–67.
[doi:10.18637/jss.v045.i03](https://doi.org/10.18637/jss.v045.i03)

von Hippel, P. T. (2016). New confidence intervals and bias comparisons
show that maximum likelihood can beat multiple imputation in small
samples. *Structural Equation Modeling, 23*(3), 422–437.
[doi:10.1080/10705511.2015.1047931](https://doi.org/10.1080/10705511.2015.1047931)

## See also

[`df_missing`](https://yelleknek.github.io/DMAR/reference/df_missing.md)
for the degrees of freedom of a FIML estimate.

Other missing data:
[`cov_fiml()`](https://yelleknek.github.io/DMAR/reference/cov_fiml.md),
[`df_missing()`](https://yelleknek.github.io/DMAR/reference/df_missing.md),
[`fraction_missing_information()`](https://yelleknek.github.io/DMAR/reference/fraction_missing_information.md),
[`missing_coverage()`](https://yelleknek.github.io/DMAR/reference/missing_coverage.md),
[`missing_patterns()`](https://yelleknek.github.io/DMAR/reference/missing_patterns.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# One coefficient from m = 5 imputations, with complete data error
# degrees of freedom of 23: B = 0.04 and a within-imputation variance of
# 0.10. floor = 0 gives Barnard and Rubin's formula itself.
df_barnard_rubin(between = 0.04, within = 0.10, m = 5, df_complete = 23,
                 floor = 0)
#>  term        df     df_unbounded gamma relative_increase total_variance nu_m  
#>  parameter_1 10.416 10.416       0.324 0.480             0.148          38.028
#>  nu_obs at_floor
#>  14.345 FALSE   
#> 
#> Degrees of freedom: Barnard and Rubin's (1999) rule for m = 5 imputations,
#> 1 / (1 / nu_m + 1 / nu_obs) with df_complete = 23, bounded below at 0

# Rubin's (1987) large-sample degrees of freedom are 38, well above the
# 23 a complete data analysis of these cases would have.
df_barnard_rubin(between = 0.04, within = 0.10, m = 5)
#>  term        df     df_unbounded gamma relative_increase total_variance nu_m  
#>  parameter_1 38.028 38.028       0.324 0.480             0.148          38.028
#>  nu_obs at_floor
#>  Inf    FALSE   
#> 
#> Degrees of freedom: Rubin's (1987) large-sample rule for m = 5 imputations,
#> (m - 1) / gamma^2 (df_complete = Inf), bounded below at 3

# Pooling by hand: five estimates of each of two coefficients and their
# squared standard errors (illustrative values), one column per
# coefficient.
est <- cbind(slope = c(0.42, 0.51, 0.38, 0.47, 0.55),
             intercept = c(1.90, 1.94, 1.88, 1.93, 1.91))
u <- cbind(slope = c(0.010, 0.011, 0.009, 0.010, 0.012),
           intercept = c(0.020, 0.021, 0.019, 0.020, 0.022))
df_barnard_rubin(between = apply(est, 2, var), within = colMeans(u),
                 m = 5, df_complete = 23)
#>  term      df     df_unbounded gamma relative_increase total_variance nu_m    
#>  slope     9.749  9.749        0.348 0.534             0.016          32.990  
#>  intercept 20.432 20.432       0.032 0.034             0.021          3800.614
#>  nu_obs at_floor
#>  13.838 FALSE   
#>  20.542 FALSE   
#> 
#> Degrees of freedom: Barnard and Rubin's (1999) rule for m = 5 imputations,
#> 1 / (1 / nu_m + 1 / nu_obs) with df_complete = 23, bounded below at 3
```
