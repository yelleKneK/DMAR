# Depression Treatment Study With a Pretest Covariate

The hypothetical three-group depression study that runs through the
analysis of covariance development of Maxwell, Delaney, and Kelley
(2027, Chapter 9, Table 9.7). Thirty depressive individuals are randomly
assigned, ten per group, to a selective serotonin reuptake inhibitor
(SSRI), a placebo, or a wait list control. The Beck Depression Inventory
(BDI) is administered before the study begins and again at its end,
giving a pretest that serves as the covariate and a posttest that serves
as the outcome.

## Usage

``` r
depression_bdi
```

## Format

A `data.frame` with 30 rows and 3 columns:

- condition:

  Factor with levels `ssri`, `placebo`, and `wait_list`: the randomly
  assigned treatment.

- bdi_pre:

  Beck Depression Inventory score before the study.

- bdi_post:

  Beck Depression Inventory score at the end of the study.

## Details

This is the worked example behind the ANCOVA contrast pages: the pretest
group means are 17, 17.7, and 17.4; the within-groups sum of squares of
the pretest is 752.5; the ANCOVA error variance is about 29; and the
covariate-adjusted posttest means are approximately 7.5, 12, and 14 for
the SSRI, placebo, and wait list groups. The examples of
[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md)
and
[`ci_sc_ancova`](https://yelleknek.github.io/DMAR/reference/ci_sc_ancova.md)
quote those summary values, and a test recomputes each of them from
these data so the printed numbers cannot drift.

The data are hypothetical, constructed for the book; higher BDI scores
mean more severe depressive symptoms, so a treatment that works pulls
the posttest down. The same numeric values ship in the book's data
companion, the AMCP package, as `chapter_9_table_7`; DMAR carries them
directly so its ANCOVA examples and tests need no package beyond this
one.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 on ANCOVA.)

## Examples

``` r
data(depression_bdi)

# The fingerprints the ANCOVA contrast pages quote.
tapply(depression_bdi$bdi_pre, depression_bdi$condition, mean)
#>      ssri   placebo wait_list 
#>      17.0      17.7      17.4 
fit <- lm(bdi_post ~ bdi_pre + condition, data = depression_bdi)
anova(fit)
#> Analysis of Variance Table
#> 
#> Response: bdi_post
#>           Df Sum Sq Mean Sq F value   Pr(>F)   
#> bdi_pre    1 336.68  336.68 11.5739 0.002174 **
#> condition  2 217.15  108.57  3.7324 0.037584 * 
#> Residuals 26 756.33   29.09                    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# The covariate-adjusted group means at the pretest grand mean.
ancova(depression_bdi, outcome = "bdi_post", treatment = "condition",
       covariates = "bdi_pre")
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution, so the lower confidence limit on partial eta squared is 0.
#> Warning: The observed F_value is below the alpha_lower critical value of the central F-distribution, so the lower confidence limit on omega squared is 0.
#>  term                         value 
#>  F_value                      3.73  
#>  df_1                         2     
#>  df_2                         26    
#>  p_value                      0.0376
#>  sum_of_squares_type          3     
#>  eta_squared_partial          0.223 
#>  eta_squared_partial_lower    0     
#>  eta_squared_partial_upper    0.421 
#>  omega_squared_partial        0.154 
#>  omega_squared_partial_lower  0     
#>  omega_squared_partial_upper  0.421 
#>  adjusted_mean[ssri]          7.54  
#>  adjusted_mean[placebo]       12    
#>  adjusted_mean[wait_list]     14    
#>  se_adjusted_mean[ssri]       1.71  
#>  se_adjusted_mean[placebo]    1.71  
#>  se_adjusted_mean[wait_list]  1.71  
#>  F_homogeneity_of_regression  1.06  
#>  df_homogeneity_of_regression 2     
#>  p_homogeneity_of_regression  0.3614
#> 
#> Confidence level: 95%
```
