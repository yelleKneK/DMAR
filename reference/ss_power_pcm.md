# Sample Size Planning for Power for Polynomial Change Models

Returns power given the sample size, or sample size given the desired
power, for the group difference in a polynomial change coefficient (a
flat-line intercept, a linear slope, a quadratic acceleration, or any
higher-order trend) in a two-group longitudinal design, following
Raudenbush and Liu (2001). The trend whose group difference is tested is
selected with `trend`; `trend = "linear"` (the default) reproduces the
straight-line case.

## Usage

``` r
ss_power_pcm(
  beta,
  tau,
  level_1_variance,
  frequency,
  duration,
  desired_power = NULL,
  N = NULL,
  alpha_level = 0.05,
  standardized = TRUE,
  directional = FALSE,
  trend = "linear"
)
```

## Arguments

- beta:

  The level two regression coefficient for the group by time interaction
  in the polynomial change coefficient selected by `trend` (linear by
  default), where the grouping variable is coded -.5 and .5 for the two
  groups. When `standardized = TRUE` (the default) this is the
  standardized change difference (see `standardized`).

- tau:

  The true between-subject variance of the individuals' change
  coefficient for the trend selected by `trend` (the variance of the
  slopes for `trend = "linear"`).

- level_1_variance:

  Level one (within-subject) error variance

- frequency:

  Number of measurements per unit of time, where the unit is the one in
  which `duration` is expressed. It need not be a whole number; for
  example, `frequency = 0.5` means one measurement every two time units.
  Together with `duration` it fixes the number of equally spaced
  measurement occasions, \\M = f \times D + 1\\, with \\f\\ the
  frequency and \\D\\ the duration.

- duration:

  Length of the study in the chosen time unit (for example, years,
  grades, or hours). Measurements are taken at \\M = f \times D + 1\\
  equally spaced occasions spanning time 0 to `duration`.

- desired_power:

  Desired power

- N:

  Total sample size (one-half in each of the two groups)

- alpha_level:

  Type I error rate

- standardized:

  The standardized change difference is the unstandardized change
  difference divided by the square root of `tau`, the between-subject
  variance of the change coefficient. `TRUE` (the default) treats `beta`
  as already standardized; `FALSE` treats it as the raw (unstandardized)
  change difference.

- directional:

  Should a one (`TRUE`) or two (`FALSE`) tailed test be performed.

- trend:

  The polynomial change coefficient whose group difference is tested,
  given either as a name (`"intercept"`, `"linear"`, `"quadratic"`,
  `"cubic"`, `"quartic"`, ...) or as a non-negative integer order (`0` =
  intercept / flat line, `1` = linear, `2` = quadratic, ...). Defaults
  to `"linear"`. The design must supply at least \\p + 1\\ measurement
  occasions to estimate a degree-\\p\\ trend.

## Value

A `data.frame` (class `dmar_tbl`) with one row per reported quantity in
a `term` / `value` layout: the per-group size per group
(`necessary_n_per_group`, or `specified_n_per_group` when `N` is
supplied) and the total (`total_N`); the achieved power
(`actual_power`); the measurement schedule (`freq`, `duration`,
`measurement_occasions`); the polynomial order of the tested change
coefficient (`polynomial_order`, 0 = intercept, 1 = linear, 2 =
quadratic, ...); the unstandardized and standardized change difference
(`unstd_coefficient`, `std_coefficient`); the level one error variance
(`l1_error_var`); the true and error variance of the change coefficient
(`true_var_of_slopes`, `error_var_of_slopes`, whose names retain
"slopes" from the linear case); the change-coefficient reliability
(`reliability`); and the noncentrality parameter of the *t* test
(`noncentral_t_parm`).

## Details

The two groups each contain \\N / 2\\ subjects measured on \\M = f
\times D + 1\\ equally spaced occasions. Each subject's degree-\\p\\
polynomial change coefficient is estimated within subject; the test
compares the two group means of that coefficient. The change coefficient
is taken in the derivative-scaled metric (\\p!\\ times the leading
coefficient of \\t^p\\), the metric in which the Raudenbush and Liu
(2001) constants apply.

The within-subject sampling variance of the estimated coefficient is
(Raudenbush & Liu, 2001, p. 392) \$\$V = \sigma^2_e\\ f^{2p}\\\frac{(M -
p - 1)!}{K_p\\(M + p)!}, \qquad \frac{1}{K_p} =
\frac{(2p)!\\(2p+1)!}{(p!)^2},\$\$ so that \\1/K_p = 1, 12, 720, 100800,
\ldots\\ for \\p = 0, 1, 2, 3, \ldots\\. This \\V\\ equals \\(p!)^2\\
times the variance of the ordinary least squares estimate of the
coefficient of \\t^p\\, reduces to \\\sigma^2_e / M\\ at \\p = 0\\ and
to \\12\\\sigma^2_e f^2 / \[M(M^2 - 1)\]\\ at \\p = 1\\. The slope
reliability is \\\tau / (\tau + V)\\ (their Equation 15), the variance
of the between-group difference is \\4(\tau + V)/N\\ (their Equation
10), and the *t* test has \\N - 2\\ degrees of freedom and noncentrality
\\\sqrt{N\\\beta^2\\\[\tau/(\tau + V)\]/4}\\ (their Equations 12 and
14). The linear case reproduces the National Youth Survey benchmark in
their Tables 1 and 2; the general-\\p\\ formula has been checked against
the exact \\(X'X)^{-1}\\ variance and against an end-to-end Monte Carlo
power study for the quadratic trend.

## References

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Kelley, K., & Rausch, J. R. (2011). Sample size planning for
longitudinal models: Accuracy in parameter estimation for polynomial
change parameters. *Psychological Methods, 16*(4), 391–405.
[doi:10.1037/a0023352](https://doi.org/10.1037/a0023352)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapters 11, 15.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Raudenbush, S. W., & Liu, X.-F. (2001). Effects of study duration,
frequency of observation, and sample size on power in studies of group
differences in polynomial change. *Psychological Methods, 6*(4),
387–401.
[doi:10.1037/1082-989X.6.4.387](https://doi.org/10.1037/1082-989X.6.4.387)

## See also

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_power_R2()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`ss_power_R2_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_R2_sensitivity.md),
[`ss_power_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_c.md),
[`ss_power_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_c_ancova.md),
[`ss_power_composite_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova.md),
[`ss_power_composite_ancova_2group()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_ancova_2group.md),
[`ss_power_composite_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_anova.md),
[`ss_power_composite_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova.md),
[`ss_power_composite_factorial_ancova_het()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_ancova_het.md),
[`ss_power_composite_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_factorial_anova.md),
[`ss_power_composite_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
[`ss_power_contrast()`](https://yelleknek.github.io/DMAR/reference/ss_power_contrast.md),
[`ss_power_equivalence_c()`](https://yelleknek.github.io/DMAR/reference/ss_power_equivalence_c.md),
[`ss_power_factorial_ancova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_ancova.md),
[`ss_power_factorial_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_factorial_anova.md),
[`ss_power_indirect_effect()`](https://yelleknek.github.io/DMAR/reference/ss_power_indirect_effect.md),
[`ss_power_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_power_mixed_effects.md),
[`ss_power_one_way_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_one_way_anova.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
[`ss_power_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef.md),
[`ss_power_reg_coef_sensitivity()`](https://yelleknek.github.io/DMAR/reference/ss_power_reg_coef_sensitivity.md),
[`ss_power_rm_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_rm_anova.md),
[`ss_power_sc()`](https://yelleknek.github.io/DMAR/reference/ss_power_sc.md),
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md),
[`ss_power_smd()`](https://yelleknek.github.io/DMAR/reference/ss_power_smd.md),
[`ss_power_split_plot_anova()`](https://yelleknek.github.io/DMAR/reference/ss_power_split_plot_anova.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# The examples reproduce the National Youth Survey illustration of
# Raudenbush and Liu (2001, p. 393). One observation per year
# (frequency = 1) over a four-year study (duration = 4) gives
# M = frequency * duration + 1 = 5 equally spaced occasions. The
# standardized slope difference is -0.40, the slope variance is
# tau = 0.003, and the level-one error variance is 0.0262, so the slope
# reliability is 0.53.

# (1) Power at a given sample size. With N = 238 (119 per group) the
#     design has power 0.61 (Raudenbush and Liu, 2001, Table 1, D = 4,
#     f = 1).
ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
             frequency = 1, duration = 4, N = 238)
#>  term                  value  
#>  specified_n_per_group 119    
#>  total_N               238    
#>  actual_power          0.612  
#>  freq                  1      
#>  duration              4      
#>  measurement_occasions 5      
#>  polynomial_order      1      
#>  unstd_coefficient     -0.0219
#>  std_coefficient       -0.4   
#>  l1_error_var          0.0262 
#>  true_var_of_slopes    0.003  
#>  error_var_of_slopes   0.00262
#>  reliability           0.534  
#>  noncentral_t_parm     2.25   

# (2) Sample size for a target power. Solving the same design for 0.80
#     power returns N = 370 (185 per group), between their Table 2 cells
#     N = 300 (power 0.71) and N = 400 (power 0.83).
ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
             frequency = 1, duration = 4, desired_power = .80)
#>  term                  value  
#>  necessary_n_per_group 185    
#>  total_N               370    
#>  actual_power          0.801  
#>  freq                  1      
#>  duration              4      
#>  measurement_occasions 5      
#>  polynomial_order      1      
#>  unstd_coefficient     -0.0219
#>  std_coefficient       -0.4   
#>  l1_error_var          0.0262 
#>  true_var_of_slopes    0.003  
#>  error_var_of_slopes   0.00262
#>  reliability           0.534  
#>  noncentral_t_parm     2.81   

# (3) Unstandardized slope. The unstandardized slope difference is
#     beta * sqrt(tau) = -0.40 * sqrt(0.003) = -0.0219. Passing it with
#     standardized = FALSE reproduces the power of 0.61 from example (1).
ss_power_pcm(beta = -.0219, tau = .003, level_1_variance = .0262,
             frequency = 1, duration = 4, N = 238, standardized = FALSE)
#>  term                  value  
#>  specified_n_per_group 119    
#>  total_N               238    
#>  actual_power          0.612  
#>  freq                  1      
#>  duration              4      
#>  measurement_occasions 5      
#>  polynomial_order      1      
#>  unstd_coefficient     -0.0219
#>  std_coefficient       -0.4   
#>  l1_error_var          0.0262 
#>  true_var_of_slopes    0.003  
#>  error_var_of_slopes   0.00262
#>  reliability           0.534  
#>  noncentral_t_parm     2.25   

# (4) Longer study, same number of occasions. Doubling the duration to
#     D = 8 while keeping M = 5 (so frequency = 0.5, one observation every
#     two years) raises the slope reliability to 0.82 and power to about
#     0.80. Spreading the same five occasions over a longer span sharply
#     increases power (Raudenbush & Liu, 2001, p. 393).
ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
             frequency = .5, duration = 8, N = 238)
#>  term                  value   
#>  specified_n_per_group 119     
#>  total_N               238     
#>  actual_power          0.795   
#>  freq                  0.5     
#>  duration              8       
#>  measurement_occasions 5       
#>  polynomial_order      1       
#>  unstd_coefficient     -0.0219 
#>  std_coefficient       -0.4    
#>  l1_error_var          0.0262  
#>  true_var_of_slopes    0.003   
#>  error_var_of_slopes   0.000655
#>  reliability           0.821   
#>  noncentral_t_parm     2.8     

# (5) More frequent sampling over a shorter span, same occasions. Halving
#     the duration to D = 2 while keeping M = 5 (so frequency = 2) drops
#     the slope reliability to 0.22 and power to about 0.31. Sampling more
#     often over a shorter study does little for power (Raudenbush & Liu,
#     2001, p. 393).
ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
             frequency = 2, duration = 2, N = 238)
#>  term                  value  
#>  specified_n_per_group 119    
#>  total_N               238    
#>  actual_power          0.305  
#>  freq                  2      
#>  duration              2      
#>  measurement_occasions 5      
#>  polynomial_order      1      
#>  unstd_coefficient     -0.0219
#>  std_coefficient       -0.4   
#>  l1_error_var          0.0262 
#>  true_var_of_slopes    0.003  
#>  error_var_of_slopes   0.0105 
#>  reliability           0.223  
#>  noncentral_t_parm     1.46   

# (6) One-sided test. A directional test of the base design places the
#     whole Type I error rate in the predicted tail, raising power from
#     0.61 to about 0.73.
ss_power_pcm(beta = -.4, tau = .003, level_1_variance = .0262,
             frequency = 1, duration = 4, N = 238, directional = TRUE)
#>  term                  value  
#>  specified_n_per_group 119    
#>  total_N               238    
#>  actual_power          0.727  
#>  freq                  1      
#>  duration              4      
#>  measurement_occasions 5      
#>  polynomial_order      1      
#>  unstd_coefficient     -0.0219
#>  std_coefficient       -0.4   
#>  l1_error_var          0.0262 
#>  true_var_of_slopes    0.003  
#>  error_var_of_slopes   0.00262
#>  reliability           0.534  
#>  noncentral_t_parm     2.25   

# (7) A higher-order trend. The same machinery plans power for the group
#     difference in any polynomial change coefficient. Here the target is
#     the quadratic trend (curvature / acceleration): with eight occasions
#     (frequency = 1, duration = 7), a between-subject quadratic-coefficient
#     variance tau = 0.002, level-one error variance 0.05, and a
#     standardized quadratic difference of 0.45, the design is planned for
#     0.80 power. A quadratic trend needs at least three occasions; a cubic
#     at least four (trend = "cubic" or trend = 3), and so on.
ss_power_pcm(beta = 0.45, tau = 0.002, level_1_variance = 0.05,
             frequency = 1, duration = 7, desired_power = .80,
             trend = "quadratic")
#>  term                  value  
#>  necessary_n_per_group 125    
#>  total_N               250    
#>  actual_power          0.801  
#>  freq                  1      
#>  duration              7      
#>  measurement_occasions 8      
#>  polynomial_order      2      
#>  unstd_coefficient     0.0201 
#>  std_coefficient       0.45   
#>  l1_error_var          0.05   
#>  true_var_of_slopes    0.002  
#>  error_var_of_slopes   0.00119
#>  reliability           0.627  
#>  noncentral_t_parm     2.82   
```
