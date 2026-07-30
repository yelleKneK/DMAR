# Sample Size for a Targeted Regression Coefficient

Determine the necessary sample size for a targeted regression
coefficient or determine the degree of power given a specified sample
size

## Usage

``` r
ss_power_reg_coef(
  rho2_Y_X = NULL,
  rho2_Y_X_without_j = NULL,
  p = NULL,
  desired_power = 0.85,
  alpha_level = 0.05,
  directional = FALSE,
  beta_j = NULL,
  sigma_X = NULL,
  sigma_Y = NULL,
  rho2_j_X_without_j = NULL,
  rho_XX = NULL,
  rho_YX = NULL,
  which_predictor = NULL,
  cohen_f2 = NULL,
  specified_N = NULL,
  print_progress = FALSE
)
```

## Arguments

- rho2_Y_X:

  Population squared multiple correlation coefficient predicting the
  dependent variable (i.e., *Y*) from the `p` predictor variables (i.e.,
  the *X* variables)

- rho2_Y_X_without_j:

  Population squared multiple correlation coefficient predicting the
  dependent variable (i.e., *Y*) from the `p`-1 predictor variables,
  where the one not used is the predictor of interest

- p:

  Number of predictor variables

- desired_power:

  Desired degree of statistical power for the test of targeted
  regression coefficient

- alpha_level:

  Type I error rate

- directional:

  Whether or not a direction or a nondirectional test is to be used
  (usually `directional=FALSE`)

- beta_j:

  Population value of the regression coefficient for the predictor of
  interest

- sigma_X:

  Population standard deviation for the predictor variable of interest

- sigma_Y:

  Population standard deviation for the outcome variable

- rho2_j_X_without_j:

  Population squared multiple correlation coefficient predicting the
  predictor variable of interest from the remaining p-1 predictor
  variables

- rho_XX:

  Population correlation matrix for the `p` predictor variables

- rho_YX:

  Population vector of correlation coefficient between the `p` predictor
  variables and the criterion variable

- which_predictor:

  Identifies the predictor of interest when `rho_XX` and `rho_YX` are
  specified

- cohen_f2:

  Cohen's (1988) definition for an effect size for a targeted regression
  coefficient: `(rho2_Y_X-rho2_Y_X_without_j)/(1-rho2_Y_X)`

- specified_N:

  Sample size for which power should be evaluated. This is the *total*
  sample size.

- print_progress:

  If the progress of the iterative procedure is printed to the screen as
  the iterations are occurring

## Value

- ss:

  Either the necessary sample size or the specified sample size,
  depending if one is interested in determining the necessary sample
  size given a desired degree of statistical power or if one is
  interested in the determining the value of statistical power given a
  specified sample size, respectively

- actual_power:

  Actual power of the situation described

- noncentral_t_parm:

  Value of the noncentral distribution for the appropriate
  *t*-distribution

- effect_size:

  Effect size for the noncentral *t*-distribution; this is the square
  root of `cohen_f2`, because `cohen_f2` is the effect size using an
  *F*-distribution

## Details

Determines the necessary sample size given a desired level of
statistical power. Alternatively, determines the statistical power for a
given a specified sample size. There are a number of ways that the
specification regarding the size of the regression coefficient can be
entered. The most basic, and often the simplest, is to specify
`rho2_Y_X` and `rho2_Y_X_without_j`. See the examples section for
several options.

Power is computed from a noncentral *t* distribution with noncentrality
\\\sqrt{N}\\f\\, which treats the predictors as fixed (their values held
constant across hypothetical replications). This is the standard
fixed-predictor power analysis; under random predictors, where the
predictor values themselves vary from sample to sample, the sample size
required for a given level of power is somewhat larger.

## References

Cohen, J. (1988). *Statistical power analysis for the behavioral
sciences* (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.
[doi:10.1037/1082-989X.8.3.305](https://doi.org/10.1037/1082-989X.8.3.305)

Kelley, K., & Maxwell, S. E. (2008). Sample size planning with
applications to multiple regression: Power and accuracy for omnibus and
targeted effects. In P. Alasuutari, L. Bickman, & J. Brannen (Eds.),
*The Sage handbook of social research methods* (pp. 166–192). Sage.

Maxwell, S. E. (2000). Sample size and multiple regression analysis.
*Psychological Methods, 5*(4), 434–458.
[doi:10.1037/1082-989X.5.4.434](https://doi.org/10.1037/1082-989X.5.4.434)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 4 on individual comparisons of means and
Chapter 6 on trend analysis.)

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
[doi:10.1146/annurev.psych.59.103006.093735](https://doi.org/10.1146/annurev.psych.59.103006.093735)

Anderson, S. F., Kelley, K., & Maxwell, S. E. (2017). Sample-size
planning for more accurate statistical power: A method adjusting sample
effect sizes for publication bias and uncertainty. *Psychological
Science, 28*(11), 1547–1562.
[doi:10.1177/0956797617723724](https://doi.org/10.1177/0956797617723724)

## See also

[`ss_aipe_reg_coef`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
[`ss_power_R2`](https://yelleknek.github.io/DMAR/reference/ss_power_R2.md),
[`conf_limits_ncf`](https://yelleknek.github.io/DMAR/reference/conf_limits_ncf.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

Other sample size for power:
[`power_fisher_exact()`](https://yelleknek.github.io/DMAR/reference/power_fisher_exact.md),
[`ss_aipe_mixed_effects()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_mixed_effects.md),
[`ss_aipe_tost_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_tost_smd.md),
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
[`ss_power_pcm()`](https://yelleknek.github.io/DMAR/reference/ss_power_pcm.md),
[`ss_power_r()`](https://yelleknek.github.io/DMAR/reference/ss_power_r.md),
[`ss_power_rc()`](https://yelleknek.github.io/DMAR/reference/ss_power_rc.md),
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
Cor.Mat <- rbind(
  c(1.00, 0.53, 0.58, 0.60, 0.46, 0.66),
  c(0.53, 1.00, 0.35, 0.07, 0.14, 0.43),
  c(0.58, 0.35, 1.00, 0.18, 0.29, 0.50),
  c(0.60, 0.07, 0.18, 1.00, 0.30, 0.26),
  c(0.46, 0.14, 0.29, 0.30, 1.00, 0.30),
  c(0.66, 0.43, 0.50, 0.26, 0.30, 1.00)
)

rho_XX <- Cor.Mat[2:6, 2:6]
rho_YX <- Cor.Mat[1, 2:6]

# Method 1
ss_power_reg_coef(rho2_Y_X = 0.7826786, rho2_Y_X_without_j = 0.7363697, p = 5,
                  alpha_level = .05, directional = FALSE, desired_power = .80)
#>  term              value
#>  necessary_N       40   
#>  actual_power      0.81 
#>  noncentral_t_parm 2.92 
#>  effect_size       0.462

# Method 2
ss_power_reg_coef(alpha_level = .05, rho_XX = rho_XX, rho_YX = rho_YX, which_predictor = 5,
                  directional = FALSE, desired_power = .80)
#>  term              value
#>  necessary_N       40   
#>  actual_power      0.81 
#>  noncentral_t_parm 2.92 
#>  effect_size       0.462

# Method 3
# Here, beta_j is the standardized regression coefficient. Had beta_j
# been the unstandardized regression coefficient, sigma_X and sigma_Y
# would have been the standard deviation for the X variable of
# interest and Y, respectively.
ss_power_reg_coef(rho2_Y_X = 0.7826786, rho2_j_X_without_j = 0.3652136, beta_j = 0.2700964,
                  p = 5, alpha_level = .05, sigma_X = 1, sigma_Y = 1, directional = FALSE,
                  desired_power = .80)
#>  term              value
#>  necessary_N       40   
#>  actual_power      0.81 
#>  noncentral_t_parm 2.92 
#>  effect_size       0.462

# Method 4
ss_power_reg_coef(alpha_level = .05, cohen_f2 = 0.2130898, p = 5,
                  directional = FALSE, desired_power = .80)
#>  term              value
#>  necessary_N       40   
#>  actual_power      0.81 
#>  noncentral_t_parm 2.92 
#>  effect_size       0.462

# Power given a specified N and squared multiple correlation coefficients.
ss_power_reg_coef(rho2_Y_X = 0.7826786, rho2_Y_X_without_j = 0.7363697, specified_N = 25,
                  p = 5, alpha_level = .05, directional = FALSE)
#>  term              value
#>  specified_N       25   
#>  actual_power      0.591
#>  noncentral_t_parm 2.31 
#>  effect_size       0.462

# Power given a specified N and effect size.
ss_power_reg_coef(alpha_level = .05, cohen_f2 = 0.2130898, p = 5, specified_N = 25,
                  directional = FALSE)
#>  term              value
#>  specified_N       25   
#>  actual_power      0.591
#>  noncentral_t_parm 2.31 
#>  effect_size       0.462

# Reproducing Maxwell's (2000, p. 445) Example
Cor.Mat.Maxwell <- rbind(
  c(1.00, 0.35, 0.20, 0.20, 0.20, 0.20),
  c(0.35, 1.00, 0.40, 0.40, 0.40, 0.40),
  c(0.20, 0.40, 1.00, 0.45, 0.45, 0.45),
  c(0.20, 0.40, 0.45, 1.00, 0.45, 0.45),
  c(0.20, 0.40, 0.45, 0.45, 1.00, 0.45),
  c(0.20, 0.40, 0.45, 0.45, 0.45, 1.00)
)

RHO.XX.Maxwell <- Cor.Mat.Maxwell[2:6, 2:6]
Rho.YX.Maxwell <- Cor.Mat.Maxwell[1, 2:6]
R2.Maxwell <- Rho.YX.Maxwell %*% solve(RHO.XX.Maxwell) %*% Rho.YX.Maxwell

RHO.XX.Maxwell.no.1 <- Cor.Mat.Maxwell[3:6, 3:6]
Rho.YX.Maxwell.no.1 <- Cor.Mat.Maxwell[1, 3:6]
R2.Maxwell.no.1 <-
  Rho.YX.Maxwell.no.1 %*% solve(RHO.XX.Maxwell.no.1) %*% Rho.YX.Maxwell.no.1

# This procedure arrives at N = 111, whereas Maxwell (2000, p. 445)
# reports N = 113. The two differ because of the noncentrality
# parameterization, not rounding: this function uses the fixed-predictor
# noncentrality sqrt(N) * f (see Details), while the tabled value rests on
# Cohen's (1988) convention. Neither is a random-predictor result; under
# random predictors, where the predictor values vary across replications,
# the sample size needed for the same power is larger still.
ss_power_reg_coef(rho2_Y_X = R2.Maxwell, rho2_Y_X_without_j = R2.Maxwell.no.1, p = 5,
                  alpha_level = .05, directional = FALSE, desired_power = .80)
#>  term              value
#>  necessary_N       111  
#>  actual_power      0.801
#>  noncentral_t_parm 2.83 
#>  effect_size       0.269
```
