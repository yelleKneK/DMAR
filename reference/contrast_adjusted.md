# Confidence Interval for a Contrast of Covariate-Adjusted Cell Means in a Factorial ANCOVA

Given a fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
[`aov`](https://rdrr.io/r/stats/aov.html) object for a factorial
analysis of covariance (one or more crossed factors plus one or more
covariates) and a numeric contrast vector over the cells of the
factorial design, `contrast_adjusted()` forms the contrast of the
covariate-adjusted cell means, \\\hat{\psi} = \sum_j c_j \\
\hat{\bar{Y}}\_j\\, where each \\\hat{\bar{Y}}\_j\\ is the model's
predicted mean for cell \\j\\ evaluated at the mean of every covariate
(the adjusted, or least-squares, cell mean). It returns the point
estimate, a *t* confidence interval on the model's residual degrees of
freedom, and the accompanying *t* statistic and two-sided *p*-value for
\\H_0\\: \psi = 0\\.

## Usage

``` r
contrast_adjusted(model, contrast, conf_level = 0.95)
```

## Arguments

- model:

  A fitted [`lm`](https://rdrr.io/r/stats/lm.html) or
  [`aov`](https://rdrr.io/r/stats/aov.html) object for a factorial
  ANCOVA: one or more crossed factors and one or more numeric covariates
  on the right-hand side of the formula.

- contrast:

  A numeric vector of contrast weights, one weight per cell of the
  factorial design (the crossing of the model's factors). Its length
  must equal the number of cells. The weights typically sum to zero.

- conf_level:

  The confidence level for the interval (default `0.95`).

## Value

A five-row `dmar_tbl` (a `data.frame` with columns `term` and `value`).
The `term` values are `"contrast"` (the point estimate \\\hat{\psi}\\ of
the contrast of adjusted cell means), `"lower_limit"` and
`"upper_limit"` (the confidence limits), `"t"` (the *t* statistic), and
`"p"` (the two-sided *p*-value). The stored `value` column is numeric at
full precision.

## Details

The adjusted cell means are the means the ANCOVA actually tests: the
predicted outcome for each combination of factor levels, holding every
covariate at its sample mean. Writing \\L\\ for the linear map that
sends the model coefficients to that contrast of adjusted means (built
by evaluating the model's design matrix at each cell with the covariates
set to their means and combining the rows with the contrast weights),
the point estimate is \\\hat{\psi} = L' \hat{\beta}\\ and its standard
error is the square root of the quadratic form \\L' \\
\mathrm{vcov}(\hat{\beta}) \\ L\\. The interval is \\\hat{\psi} \pm
t\_{1 - \alpha/2,\\ \nu}\\ \mathrm{SE}\\, with \\\nu\\ the residual
degrees of freedom of the fitted model.

Because \\L\\ is read off the fitted model's own design matrix, the
function is agnostic to how the factors are parameterized: a cell-means
parameterization (`y ~ 0 + cell + x`) and the crossed-factor
parameterization (`y ~ A * B + x`) give the same contrast estimate and
standard error, provided the contrast vector is ordered to match the
cells of the reference grid (see `Note`).

The *t* interval returned here has exact per-comparison coverage for a
single contrast chosen in advance. For a family of contrasts examined
together, adjust the critical value for multiplicity (for example the
Scheffe critical value for the full cell space,
[`cv_scheffe`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
or a Bryant–Paulson simultaneous interval,
[`ci_c_ancova_bp`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md)).

## Note

The contrast weights are matched to the cells of the reference grid,
which is the crossing of the model's factors in the order the factors
appear in the model formula, with the first factor varying fastest (the
order [`expand.grid`](https://rdrr.io/r/base/expand.grid.html) produces
over the factor levels). For a single factor this is simply the order of
its levels. When in doubt, fit the cell-means form
`y ~ 0 + cell + covariates` with `cell = interaction(A, B, ...)` and
order the weights to match `levels(cell)`.

## References

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (See Chapter 9 on designs with covariates.)

## See also

[`contrast_test`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
for contrasts of unadjusted group means in a one-way design;
[`ci_c_ancova`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md)
for a single-covariate ANCOVA contrast from summary statistics.

Other confidence intervals for effect sizes:
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ci_c()`](https://yelleknek.github.io/DMAR/reference/ci_c.md),
[`ci_c_ancova()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova.md),
[`ci_c_ancova_bp()`](https://yelleknek.github.io/DMAR/reference/ci_c_ancova_bp.md),
[`ci_correlation`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_cv()`](https://yelleknek.github.io/DMAR/reference/ci_cv.md),
[`ci_eta_squared()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared.md),
[`ci_eta_squared_generalized()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_generalized.md),
[`ci_eta_squared_partial()`](https://yelleknek.github.io/DMAR/reference/ci_eta_squared_partial.md),
[`ci_mahalanobis()`](https://yelleknek.github.io/DMAR/reference/ci_mahalanobis.md),
[`ci_omega_squared()`](https://yelleknek.github.io/DMAR/reference/ci_omega_squared.md),
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md),
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
[`plot_smd()`](https://yelleknek.github.io/DMAR/reference/plot_smd.md)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# A 2 x 2 factorial ANCOVA with one covariate.
set.seed(113)
d <- data.frame(
  A = factor(rep(c("a1", "a2"), each = 40)),
  B = factor(rep(rep(c("b1", "b2"), each = 20), 2)),
  x = rnorm(80)
)
d$y <- 5 + 2 * (d$A == "a2") + 1.5 * (d$B == "b2") +
  0.8 * d$x + rnorm(80)
fit <- lm(y ~ A * B + x, data = d)

# Cells in reference-grid order: (a1,b1), (a2,b1), (a1,b2), (a2,b2).
# Main effect of A, averaged over B: mean(a2 cells) - mean(a1 cells).
contrast_adjusted(fit, contrast = c(-0.5, 0.5, -0.5, 0.5))
#>  term        value   
#>  contrast    1.9     
#>  lower_limit 1.49    
#>  upper_limit 2.31    
#>  t           9.18    
#>  p           < 0.0001
#> 
#> Confidence level: 95%
```
