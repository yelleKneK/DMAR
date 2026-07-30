# Print a Model Summary With DMAR *p*-value Formatting

Pretty-print a model summary (the output of `summary.lm`, `summary.glm`,
or `summary` on an lme4 or lmerTest fit) with *p*-values formatted at a
fixed number of decimal places (default 4) and with a “\<
10^(-digits_p)” floor for values too small to express. The default
`print.summary.lm` / `print.summary.merMod` routes *p*-values through
`stats::format.pval`, which applies its own digit rule and switches to
scientific notation for tiny values. `print_summary()` sidesteps that by
converting the *p*-value columns to character strings up front and
printing as a data frame.

## Usage

``` r
print_summary(fit, digits_p = 4L)
```

## Arguments

- fit:

  A fitted model object with a `summary` method that returns
  coefficients via `coef(summary(fit))`, including a `Pr(...)` column.
  Tested with `lm`, `glm`,
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html), and
  [`lmerTest::lmer`](https://rdrr.io/pkg/lmerTest/man/lmer.html).

- digits_p:

  Integer number of decimal places for the *p*-value column(s). Default
  `4L`.

## Value

The model summary, invisibly and unchanged.

## Details

For a linear model, the function prints the coefficient table, the
residual standard error and degrees of freedom, the multiple and
adjusted \\R^2\\, and the omnibus *F* test and its *p*-value. For a
mixed-effects model fit through lme4 / lmerTest, the function prints the
random-effect variances (from
[`lme4::VarCorr`](https://rdrr.io/pkg/nlme/man/VarCorr.html)) and the
fixed-effect coefficient table.

The returned object is the model summary, invisibly and unchanged: the
underlying numeric *p*-values retain full precision and can still be
indexed (for example as `coef(summary(fit))[, "Pr(>|t|)"]`).

## See also

[`format_p`](https://yelleknek.github.io/DMAR/reference/format_p.md),
[`print_anova`](https://yelleknek.github.io/DMAR/reference/print_anova.md).

## Author

Ken Kelley

## Examples

``` r
fit_lm <- lm(weight ~ Time + Diet, data = ChickWeight)
print_summary(fit_lm)
#> Coefficients:
#>              Estimate Std. Error   t value Pr(>|t|)
#> (Intercept) 10.924391  3.3606567  3.250672   0.0012
#> Time         8.750492  0.2218052 39.451248 < 0.0001
#> Diet2       16.166074  4.0858416  3.956608 < 0.0001
#> Diet3       36.499407  4.0858416  8.933143 < 0.0001
#> Diet4       30.233456  4.1074850  7.360576 < 0.0001
#> 
#> Residual standard error: 35.99 on 573 degrees of freedom
#> Multiple R-squared: 0.7453,  Adjusted R-squared: 0.7435
#> F-statistic: 419.2 on 4 and 573 DF, p-value: < 0.0001

fit_lmer <- lme4::lmer(weight ~ Time + (1 | Chick), data = ChickWeight)
print_summary(fit_lmer)
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  Chick    (Intercept) 26.793  
#>  Residual             28.274  
#> 
#> Fixed effects:
#>              Estimate Std. Error   t value
#> (Intercept) 27.845104  4.3876736  6.346211
#> Time         8.726062  0.1755185 49.715925

# Underlying numeric p-values are untouched:
sm <- summary(fit_lm)
sm$coefficients[, "Pr(>|t|)"]   # full-precision doubles
#>   (Intercept)          Time         Diet2         Diet3         Diet4 
#>  1.218886e-03 1.803038e-165  8.556049e-05  5.628378e-18  6.391748e-13 
```
