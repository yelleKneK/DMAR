# Tidy / Glance Methods for R2_mixed_effects Output

Returns the broom-style summary of the marginal and conditional \\R^2\\:
`term` (`"R2_marginal"` or `"R2_conditional"`), `estimate`, and, when a
bootstrap interval was requested, `ci_lower`, `ci_upper`, and
`conf_level`. `glance` returns the same information in one row per
quantity.

## Usage

``` r
# S3 method for class 'dmar_R2_mixed_effects'
tidy(x, ...)

# S3 method for class 'dmar_R2_mixed_effects'
glance(x, ...)
```

## Arguments

- x:

  A `dmar_R2_mixed_effects` object returned by
  [`R2_mixed_effects`](https://yelleknek.github.io/DMAR/reference/R2_mixed_effects.md).

- ...:

  Unused.

## Value

A `data.frame` in broom convention.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
fit <- lme4::lmer(Reaction ~ Days + (Days | Subject),
                  data = lme4::sleepstudy)
res <- R2_mixed_effects(fit)
generics::tidy(res)
#>             term  estimate
#> 1    R2_marginal 0.2786511
#> 2 R2_conditional 0.7992199
generics::glance(res)
#>             term  estimate
#> 1    R2_marginal 0.2786511
#> 2 R2_conditional 0.7992199
```
