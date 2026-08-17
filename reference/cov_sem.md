# Model Implied Covariance Matrix From a Lavaan-Specified SEM

Given a structural equation model written in lavaan model syntax with
all of its parameters fixed to their population values, compute the
model implied population covariance matrix \\\Sigma(\theta)\\ of the
observed variables and, when the model has a mean structure, the model
implied population mean vector \\\mu(\theta)\\. This function requires
lavaan to be installed.

This is the helper that drives the *population* side of the sample size
planning workflow for SEM: it lets the user state a population model,
obtain the \\\Sigma(\theta)\\ (and \\\mu(\theta)\\) those fixed values
imply, and then pass that population to
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md),
[`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md),
[`ss_aipe_rmsea_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea_sensitivity.md),
[`ss_power_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_power_composite_sem.md),
or
[`ss_aipe_composite_sem`](https://yelleknek.github.io/DMAR/reference/ss_aipe_composite_sem.md).

## Usage

``` r
cov_sem(model)
```

## Arguments

- model:

  A single character string giving a structural equation model in lavaan
  model syntax (see
  [`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)),
  with every parameter fixed to its population value. A factor loading,
  a structural path, a variance, or a covariance is fixed by prefixing
  the numeric value to the variable with the `*` operator, for example
  `"f1 =~ 1*y1 + 0.8*y2 + 0.8*y3"` for the loadings, `"f2 ~ 0.5*f1"` for
  a structural path, and `"y1 ~~ 0.5*y1"` for a residual variance. A
  model with a mean structure (for example a latent growth curve model)
  also fixes every intercept and latent mean, for example `"t1 ~ 0*1"`
  and `"s ~ 0.3*1"`. lavaan syntax embeds the values and knows which
  variables are latent, so no separate parameter vector or list of
  latent variables is needed.

## Value

A list with components:

- `sigma_theta`:

  The model implied population covariance matrix of the observed
  variables, with rows and columns named.

- `mu_theta`:

  The model implied population mean vector of the observed variables,
  named, in the row order of `sigma_theta`. A vector of zeros when the
  model has no mean structure.

- `observed_vars`:

  Character vector of observed variable names in the row/column order of
  `sigma_theta`.

## Details

The function builds a non-fitted lavaan object from `model` with all
parameters held at the population values written into the syntax, and
reads back the model implied covariance matrix of the observed
variables. Because the object is created with `do.fit = FALSE`, no
estimation is performed and the placeholder sample covariance lavaan
needs to construct the object is never used; the returned
\\\Sigma(\theta)\\ comes entirely from the fixed parameter values. The
observed-variable names are taken from the model syntax and fix the row
and column order of the returned matrix.

## References

Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
targeted effects in structural equation modeling: Sample size planning
for narrow confidence intervals. *Psychological Methods, 16*(2),
127–148. [doi:10.1037/a0021764](https://doi.org/10.1037/a0021764)

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

## See also

[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html),
[`ss_aipe_sem_path`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path.md),
[`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md),
[`ss_aipe_rmsea_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_rmsea_sensitivity.md),
[`covmat_from_cfa`](https://yelleknek.github.io/DMAR/reference/covmat_from_cfa.md).

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Population model with all parameters fixed to their values: two factors,
# three indicators each, and a structural path f2 ~ f1 of 0.5.
pop_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f2 ~ 0.5*f1
  f1 ~~ 1*f1
  f2 ~~ 0.75*f2
  y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
  y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
"
cov_sem(pop_model)$sigma_theta
#>     y1   y2   y3  y4   y5   y6
#> y1 1.5 0.80 0.80 0.5 0.40 0.40
#> y2 0.8 1.14 0.64 0.4 0.32 0.32
#> y3 0.8 0.64 1.14 0.4 0.32 0.32
#> y4 0.5 0.40 0.40 1.5 0.80 0.80
#> y5 0.4 0.32 0.32 0.8 1.14 0.64
#> y6 0.4 0.32 0.32 0.8 0.64 1.14

# A population model with a mean structure: a linear latent growth curve
# over four waves. The intercepts and latent means are fixed too, and
# mu_theta carries the model implied means (5.0, 5.3, 5.6, 5.9).
pop_lgm <- "
  i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
  i ~~ 1*i
  s ~~ 0.2*s
  i ~~ -0.15*s
  t1 ~~ 0.5*t1; t2 ~~ 0.5*t2; t3 ~~ 0.5*t3; t4 ~~ 0.5*t4
  t1 ~ 0*1; t2 ~ 0*1; t3 ~ 0*1; t4 ~ 0*1
  i ~ 5*1
  s ~ 0.3*1
"
cov_sem(pop_lgm)$mu_theta
#>  t1  t2  t3  t4 
#> 5.0 5.3 5.6 5.9 
```
