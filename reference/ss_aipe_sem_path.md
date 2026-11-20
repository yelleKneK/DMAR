# Sample Size Planning for SEM Targeted Effects

Plan sample size for structural equation models so that the confidence
interval for the targeted model parameter is sufficiently narrow

## Usage

``` r
ss_aipe_sem_path(
  model,
  true_cov_matrix,
  width,
  which_path,
  conf_level = 0.95,
  assurance = NULL,
  detail = FALSE,
  internal = FALSE,
  ...
)
```

## Arguments

- model:

  A single character string giving the free analysis model in lavaan
  model syntax (see
  [`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)).
  The target path must carry a parameter label so it can be referred to
  by name, for example `"f2 ~ b*f1"` labels the structural path `b`.
  This is the model that would be fit to the data; its parameters are
  free, not fixed to population values

- true_cov_matrix:

  Population covariance matrix of the observed variables planned for,
  with row and column names matching the observed variables in `model`.
  It is typically obtained from a fully fixed population model via
  [`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md)

- width:

  Desired confidence interval width for the model parameter of interest

- which_path:

  The parameter label of the targeted path, given as a character string,
  for example `"b"` for the path labeled `f2 ~ b*f1` in `model`

- conf_level:

  Confidence level (i.e., 1 - Type I error rate)

- assurance:

  The assurance that the confidence interval obtained in a particular
  study will be no wider than desired (must be `NULL` or a value between
  0.50 and 1)

- detail:

  if `TRUE`, additionally print the model parameter names and the
  observed variable names (the returned table is unchanged)

- internal:

  option to output a list for internal use (for
  ss_aipe_sem_path_sensitivity)

- ...:

  Further arguments passed to
  [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) for the fit to
  `true_cov_matrix`, from which the plan, with or without `assurance`,
  is read (for example `std.lv = TRUE`). Every name given here must be
  an argument or an option of lavaan (see
  [`lavOptions`](https://rdrr.io/pkg/lavaan/man/lavOptions.html)); any
  other name, and any argument without a name, is refused before a model
  is fitted, with a condition of class `dmar_argument_invalid` that
  gives the name as it was typed

## Value

A `data.frame` (a `dmar_tbl`) with `term` and `value` columns whose rows
are `necessary_N` (the planned sample size), `path_index` (the position
of the target path among the model parameters), and `var_theta_j` (the
population sampling variance of the target path at the planned sample
size). The returned table is the same whether or not `detail = TRUE`.
When `internal = TRUE` a list is returned for use by
[`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md).

## Details

This function implements the sample size planning methods proposed in
Lai and Kelley (2011). It requires lavaan to be installed and uses
[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) to obtain the expected
information, that is the asymptotic covariance matrix of the parameter
estimates, by fitting the free analysis model to the population
covariance matrix `true_cov_matrix` at a very large sample size. The
analysis model is written in lavaan model syntax with the targeted path
given a parameter label; see
[`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html) for
the syntax and [`sem`](https://rdrr.io/pkg/lavaan/man/sem.html) for the
fitting machinery. The population covariance matrix `true_cov_matrix` is
most naturally produced by
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md) from
a fully fixed population model. An analysis model whose likelihood has
no maximum at `true_cov_matrix`, an error variance running off below
zero while the likelihood is flat or still rising (see
[`fa_bifactor`](https://yelleknek.github.io/DMAR/reference/fa_bifactor.md)),
stops the function with an error of class `dmar_no_maximum` that prints
the lavaan code of the fit, since the asymptotic covariance matrix would
describe where the optimizer stopped rather than the model.

When `assurance` is supplied, the assurance adjustment is based on a chi
square approximation to the sampling variability of the confidence
interval width and can undershoot the nominal assurance in finite
samples; use
[`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md)
to check the realized width and coverage at the planned sample size.

## References

Lai, K., & Kelley, K. (2011). Accuracy in parameter estimation for
targeted effects in structural equation modeling: Sample size planning
for narrow confidence intervals. *Psychological Methods, 16*(2),
127–148. [doi:10.1037/a0021764](https://doi.org/10.1037/a0021764)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Rosseel, Y. (2012). lavaan: An R package for structural equation
modeling. *Journal of Statistical Software, 48*(2), 1–36.
[doi:10.18637/jss.v048.i02](https://doi.org/10.18637/jss.v048.i02)

## See also

[`sem`](https://rdrr.io/pkg/lavaan/man/sem.html),
[`model.syntax`](https://rdrr.io/pkg/lavaan/man/model.syntax.html),
[`cov_sem`](https://yelleknek.github.io/DMAR/reference/cov_sem.md),
[`ss_aipe_sem_path_sensitivity`](https://yelleknek.github.io/DMAR/reference/ss_aipe_sem_path_sensitivity.md)

[`design_consequences`](https://yelleknek.github.io/DMAR/reference/design_consequences.md)
for what a chosen design delivers: power, the Type S (sign) and Type M
(exaggeration) errors of the significance filter, and the expected
confidence interval width.

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
# Population covariance from a fully fixed model (see cov_sem()).
pop_model <- "
  f1 =~ 1*y1 + 0.8*y2 + 0.8*y3
  f2 =~ 1*y4 + 0.8*y5 + 0.8*y6
  f2 ~ 0.5*f1
  f1 ~~ 1*f1
  f2 ~~ 0.75*f2
  y1 ~~ 0.5*y1; y2 ~~ 0.5*y2; y3 ~~ 0.5*y3
  y4 ~~ 0.5*y4; y5 ~~ 0.5*y5; y6 ~~ 0.5*y6
"
Sigma <- cov_sem(pop_model)$sigma_theta

# Free analysis model with the target structural path labeled "b".
analysis_model <- "
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f2 ~ b*f1
"
ss_aipe_sem_path(model = analysis_model, true_cov_matrix = Sigma,
                 width = 0.30, which_path = "b")
#>  term        value
#>  necessary_N 264  
#>  path_index  5    
#>  var_theta_j 0.006
#> 
#> Confidence level: 95%

# That sample size holds the interval to the desired width on average, so
# about half of the studies it plans return a wider one. Adding assurance
# plans for the width to be met in 90 percent of studies instead, at the
# cost of a larger sample size.
ss_aipe_sem_path(model = analysis_model, true_cov_matrix = Sigma,
                 width = 0.30, which_path = "b",
                 assurance = 0.90)
#>  term        value
#>  necessary_N 293  
#>  path_index  5    
#>  var_theta_j 0.005
#> 
#> Confidence level: 95%
```
