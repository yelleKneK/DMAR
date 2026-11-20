# Structural Relations Estimated After the Measurement Models

Fits local structural-after-measurement (SAM) models using
[`lavaan::sam()`](https://rdrr.io/pkg/lavaan/man/sam.html). Each
measurement block is estimated separately, and the structural relations
use a covariance matrix corrected for measurement error. Two-step
standard errors account for estimating the measurement parameters. This
offers an alternative to treating estimated factor scores as perfectly
measured observations. A latent SEM remains available when this
wrapper's narrower scope does not fit the study.

## Usage

``` r
sem_sam(
  model,
  data,
  estimator = "ML",
  conf_level = 0.95,
  gamma_unbiased = TRUE
)
```

## Arguments

- model:

  A single-group lavaan model string with at least two simple-structure
  factors, at least three distinct continuous indicators per factor, and
  structural regressions among those factors. Indicators cannot
  cross-load. Derived parameters, explicit equality constraints,
  observed covariates, residual covariances, feedback paths, and latent
  interactions are outside this interface's validated scope.

- data:

  Raw numeric observations on the model's indicators. Incomplete
  observations are refused; the function never drops cases
  automatically. Use a joint latent SEM with FIML for incomplete data.
  An indicator of `model` that `data` lacks is refused with an error of
  class `dmar_items_invalid` (see
  [`dmar_metadata`](https://yelleknek.github.io/DMAR/reference/dmar_metadata.md)),
  as are columns of `data` that share a name.

- estimator:

  `"ML"` (default) with normal-theory `"twostep"` errors, or `"MLR"`
  with `"twostep.robust"` errors. Both use local ML mapping and the ML
  measurement estimators; the covariance assumptions for inference
  differ. Numerical optimization versus an analytic solution does not
  determine which SE is appropriate.

- conf_level:

  Confidence level, default 0.95.

- gamma_unbiased:

  Logical, default `TRUE`. For robust inference, use lavaan's
  finite-sample unbiased moment-covariance correction. `FALSE` uses the
  empirical moment covariance (HC0 convention).

## Value

A `dmar_tbl` with parameter labels, `estimate`, `se`, `statistic`,
`p_value`, confidence limits, `conf_level`, and `stage`. Attributes
retain the `fit`, requested and actual `method`, and the covariance
limitation. `tidy()` retains the rows and `glance()` gives the method
and sample size.

## Details

This interface currently supports unweighted, single-level, complete
continuous data. It verifies the actual local mapping and two-step SE
method after fitting and refuses a substitution. These restrictions
distinguish numerical agreement with independent complete-data oracles
from evidence for other models or missing-data uncertainty calculations.
Small or moderate samples can benefit from separating measurement and
structural estimation (Rosseel & Loh, 2024), but SAM does not repair a
misspecified or unidentified measurement model.

Returned intervals concern individual, original-unit parameters. The
underlying fit's joint covariance matrix does not populate every
cross-stage covariance. Consequently this function refuses user-defined
parameters and does not report standardized structural intervals or
intervals for effects combining measurement and structural estimates.
The attached fit is retained for inspection with that limitation, not as
a complete covariance matrix for arbitrary delta-method calculations.

## References

Rosseel, Y., & Loh, W. W. (2024). A structural after measurement (SAM)
approach to structural equation modeling. *Psychological Methods,
29*(3), 561–588.
[doi:10.1037/met0000503](https://doi.org/10.1037/met0000503)

## Author

Ken Kelley <kkelley@nd.edu>

## Examples

``` r
model <- "visual =~ t1_visual_perception + t2_cubes + t4_lozenges
          verbal =~ t6_paragraph_comprehension + t7_sentence + t9_word_meaning
          verbal ~ visual"
sem_sam(model, holzinger_swineford)
#>  term                                                    
#>  visual =~ t1_visual_perception                          
#>  visual =~ t2_cubes                                      
#>  visual =~ t4_lozenges                                   
#>  verbal =~ t6_paragraph_comprehension                    
#>  verbal =~ t7_sentence                                   
#>  verbal =~ t9_word_meaning                               
#>  verbal ~ visual                                         
#>  t1_visual_perception ~~ t1_visual_perception            
#>  t2_cubes ~~ t2_cubes                                    
#>  t4_lozenges ~~ t4_lozenges                              
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension
#>  t7_sentence ~~ t7_sentence                              
#>  t9_word_meaning ~~ t9_word_meaning                      
#>  visual ~~ visual                                        
#>  verbal ~~ verbal                                        
#>  lhs                        op rhs                        estimate se   
#>  visual                     =~ t1_visual_perception       1.000    0.000
#>  visual                     =~ t2_cubes                   0.519    0.094
#>  visual                     =~ t4_lozenges                1.476    0.285
#>  verbal                     =~ t6_paragraph_comprehension 1.000    0.000
#>  verbal                     =~ t7_sentence                1.511    0.089
#>  verbal                     =~ t9_word_meaning            2.156    0.132
#>  verbal                     ~  visual                     0.264    0.058
#>  t1_visual_perception       ~~ t1_visual_perception       30.047   4.253
#>  t2_cubes                   ~~ t2_cubes                   17.039   1.674
#>  t4_lozenges                ~~ t4_lozenges                40.497   8.267
#>  t6_paragraph_comprehension ~~ t6_paragraph_comprehension 3.435    0.440
#>  t7_sentence                ~~ t7_sentence                6.658    0.946
#>  t9_word_meaning            ~~ t9_word_meaning            18.068   2.159
#>  visual                     ~~ visual                     18.854   4.688
#>  verbal                     ~~ verbal                     7.405    0.915
#>  statistic p_value  ci_lower ci_upper conf_level stage      
#>  <NA>      <NA>     1.000    1.000    0.950      measurement
#>  5.532     < 0.0001 0.335    0.702    0.950      measurement
#>  5.173     < 0.0001 0.917    2.036    0.950      measurement
#>  <NA>      <NA>     1.000    1.000    0.950      measurement
#>  16.906    < 0.0001 1.335    1.686    0.950      measurement
#>  16.391    < 0.0001 1.899    2.414    0.950      measurement
#>  4.560     < 0.0001 0.151    0.378    0.950      structural 
#>  7.064     < 0.0001 21.711   38.383   0.950      measurement
#>  10.177    < 0.0001 13.757   20.320   0.950      measurement
#>  4.899     < 0.0001 24.295   56.700   0.950      measurement
#>  7.805     < 0.0001 2.572    4.298    0.950      measurement
#>  7.038     < 0.0001 4.804    8.513    0.950      measurement
#>  8.367     < 0.0001 13.836   22.300   0.950      measurement
#>  4.021     < 0.0001 9.665    28.043   0.950      structural 
#>  8.089     < 0.0001 5.611    9.199    0.950      structural 
#> 
#> Confidence level: 95%
#> 
#> Local SAM: ML; SE: twostep; complete continuous data; gamma_unbiased = TRUE. Intervals concern individual original-unit parameters. Cross-stage covariance cells are incomplete; do not use the attached fit for arbitrary derived or standardized structural intervals.
```
