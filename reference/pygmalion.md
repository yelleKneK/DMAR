# Pygmalion in the Classroom Teacher-Expectancy Data

The teacher-expectancy data from Rosenthal and Jacobson's (1968)
*Pygmalion in the Classroom*, the study that introduced the "Pygmalion
effect": the hypothesis that a teacher's expectations can become a
self-fulfilling prophecy for a pupil's intellectual growth.
Intelligence-test scores were obtained for *N* = 310 elementary school
children in grades 1 through 6, of whom *n* = 64 were randomly
designated to their teachers as likely "intellectual bloomers" while the
remaining *n* = 246 served as controls. The data set is a classic
benchmark for the analysis of covariance (ANCOVA) and, in particular,
for ANCOVA with *heterogeneity of regression*: it is the running example
for that topic in Maxwell, Delaney, and Kelley, *Designing Experiments
and Analyzing Data: A Model Comparison Perspective* (Routledge), where
it appears as a Chapter 9 example (and as a Chapter 3 exercise).

## Usage

``` r
pygmalion
```

## Format

A data frame with 310 observations on 6 variables.

- `grade`:

  Grade in school at the start of the study, an integer from 1 to 6.

- `treatment`:

  Factor with levels `Control` (reference, *n* = 246) and `Bloomer` (*n*
  = 64). The `Bloomer` children were a randomly selected ~20% of each
  classroom whose teachers were told, on the basis of a fictitious test
  purportedly predicting intellectual blooming, that they were likely to
  show unusual gains during the year; the `Control` children were not
  singled out. In the AMCP source this variable is coded `1` = Bloomer,
  `0` = Control.

- `iq_pre`:

  Pretest total IQ, measured before the expectancy manipulation. The
  covariate in the analysis of covariance.

- `iq_4`:

  Total IQ at an intermediate follow-up assessment.

- `iq_8`:

  Total IQ at the end-of-study follow-up assessment. This is the
  dependent variable in the book's Chapter 9 analysis of covariance.

- `iq_gain`:

  Total IQ change from pretest to the end-of-study assessment, equal to
  `iq_8 - iq_pre`.

## Source

Rosenthal, R., & Jacobson, L. (1968). *Pygmalion in the classroom:
Teacher expectation and pupils' intellectual development*. Holt,
Rinehart and Winston.

Distributed with the AMCP data companion to Maxwell, Delaney, and Kelley
(see References) as `chapter_9_exercise_15`.

## Details

**The study.** Robert Rosenthal (Harvard University) and Lenore Jacobson
(principal of an elementary school in South San Francisco referred to as
"Oak School") set out to test experimentally whether teacher
expectations influence pupil achievement. At the start of the school
year all children were given a standardized test of general ability,
described to teachers as the "Harvard Test of Inflected Acquisition," a
test said to identify children poised for an intellectual growth spurt.
In reality the instrument was Flanagan's Tests of General Ability (TOGA)
and the children identified as likely "bloomers" were chosen *at
random*, about one in five per classroom. The only experimental
manipulation was the expectation planted in the teachers' minds.
Children were re-tested over the following year(s), and the question was
whether the randomly labeled bloomers would out-gain their controls in
measured IQ. Rosenthal and Jacobson reported that they did, most
strongly in the earliest grades, and interpreted the difference as
evidence that teacher expectations operate as a self-fulfilling
prophecy. The study became one of the most famous and most debated
experiments in the social sciences; subsequent critiques (e.g.,
Thorndike, 1968) questioned the reliability of the TOGA at the extremes
of the score range for the youngest children, which is itself part of
why the data are instructive for teaching careful analysis.

**Why it is a benchmark for heterogeneity of regression.** A standard
ANCOVA adjusts the group comparison for the pretest covariate under the
assumption that the regression of the outcome on the covariate has the
*same* slope in every group (homogeneity of regression). In these data
that assumption is questionable: the within-group regression of `iq_8`
on `iq_pre` is steeper for the bloomers than for the controls, so the
estimated treatment effect depends on the covariate value at which it is
evaluated. This makes the data an ideal teaching example for (a) testing
the homogeneity-of-regression assumption, (b) interpreting a
treatment-by-covariate interaction, and (c) estimating the treatment
effect, and its sampling variance, *at chosen covariate values* rather
than only at the grand mean.

**Reproducible quantities.** Fitting the separate-slopes model
`lm(iq_8 ~ iq_pre * treatment)` gives a within-group slope of
\\0.77799\\ for the controls and \\0.96894\\ for the bloomers (reported
as \\0.96895\\ in `MBESS::var.ete`, a fifth-decimal rounding
difference). The pooled within-group residual variance is \\\hat\sigma^2
= 175.3251\\ on 306 degrees of freedom, and the sample variance of the
covariate is \\348.91\\. These are exactly the inputs used in the worked
example for the variance of the estimated treatment effect at selected
covariate values under heterogeneity of regression (Li, McLouth, and
Delaney; see `MBESS::var.ete`).

**Relationship to the AMCP package.** The same numeric data ship with
the book's data companion, the AMCP package, as `chapter_9_exercise_15`
and `chapter_9_extension_exercise_3` (with `IQGain`) and
`chapter_3_exercise_22` (without it). The version here renames the
columns to DMAR's descriptive snake_case style and labels the
experimental condition as a factor; no measured value has been altered.
See `data-raw/pygmalion.R` for the construction script and its
verification checks.

## References

Rosenthal, R., & Jacobson, L. (1968). *Pygmalion in the classroom:
Teacher expectation and pupils' intellectual development*. Holt,
Rinehart and Winston.

Rosenthal, R., & Jacobson, L. (1968). Pygmalion in the classroom. *The
Urban Review, 3*(1), 16–20.
[doi:10.1007/BF02322211](https://doi.org/10.1007/BF02322211)

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge. (Heterogeneity-of-regression ANCOVA example, Chapter
9.)

Thorndike, R. L. (1968). Review of *Pygmalion in the Classroom*.
*American Educational Research Journal, 5*(4), 708–711.

## See also

[`ancova`](https://yelleknek.github.io/DMAR/reference/ancova.md) for an
ANCOVA that returns adjusted means, effect size confidence intervals,
and a homogeneity-of-regression test.

## Author

Ken Kelley

## Examples

``` r
data(pygmalion)
str(pygmalion)
#> 'data.frame':    310 obs. of  6 variables:
#>  $ grade    : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ treatment: Factor w/ 2 levels "Control","Bloomer": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ iq_pre   : int  45 75 61 84 65 72 85 100 88 57 ...
#>  $ iq_4     : int  58 62 82 86 76 94 87 70 85 84 ...
#>  $ iq_8     : int  76 85 87 86 78 94 96 95 81 93 ...
#>  $ iq_gain  : int  31 10 26 2 13 22 11 -5 -7 36 ...

# Design: pupils per condition within each grade.
table(pygmalion$treatment, pygmalion$grade)
#>          
#>            1  2  3  4  5  6
#>   Control 45 46 40 47 25 43
#>   Bloomer  7 12 13 12  9 11

# ---- Heterogeneity-of-regression ANCOVA (book Chapter 9) ----
# Separate IQ8-on-IQpre slopes for the two conditions.
fit_het <- lm(iq_8 ~ iq_pre * treatment, data = pygmalion)
coef(fit_het)
#>             (Intercept)                  iq_pre        treatmentBloomer 
#>              30.3658665               0.7779856             -14.8328265 
#> iq_pre:treatmentBloomer 
#>               0.1909591 
# Control slope = 0.778; the interaction (0.191) gives the
# steeper Bloomer slope of 0.969.

# The treatment-by-covariate interaction is the
# heterogeneity-of-regression test (1 df): compare the additive
# ANCOVA model to the separate-slopes model.
fit_add <- lm(iq_8 ~ iq_pre + treatment, data = pygmalion)
anova(fit_add, fit_het)
#> Analysis of Variance Table
#> 
#> Model 1: iq_8 ~ iq_pre + treatment
#> Model 2: iq_8 ~ iq_pre * treatment
#>   Res.Df   RSS Df Sum of Sq      F  Pr(>F)  
#> 1    307 54330                              
#> 2    306 53649  1    680.14 3.8793 0.04979 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Pooled within-group residual variance (175.3251) and the
# covariate variance (348.91), as used by MBESS::var.ete.
sum(residuals(fit_het)^2) / fit_het$df.residual
#> [1] 175.3251
var(pygmalion$iq_pre)
#> [1] 348.9097

# ---- DMAR's ANCOVA, with the homogeneity-of-regression check ----
ancova(pygmalion, outcome = "iq_8", treatment = "treatment",
       covariates = "iq_pre")
#>  term                         value   
#>  F_value                      5.38    
#>  df_1                         1       
#>  df_2                         307     
#>  p_value                      0.0210  
#>  sum_of_squares_type          3       
#>  eta_squared_partial          0.0172  
#>  eta_squared_partial_lower    0.000201
#>  eta_squared_partial_upper    0.056   
#>  omega_squared_partial        0.0139  
#>  omega_squared_partial_lower  0.000201
#>  omega_squared_partial_upper  0.056   
#>  adjusted_mean[Control]       107     
#>  adjusted_mean[Bloomer]       111     
#>  se_adjusted_mean[Control]    0.849   
#>  se_adjusted_mean[Bloomer]    1.67    
#>  F_homogeneity_of_regression  3.88    
#>  df_homogeneity_of_regression 1       
#>  p_homogeneity_of_regression  0.0498  
#> 
#> Confidence level: 95%
```
