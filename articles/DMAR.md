# Introduction to DMAR

## What Is DMAR?

DMAR is a more modern, more general, and greatly expanded reimagining of
the **MBESS** package (Kelley, 2007a, *Journal of Statistical Software*;
2007b, *Behavior Research Methods*). MBESS itself was originally framed
as *Methods for the Behavioral, Educational, and Social Sciences*, but
the methods it implements became useful far beyond that scope; DMAR
ships under a new name to make that scope expansion explicit. Its
emphasis is the methodology that the substantive disciplines depend on
most: confidence intervals for effect sizes built from noncentral
sampling distributions, sample size planning for both power and accuracy
in parameter estimation (AIPE), reliability, factor analysis, and the
critical values used in classical inference. DMAR is intended to serve
as the computational companion to Maxwell, Delaney, and Kelley’s
*Designing Experiments and Analyzing Data* (4th ed., 2027; henceforth
**MDK**), and it is suitable for both instruction and substantive
research.

Although the package’s origins are in psychology, the questions it
answers are not field-specific. Researchers in biostatistics
(clinical-trial effect size CIs and AIPE-based trial sizing),
organizational behavior and management science (intervention, training,
and program evaluation), information systems (user-experience A/B
comparisons; technology-adoption studies), education (instructional
effectiveness, longitudinal achievement gains), sociology (program
evaluation, attitudinal surveys), and the methodological literature
itself will find tools here for the same three tasks: estimating effect
sizes, attaching confidence intervals that hold their stated coverage,
and planning future studies that will not be underpowered or
underestimated.

What DMAR keeps from MBESS is the philosophy, that is, the same
insistence on *effect sizes with uncertainty quantified*, and on
planning studies for *accurate estimation* rather than only for
sufficient power. What DMAR changes is the implementation: a uniform
`snake_case` interface, tidy `data.frame` returns with stable column
schemas, native `ggplot2` visualization, dependencies that follow
current best practice, and a test suite. The function families have also
been broadened, particularly the `ss_power_*` family, which now covers
between-subjects, within-subjects, mixed, and multi-level designs across
the chapters of MDK.

## Why Effect Size Confidence Intervals?

A *p*-value tells you whether an effect is plausibly nonzero in the
direction you tested. It does not tell you how *large* the effect is,
and it does not tell you how *precisely* you have estimated it. Both of
those questions matter scientifically: a study can produce a
“significant” but trivially small effect, or a “nonsignificant” effect
whose confidence interval is so wide that it is consistent with both a
small and a large true value. Reporting an effect size *and* a
confidence interval for it addresses both shortcomings, and the practice
is now standard in most quantitative areas of psychology and education
(Kelley & Preacher, 2012).

The technical machinery underneath is less familiar than it should be.
For most standardized effect sizes, the sampling distribution is
*noncentral* (e.g., the noncentral *t* underlying Cohen’s $`d`$, the
noncentral *F* underlying $`R^2`$, the noncentral $`\chi^2`$ underlying
RMSEA), which means the correct interval is asymmetric about the point
estimate. The usual symmetric “estimate $`\pm`$ standard error” shortcut
ignores that skew, so it puts the limits in the wrong places. DMAR
computes the correct intervals from the noncentral distributions
directly. This is the common engine behind every `ci_*()` function in
the package.

## Two Philosophies of Sample Size Planning

Researchers who plan studies in DMAR can do so from either of two
complementary perspectives:

- **Power.** Functions in the `ss_power_*` family answer the question
  *“How many subjects do I need so that, if the effect in the population
  is what I expect, I have at least probability $`\pi`$ of rejecting the
  null?”* This is the familiar Cohen-style framing.
- **Accuracy in parameter estimation (AIPE).** Functions in the
  `ss_aipe_*` family answer a different question: *“How many subjects do
  I need so that the confidence interval around my effect size estimate
  is no wider than $`w`$?”* AIPE is appropriate when estimation is the
  primary goal, for example, when an effect is already known to be
  nonzero and the next study should pin down its magnitude rather than
  re-establish significance.

It is worth being precise about what *accuracy* means in this framework,
because the term is sometimes used loosely. Accuracy considers precision
and bias together: an accurate estimator is one whose sampling
distribution is concentrated *and* centered on the target parameter.
Precision alone, in contrast, can always be improved by trading bias for
precision. A trivial estimator that returns the same constant value for
every sample has zero variance and is therefore arbitrarily “precise,”
but it is biased unless the constant happens to equal the population
parameter. The AIPE framing asks for sample sizes that buy narrow
confidence intervals (a precision target) *without* trading bias to get
there. The confidence interval width target $`w`$ in `ss_aipe_*` is
shorthand for that joint accuracy goal: a CI that is informative because
it is both narrow and centered on the parameter (Kelley & Rausch, 2006;
Kelley & Maxwell, 2003; Maxwell, Kelley, & Rausch, 2008).

The two perspectives often imply different sample sizes for the same
study. An AIPE plan tends to require more subjects than a power plan
when the expected effect is moderate, because narrow intervals are more
demanding than null rejection. DMAR provides both side by side so that
the choice can be made deliberately, in light of the scientific
question.

## When Do You Reach for DMAR?

The package is organized around five families of tasks. Pick the one
that matches the question you are trying to answer.

### 1. *I Just Ran an ANOVA, t-Test, or Regression. What’s the Effect Size and a CI for It?*

Use the `ci_*` family.

``` r

# Confidence interval for omega squared (proportion of variance accounted
# for) from a one-way ANOVA with F = 11.221 on (4, 50) df, total N = 55.
ci_omega_squared(F_value = 11.221, df_effect = 4, df_error = 50, N = 55)
```

| effect  | omega_squared | lower_limit | upper_limit | F_value | df_effect | df_error | N   |
|:--------|:--------------|:------------|:------------|:--------|:----------|:---------|:----|
| overall | 0.426         | 0.226       | 0.587       | 11.2    | 4         | 50       | 55  |

``` r


# More naturally, from a fitted model:
fit <- aov(weight ~ group, data = PlantGrowth)
ci_omega_squared(fit)
```

| effect | omega_squared | lower_limit | upper_limit | F_value | df_effect | df_error | N   |
|:-------|:--------------|:------------|:------------|:--------|:----------|:---------|:----|
| group  | 0.204         | 0.0099      | 0.464       | 4.85    | 2         | 27       | 30  |

The output is a tidy `data.frame`: an effect label, the point estimate
(`omega_squared`), the lower and upper limits of the CI, and the inputs
that produced it (`F_value`, `df_effect`, `df_error`, and *N*). That
format makes it easy to feed the result directly into a table or a
forest plot
([`plot_ci()`](https://yelleknek.github.io/DMAR/reference/plot_ci.md) is
built to consume it).

CI variants exist for many other effect sizes:
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md) for
the standardized mean difference (Cohen’s $`d`$);
[`ci_smd_c()`](https://yelleknek.github.io/DMAR/reference/ci_smd_c.md)
for the standardized mean difference using the control-group SD as the
divisor;
[`ci_pvaf()`](https://yelleknek.github.io/DMAR/reference/ci_pvaf.md) for
the proportion of variance accounted for by a single predictor;
[`ci_snr()`](https://yelleknek.github.io/DMAR/reference/ci_snr.md) and
[`ci_srsnr()`](https://yelleknek.github.io/DMAR/reference/ci_srsnr.md)
for signal-to-noise ratios;
[`ci_R()`](https://yelleknek.github.io/DMAR/reference/ci_correlation.md),
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md), and
[`ci_rc()`](https://yelleknek.github.io/DMAR/reference/ci_rc.md) for
correlations and standardized regression coefficients;
[`ci_rmsea()`](https://yelleknek.github.io/DMAR/reference/ci_rmsea.md)
for the RMSEA model-fit index; and others.

### 2. *I’m Planning a Study. How Many Subjects Do I Need?*

Use `ss_power_*` for sufficient power, or `ss_aipe_*` for sufficient
precision. The same effect can be approached from either angle.

``` r

# Power-based: per-group n needed for 90% power to detect a contrast among
# 4 groups with population means (90, 92, 88, 81) and within-group
# variance 144. The contrast compares the average of the first three
# groups against the fourth.
ss_power_contrast(
  c_weights     = c(1/3, 1/3, 1/3, -1),
  mu            = c(90, 92, 88, 81),
  sigma_squared = 144,
  desired_power = 0.90
)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 26    |
| total_N               | 104   |
| actual_power          | 0.907 |
| noncentral_t_parm     | 3.31  |
| effect_size_f         | 0.325 |

``` r

# AIPE-based: n needed so a 95% CI around a standardized mean difference
# (delta) of 0.50 is no wider than 0.30 units (a tight interval).
ss_aipe_smd(delta = 0.50, conf_level = 0.95, width = 0.30)
```

| term                  | value |
|:----------------------|:------|
| necessary_n_per_group | 353   |
| supposed_smd          | 0.5   |
| width                 | 0.3   |

Confidence level: 95%

The `ss_power_*` family in v1.0.0 spans most of the standard designs
covered in MDK: one-way ANOVA (`ss_power_one_way_anova`), factorial
ANOVA (`ss_power_factorial_anova`), repeated measures ANOVA with
sphericity adjustment (`ss_power_rm_anova`), split-plot mixed ANOVA
(`ss_power_split_plot_anova`), Pearson correlation (`ss_power_r`),
two-independent-groups standardized mean difference (`ss_power_smd`),
contrasts standardized and unstandardized (`ss_power_c`, `ss_power_sc`,
`ss_power_c_ancova`, `ss_power_contrast`), regression coefficients
(`ss_power_reg_coef`), and two-level random-intercept models for
cluster-randomized trials (`ss_power_mixed_effects`). For SEM-style
designs, use
[`ss_power_sem()`](https://yelleknek.github.io/DMAR/reference/ss_power_sem.md).
AIPE counterparts exist for the correspondingly named effects.

### 3. *I Want to Test a Specific Contrast or Set of Pairwise Comparisons.*

Use
[`contrast_test()`](https://yelleknek.github.io/DMAR/reference/contrast_test.md)
for arbitrary contrasts in a one-way ANOVA, with optional
multiple-comparison adjustment.

``` r

fit <- aov(weight ~ group, data = PlantGrowth)
# All pairwise comparisons, Tukey-protected family-wise error rate:
contrast_test(fit, contrasts = "pairwise", adjust = "tukey")
```

| contrast    | estimate | se    | t     | df  | p_value | p_adjusted | ci_lower | ci_upper |
|:------------|:---------|:------|:------|:----|:--------|:-----------|:---------|:---------|
| trt1 - ctrl | -0.371   | 0.279 | -1.33 | 27  | 0.1944  | 0.3909     | -1.06    | 0.32     |
| trt2 - ctrl | 0.494    | 0.279 | 1.77  | 27  | 0.0877  | 0.1980     | -0.197   | 1.19     |
| trt2 - trt1 | 0.865    | 0.279 | 3.1   | 27  | 0.0045  | 0.0120     | 0.174    | 1.56     |

Confidence level: 95%

For the bare critical value of a particular procedure, see the `cv_*`
family: [`cv_t()`](https://yelleknek.github.io/DMAR/reference/cv_t.md),
[`cv_tukey_hsd()`](https://yelleknek.github.io/DMAR/reference/cv_tukey_hsd.md),
[`cv_scheffe()`](https://yelleknek.github.io/DMAR/reference/cv_scheffe.md),
[`cv_smm()`](https://yelleknek.github.io/DMAR/reference/cv_smm.md) (the
studentized maximum modulus),
[`cv_dunnett()`](https://yelleknek.github.io/DMAR/reference/cv_dunnett.md),
and [`cv_z()`](https://yelleknek.github.io/DMAR/reference/cv_z.md).
These return the threshold value alone, suitable for inclusion in user
code or hand-checking textbook examples.

### 4. *I Have Within-Subjects (Repeated Measures) Data.*

Start with
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md),
which combines the univariate within-subjects *F* test with a sphericity
diagnostic (Mauchly’s *W*) and the three standard
$`\varepsilon`$-corrected *p*-values (Greenhouse-Geisser, Huynh-Feldt,
lower-bound). The motivation is that an uncorrected within-subjects *F*
test can be substantially anticonservative when the sphericity
assumption fails, and
[`anova_within()`](https://yelleknek.github.io/DMAR/reference/anova_within.md)
puts the diagnostic and the correction in one place.

``` r

if (requireNamespace("nlme", quietly = TRUE)) {
  res <- anova_within(nlme::Orthodont,
                      id = "Subject", time = "age", outcome = "distance")
  res
}
```

| adjustment         | F_value | df_1 | df_2 | p_value   | epsilon |
|:-------------------|:--------|:-----|:-----|:----------|:--------|
| none               | 38      | 3    | 78   | \< 0.0001 | NA      |
| Greenhouse-Geisser | 38      | 2.63 | 68.4 | \< 0.0001 | 0.877   |
| Huynh-Feldt        | 38      | 2.95 | 76.8 | \< 0.0001 | 0.984   |
| lower_bound        | 38      | 1    | 26   | \< 0.0001 | 0.333   |

For the underlying components separately, see
[`mauchly_test()`](https://yelleknek.github.io/DMAR/reference/mauchly_test.md)
and
[`epsilon_corrections()`](https://yelleknek.github.io/DMAR/reference/epsilon_corrections.md).
For visualizing individual subject trajectories, see
[`plot_trajectories()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories.md)
and
[`plot_trajectories_fitted()`](https://yelleknek.github.io/DMAR/reference/plot_trajectories_fitted.md).

### 5. *I Want Descriptive Statistics, a Correlation Table, or Distribution-Shape Diagnostics.*

[`descriptives()`](https://yelleknek.github.io/DMAR/reference/descriptives.md)
returns a tidy summary suitable for screening:

``` r

descriptives(attitude)$descriptives
#>     variable    type  n n_missing prop_missing     mean median        sd min
#> 1     rating numeric 30         0            0 64.63333   65.5 12.172562  40
#> 2 complaints numeric 30         0            0 66.60000   65.0 13.314757  37
#> 3 privileges numeric 30         0            0 53.13333   51.5 12.235430  30
#> 4   learning numeric 30         0            0 56.36667   56.5 11.737013  34
#> 5     raises numeric 30         0            0 64.63333   63.5 10.397226  43
#> 6   critical numeric 30         0            0 74.76667   77.5  9.894908  49
#> 7    advance numeric 30         0            0 42.93333   41.0 10.288706  25
#>   max   q25   q75   skewness    kurtosis
#> 1  85 58.75 71.75 -0.3967148 -0.49460977
#> 2  90 58.50 77.00 -0.2387632 -0.38172518
#> 3  83 45.00 62.50  0.4202101 -0.04219127
#> 4  75 47.00 66.75 -0.0598894 -1.07638294
#> 5  88 58.25 71.00  0.2189518 -0.28201359
#> 6  92 69.25 80.00 -0.9596072  0.69175815
#> 7  72 35.00 47.75  0.9425594  1.07314413
```

Adding `correlations = TRUE` appends a Pearson correlation matrix. For a
publication-ready correlation table with *p*-values, confidence
intervals, and significance stars (and, optionally, HTML or LaTeX
export), use
[`correlations_test()`](https://yelleknek.github.io/DMAR/reference/correlations_test.md):

``` r

correlations_test(attitude, stars = TRUE)
#> Correlations (Pearson, 95% CI)
#> 
#>               rating        complaints    privileges    learning      raises        critical      advance       
#> ----------------------------------------------------------------------------------------------------------------
#> rating        -                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#>                                                                                                                 
#> 
#> complaints    .83***        -                                                                                   
#>               p < .0001                                                                                         
#>               [.66, .91]                                                                                        
#>               N = 30                                                                                            
#> 
#> privileges    .43*          .56**         -                                                                     
#>               p = .0189     p = .0013                                                                           
#>               [.08, .68]    [.25, .76]                                                                          
#>               N = 30        N = 30                                                                              
#> 
#> learning      .62***        .60***        .49**         -                                                       
#>               p = .0002     p = .0005     p = .0056                                                             
#>               [.34, .80]    [.30, .79]    [.16, .72]                                                            
#>               N = 30        N = 30        N = 30                                                                
#> 
#> raises        .59***        .67***        .45*          .64***        -                                         
#>               p = .0006     p < .0001     p = .0136     p = .0001                                               
#>               [.29, .78]    [.41, .83]    [.10, .69]    [.36, .81]                                              
#>               N = 30        N = 30        N = 30        N = 30                                                  
#> 
#> critical      .16           .19           .15           .12           .38*          -                           
#>               p = .4091     p = .3205     p = .4375     p = .5417     p = .0401                                 
#>               [-.22, .49]   [-.19, .51]   [-.22, .48]   [-.25, .46]   [.02, .65]                                
#>               N = 30        N = 30        N = 30        N = 30        N = 30                                    
#> 
#> advance       .16           .22           .34           .53**         .57***        .28           -             
#>               p = .4132     p = .2328     p = .0633     p = .0025     p = .0009     p = .1292                   
#>               [-.22, .49]   [-.15, .54]   [-.02, .63]   [.21, .75]    [.27, .77]    [-.09, .58]                 
#>               N = 30        N = 30        N = 30        N = 30        N = 30        N = 30                      
#> 
#> Note. * p < .05, ** p < .01, *** p < .001.
```

For interrater agreement, use
[`cohen_kappa()`](https://yelleknek.github.io/DMAR/reference/cohen_kappa.md)
(two raters),
[`fleiss_kappa()`](https://yelleknek.github.io/DMAR/reference/fleiss_kappa.md)
(multiple raters), or
[`icc()`](https://yelleknek.github.io/DMAR/reference/icc.md) (intraclass
correlation, all six classical forms). For distribution-shape statistics
in isolation, use
[`skewness()`](https://yelleknek.github.io/DMAR/reference/skewness.md)
and
[`kurtosis()`](https://yelleknek.github.io/DMAR/reference/kurtosis.md).

## What DMAR Is *Not*

DMAR deliberately overlaps with, but does not duplicate, several
specialized packages.

- For full-featured factor analysis and SEM, use `lavaan`.
- For coefficient $`\alpha`$ / $`\omega`$ with extensive options, use
  `psych::omega()`.
- For mixed-effects models, use `lme4` or `nlme`.
- For robust ANOVA, use `WRS2`.

DMAR layers tidy effect size CIs, sample size planning, and
publication-ready output over those tools, rather than replacing them.

## A Note on Argument Naming

The package uses `snake_case` throughout. Capital letters in identifiers
are reserved for statistically meaningful symbols, for example `R2` for
the squared multiple correlation, `N` for total sample size, `S` for a
covariance matrix, `Lambda` for a factor-loadings matrix, `F_value` for
an *F* statistic. Every other argument is lowercase, with the result
that examples read consistently regardless of discipline:

``` r

ci_R2(R2 = 0.25, N = 100, p = 5, random_predictors = TRUE)
ss_power_reg_coef(rho2_Y_X = 0.78, rho2_Y_X_without_j = 0.74,
                  p = 5, desired_power = 0.85)
ss_aipe_cv_sensitivity(true_cv = 0.25, estimated_cv = 0.25,
                       width = 0.10, conf_level = 0.95, G = 200)
```

The MBESS dot.case style (`conf.level`, `Random.Predictors`,
`Specified.N`) is no longer accepted. Scripts written against MBESS will
need to convert dots to underscores in argument names and lower-case any
prefix words that are not themselves statistical notation. Function
names are unchanged in nearly all cases
([`smd()`](https://yelleknek.github.io/DMAR/reference/smd.md),
[`ci_smd()`](https://yelleknek.github.io/DMAR/reference/ci_smd.md),
[`ci_R2()`](https://yelleknek.github.io/DMAR/reference/ci_R2.md),
[`ss_aipe_reg_coef()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_reg_coef.md),
etc.), and *return* shapes are unchanged or strictly richer.

If you have scripts written against MBESS or against pre-1.0 DMAR
snapshots and you see errors of the form

    Error in foo(...) : unused argument (Group.1 = ...)

the fix is mechanical: replace dots with underscores in the argument
name (`Group_1`, `n_1`, `Mean_1`, `Specified_N`, etc.) and then
lowercase the non-meaningful capitalized prefix words (`group_1`,
`mean_1`, `specified_N`, …). See `NEWS.md` for the full mapping. The
exception worth noting on the function side is `aipe_smd()`, which has
been renamed to
[`ss_aipe_smd()`](https://yelleknek.github.io/DMAR/reference/ss_aipe_smd.md)
to align it with the rest of the `ss_aipe_*` family.

The motivation for the unification is simple: the former mix of
`Random.Predictors` (dot.case for legacy MBESS arguments) and
`conf_level` (snake_case for newer arguments) was confusing in
side-by-side examples, especially in teaching, and a single convention
removes the source of friction.

## Where to Next

- For a worked walkthrough of the visualization functions and effect
  size confidence intervals on a simulated psychometric study, see
  [`vignette("effect-size-visualization", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/effect-size-visualization.md).
- For when full information maximum likelihood regression actually helps
  over [`lm()`](https://rdrr.io/r/stats/lm.html), with
  [`mlmr()`](https://yelleknek.github.io/DMAR/reference/mlmr.md) and
  [`mlmr_mv()`](https://yelleknek.github.io/DMAR/reference/mlmr_mv.md),
  see
  [`vignette("mlmr", package = "DMAR")`](https://yelleknek.github.io/DMAR/articles/mlmr.md).
- For the function reference, `?DMAR-package` lists the exported
  functions grouped by family.
- Citation: `citation("DMAR")`.
- Author web site: <https://kenkelley.org>; related publications:
  <https://kenkelley.org/publications/>.

## References

Kelley, K. (2007a). Confidence intervals for standardized effect sizes:
Theory, application, and implementation. *Journal of Statistical
Software, 20*(8), 1–24. <https://doi.org/10.18637/jss.v020.i08>

Kelley, K. (2007b). Methods for the behavioral, educational, and social
sciences: An R package. *Behavior Research Methods, 39*(4), 979–984.
<https://doi.org/10.3758/BF03192993>

Kelley, K., & Maxwell, S. E. (2003). Sample size for multiple
regression: Obtaining regression coefficients that are accurate, not
simply significant. *Psychological Methods, 8*(3), 305–321.

Kelley, K., & Preacher, K. J. (2012). On effect size. *Psychological
Methods, 17*(2), 137–152.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the
standardized mean difference: Accuracy in parameter estimation via
narrow confidence intervals. *Psychological Methods, 11*(4), 363–385.

Maxwell, S. E., Delaney, H. D., & Kelley, K. (2027). *Designing
experiments and analyzing data: A model comparison perspective* (4th
ed.). Routledge.

Maxwell, S. E., Kelley, K., & Rausch, J. R. (2008). Sample size planning
for statistical power and accuracy in parameter estimation. *Annual
Review of Psychology, 59*, 537–563.
<https://doi.org/10.1146/annurev.psych.59.103006.093735>
